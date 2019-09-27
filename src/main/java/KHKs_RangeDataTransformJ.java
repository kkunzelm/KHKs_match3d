import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.transform.Transform;
import java.awt.Button;
import java.awt.FileDialog;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Panel;
import java.awt.Point;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;

/**
 * Performs an affine transformation of range data. 
 * range data = synonym: 2.5D data, dem as in digital elevation model
 * range data = x,y,z coordinates organized in image files with x = column, y = row and z = f(x,y) = value at position of (x,y)
 * The current implementation does not extend the image size. 
 * For my proposes it was necessary to keep the size the same as the original image.
 * But take care. A rotation might easily move the data outside the image boundaries, which means that you will see nothing.
 * So far I did not take care of hidden surfaces: my data will be rotated only a little. 
 * This has to be done for a more generic version, however.
 * The file format of the matrix file is simple:
 * the data are written in four lines corresponding to the matrix rows, the columns are separated with \t (tab)
 * The files are identical with the files generated with TransformJ's matrix routine, which you can use to 
 * make the matrix files.
 * 
 * @author Karl-Heinz Kunzelmann (karl-heinz@kunzelmann.de)
 * 
 */

// This is the frontend, which calles the class KHKsRangeDataTransformJ, which does the real work.
// The code of the frontend is copied from TransformJ's TJ_Affine.java (http://www.imagescience.org/meijering/software/transformj/).
// Only a few minor modifications were made.
// Thank you, Prof. Erik Meijering, for this beautiful piece of code!

// KHK todo save matrix dialog... does not work flawless when a fully qualified path is entered in the filename menu.
//     always need a click on a directory
public class KHKs_RangeDataTransformJ implements PlugInFilter, ActionListener, WindowListener {
	
	private static String file = "";
        private static final Point pos = new Point(-1,-1);
        ImagePlus imp;

        private static final String[] schemes = {
		"nearest neighbor",
		"linear (default)",
		"bicubic (ImageJ)",
		"bicubic (Bob Dougherty)"
	};
	private static int scheme = 1;
	
        private Button browseButton, createButton;
	private TextField fileField;
        
        @Override
	public int setup(String arg, ImagePlus imp) {
            return DOES_32;
	}

            @Override
    public void run(ImageProcessor ip) {
		
            imp = WindowManager.getCurrentImage();
           
            
            GenericDialog gd = new GenericDialog("Affine Transformation 2.5D");
            gd.addStringField("Matrix file:",file,30);
            fileField = (TextField)gd.getStringFields().get(0);

            final Panel buttons = new Panel();
            GridBagLayout bgbl = new GridBagLayout();
            buttons.setLayout(bgbl);
            browseButton = new Button("    Browse    ");
            browseButton.addActionListener(this);
            createButton = new Button("     Create     ");
            createButton.addActionListener(this);
            GridBagConstraints bgbc = new GridBagConstraints();
            bgbc.anchor = GridBagConstraints.WEST;
            bgbc.insets = new Insets(0,0,0,5);
            bgbl.setConstraints(browseButton,bgbc);
            buttons.add(browseButton);
            bgbc.insets = new Insets(0,0,0,0);
            bgbl.setConstraints(createButton,bgbc);
            buttons.add(createButton);
            gd.addPanel(buttons,GridBagConstraints.WEST,new Insets(0,0,20,0));
            bgbl = (GridBagLayout)gd.getLayout();
            bgbc = bgbl.getConstraints(buttons); bgbc.gridx = 1;
            bgbl.setConstraints(buttons,bgbc);

            gd.addChoice("Interpolation scheme:",schemes,schemes[scheme]);

            if (pos.x >= 0 && pos.y >= 0) {
                    gd.centerDialog(false);
                    gd.setLocation(pos);
            } else gd.centerDialog(true);
            gd.addWindowListener(this);
            gd.showDialog();

            if (gd.wasCanceled()) return;

            file = gd.getNextString();
            scheme = gd.getNextChoiceIndex();

            long startTime = System.currentTimeMillis();

            (new modifiedTJAffine()).run(imp,file,scheme);
            
            long stopTime = System.currentTimeMillis();
            long elapsedTime = stopTime - startTime;
            System.out.println("Time to perform the target-to-source transformation: " + elapsedTime + "ms");
             
	}
	
        @Override
	public void actionPerformed(final ActionEvent e) {
		
		if (e.getSource() == browseButton) {
			final FileDialog fdg = new FileDialog(IJ.getInstance(),"Load",FileDialog.LOAD);
			fdg.setFile(""); fdg.setVisible(true);
			final String d = fdg.getDirectory();
			final String f = fdg.getFile();
			fdg.dispose();
			if (d != null && f != null) {
				String path = d + f;
				if (File.separator.equals("\\"))
					path = path.replace('\\','/');
				fileField.setText(path);
			}
		} else if (e.getSource() == createButton) {
                        
			final TJ_Matrix tjm = new TJ_Matrix();
			try { tjm.load(fileField.getText()); }
			catch (Throwable x) { }
			tjm.run("");
			final String path = tjm.saved();
			if (path != null) fileField.setText(path);
		}
	}
	
        @Override
	public void windowActivated(final WindowEvent e) { }
	
        @Override
	public void windowClosed(final WindowEvent e) {
		
		pos.x = e.getWindow().getX();
		pos.y = e.getWindow().getY();
	}
	
        @Override
	public void windowClosing(final WindowEvent e) { }
	
        @Override
	public void windowDeactivated(final WindowEvent e) { }
	
        @Override
	public void windowDeiconified(final WindowEvent e) { }
	
        @Override
	public void windowIconified(final WindowEvent e) { }
	
        @Override
	public void windowOpened(final WindowEvent e) { }

}

class modifiedTJAffine {
	
    	double[][] transformationMatrix = new double[4][4];
    
	void run(final ImagePlus imp, final String file, final int interpolationMethod) {
		
		try {
			if (file == null || file.equals(""))
				throw new IllegalArgumentException("Empty matrix file name");
			final TJ_Matrix tjm = new TJ_Matrix();
			tjm.load(file);
			final Transform tfm = tjm.get();
			
                        transformationMatrix = tfm.get();
                        KHKsRangeDataTransformJ.applyT2sTransformationRangeData(imp, transformationMatrix, interpolationMethod); 
                        System.out.println("Interpolation Method: "+ interpolationMethod);
			
		} catch (OutOfMemoryError e) {
			TJ.error("Not enough memory for this operation");
			
		} catch (UnknownError e) {
			TJ.error("Could not create output image for some reason.\nPossibly there is not enough free memory");
			
		} catch (IllegalArgumentException e) {
			TJ.error(e.getMessage());
			
		}
	}
	
}
