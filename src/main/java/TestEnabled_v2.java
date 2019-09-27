
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import static ij.plugin.filter.PlugInFilter.DOES_32;
import imagescience.transform.Transform;



/**
 * Some experiments with GUI in ImageJ
 * see also Tutorial: "Getting Started with Swing" from Sun
 * 
 * @author Karl-Heinz Kunzelmann (karl-heinz@kunzelmann.de)
 * 
 */


//
public class TestEnabled_v2 implements PlugIn {
    
    
        private static final String[] schemes = {
		"refine_clamp",
                "refine_clip",
		"refine_sparse",
		"refine_unique"
	};
	private static int scheme = 0;                  //default = clamp
    
        protected boolean refine_clamp = true;          // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
        protected double refine_clamp_par = 3.0;
        
        protected boolean refine_clip = false;          // clipping ab bestimmter Position im Abstandshistogram, 
        protected double refine_clip_par = 0.75;         // z. B. nur 75 % der niedrigeren Abstände berücksichtigt, 
                                                        // Rest der Punktepaare verworfen

        protected boolean refine_sparse = false;        // nur jeder  xte Punkte wird verwendet
        protected double refine_sparse_par = 0.5;
        
        protected boolean refine_unique = false;        // klingt als ob jeder Punkt berücksichtigt wird
   
     
        public String title1 = "";
        public String title2 = "";
        String[] titles;
    
        ImagePlus imp1;                     //Mnemo: imp = imageplus
        ImagePlus imp2;

        int[] wList;

       
	public int setup(String arg, ImagePlus imp) {
            return DOES_32;
	}

           
    public void run(String arg) {
	 // check: is there an open image?

        wList = WindowManager.getIDList();
        if (wList==null) {
                IJ.noImage();
                return;
        }

        // check: we need at least 2 images for matching

        titles = new String[wList.length];
        if (wList.length < 2) {
                IJ.error(
                        "At least two match3d images (type float or Gray32) are required");
                return;
        }


        // make a list of image titles of all open images

        for (int i=0; i<wList.length; i++) {
                ImagePlus impTemp = WindowManager.getImage(wList[i]);
                if (impTemp!=null) {
                    titles[i] = impTemp.getTitle();
                }
                else {
                    titles[i] = "";
                }
        }
        

        // user interaction to get matching parameters
        if (!showDialog()) {
            return;
        }
         
        
        // **************************************************************************


            long startTime = System.currentTimeMillis();

            ///   (new modifiedTJAffine()).run(imp,file,scheme);
            // DO SOMETHING 
            
            
            long stopTime = System.currentTimeMillis();
            long elapsedTime = stopTime - startTime;
            System.out.println("Time to perform the target-to-source transformation: " + elapsedTime + "ms");
             
	}


/**
 * Show plugin configuration dialog.
 *
 * @return <code>true</code> when user clicked OK (confirmed changes, <code>false</code>
 *         otherwise.
 */
public boolean showDialog() {
    
        GenericDialog gd = new GenericDialog("KHKs jMatch3D", IJ.getInstance());
                
        gd.addMessage("Two images of file type '32bit' or 'float' are required. \nx and y are coordinates on a rectangular grid, \nz represents the height information z = f(x,y).");
        
        String defaultItem;
        if (title1.equals(""))
                defaultItem = titles[0];
        else
                defaultItem = title1;
        gd.addChoice("Image1 (source):", titles, defaultItem);

        if (title2.equals(""))
                defaultItem = titles[0];
        else
                defaultItem = title2;
        gd.addChoice("Image2 (target):", titles, defaultItem);
        
        gd.addMessage("");
        
        gd.addChoice("ICP Point Selection Method: ",schemes,schemes[scheme]);
        

        gd.addMessage("refine_clamp:\tremoves points further away than multiple of mean/median distance.\n"
                    + "refine_clip:\tremoves percentage of points with highest distance.\n"
                    + "refine_sparse:\tremoves percentage of all points\n"
                    + "refine_unique:\tall nearest points (very slow)\n");
        gd.addMessage("");
        
        // refine_clamp = kill all points farther away than a multiple of
        // the last reported mean distance and a multiple of
        // the current median distance
        // refine_clip =    kill percentage of points with highest distance
        //                  (approximated value from non-transformed distance values)
        //  refine_sparse = kill a certain percentage of points.
        //  refine_unique = all nearest points (very slow).
        gd.addMessage("Parameters: ");
        gd.addNumericField("refine_clamp multiple of mean/med error: ", refine_clamp_par, 2);
        gd.addNumericField("refine_clip percentage (0.0 - 1.0): ", refine_clip_par, 2);
        gd.addNumericField("refine_sparce percentage (0.0 - 1.0): ", refine_sparse_par, 2);
        
        gd.addMessage("");
        gd.addMessage("Developed by Karl-Heinz Kunzelmann.\nBased on plenty of code from the Internet\n");

        
        gd.showDialog();

        if (gd.wasCanceled())
                return false;


        int index1 = gd.getNextChoiceIndex();
        title1 = titles[index1];

        int index2 = gd.getNextChoiceIndex();
        title2 = titles[index2];
        
        imp1 = WindowManager.getImage(wList[index1]);
        imp2 = WindowManager.getImage(wList[index2]);	
        

        // parameters to initialize the icp algorithm     

        scheme = gd.getNextChoiceIndex();
        
        refine_clamp_par = gd.getNextNumber();
        refine_clip_par = gd.getNextNumber();
        refine_sparse_par = gd.getNextNumber();
        
        // plausibility check
        if (refine_clip_par < 0.0 || refine_clip_par > 1.0){
            refine_clip_par = 0.75;         // default
        }
        
        if (refine_sparse_par < 0.0 || refine_sparse_par > 1.0){
            refine_sparse_par = 0.5;
        }
        
        switch(scheme){
            case 0: 
                refine_clamp = true;
                refine_clip = false;
                refine_sparse = false;
                refine_unique = false;
                break;
            case 1:
                refine_clamp = false;
                refine_clip = true;
                refine_sparse = false;
                refine_unique = false;
                break;
            case 2:
                refine_clamp = false;
                refine_clip = false;
                refine_sparse = true;
                refine_unique = false;
                break;
            case 3:
                refine_clamp = false;
                refine_clip = false;
                refine_sparse = false;
                refine_unique = true;
                break;
        }
        
        System.out.println(  "refine_clamp: " + refine_clamp + " = " + refine_clamp_par + "\n" +
                             "refine_clip: " + refine_clip + " = " + refine_clip_par +"\n" +
                             "refine_sparse: " + refine_sparse + " = " + refine_sparse_par +"\n" +
                             "refine_unique: " + refine_unique);
        return true;

        }
    
}
        
class KHKDoSomething {
	
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
        
 
