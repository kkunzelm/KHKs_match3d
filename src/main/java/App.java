
import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * File open -> ZIP -> save
 *
 */
public class App {
	
	private class SelectedFile {
		private String fileName;
		private byte[] rawData;
	}

	private List<SelectedFile> filesToZip = new ArrayList<SelectedFile>();
	
	public static void main(String[] args) {
        Runnable r = new Runnable() {
            @Override
            public void run() {
                new App().createUI();
            }
        };
        EventQueue.invokeLater(r);
    }

	/**
	 * Creates an zip archive with all entries from List<SelectedFile> filesToZip
	 */
	private void createZipFile(String pathWithFilename) {
		byte[] buffer = new byte[1024 * 4];

		pathWithFilename = pathWithFilename.toLowerCase();
		if(!pathWithFilename.endsWith(".zip")){
			pathWithFilename += ".zip";
		}
		
		try {
			String outputFileName = pathWithFilename;
			FileOutputStream fos = new FileOutputStream(outputFileName);
			ZipOutputStream zos = new ZipOutputStream(fos);

			for (SelectedFile file : filesToZip) {
				ZipEntry ze = new ZipEntry(file.fileName);
				zos.putNextEntry(ze);

				ByteArrayInputStream bis = new ByteArrayInputStream(file.rawData);

				int len;
				while ((len = bis.read(buffer)) > 0) {
					zos.write(buffer, 0, len);
				}
				bis.close();
			}
			zos.closeEntry();
			zos.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		JOptionPane.showMessageDialog(null, "Die Datei(n) wurde(n) gezipt...",
				"Info", JOptionPane.INFORMATION_MESSAGE);
	}
	
    private void createUI() {
        final JFrame frame = new JFrame();
        frame.setLayout(new BorderLayout());

        JButton btnOpen = new JButton("Datei öffnen");
        JButton btnSave = new JButton("ZIP speichern");

        btnOpen.addActionListener(new ActionListener() {
        	@Override
        	public void actionPerformed(ActionEvent ae) {
        		JFileChooser openDialog = new JFileChooser();
        		int openDialogReturn = openDialog.showOpenDialog(frame);
        		openDialog.setMultiSelectionEnabled(false);
        	        if(openDialogReturn == JFileChooser.APPROVE_OPTION)
        	        {
        	        	SelectedFile sf = new SelectedFile();
        	        	sf.fileName = openDialog.getSelectedFile().getName();
        	        	String pathToFile = openDialog.getSelectedFile().getAbsolutePath();
        	        	try {
							sf.rawData = Files.readAllBytes(Paths.get(pathToFile));
						} catch (IOException e) {
							e.printStackTrace();
						}
        	        	filesToZip.add(sf);
        	        }
        	}
        });

		btnSave.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent ae) {
				JFileChooser saveDialog = new JFileChooser();
			    FileNameExtensionFilter zipExtensionFilter =
			            new FileNameExtensionFilter("ZIP File(*.zip)", "zip");
			    saveDialog.addChoosableFileFilter(zipExtensionFilter);
			    int saveDialogReturn = saveDialog.showSaveDialog(frame);
				
			    if (saveDialogReturn == JFileChooser.APPROVE_OPTION) {
					createZipFile(saveDialog.getSelectedFile().getAbsolutePath());
				}
			}
		});

        frame.add(new JLabel("Bitte zu zippende Datei auswählen:"), BorderLayout.NORTH);
        frame.add(btnOpen, BorderLayout.CENTER);
        frame.add(btnSave, BorderLayout.SOUTH);
        frame.setTitle("ZIP Files");
        frame.pack();
        frame.setLocationRelativeTo(null); // Desktop mittig
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
}