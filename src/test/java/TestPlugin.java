import ij.IJ;
import ij.ImageJ;

import java.net.URL;

public class TestPlugin {
    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = Match3d_withFiducialMarkersAndICPv2_1.class;
        URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
        java.io.File file = new java.io.File(url.toURI());
        System.setProperty("plugins.dir", file.getAbsolutePath());

        // start ImageJ
        new ImageJ();

        // open the Clown sample
//        ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
//        image.show();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}
