/*
Format Specifications

An ASCII STL file begins with the line:

 solid name

where name is an optional string. The file continues with any number of triangles, each represented as follows:

solid name
 facet normal n1 n2 n3
   outer loop
     vertex v11 v12 v13
     vertex v21 v22 v23
     vertex v31 v32 v33
   endloop
 endfacet
 endsolid name

where n1-n3 and v11-v33 are floating point numbers in sign-mantissa'e'-sign-exponent format and concludes with:

 endsolid name

The structure of the format suggests that other possibilities exist (eg Facets with more than one 'loop' or loops with other than three vertices) but in practice, all facets are simple triangles.

White space (spaces, tabs, newlines) may be used anywhere in the file except within numbers or words. The spaces between 'facet' and 'normal' and between 'outer' and 'loop' are required.

The Facet Normal

In both ASCII and binary versions of STL, the facet normal should be a unit vector pointing outwards from the solid object. In most software this may be set to (0,0,0) and the software will automatically calculate a normal based on the order of the triangle vertices using the 'right hand rule'. 
 */

import java.io.*;
import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

// Saves an image described by an ImageProcessor object as a tab-delimited text file. */
// STL File Format
// generated very simple as triangle strips 
// for use with our match3d images
// input image has heigth data encoded as grey value z = I(x,y)
// input data should be on regular grid already 

public class SimpleSTL_Writer implements PlugInFilter {

	private static ImagePlus imp;
	private ImageProcessor ip;
	private Calibration cal;
	private double scaleX;
	private double scaleY;
	private int width;
	private int height;

	public int setup(String arg, ImagePlus imp) {

		return DOES_32;
	}

	public void run(ImageProcessor ip) {
		imp = WindowManager.getCurrentImage();

		FileInfo fi = imp.getFileInfo();
		scaleX = fi.pixelWidth;
		scaleY = fi.pixelHeight;

		this.ip = ip;

		String path = getPath("STL-File", ".stl");

		try {

			// TextEncoder file = new TextEncoder(imp.getProcessor(), cal, precision);
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
			// file.
			write(out);
			out.close();
		} catch (IOException e) {
			showErrorMessage(e);
		}
	}

	private String getPath(String type, String extension) {
		String name = imp.getTitle();
		SaveDialog sd = new SaveDialog("Save as " + type, name, extension);
		name = sd.getFileName();
		if (name == null)
			return null;
		String directory = sd.getDirectory();
		imp.startTiming();
		return directory + name;
	}

	/** Saves the image as a text file. */
	private void write(DataOutputStream out) {

		PrintWriter pw = new PrintWriter(out);
		boolean calibrated = cal != null && cal.calibrated(); // calibrated wird im MOment nicht verwendet. Ich gehe den
																// Umweg über FileInfo und scaleX/Y
		if (calibrated)
			ip.setCalibrationTable(cal.getCTable());
		else
			ip.setCalibrationTable(null);
		width = ip.getWidth();
		height = ip.getHeight();

		// IJ.showStatus("Exporting as text...");

		// Pnt3d = Point class defined at the end of this class file, similar to
		// vecmath.Vector3d etc.
		Pnt3d[] validPoints = readVerticesForObjectFile();

		pw.println("solid khk");

		for (int p = 0; p < validPoints.length; p = p + 3) {
			Pnt3d a;

			a = dotProduct(validPoints[p], validPoints[p + 1], validPoints[p + 2]);
			pw.println("facet normal " + a.x + " " + a.y + " " + a.z);
			pw.println("outer loop");
			// todo check: if scaleX and Y are different, the it could be that I have to
			// swap scaleX and scaleY here --- see comment line 260ff
			pw.println("vertex " + validPoints[p].x + " " + validPoints[p].y + " " + validPoints[p].z);
			pw.println("vertex " + validPoints[p + 1].x + " " + validPoints[p + 1].y + " " + validPoints[p + 1].z);
			pw.println("vertex " + validPoints[p + 2].x + " " + validPoints[p + 2].y + " " + validPoints[p + 2].z);
			pw.println("endloop");
			pw.println("endfacet");

		}
		pw.println("endsolid khk");
		pw.println();

		pw.close();
		IJ.showProgress(1.0);
		// IJ.showStatus("");
	}
	private void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length() > 100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured writing the file.\n \n" + msg);
	}

	private Pnt3d[] readVerticesForObjectFile() {

		// mentally arrange the lines like:

		// first line: 0 2 4 6 8 ...
		// second line: 1 3 5 7 9 ...
		// etc.
		//
		// Then make Triangles always in the same rotation direction - clockwise or
		// alternatively counterclockwise (here clockwise)
		// Triangle 1: 021 = ----
		// | /
		// |/
		// Triangle 2: 231 = /|
		// ... / |
		// ----
		//
		// translate 0,2,4..., 1,3,5... in k, k+1, l, k+1
		// this should give a triangulated surface

		// KHK todo: vector() is an obsolete or deprecated collection -> use ArrayList
		// instead!
		Vector<Pnt3d> temp = new Vector<>(); // wie ein Array nur ohne vordefinierte Länge

		// started to modify code to ArrayList, but postponed it for the moment.
		/*
		 * List<Vector3f> vertexArrayList = new ArrayList<Vector3f>(); Vector3f vertex;
		 * 
		 * vertex.x = ... ; vertex.y = ... ; vertex.z = ... ;
		 */

		/*
		 * ... just a few examples... for (Vector3f p : mesh) { x = p.getX(); y =
		 * p.getY(); z = p.getZ();
		 * 
		 * vertexArrayList.add(new Vector3f(x, y, z)); }
		 */

		float tempX, tempY, tempZ;
		for (int k = 0; k < width - 1; k++) {
			for (int l = 0; l < height - 1; l++) {
				// first triangle
				if (ip.getPixelValue(k, l) > 0 && ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k, l + 1) > 0) {
					// first vertex
					tempX = k * (float) scaleX;
					tempY = l * (float) scaleY;
					tempZ = ip.getPixelValue(k, l);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
					// second vertex
					tempX = (k + 1) * (float) scaleX;
					tempY = l * (float) scaleY;
					tempZ = ip.getPixelValue(k + 1, l);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
					// third vertex
					tempX = k * (float) scaleX;
					tempY = (l + 1) * (float) scaleY;
					tempZ = ip.getPixelValue(k, l + 1);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
				}
				// second triangle with same k,l,k+1,l+1
				if (ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k + 1, l + 1) > 0
						&& ip.getPixelValue(k, l + 1) > 0) {
					// first vertex
					tempX = (k + 1) * (float) scaleX;
					tempY = l * (float) scaleY;
					tempZ = ip.getPixelValue(k + 1, l);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
					// second vertex
					tempX = (k + 1) * (float) scaleX;
					tempY = (l + 1) * (float) scaleY;
					tempZ = ip.getPixelValue(k + 1, l + 1);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
					// third vertex
					tempX = k * (float) scaleX;
					tempY = (l + 1) * (float) scaleY;
					tempZ = ip.getPixelValue(k, l + 1);
					temp.add(new Pnt3d(tempX, tempY, tempZ));
				}
			}
		}

		Pnt3d[] ret = new Pnt3d[temp.size()];
		temp.toArray(ret); // to preserve type Pnt3d
		return ret;
	}

	private Pnt3d dotProduct(Pnt3d P1, Pnt3d P2, Pnt3d P3) {

		// for Paraview normalization not necessary

		/*
		 * dot product a = b x c -> Wikipedia.org (a_x) a = (a_y) , same for b and c
		 * (a_z)
		 * 
		 * a_x = b_y c_z - b_z c_y a_y = b_z c_x - b_x c_z a_z = b_x c_y - b_y c_x
		 */

		/*
		 * Normalizing a vector involves two steps: 1 calculate its length, then, 2
		 * divide each of its (xy or xyz) components by its length.
		 */
		/*
		 * double norm_b = Math.sqrt((b.x * b.x) + (b.y * b.y) + (b.z * b.z)); double
		 * norm_c = Math.sqrt((c.x * c.x) + (c.y * c.y) + (c.z * c.z));
		 * 
		 * // normalize b and c
		 * 
		 * b.x = b.x/norm_b; b.y = b.y/norm_b; b.z = b.z/norm_b; c.x = c.x/norm_c; c.y =
		 * c.y/norm_c; c.z = c.z/norm_c;
		 * 
		 * 
		 * The cross product of two sides of the triangle equals the surface normal. So,
		 * if V = P2 - P1 and W = P3 - P1, and N is the surface normal, then:
		 * 
		 * Nx=(Vy∗Wz)−(Vz∗Wy) Ny=(Vz∗Wx)−(Vx∗Wz) Nz=(Vx∗Wy)−(Vy∗Wx)
		 */

		// bei dieser Regel gehen die beiden Vektoren, mit denen die Normale berechnet
		// wird
		// beide vom Punkt P1 aus. Die Richtung der Normale ergibt ich nach der Rechten
		// Handregel:
		// Daumen zeigt von P1 nach P2, Zeigefinger zeigt von P1 nach P3, Normale zeigt
		// in Richtung Mittelfinger

		Pnt3d V = new Pnt3d(0, 0, 0);
		Pnt3d W = new Pnt3d(0, 0, 0);
		Pnt3d N = new Pnt3d(0, 0, 0);
		Pnt3d normOfN = new Pnt3d(0, 0, 0);

		V.x = P2.x - P1.x;
		V.y = P2.y - P1.y;
		V.z = P2.z - P1.z;

		W.x = P3.x - P1.x;
		W.y = P3.y - P1.y;
		W.z = P3.z - P1.z;

		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

		/*
		 * 
		 * Cross product vectors can vary in size depending upon the components, so in
		 * order for OpenGL to make some comparisons needed for lighing, it is necessary
		 * to "normalize" the result by making them all uniform in some way. This is
		 * done by creating a normal vector with distance of one (1) unit long. To scale
		 * the normal, it is necessary to divide each component part by the total length
		 * or distance of the cross product vector which will scale each portion
		 * appropriately.
		 * 
		 * The length of the cross product vector or its distance is calculated with the
		 * following formula: dist = SQRT( cross[x]2 + cross[y]2 + cross[z]2)
		 * 
		 * The components of the normal would then become: norm[x] = cross[x] / dist
		 * norm[y] = cross[y] / dist norm[z] = cross[z] / dist
		 * 
		 */

		float dist = (float) Math.sqrt(N.x * N.x + N.y * N.y + N.z * N.z);

		normOfN.x = N.x / dist;
		normOfN.y = N.y / dist;
		normOfN.z = N.z / dist;

		// System.out.println("normal x/y/z: "+a.x+"/"+a.y+"/"+a.z);
		return normOfN;
	}

	static class Pnt3d {
		float x;
		float y;
		float z;

		Pnt3d(float x, float y, float z) {
			this.x = x;
			this.y = y;
			this.z = z;
			// System.out.println("Point3d: " + x + " / " + y + " / " + z);
		}
	}

}