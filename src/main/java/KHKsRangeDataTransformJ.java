// KHK status 30.4.2014
// Dieses File ist eine "Studie". Ich wollte unbedingt die target-to-source Abbildung bei 2.5 D Daten verstehen.
// Es hat mich eine Woche gekostet. Ich habe Umwege gemacht, über Raytracing, baryzentrische Interpolation, ...
// Ich habe viel gelernt, z. B. die Verwendung von List, Arraylist etc. 
// Ich habe auch die sinnvolle Anwendung für Objekte (hier datastruct.Triangle3D) mit ihren eigenen Methoden schätzen gelernt!!
// Die wichtigste Erkenntnis war aber, dass ich das Problem nicht über den anschaulichen Ansatz: Skizzen und probieren 
// lösen konnte. Erst als ich mit Papier und Kugelschreiber über die Transformationen ausgehend von der Ebenengleichung
// das Problem mit Matrizen gelöst hatte, konnte ich es in 30 min implementieren und es hat auf Anhieb funktioniert.
// In dieser Version lösche ich die Entwicklungsschritte über Raytracing, Quadtrees etc.

// Analyse von TransformJ:
// - verwendet zum Einlesen der Matrix: double[r][c]
// - es gibt eine Klasse imagescience/transform/Transform.java in der die Matrix weiterverwendet wird
// - in der Klasse TJ_Matrix.java wird die Matrix gelesen/geschrieben, die gleiche Leseroutine wird in TJ_Affine.java verwendet
// - TJ_Affine frage die Parameter ab, liest die Matrix ein, ruft das Objekt "Affine" auf. Übergaben von imp, Matrix etc., mit
//   TJ_Affine wird dann auch das Ergebnis wieder dargestellt.
// - TransformJ verwendet einen Wrapper für ImagePlus (imagescience/image/Image.java)
// - nebenbei: zur Verwendung von final beim Aufruf einer Methode/Klasse (von: http://stackoverflow.com/questions/500508/why-should-i-use-the-keyword-final-on-a-method-parameter-in-java)
//  If you create an anonymous inner class in your method and use a local variable (such as a method parameter) inside that class, then the compiler forces you to make the parameter final:
//  public Iterator<Integer> createIntegerIterator(final int from, final int to)

import static java.lang.Math.rint;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import datastruct.Triangle3D;
import datastruct.Vec3;
import ij.ImagePlus;
import ij.gui.NewImage;
import ij.io.FileInfo;
import ij.process.ImageProcessor;

/**
 * Does target-to-source transformation of range data. Data must be on a regular
 * grid for this implementation.
 * 
 * @author kkunzelm
 * 
 *
 */

// This class is called by KHKs_RangeDataTransformJ
public class KHKsRangeDataTransformJ {

	public KHKsRangeDataTransformJ() {

	}

	// KHK todo: klären - Übergabe der Transformationsmatix in Zukunft als JAMA
	// matrix oder als double[][] Array
	// KHK todo: überlegen, ob analog zu TransformJ eine Option zur Erweiterung der
	// Bildgröße eingeplant werden sollte??
	public static void applyT2sTransformationRangeData(ImagePlus imp, double[][] transformationMatrix,
			int interpolationMethod) {

		boolean debug = true;

		// Prepare new image, copy dimensions and scaling from input image
		ImagePlus impTransformedImage;
		FileInfo fi = imp.getFileInfo();
		ImageProcessor ip = imp.getProcessor();
		Match3d_withFiducialMarkersAndICPv2_1.setZeroToNan(ip);// masking: zero = NaN
		int w = ip.getWidth();
		int h = ip.getHeight();

		// debugging:
		if (debug) {
			System.out.println("applyT2sTransformationRangeData: Dimensions impTransformed (w x t): " + w + " x " + h);
		}

		// Syntax: NewImage.createFloatImage(java.lang.String title, int width, int
		// height, int slices, int options)
		impTransformedImage = NewImage.createFloatImage("Transformed t2s Img " + interpolationMethod, w, h, 1,
				NewImage.FILL_BLACK);

		// transfer scale from input image to the new image
		impTransformedImage.copyScale(imp);

		ImageProcessor ipTransformedImage = impTransformedImage.getProcessor();

		double scalex = fi.pixelWidth;
		double scaley = fi.pixelHeight;

		// List triangleList = makeTriangleList(imp);

		// inverse transform of all triangles
		// List triangleList = makeInvTransformedTriangleList(imp,transformationMatrix);

		// für forward Transformation verwendete Zeile
		// List triangleList = makeTransformedTriangleList(imp,transformationMatrix);

		Vec3 origin = new Vec3();
		Vec3 direction = new Vec3();

		// sample direction
		// the direction is the same for all points.
		// it is calculated with the inverse transformation matrix.
		// only the rotation component of the transformation matrix is needed.
		// the inverse transformation is applied to the normal vector pointing in the
		// target frame of reference in z-direction

		Matrix directionAsMatrix = new Matrix(3, 1); // I make a matrix from the array to use available JAMA matrix
														// multiplication later

		directionAsMatrix.set(0, 0, 0);
		directionAsMatrix.set(1, 0, 0);
		directionAsMatrix.set(2, 0, 1.0); // pointing in z direction

		System.out.println("direction as Matrix: ");
		directionAsMatrix.print(10, 3);

		// directions should not be translated, only rotated
		Matrix transformation = new Matrix(transformationMatrix);
		Matrix rotation = transformation.getMatrix(0, 2, 0, 2);
		Matrix invRotation = rotation.inverse();

		boolean printInfos = true;
		Matrix invTM = Match3d_withFiducialMarkersAndICPv2_1.getInverseTransformationMatrix(transformationMatrix,
				printInfos);

		// the target point at position u,v is mapped with the inverse transformation to
		// position x,y of source data set
		// directionAsMatrix = invRotation.times(directionAsMatrix);

		System.out.println("direction as Matrix after inverse transformation: ");
		directionAsMatrix.print(10, 3);

		// direction where to look for triangles - after inverse transformation
		direction.x = directionAsMatrix.get(0, 0);
		direction.y = directionAsMatrix.get(1, 0);
		direction.z = directionAsMatrix.get(2, 0);

		System.out.println("direction as Vector after inverse transformation: " + direction.toString());

		Matrix normOfPlane = new Matrix(4, 1); // I make a matrix from the array to use available JAMA matrix
												// multiplication later

		normOfPlane.set(0, 0, 0);
		normOfPlane.set(1, 0, 0);
		normOfPlane.set(2, 0, 1.0); // pointing in z direction
		normOfPlane.set(3, 0, 1.0);

		System.out.println("normOfPlane as Matrix: ");
		normOfPlane.print(10, 3);

		// transform norm
		normOfPlane = transformation.times(normOfPlane);

		System.out.println("normOfPlane as Matrix after transformation: ");
		normOfPlane.print(10, 3);

		for (int v = 0; v < h; v++) { // rows of the image = y coord.
			// System.out.print(".");
			for (int u = 0; u < w; u++) { // columns of the image = x coord.

				// Ebene in Normalenform: (vec_x - vec_a) * vec_n
				// vec_x = beliebiger Vector auf Ebene
				// vec_a = Aufpunkt zur Ebene
				// vec_n = Flächennormale der Ebene

				// KHK todo zero und NaN für Differenzen und Interpolation nachvollziehen
				// KHK todo Dokumentation aus blauem Heft in Source code ergänzen

				// vector to target plane
				Matrix vecTarget = new Matrix(4, 1);

				vecTarget.set(0, 0, u * scalex);
				vecTarget.set(1, 0, v * scaley);
				vecTarget.set(2, 0, 0);
				vecTarget.set(3, 0, 1.0);

				// System.out.println("vecTarget as Matrix: ");
				// vecTarget.print(10, 3);

				// multiply vecTarget with transformed normOfPlane
				vecTarget.times(normOfPlane.transpose());

				// inversely transform vecTarget to get x,y coordinates
				vecTarget = invTM.times(vecTarget);

				double result = 0;
				switch (interpolationMethod) {
					case 0 :
						result = nearestNeighborInterpolation(imp, vecTarget);
						break;
					case 1 :
						result = bilinearInterpolation(imp, vecTarget);
						break;
					case 2 :
						result = bicubicInterpolation(imp, vecTarget);
						break;
					case 3 :
						result = bicubicInterpolationBobDoughertyStyle(imp, vecTarget);
						break;
				}

				vecTarget.set(2, 0, result);

				vecTarget = transformation.times(vecTarget);

				result = vecTarget.get(2, 0);

				ipTransformedImage.setf(u, v, (float) result);

			}
		} // Ende

		ipTransformedImage.resetMinAndMax();
		double min = ipTransformedImage.getMin();
		double max = ipTransformedImage.getMax();
		ipTransformedImage.setMinAndMax(min, max);

		impTransformedImage.updateAndDraw();
		impTransformedImage.show();

		// return impTransformedImage;
	}

	public static List makeInvTransformedTriangleList(ImagePlus imp, double[][] transformationMatrix) {

		FileInfo fi = imp.getFileInfo();
		double scaleX = fi.pixelWidth;
		double scaleY = fi.pixelHeight;

		ImageProcessor ip = imp.getProcessor();
		Match3d_withFiducialMarkersAndICPv2_1.setZeroToNan(ip);// masking: zero = NaN
		int width = ip.getWidth();
		int height = ip.getHeight();

		boolean printInfos = true;
		Matrix invTM = Match3d_withFiducialMarkersAndICPv2_1.getInverseTransformationMatrix(transformationMatrix,
				printInfos);

		// List with triangles
		List<Triangle3D> triangleList = new ArrayList<>();

		for (int k = 0; k < width - 1; k++) {
			for (int l = 0; l < height - 1; l++) {
				// first triangle
				if (ip.getPixelValue(k, l) > 0 && ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k, l + 1) > 0) { // pixels
																														// with
																														// value
																														// zero
																														// are
																														// not
																														// used

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = k * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = l * scaleY;
					B.z = ip.getPixelValue(k + 1, l);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					A = transformPoint(A, invTM);
					B = transformPoint(B, invTM);
					C = transformPoint(C, invTM);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}

				// second triangle with same k,l,k+1,l+1
				if (ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k + 1, l + 1) > 0
						&& ip.getPixelValue(k, l + 1) > 0) {

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = (k + 1) * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k + 1, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = (l + 1) * scaleY;
					B.z = ip.getPixelValue(k + 1, l + 1);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					A = transformPoint(A, invTM);
					B = transformPoint(B, invTM);
					C = transformPoint(C, invTM);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}
			}
		}
		return triangleList;
	}

	public static List makeTransformedTriangleList(ImagePlus imp, double[][] transformationMatrix) {

		FileInfo fi = imp.getFileInfo();
		double scaleX = fi.pixelWidth;
		double scaleY = fi.pixelHeight;

		ImageProcessor ip = imp.getProcessor();
		Match3d_withFiducialMarkersAndICPv2_1.setZeroToNan(ip);// masking: zero = NaN
		int width = ip.getWidth();
		int height = ip.getHeight();

		Matrix transformation = new Matrix(transformationMatrix);

		// List with triangles
		List<Triangle3D> triangleList = new ArrayList<>();

		for (int k = 0; k < width - 1; k++) {
			for (int l = 0; l < height - 1; l++) {
				// first triangle
				if (ip.getPixelValue(k, l) > 0 && ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k, l + 1) > 0) { // pixels
																														// with
																														// value
																														// zero
																														// are
																														// not
																														// used

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = k * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = l * scaleY;
					B.z = ip.getPixelValue(k + 1, l);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					A = transformPoint(A, transformation);
					B = transformPoint(B, transformation);
					C = transformPoint(C, transformation);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}

				// second triangle with same k,l,k+1,l+1
				if (ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k + 1, l + 1) > 0
						&& ip.getPixelValue(k, l + 1) > 0) {

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = (k + 1) * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k + 1, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = (l + 1) * scaleY;
					B.z = ip.getPixelValue(k + 1, l + 1);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					A = transformPoint(A, transformation);
					B = transformPoint(B, transformation);
					C = transformPoint(C, transformation);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}
			}
		}
		return triangleList;
	}

	private static Vec3 transformPoint(Vec3 vector, Matrix transformation) {

		Matrix point = new Matrix(4, 1); // point in target frame of reference
		point.set(0, 0, vector.x);
		point.set(1, 0, vector.y);
		point.set(2, 0, vector.z); // z is zero, only the grid is mapped
		point.set(3, 0, 1.0);

		point = transformation.times(point); // point in source frame of reference (after inverse transformation)

		Vec3 transformedPoint;

		transformedPoint = new Vec3(point.get(0, 0), point.get(1, 0), point.get(2, 0));

		return transformedPoint;

	}

	public static List makeTriangleList(ImagePlus imp) {

		ImageProcessor ip = imp.getProcessor();

		int width = ip.getWidth();
		int height = ip.getHeight();

		FileInfo fi = imp.getFileInfo();
		double scaleX = fi.pixelWidth;
		double scaleY = fi.pixelHeight;

		// List with triangles
		List<Triangle3D> triangleList = new ArrayList<>();

		// get range data and triangulate them
		//
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

		for (int k = 0; k < width - 1; k++) {
			for (int l = 0; l < height - 1; l++) {
				// first triangle
				if (ip.getPixelValue(k, l) > 0 && ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k, l + 1) > 0) { // pixels
																														// with
																														// value
																														// zero
																														// are
																														// not
																														// used

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = k * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = l * scaleY;
					B.z = ip.getPixelValue(k + 1, l);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}

				// second triangle with same k,l,k+1,l+1
				if (ip.getPixelValue(k + 1, l) > 0 && ip.getPixelValue(k + 1, l + 1) > 0
						&& ip.getPixelValue(k, l + 1) > 0) {

					// vertices of triangle
					Vec3 A = new Vec3();
					Vec3 B = new Vec3();
					Vec3 C = new Vec3();

					// first vertex
					A.x = (k + 1) * scaleX;
					A.y = l * scaleY;
					A.z = ip.getPixelValue(k + 1, l);

					// second vertex
					B.x = (k + 1) * scaleX;
					B.y = (l + 1) * scaleY;
					B.z = ip.getPixelValue(k + 1, l + 1);

					// third vertex
					C.x = k * scaleX;
					C.y = (l + 1) * scaleY;
					C.z = ip.getPixelValue(k, l + 1);

					// accumulate temporary vertex array
					Vec3[] vertices = new Vec3[3];

					// vertex array
					vertices[0] = A;
					vertices[1] = B;
					vertices[2] = C;

					// make triangle
					Triangle3D triangle = new Triangle3D();
					triangle.makeTriangle3D(vertices);

					// add current triangle to triangleList
					triangleList.add(triangle);
				}

			}
		}

		return triangleList;

	}

	private static double nearestNeighborInterpolation(ImagePlus imp, Matrix pointInvTransformed) {

		double interpolatedPixelValue;

		ImageProcessor ip = imp.getProcessor();
		int w = ip.getWidth();
		int h = ip.getHeight();

		FileInfo fi = imp.getFileInfo();
		double scalex = fi.pixelWidth;
		double scaley = fi.pixelHeight;

		// for the moment I use the nearest neighbor interpolation to get the z-value at
		// x,y
		// add interpolations hooks here
		int x = (int) rint(pointInvTransformed.get(0, 0) / scalex);
		int y = (int) rint(pointInvTransformed.get(1, 0) / scaley);

		// we have no meaningful data outside the image boundaries
		if (x >= 0 && x < w && y >= 0 && y < h) {
			interpolatedPixelValue = ip.getf(x, y);
		} else {
			interpolatedPixelValue = 0;
		}

		return interpolatedPixelValue;

	}

	private static double bilinearInterpolation(ImagePlus imp, Matrix pointInvTransformed) {

		double interpolatedPixelValue;

		ImageProcessor ip = imp.getProcessor();
		int w = ip.getWidth();
		int h = ip.getHeight();

		FileInfo fi = imp.getFileInfo();
		double scalex = fi.pixelWidth;
		double scaley = fi.pixelHeight;

		double x = pointInvTransformed.get(0, 0) / scalex;
		double y = pointInvTransformed.get(1, 0) / scaley;

		// bilinear is 1, bicubic would be 2
		ip.setInterpolationMethod(1);

		// we have no meaningful data outside the image boundaries
		if (x >= 0 && x < w && y >= 0 && y < h) {
			interpolatedPixelValue = ip.getInterpolatedValue(x, y);
		} else {
			interpolatedPixelValue = 0;
		}

		return interpolatedPixelValue;

	}

	private static double bicubicInterpolation(ImagePlus imp, Matrix pointInvTransformed) {

		double interpolatedPixelValue;

		ImageProcessor ip = imp.getProcessor();
		int w = ip.getWidth();
		int h = ip.getHeight();

		FileInfo fi = imp.getFileInfo();
		double scalex = fi.pixelWidth;
		double scaley = fi.pixelHeight;

		double x = pointInvTransformed.get(0, 0) / scalex;
		double y = pointInvTransformed.get(1, 0) / scaley;

		// bilinear is 1, bicubic would be 2
		ip.setInterpolationMethod(2);

		// we have no meaningful data outside the image boundaries
		if (x >= 0 && x < w && y >= 0 && y < h) {
			interpolatedPixelValue = ip.getInterpolatedValue(x, y);
		} else {
			interpolatedPixelValue = 0;
		}

		return interpolatedPixelValue;

	}

	// File: CubicFloatProcessore.java needed!!
	// Bob Dougherty 6/19/2008. Modified from FloatProcesor by Wayne Rasband.
	// http://www.optinav.com/CubicFloatResizeRotate.htm
	private static double bicubicInterpolationBobDoughertyStyle(ImagePlus imp, Matrix pointInvTransformed) {

		double interpolatedPixelValue;

		ImageProcessor ip = imp.getProcessor();

		// bicubic interpolation developed by Bob Dougherty
		CubicFloatProcessor cfp = new CubicFloatProcessor(ip.getFloatArray());

		int w = ip.getWidth();
		int h = ip.getHeight();

		FileInfo fi = imp.getFileInfo();
		double scalex = fi.pixelWidth;
		double scaley = fi.pixelHeight;

		double x = pointInvTransformed.get(0, 0) / scalex;
		double y = pointInvTransformed.get(1, 0) / scaley;

		// we have no meaningful data outside the image boundaries
		if (x >= 0 && x < w && y >= 0 && y < h) {
			interpolatedPixelValue = cfp.getInterpolatedPixel(x, y);
		} else {
			interpolatedPixelValue = 0;
		}

		return interpolatedPixelValue;

	}

	// efficient implementation of matrix inverse for affine transform, not used
	// here so far
	// http://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform

	public static double[][] GetInverse(double[][] a) {

		double s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1];
		double s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
		double s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3];
		double s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
		double s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3];
		double s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];

		double c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3];
		double c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
		double c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2];
		double c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
		double c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2];
		double c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];

		// Should check for 0 determinant
		double invdet = 1.0 / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

		double[][] b = new double[4][4];

		b[0][0] = (a[1][1] * c5 - a[1][2] * c4 + a[1][3] * c3) * invdet;
		b[0][1] = (-a[0][1] * c5 + a[0][2] * c4 - a[0][3] * c3) * invdet;
		b[0][2] = (a[3][1] * s5 - a[3][2] * s4 + a[3][3] * s3) * invdet;
		b[0][3] = (-a[2][1] * s5 + a[2][2] * s4 - a[2][3] * s3) * invdet;

		b[1][0] = (-a[1][0] * c5 + a[1][2] * c2 - a[1][3] * c1) * invdet;
		b[1][1] = (a[0][0] * c5 - a[0][2] * c2 + a[0][3] * c1) * invdet;
		b[1][2] = (-a[3][0] * s5 + a[3][2] * s2 - a[3][3] * s1) * invdet;
		b[1][3] = (a[2][0] * s5 - a[2][2] * s2 + a[2][3] * s1) * invdet;

		b[2][0] = (a[1][0] * c4 - a[1][1] * c2 + a[1][3] * c0) * invdet;
		b[2][1] = (-a[0][0] * c4 + a[0][1] * c2 - a[0][3] * c0) * invdet;
		b[2][2] = (a[3][0] * s4 - a[3][1] * s2 + a[3][3] * s0) * invdet;
		b[2][3] = (-a[2][0] * s4 + a[2][1] * s2 - a[2][3] * s0) * invdet;

		b[3][0] = (-a[1][0] * c3 + a[1][1] * c1 - a[1][2] * c0) * invdet;
		b[3][1] = (a[0][0] * c3 - a[0][1] * c1 + a[0][2] * c0) * invdet;
		b[3][2] = (-a[3][0] * s3 + a[3][1] * s1 - a[3][2] * s0) * invdet;
		b[3][3] = (a[2][0] * s3 - a[2][1] * s1 + a[2][2] * s0) * invdet;

		return b;
	}

}
