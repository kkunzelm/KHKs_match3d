import static java.lang.Math.rint;

import java.awt.*;
import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.List;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import datastruct.KDNode;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.*;
import ij.io.FileInfo;
import ij.plugin.ContrastEnhancer;
import ij.plugin.PlugIn;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import vecmath.Matrix3d;
import vecmath.Point3d;
import vecmath.Vector3d;

// Point based rigid registration
// ==============================
// 
// The following matching procedure is based on the paper of Kenichi Kanatani: Analysis of 3-D Rotation Fitting,
// IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 16, No. 5, pp. 543ff
// 
// First the centroid of two corresponding pointlists are determined and the data are translated so that the centroids
// are the origin. Then the rotation is calculated with the SVD of the correlation matrix
//
// The registration will be implemented as a two-step process
// 1. corresponding point-pairs help to find a coarse fit/inital guess for the consecutive ICP
// 2. the result will be refined with an ICP approach
// 
// This version uses the target-to-source mapping to display the transformed image

// A brief note on semantics:
// ==========================
// In the literature we find: model image, data image, source image, target image... and other terms.
// I use the following interpretation of the terms:
// the source image is the original data set which will never be changed. in our special case this refers to the
// "baseline situation" of a dental restoration in the mouth. the terms "model" image - as in the model to which
// anything else is matched -, or also image1/img1/...,  are synonymously used for source image. the data image is also
// refered to as image2, img2 etc. This image contains the data which are transformed to place them in the same position
// as the model/source image. In dentistry we call this the follow up image as it is a replica of the restoration after
// years of service in the mouth.
// The target image is the image which accepts the data image AFTER the transformation is applied.
// After transformation we are subtracting the target image from the source image for further evaluation
// 
// History: 
// ========
// 2008 the first draft was functional
// October 2012 the plugin was working as intended        
// 27th March 2014 Wolfram Gloger, a long time companion, died ... I must finish this "in memoriam wg"
// In 2014 I started with the version from 2012 and decided to clean and rewrite the code 
// 
// This individual class here is just the implementation for the initial guess based on fiducial markers
// = manually added landmark points.
// The landmark points have to be entered as a polygon of type polyroi. 
// In addition, an interpretation of the Iterative Closes Point algorithm is implemented for fine adjustments.
// The source of the ICP applet which help me a lot to understand registration is from: 
// http://www9.informatik.uni-erlangen.de:81/sfb603/Saeulen/Optimierung/Allgemein/applet2 (accessed 2008)
// another excellent source of information was the thesis of Peter Neugebauer, Feinjustierung von Tiefenbildern, 1991

/**
 * @author Prof. Dr. Karl-Heinz Kunzelmann, Ludwig Maximilians Universität,
 * Poliklinik für Zahnerhaltung und Parodontologie, 80336 München email:
 * karl-heinz@kunzelmann.de www: www.kunzelmann.de
 * @version 0.1
 */
public class Match3d_withFiducialMarkersAndICPv2_1 implements PlugIn, DialogListener {

    private static final String[] schemes = {"refine_clamp", "refine_sd", "refine_clip", "refine_sparse",
            "refine_unique"};

    private static final String[] interpolationMethod = {"nearest neigbor", "bilinear", "bicubic",
            "bicubic (Bob Dougherty's Style)" // File: CubicFloatProcessore.java needed!!
            // Bob Dougherty 6/19/2008. Modified from FloatProcesor by Wayne Rasband.
            // http://www.optinav.com/CubicFloatResizeRotate.htm
    };

    private static int interpol = 1;
    private static int scheme = 1; // default = clamp
    private static String file = "";
    // Rest der Punktepaare verworfen
    private boolean refine_clamp = false; // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer
    // Abstand, definiert
    private double refine_clamp_par = 3.0;
    private boolean refine_sd = true; // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer
    // Abstand, definiert
    private double refine_sd_par = 2.5;
    private boolean refine_clip = false; // clipping ab bestimmter Position im Abstandshistogram,
    private double refine_clip_par = 0.90; // z. B. nur 75 % der niedrigeren Abstände berücksichtigt,
    private boolean refine_sparse = false; // nur jeder xte Punkte wird verwendet
    private double refine_sparse_par = 0.50;
    private boolean refine_unique = false; // klingt als ob jeder Punkt berücksichtigt wird
    private int minimum_valid_points = 800; // Anzahl der minimal notwendigen Punkte
    private boolean refine_a_priori_landmark = false;
    private boolean refine_a_priori_confidence_interval = false;
    private boolean refine_a_priori_show_diff_slider = false;
    private boolean use_landmark_mask = false;
    private ParameterICP refineParameters;
    private String title1 = "";
    private String title2 = "";
    private String[] titles;
    private ImagePlus imp1; // Mnemo: imp = imageplus
    private ImagePlus imp2;
    private int[] wList;
    private boolean loadMatrixFileFlag = false;
    private boolean saveMatrixFileFlag = false;
    private double[][] transformationMatrix;
    private boolean debug = false;
    private ImageProcessor ip1; // Mnemo: ip = image processor
    private ImageProcessor ip2;

    private double[][] correlationMatrixK = new double[3][3];

    double[] differences;
    ImagePlus impTransformedImage;
    ImagePlus ipOriginal;

    /*##################################################################################################################
	*
	* 												MAIN METHOD
	*
	##################################################################################################################*/

    public void run(String arg) {
        //--------------------------------------------------------------------------------------------------------------
        // some tests to be sure everything is fine before we start
        //--------------------------------------------------------------------------------------------------------------

        // get list of images & check for correct image type and number
        if (checkIJImageList()) return;

        // check for correct polygon selection
        if (checkLandmarkPolygons()) return;

        // ***************** shortcut ************ ***************
        // shortcut: Matrix file loaded... can be applied directly
        if (loadMatrixFileFlag) {
            impTransformedImage = makeTarget(imp1);
            impTransformedImage.setTitle("TransformedLoadedFromMatix");
            applyGeneralTransformation3D(imp1, imp2, impTransformedImage, transformationMatrix);
            return;
        }

        // ***************** lets start the real work ***************
        //--------------------------------------------------------------------------------------------------------------
        // get landmark centroids
        //--------------------------------------------------------------------------------------------------------------

        List<Vector3d> landmarks1 = getCorrespondingPointListFromPolygonRoi(imp1, (PolygonRoi) imp1.getRoi());
        List<Vector3d> landmarks2 = getCorrespondingPointListFromPolygonRoi(imp2, (PolygonRoi) imp2.getRoi());

        double[] centroidLandmarks1, centroidLandmarks2;

        centroidLandmarks1 = computeCentroid(landmarks1);
        centroidLandmarks2 = computeCentroid(landmarks2);

        List<Vector3d> relativeLandmarks1 = relativeCord(landmarks1, centroidLandmarks1);
        List<Vector3d> relativeLandmarks2 = relativeCord(landmarks2, centroidLandmarks2);

        if (debug) {
            System.out.println("PolygonRoi Centroids Landmarks1: " + Arrays.toString(centroidLandmarks1));
            System.out.println("PolygonRoi Centroids Landmarks2: " + Arrays.toString(centroidLandmarks2));
        }

        //--------------------------------------------------------------------------------------------------------------
        // calculate landmark rotation
        //--------------------------------------------------------------------------------------------------------------

        correlationMatrixK = correlationMatrixK(relativeLandmarks1, relativeLandmarks2);

        Matrix rotation = svd(correlationMatrixK); // singular value decomposition

        // for debugging
        if (debug = true) {
            System.out.println(" ");
            System.out.println("Original Rotationsmatrix wie nach svd geliefert: ");
            System.out.println(rotation.get(0, 0) + "\t" + rotation.get(0, 1) + "\t" + rotation.get(0, 2));
            System.out.println(rotation.get(1, 0) + "\t" + rotation.get(1, 1) + "\t" + rotation.get(1, 2));
            System.out.println(rotation.get(2, 0) + "\t" + rotation.get(2, 1) + "\t" + rotation.get(2, 2));
            System.out.println(" ");
        }

        Matrix invRotation = rotation.inverse();

        double[][] rotationMatrixR = rotation.getArrayCopy();
        double[][] invRotationMatrixInvR = invRotation.getArrayCopy();

        // for debugging
        if (debug = true) {
            System.out.print("svd() RotationMatrix K = ");
            rotation.print(9, 6);
            System.out.print("svd ()inverseRotationMatrix K^-1 = ");
            invRotation.print(9, 6);

            // for debugging
            Matrix rotTimeInv = rotation.times(invRotation);
            System.out.println("Rotation times Inverted Rotation (should be Identity I): ");
            rotTimeInv.print(9, 6);
        }

        //--------------------------------------------------------------------------------------------------------------
        // calculate landmark translation
        //--------------------------------------------------------------------------------------------------------------

        double[] translation = computeTranslation(centroidLandmarks1, centroidLandmarks2, rotationMatrixR);

        //--------------------------------------------------------------------------------------------------------------
        // calculate transformation matrix
        //--------------------------------------------------------------------------------------------------------------

        double[][] transformationMatrix = populateTransformationMatrix(rotationMatrixR, translation);
        if (debug = true) {
            System.out.println("TransformationMatrix in homogeneous coordinates: " + transformationMatrix[0][0] + "  "
                    + transformationMatrix[0][1] + "  " + transformationMatrix[0][2] + "  "
                    + transformationMatrix[0][3]);
            System.out.println("TransformationMatrix in homogeneous coordinates: " + transformationMatrix[1][0] + "  "
                    + transformationMatrix[1][1] + "  " + transformationMatrix[1][2] + "  "
                    + transformationMatrix[1][3]);
            System.out.println("TransformationMatrix in homogeneous coordinates: " + transformationMatrix[2][0] + "  "
                    + transformationMatrix[2][1] + "  " + transformationMatrix[2][2] + "  "
                    + transformationMatrix[2][3]);
            System.out.println("TransformationMatrix in homogeneous coordinates: " + transformationMatrix[3][0] + "  "
                    + transformationMatrix[3][1] + "  " + transformationMatrix[3][2] + "  "
                    + transformationMatrix[3][3]);
        }

        /*
         * // for debugging - testing my matrix multiplication double[] cLA = new double[3];
         * cLA[0] = centroidLandmarks2[0]; cLA[1] = centroidLandmarks2[1]; cLA[2] = centroidLandmarks2[2];
         *
         * Matrix cL = new Matrix(cLA, 1);
         *
         * Matrix c = rotation.times(cL.transpose()); // R*centroid2
         *
         * System.out.println("Matrix c = rotation.times(cL.transpose()): ");
         * c.print(10,2);
         *
         * double[] translationM = new double[3]; //translation calculated based on matrix operations
         * translationM[0] = centroidLandmarks1[0] - c.get(0,0);
         * translationM[1] = centroidLandmarks1[1] - c.get(1,0); translationM[2] = centroidLandmarks1[2] - c.get(2,0);
         * System.out.println("Translation von Matrix Mult: " + translationM[0] + "  " + translationM[1] + "  "
         * + translationM[2]);
         *
         * // end: for debugging - testing my matrix multiplication
         */

        // make the target image for the translation hier: source to target
        // ImagePlus impForwardTransformedImage = makeTarget(imp1);
        // impForwardTransformedImage.setTitle("ForwardTransformed");
        // applyGeneralForwardTransformation3D(imp2, impForwardTransformedImage, transformationMatrix);

        /*
         * ImagePlus impTransformedImage = makeTarget(imp1);
         * impTransformedImage.setTitle("Transformed Landmark based");
         * applyGeneralTransformation3D(imp1, imp2, impTransformedImage, transformationMatrix);
         */

        // ***** KHK Grobjustierung fertig, Vorbereitung für ICP
        //--------------------------------------------------------------------------------------------------------------
        // prepare for ICP
        //--------------------------------------------------------------------------------------------------------------

        Point3d[] vectorArrayImg1;
        vectorArrayImg1 = getCorrespondingPointListFromImageDataAsArray(imp1);

        Point3d[] vectorArrayImg2;
        vectorArrayImg2 = getCorrespondingPointListFromImageDataAsArray(imp2);

        Vector3d trans = new Vector3d();
        trans.x = translation[0];
        trans.y = translation[1];
        trans.z = translation[2];

        Matrix3d rotMatrix = new Matrix3d(rotation.get(0, 0), rotation.get(0, 1), rotation.get(0, 2),
                rotation.get(1, 0), rotation.get(1, 1), rotation.get(1, 2), rotation.get(2, 0), rotation.get(2, 1),
                rotation.get(2, 2));

        // Build Tree
        KDNode modelTree = new KDNode(vectorArrayImg1);
        // bounding box
        double x0, x1, y0, y1, z0, z1;
        x0 = y0 = z0 = Double.POSITIVE_INFINITY;
        x1 = y1 = z1 = Double.NEGATIVE_INFINITY;

        for (Point3d vectorArrayImg11 : vectorArrayImg1) {
            if (vectorArrayImg11.x < x0) {
                x0 = vectorArrayImg11.x;
            }
            if (vectorArrayImg11.x > x1) {
                x1 = vectorArrayImg11.x;
            }
            if (vectorArrayImg11.y < y0) {
                y0 = vectorArrayImg11.y;
            }
            if (vectorArrayImg11.y > y1) {
                y1 = vectorArrayImg11.y;
            }
            if (vectorArrayImg11.z < z0) {
                z0 = vectorArrayImg11.z;
            }
            if (vectorArrayImg11.z > z1) {
                z1 = vectorArrayImg11.z;
            }
        }
        System.out.println("Building model tree...");

        if (debug = true) {
            System.out.println("Model points boundaries:");
            System.out.println("  " + x0);
            System.out.println("  " + x1);
            System.out.println("  " + y0);
            System.out.println("  " + y1);
            System.out.println("  " + z0);
            System.out.println("  " + z1);
        }

        modelTree.build(0, vectorArrayImg1.length - 1, x0, x1, y0, y1, z0, z1);

        Point3d modelCorner0 = new Point3d(x0, y0, z0);
        Point3d modelCorner1 = new Point3d(x1, y1, z1);

        ICPAlgorithm2014 icp = new ICPAlgorithm2014();

        //--------------------------------------------------------------------------------------------------------------
        // run ICP
        //--------------------------------------------------------------------------------------------------------------
        // initialization of ICP algorithm
        icp.init(vectorArrayImg1, vectorArrayImg2, modelTree, modelCorner0, modelCorner1, rotMatrix, trans,
                refineParameters, centroidLandmarks1, centroidLandmarks2);

        // return value from ICP
        double[][] tMat = icp.runICP();

        // write matrix to file
        if (saveMatrixFileFlag) {
            new WriteMatrix().run(tMat);
        }

        //--------------------------------------------------------------------------------------------------------------
        // apply transformation to new image
        //--------------------------------------------------------------------------------------------------------------
        // make the target image for the translation hier: source to target
        // ImagePlus impForwardTransformedImageICP = makeTarget(imp1);
        // impForwardTransformedImageICP.setTitle("ForwardTransformedPostICP");
        // applyGeneralForwardTransformation3D(imp2, impForwardTransformedImageICP, tMat);

        // apply transformation matrix
        impTransformedImage = makeTarget(imp1);
        impTransformedImage.setTitle("TransformedPostICP");
        applyGeneralTransformation3D(imp1, imp2, impTransformedImage, tMat); // tMat = Transformationsmatrix

        // ********************************************
        // Dieses Programm ist eine coole Implementierung der Projektion von irregulären xyz-Koordinaten und deren
        // Interpolation auf ein reguläres Gitter.
        //
        // XYZ2DEM_ImporterHack hack0 = new XYZ2DEM_ImporterHack();
        // hack0.khkDisplayXYZ(vectorListImg2);
    }

	/*##################################################################################################################
	*
	* 												LANDMARKS
	*
	##################################################################################################################*/

    /**
     * Returns ArrayList of Vector3d
     *
     * @param imp
     * @param polyRoi
     * @return
     */
    private List<Vector3d> getCorrespondingPointListFromPolygonRoi(ImagePlus imp, PolygonRoi polyRoi) {
        List<Vector3d> vectorList;
        vectorList = new ArrayList<>();

        ImageProcessor ip = imp.getProcessor();
        FileInfo fi;

        fi = imp.getFileInfo();

        if (fi.fileType != FileInfo.GRAY32_FLOAT) {
            IJ.error("Error Image1 is not type 'float' - will try to convert");
            ip = ip.convertToFloat();
        }

        int[] roiX;
        int[] roiY;
        int roiLength = polyRoi.getNCoordinates();

        // later I need the absolute coordinates.
        // the method getXCoordinates or getYcoordinates returns relative coordinates,
		// relative to the bouncing box of the ROI.

        // Rectangle is the bouncing box
        Rectangle rectRoiSource = ip.getRoi();

        // Source Image: Creating array for SVD
        roiX = polyRoi.getXCoordinates();
        roiY = polyRoi.getYCoordinates();

        Polygon poly = polyRoi.getPolygon();

        System.out.println("***************************************************");
        System.out.println("PolygonRoi-Koordinaten: ");

        for (int i = 0; i < roiLength; i++) {
            Vector3d vert = new Vector3d();
            vert.x = (float) fi.pixelWidth * (roiX[i] + rectRoiSource.x);
            vert.y = (float) fi.pixelHeight * (roiY[i] + rectRoiSource.y);
            vert.z = ip.getf((roiX[i] + rectRoiSource.x), (roiY[i] + rectRoiSource.y));

            // compare pixels in 5x5 mask around the selected landmark
            if (use_landmark_mask) {
                System.out.println("use_landmark_mask");
                final int maskSizeX = 5;
                final int maskSizeY = 5;

                int ctr = 0;
                float mean = 0;

                for (int j = -2; j <= 2; j++) {
                    for (int k = -2; k <= 2; k++) {
                        mean += ip.getf(poly.xpoints[i] + j, poly.ypoints[i] + k);

                        // find pixel with highest z-value in mask and choose as new landmark
                        if (ip.getf(poly.xpoints[i] + j, poly.ypoints[i] + k) > vert.z) {
                            vert.x = poly.xpoints[i] + j;
                            vert.y = poly.ypoints[i] + k;
                            vert.z = ip.getf(poly.xpoints[i] + j, poly.ypoints[i] + k);
                        }

                        System.out.println("Poly: " + i + " X: " + poly.xpoints[i] + j + " Y: " + poly.ypoints[i] + k
                                + " Z: " + ip.getf(poly.xpoints[i] + j, poly.ypoints[i] + k));
                    }
                }
            }

            vectorList.add(vert);
//			System.out.println("New" + vert.x + "\t" + vert.y + "\t" + vert.z);
        }
        System.out.println("***************************************************");
        return vectorList;
    }

	private double[] computeCentroid(List<Vector3d> mesh) {
		double x = 0, y = 0, z = 0;

		for (Vector3d p : mesh) {
			x = x + p.getX();
			y = y + p.getY();
			z = z + p.getZ();
		}

		x = x / mesh.size();
		y = y / mesh.size();
		z = z / mesh.size();

		return new double[]{x, y, z};
	}

	private List<Vector3d> relativeCord(List<Vector3d> mesh, double[] centroid) {
		List<Vector3d> relative = new ArrayList<>();
		double x, y, z;

		for (Vector3d p : mesh) {
			x = p.getX() - centroid[0];
			y = p.getY() - centroid[1];
			z = p.getZ() - centroid[2];

			relative.add(new Vector3d(x, y, z));
		}

		return relative;
	}

    /**
     * Accumulates the Correlation Matrix K which is needed for Singular Value
     * Decomposition
     *
     * @param relativeLandmarks1 Landmark list number 1
     * @param relativeLandmarks2 Landmark list number 2
     * @return correlationMatrix in Kanatani's paper called K
     */
    private double[][] correlationMatrixK(List<Vector3d> relativeLandmarks1, List<Vector3d> relativeLandmarks2) {

        /*
         * this method implements formula (3) of Kanatani the weight is set to "1". It is a scaling factor only.
         * Formula (3):
         * K = SUM(vectorPoint1 * vectorPoint'1^T) --> ^T means transposed
         *
         * in detail:
         *
         * Each point consists of 3 coordinates x, y, z Each point is treated as a vector. The coordinates of these
         * vectors are arranged vertically:
         * Point1 = ( x1 ) ( y1 ) ( z1 )
         * Point'1^T = (x'1, y'1, z'1)
         *
         * The mark "'" means the corresponding point in the second image
         *
         * The multiplication of the two 3x1 vectors results in a 3x3 matrix = correlationMatrixK.
         * The matrices off all points are summarized and result in K
         *
         * KHK Added: 19.1.13 from: http://nghiaho.com/?page_id=671 Pay close attention to the transpose symbol. It’s
         * doing a multiplication between 2 matrices where the dimensions effectively are, 3×1 and 1×3, respectively.
         * The ordering of the multiplication is also important, doing it the other way will find a rotation from B to A
         * instead.
         */

        double[][] modelPoint = new double[3][3]; // model = baseline, unchanged, reference
        double[][] dataPoint = new double[3][3]; // data = follow-up
        double[][] correlationMatrixK_temp = new double[3][3];

        for (int x = 0; x < relativeLandmarks1.size(); x++) {

            modelPoint[0][0] = relativeLandmarks1.get(x).getX();
            modelPoint[1][0] = relativeLandmarks1.get(x).getY();
            modelPoint[2][0] = relativeLandmarks1.get(x).getZ();
            dataPoint[0][0] = relativeLandmarks2.get(x).getX();
            dataPoint[0][1] = relativeLandmarks2.get(x).getY();
            dataPoint[0][2] = relativeLandmarks2.get(x).getZ();

            System.out.println("Accumulation of the correlationMatrixK - step no: x = " + x);

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        correlationMatrixK_temp[i][k] += modelPoint[i][j] * dataPoint[j][k];
                    }
                }
            }

        }

        return correlationMatrixK_temp;
    }

    /**
     * Singular Value Decomposition of the correlation matrix.
     *
     * @param correlationMatrixK
     * @return rotation rotation matrix 3x3
     */
    private Matrix svd(double[][] correlationMatrixK) {

        /*
         * To perform the singular value decomposition svd the package Jama is needed As soon as we have the correlation
         * matrix, the next step in the algorithm as described by Kanatani is to decompose the correlation matrix
         * (see formulae 12 and 13)
         *
         * It is important to mention that Kanatani writes: K = VAU^T In Jama there is a svd.jama file which uses a
         * slightly different syntax A = USV^T Please concentrate not to mix the meaning of U and V Kanatini's U will be
         * the V of SVD to stay in the logic of Jama we write for formula (13)
         *
         * (1 0 0 ) R = U * (0 1 0 )*V^T formula 13 later called intermedResult (0 0 det(UV^T) )
         *
         * Just to mention a few details:
         *
         * Kanatani writes in formula (4) that to minimize the distance between the two pointsets for the given
         * correlation matrix (formula 3), the trace(R^T*K) musst be maximized in formula (13) he states - and proves
         * too - when the trace(R^T*K) is maximized
         *
         * What we do is, we use K make a SVD and get U and V Using U and V we can calculate R with formula 13 R is our
         * final result for the rotation
         *
         * to calculate R we need the intermediate result for det(VU^T) - Matrix.de()
         * then we do matrix multiplication
         *
         * ( 1 0 0 ) 1) U * ( 0 1 0 ) ( 0 0 det(VU^T))
         *
         * 2) |_____________________|
         *
         * = C * V^T |___________|
         *
         * = R
         *
         */

        this.correlationMatrixK = correlationMatrixK;

        Matrix A = new Matrix(this.correlationMatrixK);
        System.out.print("correlationMatrixK as Jama matrix object: ");
        A.print(9, 6);

        // compute the singular value decomposition
        System.out.println("A = U S V^T");
        System.out.println();

        SingularValueDecomposition s = A.svd();

        System.out.print("U = ");
        Matrix U = s.getU();
        U.print(9, 6);
        System.out.print("Sigma = ");
        Matrix S = s.getS();
        S.print(9, 6);
        System.out.print("V = ");
        Matrix V = s.getV();

        System.out.println("det(VU^T): " + s.getV().times(s.getU().transpose()).det());

        Matrix intermedResult = new Matrix(
                new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, s.getV().times(s.getU().transpose()).det()}});

        return s.getU().times(intermedResult).times(s.getV().transpose());

    }

    /**
     * Calculate the Translation between the two datasets
     *
     * @param centroid1       Source image center
     * @param centroid2       Target image center
     * @param rotationMatrixR
     * @return translation Translation vector
     */
    private double[] computeTranslation(double[] centroid1, double[] centroid2, double[][] rotationMatrixR) {

        double[] translation = new double[3];

        double[] rotCentroid2 = new double[3];

        // Calculate the translation.
        // Important: to calculate the translation, the centroid of the target image, which will be transformed,
		// has to be rotated!!!
        // T = Sc - Tc*R
        // Sc = Source image center, Tc = Target image center

        // Lets assume now that we have a matrix R. When P is multiplied by R it transforms P to PT. Considering what we
		// know about matrix multiplication lets see how we can re-write a point-matrix multiplication and isolate the
        // computation of each of the transformed point coordinates:
        // PT.x=P.x∗R00+P.y∗R10+P.z∗R20
        // PT.y=P.x∗R01+P.y∗R11+P.z∗R21
        // PT.z=P.x∗R02+P.y∗R12+P.z∗R22

        rotCentroid2[0] = (float) (rotationMatrixR[0][0] * centroid2[0] + rotationMatrixR[0][1] * centroid2[1]
                + rotationMatrixR[0][2] * centroid2[2]);
        rotCentroid2[1] = (float) (rotationMatrixR[1][0] * centroid2[0] + rotationMatrixR[1][1] * centroid2[1]
                + rotationMatrixR[1][2] * centroid2[2]);
        rotCentroid2[2] = (float) (rotationMatrixR[2][0] * centroid2[0] + rotationMatrixR[2][1] * centroid2[1]
                + rotationMatrixR[2][2] * centroid2[2]);

        translation[0] = centroid1[0] - rotCentroid2[0];
        translation[1] = centroid1[1] - rotCentroid2[1];
        translation[2] = centroid1[2] - rotCentroid2[2];

        System.out.println("Centroid1: " + centroid1[0] + " " + centroid1[1] + " " + centroid1[2]);
        System.out.println("Centroid2: " + centroid2[0] + " " + centroid2[1] + " " + centroid2[2]);
        System.out.println("Translation: " + translation[0] + ", " + translation[1] + ", " + translation[2]);

        return translation;

    }

    private double[][] populateTransformationMatrix(double[][] rotationMatrixR, double[] translation) {

        double[][] transformationMatrix = new double[4][4];

        for (int k = 0; k < 3; k++) {
            System.arraycopy(rotationMatrixR[k], 0, transformationMatrix[k], 0, 3);
        }
        transformationMatrix[3][0] = 0;
        transformationMatrix[3][1] = 0;
        transformationMatrix[3][2] = 0;

        transformationMatrix[0][3] = translation[0];
        transformationMatrix[1][3] = translation[1];
        transformationMatrix[2][3] = translation[2];
        transformationMatrix[3][3] = 1.0d;

        if (debug = true) {
            System.out.println("Transformation Matrix formated for TransformJ: ");
            System.out.println(transformationMatrix[0][0] + "\t" + transformationMatrix[0][1] + "\t"
                    + transformationMatrix[0][2] + "\t" + transformationMatrix[0][3]);
            System.out.println(transformationMatrix[1][0] + "\t" + transformationMatrix[1][1] + "\t"
                    + transformationMatrix[1][2] + "\t" + transformationMatrix[1][3]);
            System.out.println(transformationMatrix[2][0] + "\t" + transformationMatrix[2][1] + "\t"
                    + transformationMatrix[2][2] + "\t" + transformationMatrix[2][3]);
            System.out.println(transformationMatrix[3][0] + "\t" + transformationMatrix[3][1] + "\t"
                    + transformationMatrix[3][2] + "\t" + transformationMatrix[3][3]);
        }
        return transformationMatrix;
    }

    /*##################################################################################################################
	*
	* 											3D TRANSFORMATION
	*
	##################################################################################################################*/

    /**
     * DESCRIPTION: Takes two input datasets and the location parameters as a transformation matrix in homogeneous
	 * coordinates and computes the differences per pixel between the rotated and translated patterns. The difference is
     * displayed in a new image.
     * <p>
     * A target-to-source-mapping is applied.
     *
     * @param imp1
     * @param imp2
     * @param differenceImage
     * @param transformationMatrix
     */
    private void applyGeneralTransformation3D(ImagePlus imp1, ImagePlus imp2, ImagePlus differenceImage,
                                              double[][] transformationMatrix) {

        // target-to-source mapping as described for 2d in Burger and Burge (http://imagingbook.com/)

        // impTransformedImage is a blank image plus template same size as the "target" image (here ip1)
        ImageProcessor ipTransformedImage = differenceImage.getProcessor();

        FileInfo fi = imp1.getFileInfo();
        this.ip1 = imp1.getProcessor();
        setZeroToNan(ip1);// masking: zero = NaN
        int w = this.ip1.getWidth();
        int h = this.ip1.getHeight();

        this.ip2 = imp2.getProcessor();
        setZeroToNan(ip2);// masking: zero = NaN

        double tempz = 0.0;
        float diff;
        boolean printInfos = false; // used to switch on/off verbose printout during processing

        double scalex = fi.pixelWidth;
        double scaley = fi.pixelHeight;

        differences = new double[(w * h)];
        int ctr = 0;
        double mean = 0;

        for (int v = 0; v < h; v++) { // rows of the image
            for (int u = 0; u < w; u++) { // columns of the image

                // set(int i,int j,double s) - first index row, second column!
                Matrix pointTarget = new Matrix(4, 1);
                // I make a matrix from the array to use available JAMA matrix multiplication later
                Matrix pointInvTransformed; // = new Matrix(4,1) for the point after the inverse transformation;

                // see Peter Neugebauer Match3d: Feinjustierung von Tiefenbildern zur Vermessung von kleinen
				// Verfomungen, Diplomarbeit 1991.
                // The baseline image is mapped with the inverse transformation matrix close to the follow-up image
				// (= source): f(u,x) -> f'(x,y)
                // At x,y there exists usually no value for the source data therefore this value has to be interpolated.
                // The difference is calculated and mapped using the transformation from x,y back to u,v.
                // In our case we do not really need to apply the transformation as we just have to use u,v which is the
				// same after all. The difference is now in the target image (= differenceImage) which has the same
				// coordinates as the baseline image.

                pointTarget.set(0, 0, u * scalex);
                pointTarget.set(1, 0, v * scaley);
                pointTarget.set(2, 0, this.ip1.getf(u, v));
                pointTarget.set(3, 0, 1.0);

                // from JAMA API: matrixA.times(Matrix B) Linear algebraic matrix multiplication, A * B

                // I calculate the transformation matrix myself
				// Transformation T and inverse Transformation T^-1 in homogeneous coordinates:
                // R = the rotations part (upper left 3x3 matrix)
                // t = the translation part (upper right 3x1 matrix)
                //
                // T = ( R t ) R^-1 = R^T, t^-1 = -R^T*t T^-1 = ( R^T -R^T*t)
                // ( 0 1 ) ( 0 1 )
                //
                // the method getInverseTransformationMatrix should work correctly
                // T and T^-1 tested with TransformJ reveals the expected results and T^-1*T = I
                Matrix invTM = getInverseTransformationMatrix(transformationMatrix, printInfos);

                // the target point at position u,v is mapped with the inverse transformation to position x,y
                pointInvTransformed = invTM.times(pointTarget);

                switch (interpol) {
                    case 0:
                        tempz = nearestNeighborInterpolation(imp2, pointInvTransformed);
                        break;
                    case 1:
                        tempz = bilinearInterpolation(imp2, pointInvTransformed);
                        break;
                    case 2:
                        tempz = bicubicInterpolation(imp2, pointInvTransformed);
                        break;
                    case 3:
                        tempz = bicubicInterpolationBobDoughertyStyle(imp2, pointInvTransformed);
                        break;
                }

                diff = (float) (pointInvTransformed.get(2, 0) - tempz); // difference

                if (!Float.isNaN(diff)) {
                    differences[ctr++] = diff;
                    mean += diff;
                }

                ipTransformedImage.setf(u, v, diff);

                printInfos = false; // need to know the details only once
            }
        }
        // KHK todo move to separate method
        // update the rotated image

        // Tipp: man kann mit Threshold im 32 bit float Bild das gleiche machen, wie bei Wolfram mit clip-Histogram
        // Hintergrund wird zu NaN (anklicken).
        ipTransformedImage.resetMinAndMax();
        ContrastEnhancer ce = new ContrastEnhancer(); // KHK todo check again whether this is what I expect it to be
        ce.stretchHistogram(ipTransformedImage, 1.0); // ist eine Art Histogram Threshold. Es werden die obersten und
        // untersten Pixel entfernt. unklar ist, ob 1% von oben und 1% von unten oder 1%/2 von oben und 1%/2 von unten??
        double min = ipTransformedImage.getMin();
        double max = ipTransformedImage.getMax();
        ipTransformedImage.setMinAndMax(min, max);
        differenceImage.updateAndDraw();
        differenceImage.show();

        DiffImageStats stats = new DiffImageStats(mean, ctr, differences);
        if (refine_a_priori_confidence_interval) {
            // duplicate original transformation for reference values in image subtraction
            ipOriginal = impTransformedImage.duplicate();

            if(refine_a_priori_show_diff_slider)
                showDiffSlider(stats.getMinDiff(), stats.getMaxDiff(), stats.getHi());
            else
                subtractImages(stats.getHi());
        }

        if (refine_a_priori_show_diff_slider && !refine_a_priori_confidence_interval) {
            // duplicate original transformation for reference values in image subtraction
            ipOriginal = impTransformedImage.duplicate();

            showDiffSlider(stats.getMinDiff(), stats.getMaxDiff(), 0);
        }
    }

    public static Matrix getInverseTransformationMatrix(double[][] transformationMatrix, boolean printInfos) {
		/*
		 * Later to get homogeneous coordinates:
		 *
		 * R is extended to a 4 x 4 matrix (r.. stands for rotation element ..)
		 * r11 r12 r13 0 r21 r22 r23 0 R = r31 r32 r33 0 0 0 0 1
		 *
		 * T is extended to a 4 x 1 matrix (t. stands for translation element .)
		 * t1 T = t2 t3 1
		 */

		// KHK ich habe mit TransformJ kontrolliert: das Ergebnis sollte stimmen.
		//
		// die übliche Transformationsmatrix T
		//
		// (R t)
		// ( ) p = t*R*p (erst R um Umsprung, dann t)
		// (0 1)
		//
		// entspricht einer Kombination aus R und t:
		// with R being a rotation and t being a translation is a combined transformation
		// man muss die Rotation R und Translation t in der TransformationsMatrix T separat betrachten.
		// aufgrund der Orthogonalität der Rotationsachsen gilt: R^T = R^-1
		// d.h. wenn man die Rotationsmatrix transformiert, hat man die Inverse
		// für die Translation t gilt: inv(t) = -t
		// die Rotationsmatrixe kann man einfach invertieren
		// ABER: für Koordinaten gilt: p' = Rp + t
		// inverse Transformation zur Rückgewinnung von p:
		// p' = Rp + t entspricht p'-t = Rp entspricht R^-1(p'-t) = R^-1*R*p
		// entspricht wegen R^-1*R = I: p = R^-1(p'-t)
		// in der homogenen TransformationsMatrix hat das die Form:
		// (R t)^-1 (R^T -R^T*t)
		// (0 1) gleich (0 1 )
		//
		// ^-1 = Inverse
		// ^T = Transponierte

		// nebenbei: es gilt immer... erst in Ursprung verschieben, dann rotieren, dann weiter translatieren
		// t2*R*t1 ... bei Spaltenvektoren/Matrizen von rechts nach links abgearbeitet.
		// ABER: die Umkehrung von Matrixverkettungen ist
		// What is the inverse of a sequence of transformations?
		// (M1M2...Mn)-­‐1 = Mn-­‐1Mn-­‐1-­‐1...M1-­‐1
		// Inverse of a sequence of transformations is the composition of the inverses of each transformation in reverse
		// order.
		// What will our sequence look like?
		// (T^‐1RT)^-1 = T^‐1R^‐1T
		// We still translate to the origin first, then translate back at the end!

		Matrix invTM = new Matrix(4, 4);
		// R^T part (R transposed)
		invTM.set(0, 0, transformationMatrix[0][0]);
		invTM.set(0, 1, transformationMatrix[1][0]);
		invTM.set(0, 2, transformationMatrix[2][0]);
		invTM.set(1, 0, transformationMatrix[0][1]);
		invTM.set(1, 1, transformationMatrix[1][1]);
		invTM.set(1, 2, transformationMatrix[2][1]);
		invTM.set(2, 0, transformationMatrix[0][2]);
		invTM.set(2, 1, transformationMatrix[1][2]);
		invTM.set(2, 2, transformationMatrix[2][2]);

		// invariant part
		invTM.set(3, 0, 0.0);
		invTM.set(3, 1, 0.0);
		invTM.set(3, 2, 0.0);
		invTM.set(3, 3, 1.0);

		// translation part -R^T*t
		Matrix invRot = new Matrix(3, 3);
		// R^T part (R transposed)
		invRot.set(0, 0, transformationMatrix[0][0]);
		invRot.set(0, 1, transformationMatrix[1][0]);
		invRot.set(0, 2, transformationMatrix[2][0]);
		invRot.set(1, 0, transformationMatrix[0][1]);
		invRot.set(1, 1, transformationMatrix[1][1]);
		invRot.set(1, 2, transformationMatrix[2][1]);
		invRot.set(2, 0, transformationMatrix[0][2]);
		invRot.set(2, 1, transformationMatrix[1][2]);
		invRot.set(2, 2, transformationMatrix[2][2]);

		Matrix transL = new Matrix(3, 1);
		transL.set(0, 0, transformationMatrix[0][3]);
		transL.set(1, 0, transformationMatrix[1][3]);
		transL.set(2, 0, transformationMatrix[2][3]);

		Matrix modTransL;
		modTransL = invRot.times(transL);

		invTM.set(0, 3, (-1 * modTransL.get(0, 0)));
		invTM.set(1, 3, (-1 * modTransL.get(1, 0)));
		invTM.set(2, 3, (-1 * modTransL.get(2, 0)));

		// Kontrolle: invTM*T = T*invTM = I;
		Matrix transformation = new Matrix(transformationMatrix);
		Matrix identity = invTM.times(transformation);

		if (printInfos) {
			System.out.println("*************************************************");
			System.out.println("* Kontrolle von getInverseTransformationMatrix:");
			System.out.println("* ");
			System.out.println("* OriginalTransformationsMatrix: ");
			System.out.println("* ");
			System.out.println(transformationMatrix[0][0] + "\t" + transformationMatrix[0][1] + "\t"
					+ transformationMatrix[0][2] + "\t" + transformationMatrix[0][3]);
			System.out.println(transformationMatrix[1][0] + "\t" + transformationMatrix[1][1] + "\t"
					+ transformationMatrix[1][2] + "\t" + transformationMatrix[1][3]);
			System.out.println(transformationMatrix[2][0] + "\t" + transformationMatrix[2][1] + "\t"
					+ transformationMatrix[2][2] + "\t" + transformationMatrix[2][3]);
			System.out.println(transformationMatrix[3][0] + "\t" + transformationMatrix[3][1] + "\t"
					+ transformationMatrix[3][2] + "\t" + transformationMatrix[3][3]);
			System.out.println("* ");
			System.out.println("* inverse Rotationsmatrix: ");
			System.out.println("* ");
			System.out.println(invRot.get(0, 0) + "\t" + invRot.get(0, 1) + "\t" + invRot.get(0, 2));
			System.out.println(invRot.get(1, 0) + "\t" + invRot.get(1, 1) + "\t" + invRot.get(1, 2));
			System.out.println(invRot.get(2, 0) + "\t" + invRot.get(2, 1) + "\t" + invRot.get(2, 2));
			System.out.println("* ");
			System.out.println("* Translation");
			System.out.println("* ");
			System.out.println(transL.get(0, 0) + "\t" + transL.get(1, 0) + "\t" + transL.get(2, 0));
			System.out.println("* ");
			System.out.println("* -R^T * Translation");
			System.out.println("* ");
			System.out.println(
					(-1 * modTransL.get(0, 0)) + "\t" + (-1 * modTransL.get(1, 0)) + "\t" + (-1 * modTransL.get(2, 0)));
			System.out.println("* ");
			System.out.println("* inverse TransformationsMatrix: ");
			System.out.println("* ");
			System.out.println(
					invTM.get(0, 0) + "\t" + invTM.get(0, 1) + "\t" + invTM.get(0, 2) + "\t" + invTM.get(0, 3));
			System.out.println(
					invTM.get(1, 0) + "\t" + invTM.get(1, 1) + "\t" + invTM.get(1, 2) + "\t" + invTM.get(1, 3));
			System.out.println(
					invTM.get(2, 0) + "\t" + invTM.get(2, 1) + "\t" + invTM.get(2, 2) + "\t" + invTM.get(2, 3));
			System.out.println(
					invTM.get(3, 0) + "\t" + invTM.get(3, 1) + "\t" + invTM.get(3, 2) + "\t" + invTM.get(3, 3));
			System.out.println("* ");
			System.out.println("* Identiy ??: ");
			System.out.println(identity.get(0, 0) + "\t" + identity.get(0, 1) + "\t" + identity.get(0, 2) + "\t"
					+ identity.get(0, 3));
			System.out.println(identity.get(1, 0) + "\t" + identity.get(1, 1) + "\t" + identity.get(1, 2) + "\t"
					+ identity.get(1, 3));
			System.out.println(identity.get(2, 0) + "\t" + identity.get(2, 1) + "\t" + identity.get(2, 2) + "\t"
					+ identity.get(2, 3));
			System.out.println(identity.get(3, 0) + "\t" + identity.get(3, 1) + "\t" + identity.get(3, 2) + "\t"
					+ identity.get(3, 3));
			System.out.println("* ");
			System.out.println("*************************************************");

		}

		return invTM;
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

	private static double bicubicInterpolationBobDoughertyStyle(ImagePlus imp, Matrix pointInvTransformed) {
		// File: CubicFloatProcessore.java needed!!
		// Bob Dougherty 6/19/2008. Modified from FloatProcesor by Wayne Rasband.
		// http://www.optinav.com/CubicFloatResizeRotate.htm

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

    private double nearestNeighborInterpolation(ImagePlus imp2, Matrix pointInvTransformed) {

        double interpolatedPixelValue;

        this.ip2 = imp2.getProcessor();
        int w = this.ip2.getWidth();
        int h = this.ip2.getHeight();

        FileInfo fi = imp2.getFileInfo();
        double scalex = fi.pixelWidth;
        double scaley = fi.pixelHeight;

        // for the moment I use the nearest neighbor interpolation to get the z-value at x,y
        // add interpolations hooks here
        int x = (int) rint(pointInvTransformed.get(0, 0) / scalex);
        int y = (int) rint(pointInvTransformed.get(1, 0) / scaley);

        // we have no meaningful data outside the image boundaries
        if (x >= 0 && x < w && y >= 0 && y < h) {
            interpolatedPixelValue = ip2.getf(x, y);
        } else {
            interpolatedPixelValue = 0;
        }

        return interpolatedPixelValue;

    }

    /*##################################################################################################################
	*
	* 												DIALOGS
	*
	##################################################################################################################*/

	/**
	 * Show plugin configuration dialog.
	 *
	 * @return <code>true</code> when user clicked OK (confirmed changes,
	 * <code>false</code> otherwise.
	 */

	// KHK todo: I need a lot of class variables ---- try to encapsulate better!
	private boolean showDialog() {

		GenericDialog gd = new GenericDialog("KHKs jMatch3D", IJ.getInstance());

		gd.addMessage(
				"Two images of file type '32bit' or 'float' are required. \n" +
						"x and y are coordinates on a rectangular grid, \n" +
						"z represents the height information z = f(x,y).");

		String defaultItem;
		if (title1.equals(""))
			defaultItem = titles[0];
		else
			defaultItem = title1;
		gd.addChoice("Image1 (source):", titles, defaultItem);

		if (title2.equals(""))
			defaultItem = titles[1];
		else
			defaultItem = title2;
		gd.addChoice("Image2 (target):", titles, defaultItem);

		gd.addMessage("");

		gd.addChoice("ICP Point Selection Method: ", schemes, schemes[scheme]);

		gd.addMessage("refine_clamp:  removes points further away than the multiple of mean or median distance.\n"
				+ "               The smaller of the two values, mean or median, is used for this condition!\n"
				+ "refine_sd:     removes points further away then multiple of std dev distance.\n"
				+ "refine_clip:   removes percentage of points with highest distance.\n"
				+ "refine_sparse: removes percentage of all points\n"
				+ "refine_unique: all nearest points (very slow)");
		gd.addMessage("");

		// refine_clamp = kill all points farther away than a multiple of the last reported mean distance and a multiple
		// of the current median distance
		// refine_clip = kill percentage of points with highest distance
		// (approximated value from non-transformed distance values)
		// refine_sparse = kill a certain percentage of points.
		// refine_unique = all nearest points (very slow).
		gd.addMessage("Parameters: ");
		gd.addNumericField("refine_clamp: multiple of mean/med error: ", refine_clamp_par, 2);
		gd.addNumericField("refine_sd: multiple of sd: ", refine_sd_par, 2);
		gd.addNumericField("refine_clip: percentage (0.0 - 1.0): ", refine_clip_par, 2);
		gd.addNumericField("refine_sparce: percentage (0.0 - 1.0): ", refine_sparse_par, 2);
		gd.addNumericField("minimum valid points: ", minimum_valid_points, 0);
		gd.addCheckbox("a-priori: use matching landmarks as z-value truth", refine_a_priori_landmark);
		gd.addCheckbox("a-priori: use 95% approximate confidence interval", refine_a_priori_confidence_interval);
        gd.addCheckbox("a-priori: show diff-slider & histogram ", refine_a_priori_show_diff_slider);
		gd.addCheckbox("landmark mask: use coordinates with highest z-values in 5x5 mask around selected pixels",
                use_landmark_mask);
		gd.addMessage("");

		gd.addCheckbox("Check to save resulting matrix: ", saveMatrixFileFlag);
		gd.addMessage("");
		gd.addChoice("Interpolation Method (Difference Image): ", interpolationMethod, interpolationMethod[interpol]);
		gd.addMessage("");

		// KHK todo... load matrix hat noch Fehler... Abfragen ergänzen
		// KHK todo... set working directory ergänzen

		gd.addStringField("Matrix file:", file, 30);

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
		refine_sd_par = gd.getNextNumber();
		refine_clip_par = gd.getNextNumber();
		refine_sparse_par = gd.getNextNumber();
		minimum_valid_points = (int) gd.getNextNumber();
		refine_a_priori_landmark = gd.getNextBoolean();
		refine_a_priori_confidence_interval = gd.getNextBoolean();
		refine_a_priori_show_diff_slider = gd.getNextBoolean();
		use_landmark_mask = gd.getNextBoolean();
		saveMatrixFileFlag = gd.getNextBoolean();

		interpol = gd.getNextChoiceIndex();

		file = gd.getNextString();

		if (file == null || !file.equals("")) {
			loadMatrixFileFlag = true;
			transformationMatrix = (new ReadMatrix()).run(file);
		}

		if (index1 == index2) {
			IJ.error("Information:\n\nImage 1 and 2 are identical. Will continue. \nBut take care with the result!");
		}
		// plausibility check
		if (refine_clip_par < 0.0 || refine_clip_par > 1.0) {
			refine_clip_par = 0.90; // default
		}

		if (refine_sparse_par < 0.0 || refine_sparse_par > 1.0) {
			refine_sparse_par = 0.5;
		}

		if (minimum_valid_points < 3) {
			minimum_valid_points = 3;
		}

		switch (scheme) {
			case 0:
				refine_clamp = true;
				refine_sd = false;
				refine_clip = false;
				refine_sparse = false;
				refine_unique = false;
				break;
			case 1:
				refine_clamp = false;
				refine_sd = true;
				refine_clip = false;
				refine_sparse = false;
				refine_unique = false;
				break;
			case 2:
				refine_clamp = false;
				refine_sd = false;
				refine_clip = true;
				refine_sparse = false;
				refine_unique = false;
				break;
			case 3:
				refine_clamp = false;
				refine_sd = false;
				refine_clip = false;
				refine_sparse = true;
				refine_unique = false;
				break;
			case 4:
				refine_clamp = false;
				refine_sd = false;
				refine_clip = false;
				refine_sparse = false;
				refine_unique = true;
				break;
		}

		refineParameters = new ParameterICP(refine_clamp, refine_clamp_par, refine_sd, refine_sd_par, refine_clip,
				refine_clip_par, refine_sparse, refine_sparse_par, refine_unique, minimum_valid_points,
				refine_a_priori_landmark, refine_a_priori_confidence_interval, refine_a_priori_show_diff_slider,
                use_landmark_mask);
		System.out.println(refineParameters.toString());

		return true;

	}

	private void showDiffSlider(double min, double max, double defaultValue) {
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("KHKs Match3D - Histogram");
		gd.addSlider("Diff Slider", min, max, defaultValue, 0.1);
		gd.addDialogListener(this);

        if (refine_a_priori_show_diff_slider)
		    gd.showDialog();
	}

	/**
	 * Listener to modifications of the input fields of the dialog.
	 * Here the parameters should be read from the input dialog.
	 *
	 * @param gd The GenericDialog that the input belongs to
	 * @param e  The input event
	 * @return whether the input is valid and the filter may be run with these parameters
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		double zValue = gd.getNextNumber();
		subtractImages(zValue);

		impTransformedImage.updateAndRepaintWindow();
		return !gd.invalidNumber();
	}

	private void subtractImages(double zValue){
        ImageProcessor ipTransformedImage = impTransformedImage.getProcessor();
        ImageProcessor ipTransformedImageOriginal = ipOriginal.getProcessor();

        for (int i = 0; i < ipTransformedImage.getHeight(); i++) {
            for (int j = 0; j < ipTransformedImage.getWidth(); j++) {
                ipTransformedImage.setf(j, i, (float) (ipTransformedImageOriginal.getf(j, i) - zValue));
            }
        }
    }

    /*##################################################################################################################
	*
	* 										IMAGEJ HELPER METHODS
	*
	##################################################################################################################*/

	// KHK todo auslagern in utility class
	public static ImageProcessor setZeroToNan(ImageProcessor ip) {

		int width = ip.getWidth();
		int height = ip.getHeight();
		int length = width * height;

		// define an array which referes to the pixels of the image

		float[] arrayOfImagePixels = (float[]) ip.getPixels();

		for (int a = 0; a < length; a++) {
			if (arrayOfImagePixels[a] == 0.0) {
				arrayOfImagePixels[a] = Float.NaN;
			}
		}
		return ip;
	}

	private ImagePlus makeTarget(ImagePlus imp1) {
		// ip1 = baseline image which will not be transformed but we use it as a template for the transformed image.

		// I need a blank image to get the rotated img1 as we want to subtract it later from the basline image (= img1)
		// the size and pixelDimensions are identical with img1
		this.ip1 = imp1.getProcessor();

		int w = this.ip1.getWidth();
		int h = this.ip1.getHeight();

		// debugging:

		System.out.println("2014 Dimensions Image1 w x t: " + w + " x " + h);

		// Syntax: NewImage.createFloatImage(java.lang.String title, int width, int height, int slices, int options)
		ImagePlus impTransformedImage = NewImage.createFloatImage("Result-rotatedImg1", w, h, 1, NewImage.FILL_BLACK);

		// need to apply the unit scale to the new image
		impTransformedImage.copyScale(imp1);

		return impTransformedImage;

	}

    private boolean checkIJImageList() {
        // check: is there an open image?
        wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return true;
        }

        // check: we need at least 2 images for matching
        titles = new String[wList.length];
        if (wList.length < 2) {
            IJ.error("At least two match3d images (type float or Gray32) are required");
            return true;
        }

        // make a list of image titles of all open images
        for (int i = 0; i < wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null) {
                titles[i] = imp.getTitle();
            } else {
                titles[i] = "";
            }
        }

        for (int i = 0; i < wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null) {
                ImageConverter imageConverter = new ImageConverter(imp);
                imageConverter.convertToGray32();
            }
        }

        // user interaction to get matching parameters
        if (!showDialog()) {
            return true;
        }
        return false;
    }

    private boolean checkLandmarkPolygons() {
        // in this case a polygon is used to mark corresponding points we have to check that there are polygon ROIs
        // marked in both images
        // in medicine corresponding points are sometimes called fiducial marks or ficudial markers
        // I have also read the term "landmarks" in some papers

        // convention in my code:
        PolygonRoi polyRoiSource = (PolygonRoi) imp1.getRoi();
        if (polyRoiSource == null) {
            IJ.error("Error SourceImage: PlugIn requires Polygon ROI!");
            return true;
        }

        if (polyRoiSource.getType() != Roi.POLYGON) {
            IJ.error("Error SourceImage: PlugIn requires Polygon ROI!!");
            return true;
        }

        // the "source" image is usually indentified with the index number 1
        // the "target" image is usually indentified with the index number 2
        PolygonRoi polyRoiTarget = (PolygonRoi) imp2.getRoi();
        if (polyRoiTarget == null) {
            IJ.error("Error TargetImage: PlugIn requires Polygon ROI!");
            return true;
        }

        if (polyRoiTarget.getType() != Roi.POLYGON) {
            IJ.error("Error TargetImage: PlugIn requires Polygon ROI!!");
            return true;
        }

        // the number of corresponding landmarks has to be the same in both images
        int lengthSource = polyRoiSource.getNCoordinates();
        int lengthTarget = polyRoiTarget.getNCoordinates();

        if (lengthSource != lengthTarget) {
            IJ.error("Error: the number of points between the source and target polylines differs");
            return true;
        }
        return false;
    }

    private Point3d[] getCorrespondingPointListFromImageDataAsArray(ImagePlus imp) {

        int count = 0;

        ImageProcessor ip = imp.getProcessor();

        FileInfo fi;
        fi = imp.getFileInfo();
        System.out.println("Pixel-Width: " + fi.pixelWidth);

        Point3d[] vectorArray = new Point3d[(ip.getHeight() * ip.getWidth())];

        for (int i = 0; i < (ip.getHeight() * ip.getWidth()); i++) {
            vectorArray[i] = new Point3d();
        }

        for (int i = 0; i < ip.getHeight(); i++) {
            for (int j = 0; j < ip.getWidth(); j++) {

                if (ip.getf(j, i) != 0.0) {
                    vectorArray[count].z = ip.getf(j, i);
                    vectorArray[count].x = j * fi.pixelWidth;
                    vectorArray[count].y = i * fi.pixelHeight;
                    count = count + 1;
                }

                if (Float.isNaN(ip.getf(j, i))) {
                    count = count - 1;
                }
            }
        }

        /*
         * for debugging ImageProcessor ipnew;
         *
         * ipnew = new FloatProcessor(ip.getWidth(), ip.getHeight());
         *
         *
         * // Adjust brightness and contrast: ip.resetMinAndMax();
         *
         * for (int i=0; i < ip.getHeight(); i++){ for(int j= 0; j < ip.getWidth();
         * j++){ ipnew.putPixelValue(j,i,vectorArray[ip.getWidth()*i+j].z) ; //} //
         * System.out.println("    x,y und Pos: "+j + ", " + i + "und" +
         * (ip.getWidth()*i+j)); //System.out.println("   i*ip.getHeight()+j"+
         * (i*ip.getHeight()+j)); } //System.out.println("i: "+i); }
         *
         * // Show image: //new ImagePlus("XYZ_Import", ip).show(); new
         * ImagePlus("debugKHK", ipnew).show(); //new ImagePlus(fileTIF, ipnew).show();
         *
         */
        // remove zeros
        Point3d[] trimmedVectorArray = new Point3d[count];
        System.arraycopy(vectorArray, 0, trimmedVectorArray, 0, count);

        System.out.println("vectorArray-Länge: " + vectorArray.length);
        System.out.println("Active now: trimmedVectorArray: " + trimmedVectorArray.length);
        // System.out.println(Arrays.toString(trimmedVectorArray));

        return trimmedVectorArray;
    }

    // added for future reference: in case we want to read wavefront.obj files directly.
    // copied from ICPDemoApplet.java

    protected Point3d[] readVerticesFromObjectFile(URL url, double sparse) {
        double sparse_accum = 0.0;
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
            Vector<Point3d> temp = new Vector<>();
            Point3d[] ret;
            String line;
            double x = 0, y = 0, z = 0;
            while ((line = in.readLine()) != null) {
                // attempt to parse line
                // //no StringTokenizers, too much overhead
                line = line.trim(); // trim whitespace
                if (line.startsWith("v ")) { // vector
                    StringTokenizer tokens = new StringTokenizer(line);
                    tokens.nextToken(); // skip "v "
                    try {
                        x = Double.parseDouble(tokens.nextToken());
                        y = Double.parseDouble(tokens.nextToken());
                        z = Double.parseDouble(tokens.nextToken());
                    } catch (NumberFormatException e) {
                        // now what?
                        e.printStackTrace();
                        System.exit(1);
                    }
                    sparse_accum += (1.0 - sparse);
                    if (sparse_accum < 1.0) {
                        temp.add(new Point3d(x, y, z));
                    } else {
                        sparse_accum -= 1.0;
                    }
                }
            }
            in.close(); // ?
            ret = new Point3d[temp.size()]; // correct size
            temp.toArray(ret); // to preserve type Point3d[]
            return ret;
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        return null; // should never get here
    }

}