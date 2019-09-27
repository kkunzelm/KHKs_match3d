import Jama.*;
import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.plugin.*;
import ij.process.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import static java.lang.Math.*;
import java.net.URL;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import javax.swing.JFileChooser;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.vecmath.*;



// Point based rigid registration
// ==============================
// 
// The following matching procedure is based on the paper of 
// Kenichi Kanatani: Analysis of 3-D Rotation Fitting, 
// IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 16, No. 5, pp. 543ff
// 
// First the centroid of two corresponding pointlists are determined and the data are translated 
// so that the centroids are the origin. 
// Then the rotation is calculated with the SVD of the correlation matrix 
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
// the source image is the original data set which will never be changed. 
// in our special case this refers to the "baseline situation" of a dental restoration in the mouth.
// the terms "model" image - as in the model to which anything else is matched -, or also image1/img1/...,  are synonymously used for source image.
// the data image is also refered to as image2, img2 etc. This image contains the data which are transformed to place them 
// in the same position as the model/source image. In dentistry we call this the follow up image as it is a replica of the 
// restoration after years of service in the mouth.
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
 *
 * @author Prof. Dr. Karl-Heinz Kunzelmann, Ludwig Maximilians Universität, Poliklinik für Zahnerhaltung und Parodontologie, 80336 München
 *         email: karl-heinz@kunzelmann.de
 *         www:   www.kunzelmann.de
 * @version 0.1 
 * 
 */
public class Match3d_withFiducialMarkersAndICPv2_1 implements PlugIn {
    
    private static final String[] schemes = {
            "refine_clamp",
            "refine_sd",
            "refine_clip",
            "refine_sparse",
            "refine_unique"
    };
    
    private static final String[] interpolationMethod = {
            "nearest neigbor",
            "bilinear",
            "bicubic",
            "bicubic (Bob Dougherty's Style)"       // File: CubicFloatProcessore.java needed!!
                                                    // Bob Dougherty 6/19/2008.  Modified from FloatProcesor by Wayne Rasband.
                                                    // http://www.optinav.com/CubicFloatResizeRotate.htm
    };
    
    private static int interpol = 1;
    private static int scheme = 1;                  //default = clamp

    protected boolean refine_clamp = false;          // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
    protected double refine_clamp_par = 3.0;
    
    protected boolean refine_sd = true;          // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
    protected double refine_sd_par = 2.5;

    protected boolean refine_clip = false;          // clipping ab bestimmter Position im Abstandshistogram, 
    protected double refine_clip_par = 0.90;         // z. B. nur 75 % der niedrigeren Abstände berücksichtigt, 
                                                    // Rest der Punktepaare verworfen

    protected boolean refine_sparse = false;        // nur jeder  xte Punkte wird verwendet
    protected double refine_sparse_par = 0.50;

    protected boolean refine_unique = false;        // klingt als ob jeder Punkt berücksichtigt wird
    
    protected int minimum_valid_points = 800;       // Anzahl der minimal notwendigen Punkte

    public ParameterICP refineParameters;

    public String title1 = "";
    public String title2 = "";
    
    String[] titles;

    ImagePlus imp1;                     //Mnemo: imp = imageplus
    ImagePlus imp2;

    int[] wList;

    private static String file = "";
    private boolean loadMatrixFileFlag = false; 
    private boolean saveMatrixFileFlag = false;
    double[][] transformationMatrix;
    
    boolean debug = false;
    
    PolygonRoi polyRoiSource;           // convention in my code:
                                        // the "source" image is usually indentified with the index number 1
                                        // the "target" image is usually indentified with the index number 2   
    PolygonRoi polyRoiTarget;

    int lengthSource;
    int lengthTarget;
     
    ImageProcessor ip1;                 // Mnemo: ip = image processor
    ImageProcessor ip2;
 
    FileInfo fi1;
    FileInfo fi2;
    FileInfo fiTransformedImage = new FileInfo();      
   
    // neu ab 2014
    List<Vector3d> mesh1 = new ArrayList<Vector3d>();
    List<Vector3d> mesh2 = new ArrayList<Vector3d>();
    
    List<Vector3d> landmarks1 = new ArrayList<Vector3d>();
    List<Vector3d> landmarks2 = new ArrayList<Vector3d>();
    
    double[][] correlationMatrixK = new double[3][3];
    double[][] rotationMatrixR = new double[3][3];
    double[][] invRotationMatrixInvR = new double[3][3];

    /* Later to get homogeneous coordinates:
     *
     * R is extended to a 4 x 4 matrix
     * (r.. stands for rotation element ..)
     *  
     *      r11 r12 r13 0
     *      r21 r22 r23 0
     * R =  r31 r32 r33 0
     *      0   0   0   1
     * 
     * 
     * 
     * T is extended to a 4 x 1 matrix
     * (t. stands for translation element .)
     * 
     *      t1
     * T =  t2
     *      t3
     *       1 
     */
    
   
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
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp!=null) {
                titles[i] = imp.getTitle();
            }
            else {
                titles[i] = "";
            }
    }
    
     // user interaction to get matching parameters
        if (!showDialog()) {
            return;
        }
    
    

    //***********************************************************
    // some tests to be sure everything is fine before we start

    // in this case a polygon is used to mark corresponding points
    // we have to check that there are polygon ROIs marked in both images
    // in medicine corresponding points are sometimes called fiducial marks or ficudial markers
    // I have also read the term "landmarks" in some papers

    polyRoiSource = (PolygonRoi) imp1.getRoi();
    if(polyRoiSource == null) {
            IJ.error("Error SourceImage: PlugIn requires Polygon ROI!");
            return;
    } 


    if(polyRoiSource.getType() != Roi.POLYGON) {
            IJ.error("Error SourceImage: PlugIn requires Polygon ROI!!");
            return;
    } 

    polyRoiTarget = (PolygonRoi) imp2.getRoi();
    if(polyRoiTarget == null) {
            IJ.error("Error TargetImage: PlugIn requires Polygon ROI!");
            return;
    } 

    if(polyRoiTarget.getType() != Roi.POLYGON) {
            IJ.error("Error TargetImage: PlugIn requires Polygon ROI!!");
            return;
    } 

    // the number of corresponding landmarks has to be the same in both images
    lengthSource = polyRoiSource.getNCoordinates();
    lengthTarget = polyRoiTarget.getNCoordinates();

    if(lengthSource != lengthTarget) {
            IJ.error("Error: the number of points between the source and target polylines differs");
            return;
    } 
    
    // ***************** shortcut ************ ***************
    // shortcut: Matrix file loaded... can be applied directly
    if (loadMatrixFileFlag == true){
        ImagePlus impTransformedImage = makeTarget(imp1);
        impTransformedImage.setTitle("TransformedLoadedFromMatix");
        applyGeneralTransformation3D(imp1, imp2, impTransformedImage, transformationMatrix);  
        return;        
    }
    
    // ***************** lets start the real work ***************
    // 


    landmarks1 = getCorrespondingPointListFromPolygonRoi(imp1,(PolygonRoi)imp1.getRoi());
    landmarks2 = getCorrespondingPointListFromPolygonRoi(imp2,(PolygonRoi)imp2.getRoi());

    double[] centroidLandmarks1, centroidLandmarks2;

    centroidLandmarks1 = computeCentroid(landmarks1);
    centroidLandmarks2 = computeCentroid(landmarks2);

    List<Vector3d> relativeLandmarks1 = relativeCord(landmarks1, centroidLandmarks1);
    List<Vector3d> relativeLandmarks2 = relativeCord(landmarks2, centroidLandmarks2);

    if(debug == true){
        System.out.println("PolygonRoi Centroids Landmarks1: " + Arrays.toString(centroidLandmarks1));
        System.out.println("PolygonRoi Centroids Landmarks2: " + Arrays.toString(centroidLandmarks2));               
    }

    correlationMatrixK = correlationMatrixK(relativeLandmarks1, relativeLandmarks2);

    Matrix rotation = svd(correlationMatrixK);                  // singular value decomposition
    
     // for debugging
    if(debug = true){   
        System.out.println(" ");   
        System.out.println("Original Rotationsmatrix wie nach svd geliefert: ");
        System.out.println(rotation.get(0,0)+"\t"+rotation.get(0,1)+"\t"+rotation.get(0,2));
        System.out.println(rotation.get(1,0)+"\t"+rotation.get(1,1)+"\t"+rotation.get(1,2));
        System.out.println(rotation.get(2,0)+"\t"+rotation.get(2,1)+"\t"+rotation.get(2,2));
        System.out.println(" ");
    }
    
    Matrix invRotation = rotation.inverse();

    rotationMatrixR = rotation.getArrayCopy();
    invRotationMatrixInvR = invRotation.getArrayCopy();

    // for debugging
    if(debug = true){
        System.out.print("svd() RotationMatrix K = ");
        rotation.print(9, 6);
        System.out.print("svd ()inverseRotationMatrix K^-1 = ");
        invRotation.print(9, 6);

        // for debugging
        Matrix rotTimeInv = rotation.times(invRotation);
        System.out.println("Rotation times Inverted Rotation (should be Identity I): ");
        rotTimeInv.print(9, 6);
    }
   
   double[] translation = computeTranslation(centroidLandmarks1,centroidLandmarks2,rotationMatrixR);
    
   double[][] transformationMatrix = populateTransformationMatrix(rotationMatrixR,translation);
   if (debug = true){
        System.out.println("TransformationMatrix in homogeneous coordinates: "+ transformationMatrix[0][0] + "  " + transformationMatrix[0][1] + "  " + transformationMatrix[0][2] + "  " + transformationMatrix[0][3]);
        System.out.println("TransformationMatrix in homogeneous coordinates: "+ transformationMatrix[1][0] + "  " + transformationMatrix[1][1] + "  " + transformationMatrix[1][2] + "  " + transformationMatrix[1][3]);
        System.out.println("TransformationMatrix in homogeneous coordinates: "+ transformationMatrix[2][0] + "  " + transformationMatrix[2][1] + "  " + transformationMatrix[2][2] + "  " + transformationMatrix[2][3]);
        System.out.println("TransformationMatrix in homogeneous coordinates: "+ transformationMatrix[3][0] + "  " + transformationMatrix[3][1] + "  " + transformationMatrix[3][2] + "  " + transformationMatrix[3][3]);
   }
    
    /*
    // for debugging - testing my matrix multiplication
    double[] cLA = new double[3];
    cLA[0] = centroidLandmarks2[0];
    cLA[1] = centroidLandmarks2[1];
    cLA[2] = centroidLandmarks2[2];

    Matrix cL = new Matrix(cLA, 1); 
    
    Matrix c = rotation.times(cL.transpose()); // R*centroid2

    System.out.println("Matrix c = rotation.times(cL.transpose()): ");
    c.print(10,2);
    
    double[] translationM = new double[3]; //translation calculated based on matrix operations
    translationM[0]  =   centroidLandmarks1[0] -  c.get(0,0);
    translationM[1]  =   centroidLandmarks1[1] -  c.get(1,0);    
    translationM[2]  =   centroidLandmarks1[2] -  c.get(2,0);
    System.out.println("Translation von Matrix Mult: " + translationM[0] + "  " + translationM[1] + "  " + translationM[2]);
    
    // end: for debugging - testing my matrix multiplication
    */
   
    // make the target image for the translation hier: source to target
    //    ImagePlus impForwardTransformedImage = makeTarget(imp1);
    //    impForwardTransformedImage.setTitle("ForwardTransformed");
    //    applyGeneralForwardTransformation3D(imp2, impForwardTransformedImage, transformationMatrix); 
  
    /*
    ImagePlus impTransformedImage = makeTarget(imp1);
    impTransformedImage.setTitle("Transformed Landmark based");
    applyGeneralTransformation3D(imp1, imp2, impTransformedImage, transformationMatrix); 
    */
  
   // ***** KHK Grobjustierung fertig, Vorbereitung für ICP

    Point3d[] vectorArrayImg1; 
    vectorArrayImg1 = getCorrespondingPointListFromImageDataAsArray(imp1);
  
    Point3d[] vectorArrayImg2; 
    vectorArrayImg2 = getCorrespondingPointListFromImageDataAsArray(imp2);  
    
    Vector3d trans = new Vector3d();
    trans.x = translation[0];
    trans.y = translation[1];
    trans.z = translation[2];   
    
    Matrix3d rotMatrix = new Matrix3d(rotation.get(0,0), rotation.get(0,1), rotation.get(0,2),
                                      rotation.get(1,0), rotation.get(1,1), rotation.get(1,2),
                                      rotation.get(2,0), rotation.get(2,1), rotation.get(2,2));
    

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

    if (debug = true){
        System.out.println("Model points boundaries:");
        System.out.println("  "+x0);
        System.out.println("  "+x1);
        System.out.println("  "+y0);
        System.out.println("  "+y1);
        System.out.println("  "+z0);
        System.out.println("  "+z1);
    }
    
    // measure the time to build the kd-tree
    
    // long startTime = System.currentTimeMillis();
    
    modelTree.build(0, vectorArrayImg1.length-1, x0, x1, y0, y1, z0, z1);
    
    // long stopTime = System.currentTimeMillis();
    // long elapsedTime = stopTime - startTime;
    // System.out.println("Time for kd-Tree building: " + elapsedTime + "ms");
      
    Point3d modelCorner0 = new Point3d(x0, y0, z0);
    Point3d modelCorner1 = new Point3d(x1, y1, z1);
    
    ICPAlgorithm2014 icp = new ICPAlgorithm2014(); 
    
    // initialization of ICP algorithm 
     icp.init(vectorArrayImg1, vectorArrayImg2, modelTree, modelCorner0, modelCorner1, rotMatrix, trans, refineParameters);
    
    // return value from ICP
    double[][] tMat = icp.runICP();
    
    // write matrix to file
    if (saveMatrixFileFlag == true){
        new WriteMatrix().run(tMat);        
    }

    
   // make the target image for the translation hier: source to target
   //ImagePlus impForwardTransformedImageICP = makeTarget(imp1);
   //impForwardTransformedImageICP.setTitle("ForwardTransformedPostICP");
   //applyGeneralForwardTransformation3D(imp2, impForwardTransformedImageICP, tMat); 
  
    // apply transformation matrix
    ImagePlus impTransformedImageICP = makeTarget(imp1);
    impTransformedImageICP.setTitle("TransformedPostICP");
    applyGeneralTransformation3D(imp1, imp2, impTransformedImageICP, tMat);    //tMat = Transformationsmatrix

    
    //********************************************
    // Dieses Programm ist eine coole Implementierung 
    // der Projektion von irregulären xyz-Koordinaten und deren 
    // Interpolation auf ein reguläres Gitter.
    //
    //XYZ2DEM_ImporterHack hack0 = new XYZ2DEM_ImporterHack();
    //hack0.khkDisplayXYZ(vectorListImg2); 
        

}
    

/**
 * Show plugin configuration dialog.
 *
 * @return <code>true</code> when user clicked OK (confirmed changes, <code>false</code>
 *         otherwise.
 */

// KHK todo: I need a lot of class variables ---- try to encapsulate better!
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
                defaultItem = titles[1];
        else
                defaultItem = title2;
        gd.addChoice("Image2 (target):", titles, defaultItem);
        
        gd.addMessage("");
        
        gd.addChoice("ICP Point Selection Method: ",schemes,schemes[scheme]);
        

        gd.addMessage("refine_clamp:  removes points further away than the multiple of mean or median distance.\n"
                    + "               The smaller of the two values, mean or median, is used for this condition!\n"
                    + "refine_sd:     removes points further away then multiple of std dev distance.\n"
                    + "refine_clip:   removes percentage of points with highest distance.\n"
                    + "refine_sparse: removes percentage of all points\n"
                    + "refine_unique: all nearest points (very slow)");
        gd.addMessage("");
        
        // refine_clamp = kill all points farther away than a multiple of
        // the last reported mean distance and a multiple of
        // the current median distance
        // refine_clip =    kill percentage of points with highest distance
        //                  (approximated value from non-transformed distance values)
        //  refine_sparse = kill a certain percentage of points.
        //  refine_unique = all nearest points (very slow).
        gd.addMessage("Parameters: ");
        gd.addNumericField("refine_clamp: multiple of mean/med error: ", refine_clamp_par, 2);
        gd.addNumericField("refine_sd: multiple of sd: ", refine_sd_par, 2);
        gd.addNumericField("refine_clip: percentage (0.0 - 1.0): ", refine_clip_par, 2);
        gd.addNumericField("refine_sparce: percentage (0.0 - 1.0): ", refine_sparse_par, 2);
        gd.addNumericField("minimum valid points: ", minimum_valid_points, 0);
        gd.addMessage("");

        gd.addCheckbox("Check to save resulting matrix: ", saveMatrixFileFlag);
        gd.addMessage("");
        gd.addChoice("Interpolation Method (Difference Image): ",interpolationMethod,interpolationMethod[interpol]);
        gd.addMessage("");
        
        // KHK todo... load matrix hat noch Fehler... Abfragen ergänzen
        // KHK todo... set working directory ergänzen
        
        gd.addStringField("Matrix file:",file,30);
                 
        
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
        minimum_valid_points = (int)gd.getNextNumber();
        saveMatrixFileFlag = gd.getNextBoolean();
        
        interpol = gd.getNextChoiceIndex();
        
        file = gd.getNextString();
        

        if (file == null || file.equals("") == false){ 
            loadMatrixFileFlag = true;
            transformationMatrix = (new ReadMatrix()).run(file);
        }
                
        if (index1 == index2){
            IJ.error("Information:\n\nImage 1 and 2 are identical. Will continue. \nBut take care with the result!");
        }
        // plausibility check
        if (refine_clip_par < 0.0 || refine_clip_par > 1.0){
            refine_clip_par = 0.90;         // default
        }
        
        if (refine_sparse_par < 0.0 || refine_sparse_par > 1.0){
            refine_sparse_par = 0.5;
        }
        
        if (minimum_valid_points < 3){
            minimum_valid_points = 3;
        }
        
        switch(scheme){
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
        
        refineParameters = new ParameterICP(refine_clamp, refine_clamp_par, refine_sd, refine_sd_par, refine_clip, refine_clip_par, refine_sparse, refine_sparse_par, refine_unique, minimum_valid_points);
        System.out.println(refineParameters.toString());
        
        return true;

        }
 

public Point3d[] getCorrespondingPointListFromImageDataAsArray(ImagePlus imp) {

        int count = 0;
        
        ImageProcessor ip = imp.getProcessor();
        
     
        FileInfo fi;
        fi = imp.getFileInfo();
        System.out.println("Pixel-Width: "+ fi.pixelWidth);
       
        Point3d vectorArray[] = new Point3d[(ip.getHeight()*ip.getWidth())];
              
        for(int i = 0; i < (ip.getHeight()*ip.getWidth()); i++) {
            vectorArray[i] = new Point3d();
        }


        for (int i=0; i < ip.getHeight(); i++){
            for(int j= 0; j < ip.getWidth(); j++){

                 //   vectorArray[ip.getWidth()*i+j].z = ip.getf(j,i);
                 //   vectorArray[ip.getWidth()*i+j].x = j*fi.pixelWidth;
                 //   vectorArray[ip.getWidth()*i+j].y = i*fi.pixelHeight;
                if (ip.getf(j,i) != 0.0 ){ 
                    

                        vectorArray[count].z = ip.getf(j,i);
                        vectorArray[count].x = j*fi.pixelWidth;
                        vectorArray[count].y = i*fi.pixelHeight;
                        count = count +1;
                    
                }
                
                if (Float.isNaN(ip.getf(j,i))){
                    count = count -1;
                }
                
                /*
                if (i % 20 == 0 && j % 20 == 0){
                    System.out.println("i: " + i + " j: " + j + " z-Wert: " + ip.getf(j,i)+ " i*Breite+j: " + (ip.getWidth()*i+j)
                    + " count: " + count);
                }
                */

                   // System.out.println("    x,y und Pos: "+j + ", " + i + "und" + (ip.getWidth()*i+j)); 
                    //System.out.println("   i*ip.getHeight()+j"+ (i*ip.getHeight()+j));
            }
            //System.out.println("i: "+i);
      } 
        
        /* for debugging
        ImageProcessor ipnew;
                                
                ipnew = new FloatProcessor(ip.getWidth(), ip.getHeight());
                

                // Adjust brightness and contrast:
                ip.resetMinAndMax();
                
                for (int i=0; i < ip.getHeight(); i++){
                    for(int j= 0; j < ip.getWidth(); j++){
                            ipnew.putPixelValue(j,i,vectorArray[ip.getWidth()*i+j].z) ;
                        //}
                           // System.out.println("    x,y und Pos: "+j + ", " + i + "und" + (ip.getWidth()*i+j)); 
                            //System.out.println("   i*ip.getHeight()+j"+ (i*ip.getHeight()+j));
                    }
                    //System.out.println("i: "+i);
                } 
                
                // Show image:
                //new ImagePlus("XYZ_Import", ip).show();
                new ImagePlus("debugKHK", ipnew).show();
                //new ImagePlus(fileTIF, ipnew).show();
                
                */
        // remove zeros
       Point3d trimmedVectorArray[] = new Point3d[count];
       System.arraycopy(vectorArray, 0 , trimmedVectorArray, 0, count);
       
       System.out.println("vectorArray-Länge: "+vectorArray.length);
       System.out.println("Active now: trimmedVectorArray: "+trimmedVectorArray.length);
      // System.out.println(Arrays.toString(trimmedVectorArray));
       
       return trimmedVectorArray;
    }

public void getCorrespondingPointBruteForce(ImagePlus imp1, ImagePlus imp2) {
        // noch nicht fertig
              
        this.ip1 = imp1.getProcessor();
        this.ip2 = imp2.getProcessor();
        FileInfo fi;
        fi = imp1.getFileInfo();
        System.out.println("Pixel-Width: "+ fi.pixelWidth);
       
        Point3d vectorArray1[] = new Point3d[(ip1.getHeight()*ip1.getWidth())];
        Point3d vectorArray2[] = new Point3d[(ip2.getHeight()*ip2.getWidth())];
              
        
        for(int i = 0; i < (ip1.getHeight()*ip1.getWidth()); i++) {
            vectorArray1[i] = new Point3d();
        }

        for(int i = 0; i < (ip2.getHeight()*ip2.getWidth()); i++) {
            vectorArray2[i] = new Point3d();
        }

        int i,j;
        double minDistance = Double.POSITIVE_INFINITY;
        for (i = 0; i < vectorArray1.length; i++) {
            for (j = 0; j < vectorArray2.length; j++) {
                //System.out.println("j: "+j);
                double dist = (vectorArray1[i].x - vectorArray2[j].x) * (vectorArray1[i].x - vectorArray2[j].x) + (vectorArray1[i].y - vectorArray2[j].y) * (vectorArray1[i].y - vectorArray2[j].y) + (vectorArray1[i].z - vectorArray2[j].z) * (vectorArray1[i].z - vectorArray2[j].z);
                if (minDistance > dist){
                    minDistance = dist;
 
                }
                                   
            }
            System.out.println("minDist(i,j) is: "+ minDistance + " at " + i + "/" + j);
        }
         System.out.println();
            System.out.println("Ready!");
           
    }


/**
     * Returns ArrayList of Vector3d 
     * 
     * @param imp
     * @param polyRoi
     * @return
     */
public List<Vector3d> getCorrespondingPointListFromPolygonRoi(ImagePlus imp, PolygonRoi polyRoi) {
        List<Vector3d> vectorList;
        vectorList = new ArrayList<Vector3d>();
        
        ImageProcessor ip = imp.getProcessor();
        FileInfo fi;
        
        fi = imp.getFileInfo();
        
        if (fi.fileType !=FileInfo.GRAY32_FLOAT) {
                    IJ.error("Error Image1 is not type 'float' - will try to convert");
                    ip = ip.convertToFloat();
                }
        
        int[] roiX;
        int[] roiY;
        int roiLength = polyRoi.getNCoordinates();
        
        // later I need the absolute coordinates. 
        // the method getXCoordinates or getYcoordinates
        // returns relative coordinates, relative to the bouncing box of the ROI.

        // Rectangle is the bouncing box
        Rectangle rectRoiSource = ip.getRoi();

        // Source Image: 
        // Creating array for SVD
        // 

        roiX = polyRoi.getXCoordinates();
        roiY = polyRoi.getYCoordinates();
        
        System.out.println("***************************************************");
        System.out.println("PolygonRoi-Koordinaten: ");
        
        for(int i=0; i < roiLength; i++) {
            Vector3d vert = new Vector3d();
            vert.x = (float)fi.pixelWidth*(roiX[i]+ rectRoiSource.x);
            vert.y = (float)fi.pixelHeight*(roiY[i]+  rectRoiSource.y);
            vert.z = ip.getf((roiX[i]+rectRoiSource.x),(roiY[i]+rectRoiSource.y));   
            vectorList.add(vert); 
            System.out.println(vert.x+"\t"+vert.y+"\t"+vert.z);
        }
        System.out.println("***************************************************");
        return vectorList;
    }

public List<Vector3d> getCorrespondingPointListFromImageData(ImagePlus imp) {
        List<Vector3d> vectorList;
        vectorList = new ArrayList<Vector3d>();
        
        int step = 10; // quick and dirty data reduction
        
        ImageProcessor ip = imp.getProcessor();
        FileInfo fi;
        
        fi = imp.getFileInfo();
        
        if (fi.fileType !=FileInfo.GRAY32_FLOAT) {
                    IJ.error("Error Image1 is not type 'float' - will try to convert");
                    ip = ip.convertToFloat();
                }
        
       
          for (int i=0; i < ip.getHeight(); i = i + step){
            for(int j= 0; j < ip.getWidth(); j = j + step){
                Vector3d vert = new Vector3d();
                
                vert.x = (float)(j*fi.pixelWidth);
                vert.y = (float)(i*fi.pixelHeight);
                if (ip.getf(j,i) > 0.001){              // quick and dirty: remove negative numbers, zero and those close to zero... think for a more generic solution
                    vert.z  = ip.getf(j,i);
                }
                
                vectorList.add(vert); 

            }
      } 
        
        return vectorList;
    }

    /** Accumulates the Correlation Matrix K which is needed for Singular Value Decomposition
     *
     * @param relativeLandmarks1 Landmark list number 1 
     * @param relativeLandmarks2 Landmark list number 2 
     * @return correlationMatrix in Kanatani's paper called K
     * 
     */
public double[][] correlationMatrixK(List<Vector3d> relativeLandmarks1, List<Vector3d> relativeLandmarks2){
        
        /* this method implements formula (3) of Kanatani
         *  the weight is set to "1". It is a scaling factor only.
         *  Formula (3):
         *
         *      K = SUM(vectorPoint1 * vectorPoint'1^T)  --> ^T means transposed
         *
         *  in detail:
         *
         *  Each point consists of 3 coordinates x, y, z
         *  Each point is treated as a vector.
         *  The coordinates of these vectors are arranged vertically
         *           ( x1 )
         *  Point1 = ( y1 )
         *           ( z1 )
         * 
         *  Point'1^T = (x'1, y'1, z'1) 
         * 
         *  The mark "'" means the corresponding point in the second image
         * 
         *  The multiplication of the two 3x1 vectors results 
         *  in a 3x3 matrix = correlationMatrixK.
         *  The matrices off all points are summarized and result in K
         * 
         * KHK Added: 19.1.13
         * from: http://nghiaho.com/?page_id=671
         * Pay close attention to the transpose symbol. It’s doing a multiplication 
         * between 2 matrices where the dimensions effectively are, 3×1 and 1×3, 
         * respectively. The ordering of the multiplication is also important, 
         * doing it the other way will find a rotation from B to A instead.
         * 
         * 
        */
        
        double[][] modelPoint = new double[3][3];               // model = baseline, unchanged, reference
        double[][] dataPoint = new double[3][3];                // data = follow-up
        double[][] correlationMatrixK_temp = new double[3][3];

        
        for (int x = 0; x < relativeLandmarks1.size(); x++){
            
            modelPoint[0][0] = relativeLandmarks1.get(x).getX();
            modelPoint[1][0] = relativeLandmarks1.get(x).getY();
            modelPoint[2][0] = relativeLandmarks1.get(x).getZ();
            dataPoint[0][0] = relativeLandmarks2.get(x).getX();
            dataPoint[0][1] = relativeLandmarks2.get(x).getY();
            dataPoint[0][2] = relativeLandmarks2.get(x).getZ();
           
            System.out.println("Accumulation of the correlationMatrixK - step no: x = "+x);
        
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        // System.out.println("ijk: "+i+"/"+j+"/"+k+" -> Matrixelement ik before "+ correlationMatrixK[i][k]);
                        correlationMatrixK_temp[i][k] += modelPoint[i][j] * dataPoint[j][k];
                        // System.out.println("ijk: "+i+"/"+j+"/"+k+" -> Matrixelement ik after "+ correlationMatrixK[i][k]);
                    }
                }
            }
            
        }
        
        return correlationMatrixK_temp;
    }
    
    
    /** Singular Value Decomposition of the correlation matrix. 
     *
     * @param correlationMatrixK
     * @return rotation rotation matrix 3x3
     */
public Matrix svd(double[][] correlationMatrixK){

      /* To perform the singular value decomposition svd the package Jama is needed
       * As soon as we have the correlation matrix, the next step in the 
       * algorithm as described by Kanatani is to decompose the correlation matrix
       * (see formulae 12 and 13) 
       * 
       * It is important to mention that Kanatani writes: K = VAU^T
       * In Jama there is a svd.jama file which uses a slightly different syntax A = USV^T
       * Please concentrate not to mix the meaning of U and V
       * Kanatini's U will be the V of SVD
       * to stay in the logic of Jama we write for formula (13)
       * 
       *            (1 0    0       )
       *    R = U * (0 1    0       )*V^T           formula 13 later called intermedResult
       *            (0 0 det(UV^T)  )
       * 
       * Just to mention a few details:
       * 
       * Kanatani writes in formula (4) that to minimize the distance between the two pointsets
       * for the given correlation matrix (formula 3), the trace(R^T*K) musst be maximized
       * in formula (13) he states - and proves too - when the trace(R^T*K) is maximized 
       * 
       * What we do is, we use K make a SVD and get U and V
       * Using U and V we can calculate R with formula 13 
       * R is our final result for the rotation
       * 
       * to calculate R we need the intermediate result for det(VU^T) - Matrix.de()
       * then we do matrix multiplication
       *           
       *        ( 1  0  0     )
       * 1) U * ( 0  1  0     )
       *        ( 0  0  det(VU^T))
       * 
       * 2) |_____________________|
       * 
       *            = C * V^T
       *           |___________|
       * 
       *                = R
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
      
      
      // "det(VU^T)" will be calculated
      
      // s.getU().transpose();
      // s.getV().times(s.getU().transpose());
      // s.getV().times(s.getU().transpose()).det();  // this is "det(VU^T)"
      
      System.out.println("det(VU^T): "+s.getV().times(s.getU().transpose()).det());
      
      Matrix intermedResult = new Matrix(new double[][] { { 1,  0, 0 },
                                             { 0,  1, 0 },
                                             { 0,  0, s.getV().times(s.getU().transpose()).det() } });
      
      // s.getU().times(intermedResult1);  // is matrix C as above mentioned
      // s.getV().transpose();             // is matrix V^T
      // s.getU().times(intermedResult1).times(s.getV().transpose());  // is R
      
      
      Matrix rotation = s.getU().times(intermedResult).times(s.getV().transpose());

      return rotation;
      
   }

    /** Calculate the Translation between the two datasets
     *
     * @param centroid1 Source image center
     * @param centroid2 Target image center
     * @param rotationMatrix R
     * @return translation Translation vector 
     */
private double[] computeTranslation(double[] centroid1, double[] centroid2, double[][] rotationMatrixR){
        
        double[] translation = new double[3];
        
        double[] rotCentroid2 = new double[3];
        
       
        // Calculate the translation.
        // Important: to calculate the translation, the centroid of the target image, which will be transformed, has to be rotated!!!
        // T = Sc - Tc*R
        // Sc = Source image center, Tc = Target image center
        
        //Lets assume now that we have a matrix R. When P is multiplied by R it transforms P to PT. Considering what we know about matrix multiplication lets see how we can re-write a point-matrix multiplication and isolate the computation of each of the transformed point coordinates:
        // PT.x=P.x∗R00+P.y∗R10+P.z∗R20
        // PT.y=P.x∗R01+P.y∗R11+P.z∗R21
        // PT.z=P.x∗R02+P.y∗R12+P.z∗R22
       
       
       
       rotCentroid2[0] = (float) (rotationMatrixR[0][0]*centroid2[0] + rotationMatrixR[0][1]*centroid2[1] + rotationMatrixR[0][2]*centroid2[2]);
       rotCentroid2[1] = (float) (rotationMatrixR[1][0]*centroid2[0] + rotationMatrixR[1][1]*centroid2[1] + rotationMatrixR[1][2]*centroid2[2]);
       rotCentroid2[2] = (float) (rotationMatrixR[2][0]*centroid2[0] + rotationMatrixR[2][1]*centroid2[1] + rotationMatrixR[2][2]*centroid2[2]);          

       translation[0] = centroid1[0] -  rotCentroid2[0];
       translation[1] = centroid1[1] -  rotCentroid2[1];
       translation[2] = centroid1[2] -  rotCentroid2[2];
      
       System.out.println("Centroid1: " + centroid1[0] + " " + centroid1[1] + " " + centroid1[2]);
       System.out.println("Centroid2: " + centroid2[0] + " " + centroid2[1] + " " + centroid2[2]);
       // System.out.println("diff Centroids: " + (centroid2[0]-centroid1[0]) + "  " + (centroid2[1]-centroid1[1]) + "   " + (centroid2[2]-centroid1[2]));
       // System.out.println("rotCentroid2: " + rotCentroid2[0] + " " + rotCentroid2[1] + " " + rotCentroid2[2]);       
       System.out.println("Translation: "+ translation[0] + ", " + translation[1] + ", " + translation[2]);

       return translation;
        
    }
    
private double[][] populateTransformationMatrix(double[][] rotationMatrixR, double[] translation){
        
        double[][]transformationMatrix = new double[4][4];
       
        for (int k=0; k<3; k++){
            System.arraycopy(rotationMatrixR[k], 0, transformationMatrix[k], 0, 3);
        }
        transformationMatrix[3][0]=0;
        transformationMatrix[3][1]=0;
        transformationMatrix[3][2]=0;
        
        transformationMatrix[0][3]=translation[0];
        transformationMatrix[1][3]=translation[1];
        transformationMatrix[2][3]=translation[2];
        transformationMatrix[3][3]=1.0d;
        
        if (debug = true){
            System.out.println("Transformation Matrix formated for TransformJ: ");
            System.out.println(transformationMatrix[0][0]+"\t"+ transformationMatrix[0][1]+"\t"+ transformationMatrix[0][2]+"\t"+ transformationMatrix[0][3]);
            System.out.println(transformationMatrix[1][0]+"\t"+ transformationMatrix[1][1]+"\t"+ transformationMatrix[1][2]+"\t"+ transformationMatrix[1][3]);
            System.out.println(transformationMatrix[2][0]+"\t"+ transformationMatrix[2][1]+"\t"+ transformationMatrix[2][2]+"\t"+ transformationMatrix[2][3]);        
            System.out.println(transformationMatrix[3][0]+"\t"+ transformationMatrix[3][1]+"\t"+ transformationMatrix[3][2]+"\t"+ transformationMatrix[3][3]);
        }
        return transformationMatrix;        
    } 
        
public ImagePlus makeTarget(ImagePlus imp1){
        // ip1 = baseline image which will not be transformed but we use it as a template for the transformed image.
           
        // I need a blank image to get the rotated img1
        // as we want to subtract it later from the basline image (= img1) 
        // the size and pixelDimensions are identical with img1
        this.ip1 = imp1.getProcessor();

        int w = this.ip1.getWidth();
        int h = this.ip1.getHeight();
 
        // debugging:
        
        System.out.println("2014 Dimensions Image1 w x t: "+w+" x "+h);
        
        //Syntax: NewImage.createFloatImage(java.lang.String title, int width, int height, int slices, int options) 
        ImagePlus impTransformedImage = NewImage.createFloatImage("Result-rotatedImg1",w,h,1,NewImage.FILL_BLACK);
                
        // need to apply the unit scale to the new image

        impTransformedImage.copyScale(imp1);

        
        return impTransformedImage;

}
    
public void applyGeneralForwardTransformation3D(ImagePlus imp2, ImagePlus impForwardTransformedImage, double[][] transformationMatrix){

        ImageProcessor ipTransformedImage = impForwardTransformedImage.getProcessor();
        FileInfo fi = imp2.getFileInfo();
        this.ip2 = imp2.getProcessor();

        int w = this.ip2.getWidth();
        int h = this.ip2.getHeight();
        
        Matrix transformation = new Matrix(transformationMatrix);
         
        for (int v=0; v<h; v++){
            for (int u=0; u<w; u++){
                
                
                double scale = fi.pixelWidth; // pixelHeigth and pixelWidth must be the same
                //set(int i,int j,double s) - Achtung: Zeile dann Spalte
                Matrix point = new Matrix(4,1);
                Matrix pointTransformed;      // = new Matrix(4,1);
                point.set(0,0,u*scale);
                point.set(1,0,v*scale);
                point.set(2,0,ip2.getf(u,v));
                point.set(3,0,1.0);
                

                
                // times(Matrix B) Linear algebraic matrix multiplication, A * B
                
                pointTransformed = transformation.times(point);
  
                int x = (int) rint(pointTransformed.get(0,0)/scale);
                int y = (int) rint(pointTransformed.get(1,0)/scale);
                                
                if(x > 0 && x < w && y > 0 && y < h){
                
                 ipTransformedImage.setf(x,y,(float) pointTransformed.get(2,0));
                }
                
                /*
                if (v % 20 == 1 && u % 20 == 1){
                    System.out.println("u, v, z', x,y,z vor nach Transformation: ");
                    System.out.println(u + ", " + v + ", " + ip2.getf(u,v) + " = " + x + "/" + pointTransformed.get(0,0) + " " + y + "/" 
                           +pointTransformed.get(1,0) + " " + pointTransformed.get(2,0));
                }
                */
                
                
           }
        }
        // update the rotated image
        ipTransformedImage.resetMinAndMax();
        double min = ipTransformedImage.getMin();
        double max = ipTransformedImage.getMax();
        ipTransformedImage.setMinAndMax(min, max);
        impForwardTransformedImage.updateAndDraw();
        impForwardTransformedImage.show();
        //System.out.println ("TransformationMatrix: "+ Arrays.deepToString(transformationMatrix));
}


/**
     * DESCRIPTION:                                                        
     * Takes two input datasets and the location parameters as a transformation matrix 
     * in homogeneous coordinates and computes the differences per pixel between the rotated and translated patterns. 
     * The difference is displayed in a new image.
     * <p>
     * A target-to-source-mapping is applied.
     * 
     * @param imp1
     * @param imp2
     * @param differenceImage
     * @param transformationMatrix
     *
     */
public void applyGeneralTransformation3D(ImagePlus imp1, ImagePlus imp2, ImagePlus differenceImage, double[][] transformationMatrix){

      
        // target-to-source mapping as described for 2d in Burger and Burge (http://imagingbook.com/)

        // impTransformedImage is a blank image plus template same size as the "target" image (here ip1)
        ImageProcessor ipTransformedImage = differenceImage.getProcessor();
        
        FileInfo fi = imp1.getFileInfo();
        this.ip1 = imp1.getProcessor(); 
        this.ip1 = setZeroToNan(ip1);       // masking: zero = NaN
        int w = this.ip1.getWidth();
        int h = this.ip1.getHeight();
        
        this.ip2 = imp2.getProcessor();
        this.ip2 = setZeroToNan(ip2);       // masking: zero = NaN

        double tempz = 0.0;
        float diff;
        boolean printInfos = false;          // used to switch on/off verbose printout during processing
        
        double scalex = fi.pixelWidth; 
        double scaley = fi.pixelHeight;

        for (int v=0; v<h; v++){            // rows of the image
            for (int u=0; u<w; u++){        // columns of the image
                
                //set(int i,int j,double s) - first index row, second column!
                Matrix pointTarget = new Matrix(4,1);  // I make a matrix from the array to use available JAMA matrix multiplication later
                Matrix pointInvTransformed;         // = new Matrix(4,1) for the point after the inverse transformation;
                
                // see Peter Neugebauer Match3d: Feinjustierung von Tiefenbildern zur Vermessung von kleinen Verfomungen, Diplomarbeit 1991.
                // The baseline image is mapped with the inverse transformation matrix close to the follow-up image (= source): f(u,x) -> f'(x,y)
                // At x,y there exists usually no value for the source data, therefore this value has to be interpolated.
                // The difference is calculated and mapped using the transformation from x,y back to u,v.
                // In our case we do not really need to apply the transformation as we just have to use u,v which is the same after all.
                // The difference is now in the target image (= differenceImage) which has the same coordinates as the baseline image.
                
                pointTarget.set(0,0,u*scalex);
                pointTarget.set(1,0,v*scaley);
                pointTarget.set(2,0,this.ip1.getf(u,v));           
                pointTarget.set(3,0,1.0);
               
                // from JAMA API: matrixA.times(Matrix B) Linear algebraic matrix multiplication, A * B
       
                // I calculate the transformation matrix myself
                // Transformation T and inverse Transformation T^-1 in homogeneous coordinates:
                // R = the rotations part (upper left 3x3 matrix)
                // t = the translation part (upper right 3x1 matrix)
                // 
                // T = ( R   t )        R^-1 = R^T, t^-1 = -R^T*t      T^-1 = ( R^T  -R^T*t)
                //     ( 0   1 )                                              (  0      1  )
                //
                // the method getInverseTransformationMatrix should work correctly
                // T and T^-1 tested with TransformJ reveals the expected results and T^-1*T = I 
                Matrix invTM = getInverseTransformationMatrix(transformationMatrix, printInfos);
                
                // the target point at position u,v is mapped with the inverse transformation to position x,y
                pointInvTransformed = invTM.times(pointTarget);

                switch(interpol){
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

                diff = (float) (pointInvTransformed.get(2, 0) - tempz);       //difference

                ipTransformedImage.setf(u,v,diff);
                
                printInfos = false;     // need to know the details only once
           }
        }
        // KHK todo move to separate method
        // update the rotated image
        
        // Tipp: man kann mit Threshold im 32 bit float Bild das gleiche machen, wie bei Wolfram mit clip-Histogram 
        //       Hintergrund wird zu NaN (anklicken).
        ipTransformedImage.resetMinAndMax();
        ContrastEnhancer ce = new ContrastEnhancer();  // KHK todo check again whether this is what I expect it to be
        ce.stretchHistogram(ipTransformedImage, 1.0);  // ist eine Art Histogram Threshold. Es werden die obersten und untersten Pixel entfernt. 
                                                       // unklar ist, ob 1% von oben und 1% von unten oder 1%/2 von oben und 1%/2 von unten???
        double min = ipTransformedImage.getMin();
        double max = ipTransformedImage.getMax();
        ipTransformedImage.setMinAndMax(min, max);
        differenceImage.updateAndDraw();
        differenceImage.show();
        //System.out.println ("TransformationMatrix: "+ Arrays.deepToString(transformationMatrix));
}
 

private double nearestNeighborInterpolation(ImagePlus imp2, Matrix pointInvTransformed){
    
    double interpolatedPixelValue;
    
    this.ip2 = imp2.getProcessor();
    int w = this.ip2.getWidth();
    int h = this.ip2.getHeight();
    
    FileInfo fi = imp2.getFileInfo();
    double scalex = fi.pixelWidth; 
    double scaley = fi.pixelHeight;
    
    // for the moment I use the nearest neighbor interpolation to get the z-value at x,y
    // add interpolations hooks here
    int x = (int) rint(pointInvTransformed.get(0,0)/scalex);
    int y = (int) rint(pointInvTransformed.get(1,0)/scaley);

    // we have no meaningful data outside the image boundaries
    if (x >= 0 && x < w && y >= 0 && y < h){
        interpolatedPixelValue = ip2.getf(x,y);     
    }
    else {
        interpolatedPixelValue = 0;
    }
        
    return interpolatedPixelValue;
    
}


private static double bilinearInterpolation(ImagePlus imp, Matrix pointInvTransformed){
    
    double interpolatedPixelValue;
    
    ImageProcessor ip = imp.getProcessor();
    int w = ip.getWidth();
    int h = ip.getHeight();
    
    FileInfo fi = imp.getFileInfo();
    double scalex = fi.pixelWidth; 
    double scaley = fi.pixelHeight;
    
    double x = pointInvTransformed.get(0,0)/scalex;
    double y = pointInvTransformed.get(1,0)/scaley;
    
    // bilinear is 1, bicubic would be 2
    ip.setInterpolationMethod(1);
     
    // we have no meaningful data outside the image boundaries
    if (x >= 0 && x < w && y >= 0 && y < h){
        interpolatedPixelValue =  ip.getInterpolatedValue(x, y);     
    }
    else {
        interpolatedPixelValue = 0;
    }
    
    return interpolatedPixelValue;
    
}

private static double bicubicInterpolation(ImagePlus imp, Matrix pointInvTransformed){
    
    double interpolatedPixelValue;
    
    ImageProcessor ip = imp.getProcessor();
    int w = ip.getWidth();
    int h = ip.getHeight();
    
    FileInfo fi = imp.getFileInfo();
    double scalex = fi.pixelWidth; 
    double scaley = fi.pixelHeight;
    
    double x = pointInvTransformed.get(0,0)/scalex;
    double y = pointInvTransformed.get(1,0)/scaley;
    
    // bilinear is 1, bicubic would be 2
    ip.setInterpolationMethod(2);

    // we have no meaningful data outside the image boundaries
    if (x >= 0 && x < w && y >= 0 && y < h){
        interpolatedPixelValue = ip.getInterpolatedValue(x, y);     
    }
    else {
        interpolatedPixelValue = 0;
    }
        
    return interpolatedPixelValue;
    
}

// File: CubicFloatProcessore.java needed!!
// Bob Dougherty 6/19/2008.  Modified from FloatProcesor by Wayne Rasband.
// http://www.optinav.com/CubicFloatResizeRotate.htm
private static double bicubicInterpolationBobDoughertyStyle(ImagePlus imp, Matrix pointInvTransformed){
    
    double interpolatedPixelValue;
    
    ImageProcessor ip = imp.getProcessor();
  
    // bicubic interpolation developed by Bob Dougherty
    CubicFloatProcessor cfp = new CubicFloatProcessor(ip.getFloatArray());
    
    int w = ip.getWidth();
    int h = ip.getHeight();
    
    FileInfo fi = imp.getFileInfo();
    double scalex = fi.pixelWidth; 
    double scaley = fi.pixelHeight;
    
    double x = pointInvTransformed.get(0,0)/scalex;
    double y = pointInvTransformed.get(1,0)/scaley;
     
    // we have no meaningful data outside the image boundaries
    if (x >= 0 && x < w && y >= 0 && y < h){
        interpolatedPixelValue =  cfp.getInterpolatedPixel(x, y);     
    }
    else {
        interpolatedPixelValue = 0;
    }
        
    return interpolatedPixelValue;
    
}
    

// not used at the moment
private double bilinearInterpolationKH(ImagePlus imp2, Matrix pointInvTransformed){
        
    double interpolatedPixelValue;
    
    this.ip2 = imp2.getProcessor();
    int w = this.ip2.getWidth();
    int h = this.ip2.getHeight();
    
    FileInfo fi = imp2.getFileInfo();
    double scalex = fi.pixelWidth; 
    double scaley = fi.pixelHeight;
    
    double x = pointInvTransformed.get(0,0)/scalex;
    double y = pointInvTransformed.get(1,0)/scaley;
    
    int     x0, y0, x1, y1;
    
    double p,q;


    // Outside the original image
    if (x < 0)
      x = 0;
    else if (x >= h-1)
      x = h - 2;
    if (y < 0)
      y = 0;
    else if (y >= w-1)
      y = w - 2;

    // Bilinear interpolation
    x0 = (int)Math.floor(x);
    y0 = (int)Math.floor(y);
    x1 = x0 + 1;
    y1 = y0 + 1;
    p = ((x1-x)*ip2.getf(x0,y0) + (x-x0)*ip2.getf(x1,y0));
    q = ((x1-x)*ip2.getf(x0,y1) + (x-x0)*ip2.getf(x1,y1));
    interpolatedPixelValue = (y1-y)*p + (y-y0)*q;

    return interpolatedPixelValue;
}

// KHK todo auslagern in utility class 
public static ImageProcessor setZeroToNan(ImageProcessor ip) {
	
    int width = ip.getWidth();
    int height = ip.getHeight();
    int length = width * height;

    // define an array which referes to the pixels of the image

    float[] arrayOfImagePixels = (float[])ip.getPixels();


    for (int a=0; a < length; a++) {
        if (arrayOfImagePixels[a] == 0.0){
            arrayOfImagePixels[a] = Float.NaN;
        }  
    }
return ip;
} 

// KHK todo auslagern in utility class 
public static ImageProcessor setNanToZero(ImageProcessor ip) {

    int width = ip.getWidth();
    int height = ip.getHeight();
    int length = width * height;

    // define an array which referes to the pixels of the image

    float[] arrayOfImagePixels = (float[])ip.getPixels();


    for (int a=0; a < length; a++) {

        // if (arrayOfImagePixels[a] == 0.0)
        if (Float.isNaN(arrayOfImagePixels[a])){
            arrayOfImagePixels[a] = 0;
        }
     } 
    return ip;
}

public float[] adjustHistogram(ImageProcessor ip){
    float[] minMaxHistogram = new float[2];
    
    // KHK todo 
    
    return minMaxHistogram;
}

public static Matrix getInverseTransformationMatrix(double[][] transformationMatrix, boolean printInfos){
    
    // KHK ich habe mit TransformJ kontrolliert: das Ergebnis sollte stimmen.
    //
    // die übliche Transformationsmatrix T 
    //
    // (R t)
    // (   ) p  = t*R*p    (erst R um Umsprung, dann t)  
    // (0 1)
    //
    //entspricht einer Kombination aus R und t:
    // with R being a rotation and t being a translation is a combined transformation
    // man muss die Rotation R und Translation t in der TransformationsMatrix T
    // separat betrachten.
    // aufgrund der Orthogonalität der Rotationsachsen gilt: R^T = R^-1
    // d.h. wenn man die Rotationsmatrix transformiert, hat man die Inverse
    // für die Translation t gilt: 
    // inv(t) = -t
    // die Rotationsmatrixe kann man einfach invertieren
    // ABER: für Koordinaten gilt
    // p' = Rp + t
    // inverse Transformation zur Rückgewinnung von p:
    // p' = Rp + t entspricht p'-t = Rp  entspricht R^-1(p'-t) = R^-1*R*p
    // entspricht wegen R^-1*R = I: p = R^-1(p'-t)
    // in der homogenen TransformationsMatrix hat das die Form:
    // (R  t)^-1          (R^T  -R^T*t)
    // (0  1)     gleich  (0      1   )
    //
    // ^-1 = Inverse
    // ^T = Transponierte
    
    // nebenbei: es gilt immer... erst in Ursprung verschieben, dann rotieren, dann weiter translatieren
    // t2*R*t1 ... bei Spaltenvektoren/Matrizen von rechts nach links abgearbeitet.
    // ABER: die Umkehrung von Matrixverkettungen ist
    // What is the inverse of a sequence of transformations? 
    // (M1M2...Mn)-­‐1 = Mn-­‐1Mn-­‐1-­‐1...M1-­‐1 
    // Inverse of a sequence of transformations is the composition of the inverses of each transformation in reverse order.
    // What will our sequence look like? 
    //(T^‐1RT)^-1 = T^‐1R^‐1T
    // We still translate to the origin first, then translate back at the end! 

     
    Matrix invTM = new Matrix(4,4);
    // R^T part (R transposed)
    invTM.set(0,0,transformationMatrix[0][0]);
    invTM.set(0,1,transformationMatrix[1][0]);
    invTM.set(0,2,transformationMatrix[2][0]);
    invTM.set(1,0,transformationMatrix[0][1]);
    invTM.set(1,1,transformationMatrix[1][1]);
    invTM.set(1,2,transformationMatrix[2][1]);
    invTM.set(2,0,transformationMatrix[0][2]);
    invTM.set(2,1,transformationMatrix[1][2]);
    invTM.set(2,2,transformationMatrix[2][2]);
    
    // invariant part
    invTM.set(3,0,0.0);
    invTM.set(3,1,0.0);
    invTM.set(3,2,0.0);
    invTM.set(3,3,1.0);
    
    // translation part -R^T*t
    Matrix invRot = new Matrix(3,3);
        // R^T part (R transposed)
    invRot.set(0,0,transformationMatrix[0][0]);
    invRot.set(0,1,transformationMatrix[1][0]);
    invRot.set(0,2,transformationMatrix[2][0]);
    invRot.set(1,0,transformationMatrix[0][1]);
    invRot.set(1,1,transformationMatrix[1][1]);
    invRot.set(1,2,transformationMatrix[2][1]);
    invRot.set(2,0,transformationMatrix[0][2]);
    invRot.set(2,1,transformationMatrix[1][2]);
    invRot.set(2,2,transformationMatrix[2][2]);
    
    Matrix transL = new Matrix(3,1);
    transL.set(0,0,transformationMatrix[0][3]);
    transL.set(1,0,transformationMatrix[1][3]);
    transL.set(2,0,transformationMatrix[2][3]);   
    
    Matrix modTransL;
    modTransL = invRot.times(transL);

    invTM.set(0,3,(-1*modTransL.get(0,0)));
    invTM.set(1,3,(-1*modTransL.get(1,0)));
    invTM.set(2,3,(-1*modTransL.get(2,0)));
    
            // Kontrolle: invTM*T = T*invTM = I;
    Matrix transformation = new Matrix(transformationMatrix);
    Matrix identity = invTM.times(transformation);
    
    
    if (printInfos == true){
        System.out.println("*************************************************");
        System.out.println("* Kontrolle von getInverseTransformationMatrix:");
        System.out.println("* ");
        System.out.println("* OriginalTransformationsMatrix: ");
        System.out.println("* ");                
        System.out.println(transformationMatrix[0][0]+"\t"+ transformationMatrix[0][1]+"\t"+ transformationMatrix[0][2]+"\t"+ transformationMatrix[0][3]);
        System.out.println(transformationMatrix[1][0]+"\t"+ transformationMatrix[1][1]+"\t"+ transformationMatrix[1][2]+"\t"+ transformationMatrix[1][3]);
        System.out.println(transformationMatrix[2][0]+"\t"+ transformationMatrix[2][1]+"\t"+ transformationMatrix[2][2]+"\t"+ transformationMatrix[2][3]);        
        System.out.println(transformationMatrix[3][0]+"\t"+ transformationMatrix[3][1]+"\t"+ transformationMatrix[3][2]+"\t"+ transformationMatrix[3][3]);
        System.out.println("* ");       
        System.out.println("* inverse Rotationsmatrix: ");
        System.out.println("* ");
        System.out.println(invRot.get(0,0)+"\t"+invRot.get(0,1)+"\t"+invRot.get(0,2));
        System.out.println(invRot.get(1,0)+"\t"+invRot.get(1,1)+"\t"+invRot.get(1,2));
        System.out.println(invRot.get(2,0)+"\t"+invRot.get(2,1)+"\t"+invRot.get(2,2));
        System.out.println("* ");
        System.out.println("* Translation");  
        System.out.println("* ");
        System.out.println(transL.get(0,0)+"\t"+transL.get(1,0)+"\t"+transL.get(2,0));        
        System.out.println("* ");
        System.out.println("* -R^T * Translation");  
        System.out.println("* ");
        System.out.println((-1*modTransL.get(0,0))+"\t"+(-1*modTransL.get(1,0))+"\t"+(-1*modTransL.get(2,0)));
        System.out.println("* ");  
        System.out.println("* inverse TransformationsMatrix: ");
        System.out.println("* ");                
        System.out.println(invTM.get(0,0)+"\t"+ invTM.get(0,1)+"\t"+ invTM.get(0,2)+"\t"+ invTM.get(0,3));
        System.out.println(invTM.get(1,0)+"\t"+ invTM.get(1,1)+"\t"+ invTM.get(1,2)+"\t"+ invTM.get(1,3));
        System.out.println(invTM.get(2,0)+"\t"+ invTM.get(2,1)+"\t"+ invTM.get(2,2)+"\t"+ invTM.get(2,3));        
        System.out.println(invTM.get(3,0)+"\t"+ invTM.get(3,1)+"\t"+ invTM.get(3,2)+"\t"+ invTM.get(3,3));
        System.out.println("* ");  
        System.out.println("* Identiy ??: ");
        System.out.println(identity.get(0,0)+"\t"+ identity.get(0,1)+"\t"+ identity.get(0,2)+"\t"+ identity.get(0,3));
        System.out.println(identity.get(1,0)+"\t"+ identity.get(1,1)+"\t"+ identity.get(1,2)+"\t"+ identity.get(1,3));
        System.out.println(identity.get(2,0)+"\t"+ identity.get(2,1)+"\t"+ identity.get(2,2)+"\t"+ identity.get(2,3));        
        System.out.println(identity.get(3,0)+"\t"+ identity.get(3,1)+"\t"+ identity.get(3,2)+"\t"+ identity.get(3,3));
        System.out.println("* ");  
        System.out.println("*************************************************");

    }
    
    
    
    return invTM;
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
        List<Vector3d> relative = new ArrayList<Vector3d>();
        double x, y, z;

        for (Vector3d p : mesh) {
            x = p.getX() - centroid[0];
            y = p.getY() - centroid[1];
            z = p.getZ() - centroid[2];

            relative.add(new Vector3d(x, y, z));
        }

        return relative;
    }

// added for future reference: in case we want to read wavefront.obj files directly.
// copied from ICPDemoApplet.java

protected Point3d[] readVerticesFromObjectFile(URL url, double sparse) {
	double sparse_accum = 0.0;
	try {
	    BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
	    Vector temp = new Vector();
	    Point3d ret[];
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

// KHK why use an inner class instead of a method?
// An instance of InnerClass can exist only within an instance of OuterClass 
// and has direct access to the methods and fields of its enclosing instance.
// from Wikipedia:
// Local inner classes are often used in Java to define callbacks for GUI code. 
// Components can then share an object that implements an event handling interface 
// or extends an abstract adapter class, containing the code to be executed when a given event is triggered.
// This avoids a large monolithic actionPerformed(ActionEvent) method with 
//multiple if-else branches to identify the source of the event. 
// This type of code is often considered messy and the inner class variations are considered to be better in all regards.

class WriteMatrix{
    
    void run(final double[][] transformationMatrix) {
        
    // basic idea from TransformJ and imagescience but modified for my purposes
       
                
        String filename = "";
        
        final String prefix = "";
        final String delim = "\t";
        final String postfix = "\n";
        
        StringBuilder sb = new StringBuilder();
        for (int r=0; r<4; ++r) {
            sb.append(prefix);
            for (int c=0; c<4; ++c) {
                // formatiert noch nicht richtig
                sb.append(d2s(transformationMatrix[r][c]));
                //sb.append(fmt.d2s(transformationMatrix[r][c]));
                if (c < 3) sb.append(delim);
            }
            sb.append(postfix);
        }        
        
        // Set System Look&Feel
        try {
        UIManager.setLookAndFeel(
            UIManager.getSystemLookAndFeelClassName());
        } 
        catch (UnsupportedLookAndFeelException e) {
           // handle exception
        }
        catch (ClassNotFoundException e) {
           // handle exception
        }
        catch (InstantiationException e) {
           // handle exception
        }
        catch (IllegalAccessException e) {
           // handle exception
        }
        JFileChooser saveDialog = new JFileChooser();
        
        FileNameExtensionFilter zipExtensionFilter = new FileNameExtensionFilter("ZIP File(*.zip)", "zip");
        saveDialog.addChoosableFileFilter(zipExtensionFilter);
        int saveDialogReturn = saveDialog.showSaveDialog(null);

        if (saveDialogReturn == JFileChooser.APPROVE_OPTION) {
            filename = saveDialog.getSelectedFile().getAbsolutePath();
        }
       
        try {
                final BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
                bw.write(sb.toString());
                bw.close();
	} catch (IOException e) {
		throw new IllegalArgumentException("Error writing to "+filename);
	}
    }
    
    // double to formated string (from imagescience/utility/Formatter.java)
     
    /** Returns a {@code String} representation of a {@code double} value.
        @param d the {@code double} value to be represented.
        @return a new {@code String} object containing a string representation of {@code d}. The maximum number of decimals used in representing {@code d} can be specified with method {@link #decs(int)}. The value of {@code d} is rounded to the specified maximum number of decimals. The returned string will contain less than the maximum number of decimals if {@code d} can be represented exactly that way. In particular, if {@code d} is equal to an integer value, the returned string represents that integer value, without decimals and preceding decimal separator symbol. The string returned when {@code Double.isNaN(d)} yields {@code true} can be specified with method {@link #nan(String)}. Similarly, the string returned when {@code Double.isInfinite(d)} yields {@code true} can be specified with method {@link #inf(String)}. The returned string is "0" if the absolute value of {@code d} is less than the limit set with method {@link #chop(double)}.
	*/
    String d2s(final double d) {            
		
                final DecimalFormat edf = new DecimalFormat("0.#E0",new DecimalFormatSymbols(Locale.US));
                final DecimalFormat ndf = new DecimalFormat("0.#",new DecimalFormatSymbols(Locale.US)); 
                ndf.setMaximumFractionDigits(8);
                
                double limit = 1.0E-12;
                double dflobo = 1.0E-10; // war original 0.1;
                
                String nan = "NaN";
                String inf = "Inf";
                
                if (Double.isNaN(d)) return nan;
		if (Double.isInfinite(d)) return (d < 0) ? ("-"+inf) : ("+"+inf);       // Abfrage pos oder neg infinity
                
		final long dr = Math.round(d);
		if (dr == d) return String.valueOf(dr);
		                
                // Conditional operator is also known as the ternary operator. 
                // This operator consists of three operands and is used to evaluate Boolean expressions. 
                // The goal of the operator is to decide which value should be assigned to the variable. 
                // The operator is written as:

                // variable x = (expression) ? value if true : value if false
                
		final double da = (d < 0) ? -d : d;
		if (da < limit) return "0";
		
		String ds;
                if (da < dflobo || da > 10000000) ds = edf.format(d);
                 else ds = ndf.format(d);

		return ds;
	}
}
        
// probably from package imagescience.transform; -> Transform.java
class ReadMatrix {
	
    	double[][] transformationMatrix = new double[4][4];
    
	double[][] run(final String file) {
		
		try {
			if (file == null || file.equals(""))
				throw new IllegalArgumentException("Empty matrix file name");

			transformationMatrix = load(file);
                        
		} catch (OutOfMemoryError e) {
			TJ.error("Not enough memory for this operation");
			
		} catch (UnknownError e) {
			TJ.error("Could not create output image for some reason.\nPossibly there is not enough free memory");
			
		} catch (IllegalArgumentException e) {
			TJ.error(e.getMessage());
			
		}
                return transformationMatrix;
	}
        
        //from TransformJ TJ_Matrix.java
        double[][] load(final String file) {
		
		// Read lines:
		final Vector lines = new Vector(10,10);     // 10 lines first, will be incremented for 10 lines, if full
		String line = null;
		try {
			final BufferedReader br = new BufferedReader(new FileReader(file));
			line = br.readLine();
			while (line != null) {
				line = line.trim();
				if (!line.equals(""))
					lines.add(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException("Unable to find "+file);
		} catch (IOException e) {
			throw new IllegalArgumentException("Error reading from "+file);
		}
		
		// Convert lines:
		if (lines.size() != 4)
			throw new IllegalArgumentException("File "+file+" does not contain a 4 x 4 matrix");
		String delim = "\t";
		line = (String)lines.get(0);
		if (line.contains(",")) delim = ",";
		else if (line.contains(" ")) delim = " ";
		final double[][] matrix = new double[4][4];
		for (int r=0; r<4; ++r) {
			line = (String)lines.get(r);
			final StringTokenizer st = new StringTokenizer(line,delim);
			if (st.countTokens() != 4)
				throw new IllegalArgumentException("File "+file+" does not contain a 4 x 4 matrix");
			for (int c=0; c<4; ++c) {
				try {
					matrix[r][c] = Double.parseDouble(st.nextToken());
				} catch (NumberFormatException e) {
					throw new IllegalArgumentException("Error reading element ("+r+","+c+") in "+file);
				}
			}
		}
		
		// Store matrix:
		return matrix;
	}
}
