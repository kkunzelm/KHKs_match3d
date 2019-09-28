import Jama.*;
import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.plugin.*;
import ij.process.*;
import java.awt.*;
import static java.lang.Math.*;
import java.util.*;
import java.util.List;
import vecmath.*;

// Point based rigid registration
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
// 
// History: 
// ========
// 2008 the first draft was functional
// October 2012 the plugin was working as intended        
// April 2014 Wolfram died ... I must finish this "in memoriam wg"
// In 2014 I started with the version from 2012 and decided to clean and rewrite the code 
// 
// This individual class here is just the implementation for the initial guess based on fiducial markers or manually added landmark points.
// The landmark points have to be entered as a polygon of type polyroi.


/**
 *
 * @author Prof. Dr. Karl-Heinz Kunzelmann, Ludwig Maximilians Universität, Poliklinik für Zahnerhaltung und Parodontologie, 80336 München
 *         email: karl-heinz@kunzelmann.de
 *         www:   www.kunzelmann.de
 * @version 0.1 
 * 
 */
public class Match3d_withFiducialMarkersOnly implements PlugIn {
    
    private static String title1 = "";
    private static String title2 = "";
    
    boolean debug = false;
    
    PolygonRoi polyRoiSource;
    PolygonRoi polyRoiTarget;

    int lengthSource;
    int lengthTarget;

     
    ImagePlus imp1;                     //Mnemo: Imp = imageplus
    ImagePlus imp2;

    ImageProcessor ip1;                 // Mnemo: Ip = image processor
    ImageProcessor ip2;
 
    FileInfo fi1;
    FileInfo fi2;
    FileInfo fiTransformedImage = new FileInfo();      
   
    // neu ab 2014
    List<Vector3f> mesh1 = new ArrayList<Vector3f>();
    List<Vector3f> mesh2 = new ArrayList<Vector3f>();
    
    List<Vector3f> landmarks1 = new ArrayList<Vector3f>();
    List<Vector3f> landmarks2 = new ArrayList<Vector3f>();
    
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
     * T is extended to a 1 x 4 matrix
     * (t. stands for translation element .)
     * 
     *      t1
     * T =  t2
     *      t3
     *       1 
     */
    
   
public void run(String arg) {
        
                // check: is there an open image?
        
		int[] wList = WindowManager.getIDList();
		if (wList==null) {
			IJ.noImage();
			return;
		}
		
                // check: we need at least 2 images for matching
                
		String[] titles = new String[wList.length];
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
		GenericDialog gd = new GenericDialog("Select Images for Matching", IJ.getInstance());
                
                gd.addMessage("Two images of file type '32bit' or 'float', which consist of height data, are required. \nx and y are coordinates on a rectangular grid, \nz represents the height information.");
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
		
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		int index1 = gd.getNextChoiceIndex();
		title1 = titles[index1];
		
		int index2 = gd.getNextChoiceIndex();
				
		title2 = titles[index2];
		imp1 = WindowManager.getImage(wList[index1]);
		imp2 = WindowManager.getImage(wList[index2]);
		
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
                
                lengthSource = polyRoiSource.getNCoordinates();
                lengthTarget = polyRoiTarget.getNCoordinates();
                
                if(lengthSource != lengthTarget) {
			IJ.error("Error: the number of points between the source and target polylines differs");
			return;
		} 
                
                // we need the scale infos for the x and y coordinates
                
                fi1 = imp1.getFileInfo();
                fi2 = imp2.getFileInfo();

                if(fi1.pixelWidth != fi2.pixelWidth) {
			IJ.error("Error: the pixel width between the source and target polylines differs");
			return;
		}      
                
                if(fi1.pixelHeight != fi2.pixelHeight) {
			IJ.error("Error: the pixel height between the source and target polylines differs");
			return;
		}   
                
                // todo: this could be fixed later, too
                
                if(fi1.pixelHeight != fi1.pixelWidth) {
			IJ.error("Error: the pixel height is different to pixel width (currently just square pixels are supported");
			return;
		}
                
                // ***************** lets start the real work ***************
                // 
                
               
                landmarks1 = getCorrespondingPointListFromPolygonRoi(imp1,(PolygonRoi)imp1.getRoi());
                landmarks2 = getCorrespondingPointListFromPolygonRoi(imp2,(PolygonRoi)imp2.getRoi());
                
                float[] centroidLandmarks1, centroidLandmarks2;

                centroidLandmarks1 = computeCentroid(landmarks1);
                centroidLandmarks2 = computeCentroid(landmarks2);
                
                List<Vector3f> relativeLandmarks1 = relativeCord(landmarks1, centroidLandmarks1);
                List<Vector3f> relativeLandmarks2 = relativeCord(landmarks2, centroidLandmarks2);
                
                if(debug == true){
                    System.out.println("2014 PolygonRoi Centroids Landmarks1: " + Arrays.toString(centroidLandmarks1));
                    System.out.println("2014 PolygonRoi Centroids Landmarks2: " + Arrays.toString(centroidLandmarks2));               
                }
                
                correlationMatrixK = calculateCorrelationMatrixK(relativeLandmarks1, relativeLandmarks2);
                
                Matrix rotation = svd(correlationMatrixK);                  // singular value decomposition
                Matrix invRotation = rotation.inverse();
    
                rotationMatrixR = rotation.getArrayCopy();
                invRotationMatrixInvR = invRotation.getArrayCopy();

                // for debugging
                if(true){
                    System.out.print("svd() RotationMatrix K = ");
                    rotation.print(9, 6);
                    System.out.print("svd ()inverseRotationMatrix K^-1 = ");
                    invRotation.print(9, 6);

                    // for debugging
                    Matrix rotTimeInv = rotation.times(invRotation);
                    System.out.println("Rotation times Inverted Rotation (should be Identity I): ");
                    rotTimeInv.print(9, 6);
                }
                
               float[] translation = calculateTranslation(centroidLandmarks1,centroidLandmarks2,rotationMatrixR);
               double[][] transformationMatrix = populateTransformationMatrix(rotationMatrixR,translation);
               System.out.println("TransformationMatrix in homogeneous coordinates: "+ Arrays.toString(transformationMatrix));
               
               // make the target image for the translation
               ImagePlus impTransformedImage = makeTarget(imp1);

              
               applyGeneralTransformation3D(imp2, impTransformedImage, transformationMatrix); 
               // display the image
               impTransformedImage.show();               
                
	}
    
/**
     *
     * @param imp
     * @param polyRoi
     * @return
     */
public List<Vector3f> getCorrespondingPointListFromPolygonRoi(ImagePlus imp, PolygonRoi polyRoi) {
        List<Vector3f> vectorList;
        vectorList = new ArrayList<Vector3f>();
        
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
        
        for(int i=0; i < roiLength; i++) {
            Vector3f vert = new Vector3f();
            vert.x = (float)fi.pixelWidth*(roiX[i]+ rectRoiSource.x);
            vert.y = (float)fi.pixelHeight*(roiY[i]+  rectRoiSource.y);
            vert.z = ip.getf((roiX[i]+rectRoiSource.x),(roiY[i]+rectRoiSource.y));   
            vectorList.add(vert); 
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
public double[][] calculateCorrelationMatrixK(List<Vector3f> relativeLandmarks1, List<Vector3f> relativeLandmarks2){
        
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
        
        double[][] sourcePoint = new double[3][3];
        double[][] targetPoint = new double[3][3];
        double[][] correlationMatrixK_temp = new double[3][3];

        
        for (int x = 0; x < relativeLandmarks1.size(); x++){
            
            sourcePoint[0][0] = relativeLandmarks1.get(x).getX();
            sourcePoint[1][0] = relativeLandmarks1.get(x).getY();
            sourcePoint[2][0] = relativeLandmarks1.get(x).getZ();
            targetPoint[0][0] = relativeLandmarks2.get(x).getX();
            targetPoint[0][1] = relativeLandmarks2.get(x).getY();
            targetPoint[0][2] = relativeLandmarks2.get(x).getZ();
           
            System.out.println("Accumulation of the correlationMatrixK - step no: x = "+x);
        
            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    for (int k = 0; k < 3; k++){
                        // System.out.println("ijk: "+i+"/"+j+"/"+k+" -> Matrixelement ik before "+ correlationMatrixK[i][k]);
                        correlationMatrixK_temp[i][k] += sourcePoint[i][j] * targetPoint[j][k];
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
     * @param rotationMatrixR
     * @return translation Translation vector 
     */
private float[] calculateTranslation(float[] centroid1, float[] centroid2, double[][] rotationMatrixR){
        
        float[] translation = new float[3];
        float[] rotCentroid2 = new float[3];
        
       
        // Calculate the translation.
        // Important: to calculate the translation, the centroid of the target image, which will be transformed, has to be rotated!!!
        // T = Sc - Tc*R
        // Sc = Source image center, Tc = Target image center

       
       rotCentroid2[0] = (float) (rotationMatrixR[0][0]*centroid2[0] + rotationMatrixR[0][1]*centroid2[1] + rotationMatrixR[0][2]*centroid2[2]);
       rotCentroid2[1] = (float) (rotationMatrixR[1][0]*centroid2[0] + rotationMatrixR[1][1]*centroid2[1] + rotationMatrixR[1][2]*centroid2[2]);
       rotCentroid2[2] = (float) (rotationMatrixR[2][0]*centroid2[0] + rotationMatrixR[2][1]*centroid2[1] + rotationMatrixR[2][2]*centroid2[2]);          

       translation[0] = centroid1[0] - rotCentroid2[0];
       translation[1] = centroid1[1] - rotCentroid2[1];
       translation[2] = centroid1[2] - rotCentroid2[2];

       System.out.println("Translation: "+ translation[0] + ", " + translation[1] + ", " + translation[2]);

       return translation;
        
    }
    
private double[][] populateTransformationMatrix(double[][] rotationMatrixR, float[] translation){
        
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
        impTransformedImage.setTitle("Transformed");
        
        return impTransformedImage;

}
        
public void applyGeneralTransformation3D(ImagePlus imp2, ImagePlus impTransformedImage ,double[][] transformationMatrix) {
 
        ImageProcessor ipTransformedImage = impTransformedImage.getProcessor();

        FileInfo fi = imp2.getFileInfo();
        this.ip2 = imp2.getProcessor();

        int w = this.ip2.getWidth();
        int h = this.ip2.getHeight();
        
        float tempz;
    
        // Nomenklatur:
        // ich habe zwei Bilder imp1 (meist: Baseline) und imp2 (meist: Follow-up)
        // beide unterscheiden sich um R, T
        // damit ich die Differenz berechnen kann, muss ich ein neues Bild anlegen (target)
        // dieses Bild (ipTransformedImage) is das zurückrotierte Bild imp2
        // das zurückrotierte Bild nennt man auch target
        // imp2 dient als source für target
        
        // target-to-source mapping
        // my target to source mapping uses:
        // 
        // x-y_grid_point_in_source = R^-1*T^-1*u-v_grid_point_in_target (R, T are 4x4 matrizes, T^-1 is done first)
        // then I collect the z value in the source image the point(x,y,z)
        // is projected into the target image with
        // u-v_point(u,v,z)_in_target = (R^-1*T^-1)^-1*x-y_point(x,y,z)_in_source
        //                            = T(4x4)*R(4x4)*x-y_point(x,y,z)_in_source
        // 
        // Kanatani uses first translation then rotation
        
         
        for (int v=0; v<h; v++){
            for (int u=0; u<w; u++){
                
                /* in the rot_target image we map the coordinates u,v to the 
                 * coordinates x,y in the rot_source image using the inverseRotationMatrix
                 * there we determine the z-value of P by interpolation
                 * afterwards we map P back to the rot_target image using the rotationMatrix
                 * 
                 * In addition we have to consider the translation, too.
                 * 
                 * The right sequence/order caused me a lot of trial/error/thinking but finally I got it right
                 */ 

              /*   KHK: Stand 11.39 h
                
                         statt pt.x double besser float pt[] = new float[3] verwenden
                                 bedenken: u und v sind integer Werte, stellen die Gitterpunkte dar,
                                 müssen in fileInfo Einheiten transformiert werden.Einheiten
                */
                
                // Punkt als Vektor in homogenen Koordinaten angegeben.
                
                double scale = fi.pixelWidth; // pixelHeigth and pixelWidth must be the same
                //set(int i,int j,double s) - Achtung: Zeile dann Spalte
                Matrix point = new Matrix(4,1);
                Matrix pointTransformed;      // = new Matrix(4,1);
                point.set(0,0,u * scale);
                point.set(1,0,v * scale);
                point.set(2,0,0.0);
                point.set(3,0,1.0);
                
                Matrix transformation = new Matrix(transformationMatrix);
                Matrix invTransformation = transformation.inverse();
                
                // times(Matrix B) Linear algebraic matrix multiplication, A * B
                
                pointTransformed = invTransformation.times(point);
                
                // get matrix single element: get(int i, int j)    
                // ACHTUNG ganz genau über Position in int und Koordinaten in unit nachdenken.
                
                int xPositionInPixels = (int) rint(pointTransformed.get(0,0)/scale);
                int yPositionInPixels = (int) rint(pointTransformed.get(1,0)/scale);
                
                if    (xPositionInPixels > 0 && xPositionInPixels < w 
                    && yPositionInPixels > 0 && yPositionInPixels < h){
                        pointTransformed.set(2,0,(double)this.ip2.getf(xPositionInPixels,yPositionInPixels));   // nearest neighbour interpol. todo ... can be improved later
                }
               

                if (pointTransformed.get(2,0) != 0.0){                      // z-value = 0.0 can be used as a mask. 0.0 will not be mapped
                                                        // z-Werte mit Null werden nicht gemappt
                           
                    point=transformation.times(pointTransformed);
                    
                tempz = (float) point.get(2,0);
                    
                }
                else {
                tempz = 0;
                }
                ipTransformedImage.setf(u,v,tempz);
           }
        }
        
        
        // update the rotated image
        ipTransformedImage.resetMinAndMax();
        double min = ipTransformedImage.getMin();
        double max = ipTransformedImage.getMax();
        ipTransformedImage.setMinAndMax(min, max);
        impTransformedImage.updateAndDraw();
        System.out.println ("Min: " + min + " Max:" + max);
        
  
  }
   
private float[] computeCentroid(List<Vector3f> mesh) {
        float x = 0, y = 0, z = 0;

        for (Vector3f p : mesh) {
            x = x + p.getX();
            y = y + p.getY();
            z = z + p.getZ();
        }

        x = x / mesh.size();
        y = y / mesh.size();
        z = z / mesh.size();

        return new float[]{x, y, z};
    }

private List<Vector3f> relativeCord(List<Vector3f> mesh, float[] centroid) {
        List<Vector3f> relative = new ArrayList<Vector3f>();
        float x, y, z;

        for (Vector3f p : mesh) {
            x = p.getX() - centroid[0];
            y = p.getY() - centroid[1];
            z = p.getZ() - centroid[2];

            relative.add(new Vector3f(x, y, z));
        }

        return relative;
    }
   
}
 

