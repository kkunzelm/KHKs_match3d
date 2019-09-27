import Jama.*;
import java.util.Arrays;
import javax.vecmath.*;

//todo Root Mean Square ... durchdenken.

// KDNode.java muss im gleichen Verzeichnis sein!

// hier wurde die rechenzeitintensive Berechnung der nächsten Nachbarn durch einen Suchbaum optimiert
// siehe z. B. file:///KHKsData/usr2/Recherchen/3d-matching-und-microCT-news/ICP-Erlangen-Java/Kd-tree.htm 
// Quelle: http://www9.informatik.uni-erlangen.de:81/sfb603/Saeulen/Optimierung/Allgemein/applet2 - Zugriff 2008 
//         (2014 nicht mehr zugängig)


public class ICPAlgorithm2014 {
  
    boolean debug = true;
    
    protected Matrix3d rotation;
    protected Vector3d translation;
   
    protected Point3d modelPoints[];
    protected Point3d dataPoints[];

    protected KDNode modelTree;
    
    protected Point3d modelCorner0;
    protected Point3d modelCorner1;

    protected Point3d workPoints[];
    protected int corresp[];
    protected int dist_order[];
    protected int model_unique[];


    public static final double ERROR_BOUND = 0.0001;
    public static final double ERROR_DIFF = 0.00001; // 0.1% change

    protected boolean firstRun;
    protected int minimum_valid_points;     // Mindestanzahl gültiger Pixel
    protected boolean refine;
    protected boolean refine_unique;        // klingt als ob jeder Punkt berücksichtigt wird
    
    protected boolean refine_clip;          // clipping ab bestimmter Position im Abstandshistogram, 
                                            // z. B. nur 75 % der niedrigeren Abstände berücksichtigt, 
                                            // Rest der Punktepaare verworfen
    protected double refine_clip_par;
    protected boolean refine_clamp;         // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
    protected double refine_clamp_par;
    protected boolean refine_sd;            // Punkte mit Abstand von sd_par mal stdDev + mean werden ausgewählt
    protected double refine_sd_par;
    protected boolean refine_sparse;        // nur jeder  xte Punkte wird verwendet
    protected double refine_sparse_par;


    protected int steps;  
    
    //      call this class like this: 

    //      ICPAlgorithm2012 instanceOfICP;
    //      instanceOfICP = new ICPAlgorithm2012(modelPoints, dataPoints, rotationMatrix, translationVector);
    //      instanceOfICP.runICP();
    
    // Constructor 
    public ICPAlgorithm2014() {
        firstRun = true;
    }
    
    
    // public void init(Point3d m[], Point3d d[], KDNode mT, KDNode dT, Point3d mC0, Point3d mC1) {
    public void init(Point3d m[], 
                     Point3d d[], 
                     KDNode mT, 
                     Point3d mC0, Point3d mC1, 
                     Matrix3d rotationAfterPrealignment, Vector3d translationAfterPrealignment,
                     ParameterICP refineParameters) {
        
        // assign the important values
        // 3D data
	modelPoints = m;
	dataPoints = d;

        modelTree = mT;
	modelCorner0 = mC0;
	modelCorner1 = mC1;
        
        // results from landmark based prealignment
        rotation = rotationAfterPrealignment;
        translation = translationAfterPrealignment;
        
        
        // diese Parameter könnte man in der Maske, in der die Bilder fürs Matching aufgerufen werden
        // setzen/editierbar machen.
        
        refine = true;
	refine_unique = refineParameters.getRefineUnique();             // nearest points are used
	refine_clip = refineParameters.getRefineClip();                 //true;
	refine_clip_par = refineParameters.getRefineClipPar();   //  = 0.75;             // the value is arbitrary here... 0.75 was used for one special example only
	refine_clamp = refineParameters.getRefineClamp();               // statistically determined outlier test
        refine_clamp_par = refineParameters.getRefineClampPar(); // multiples of standard deviations which are clipped/clamped
	refine_sd = refineParameters.getRefineSd();
        refine_sd_par = refineParameters.getRefineSdPar();
        refine_sparse = refineParameters.getRefineSparse();
	refine_sparse_par = refineParameters.getRefineSparsePar();            // data reduction, here 50 %
        minimum_valid_points = refineParameters.getMinValidPoints();
        
        
        // some preparations
	dist_order = new int[d.length];
	if (refine_unique) {
                  model_unique = new int[m.length];
             }
	else {
                  model_unique = null;
             }
         
    }
 
    
        // wichtigste methode... der eigentliche Kern des ICP Matching
    public double[][] runICP() {
	Point3d point;

	double x0, x1, y0, y1, z0, z1;

	boolean run_verbose = false;

	int dist_order_count;
	int dist_order_write;

	// start calculating

	if (run_verbose) System.out.println("Calculating...");

	Matrix3d rotMatrix;// = new Matrix3d();
	Vector3d transVector;

        // KHK for debugging:
        //rotMatrix.setIdentity();

        rotMatrix = rotation;
        transVector = translation;
        
	if (run_verbose) System.out.println("Rotation Matrix:" + rotMatrix);
	if (run_verbose) System.out.println(rotMatrix.toString());
	if (run_verbose) System.out.println("Translation Vector:"+transVector);
	
	workPoints = new Point3d[dataPoints.length];
	corresp = new int[workPoints.length];
	Point3d dataCenter = new Point3d();
	Point3d modelCenter = new Point3d();

	GMatrix rMatrix = new GMatrix(3, 3);            // GMatrix = general matrix 
	GMatrix addMatrix = new GMatrix(3, 3);

	// rotation matrix, as two-dim. array
	double r_array[][] = {{0.0, 0.0, 0.0},
			      {0.0, 0.0, 0.0},
			      {0.0, 0.0, 0.0}};

	double error;
	double last_error = Double.POSITIVE_INFINITY;
        double last_stdDev = Double.POSITIVE_INFINITY;
        
	/* as far as I understand the algorithm, the procedure is 
	 * to find the corresponding point-pairs via modelPoints and workPoints (workPoints is a copy of dataPoints)
	 * based on the workPoints to which the rotation and translation is applied the 
	 * next set of corresponding point-pairs is determined in the next iteration
	 * the rotation and translation, however, is always determined in total using the modelPoints and dataPoints.
	 */
  

	// ******************************************
	// - apply transform to data points, creating new set
	//*******************************************
	for (int i=0; i<dataPoints.length; i++) {
	    workPoints[i] = new Point3d(dataPoints[i]);
            // hier: erst Rotation, dann Translation...
	    rotMatrix.transform(workPoints[i]);
	    workPoints[i].add(transVector);
            // System.out.println("DataPoints und Workpoints:"+dataPoints[i]+" "+workPoints[i]);
	}

	if (run_verbose) System.out.println("Finished initial transform");

	boolean finished = false;

	//steps=0;

	// *******************************************************************************************************************************

	while(!finished) {

	    steps++;

	    // - find corresponding pairs into int list
	    for (int i=0; i<workPoints.length; i++) {
		point = workPoints[i];
                // KHK den Teil mit x0, x1 verstehe ich nicht. 
                // So wie ich das verstehe, müssen alle Werte gleich Null werden und werden sie auch. 
                // Soll das sicherstellen, dass es keine Werte außerhalb der Bounding-Box gibt?
                
		// x0..z1 are the distances OUTSIDE the bounding box
		x0 = modelCorner0.x - point.x; if (x0 < 0.0) { x0 = 0.0; }
		x1 = point.x - modelCorner1.x; if (x1 < 0.0) { x1 = 0.0; }
		y0 = modelCorner0.y - point.y; if (y0 < 0.0) { y0 = 0.0; }
		y1 = point.y - modelCorner1.y; if (y1 < 0.0) { y1 = 0.0; }
		z0 = modelCorner0.z - point.z; if (z0 < 0.0) { z0 = 0.0; }
		z1 = point.z - modelCorner1.z; if (z1 < 0.0) { z1 = 0.0; }

		corresp[i] = modelTree.findNearest(point,
						   Double.POSITIVE_INFINITY, 
						   (x0 * x0) + (x1 * x1),
						   (y0 * y0) + (y1 * y1),
						   (z0 * z0) + (z1 * z1));
                
	    }

	    if (run_verbose) {
                System.out.println("Computed nearest points");
            }

	    // dist_order lists the indices of the workPoints (dataPoints)
	    // in increasing distance to their corresponding points.
	    // dist_order_count describes how many "used" data/work points
	    // are actually in the list. All beyond this count have no
	    // corresponding model point.
	    // dist_order_write is a writing index
            
	    if ((dist_order == null) || (dist_order.length < workPoints.length)) {
		dist_order = new int[workPoints.length];
		System.out.println("Warning: Redoing dist_order");
	    }
		
	    // init ("unsorted");
	    dist_order_write = 0;
	    for (int i=0; i<dist_order.length; i++) {
		if (corresp[i] >= 0) {
		    dist_order[dist_order_write++] = i;
		}
	    }
	    dist_order_count = dist_order_write;
		
 
            // KHK todo: clip_gradient implementieren
	    if (refine) {

		// REFINE HERE
		// Here, the set of corresponding points can be further trimmed
		// down to improve convergence criteria, such as rejecting
		// pairs with too high a distance (multiple of average distance,
		// or a certain percentage of points
		
		if (refine_sparse) {
		    
		    // kill a certain percentage of points.
		    // This works WITHOUT sorting!
		    double sparse_accum = 0.0;
		    dist_order_write = 0;
		    for (int i=0; i<dist_order_count; i++) {
			sparse_accum += refine_sparse_par;
			if (sparse_accum < 1.0) {
			    // write back work point index into list
			    dist_order[dist_order_write++] = dist_order[i];
			} else {
			    sparse_accum -= 1.0;
			}
		    }
		    
		    dist_order_count = dist_order_write; // update
		
		}

		// bring in sorted order
		// KHK das war der Originalaufruf: sort(0, dist_order_count-1);
                sortKH(0, dist_order_count-1);      // Aufruf meiner Variante
                
		//sort(0, 5);
		
		// Code for printing all distances, used with no starting
		// transformation
		//for (int i=0; i<dist_order.length; i++) {
		//	System.out.println(sort_dist(i));
		//}
		//System.exit(0);
		
		if (refine_clip) {
		    // kill percentage of points with highest distance
		    // The two faces have about 75% overlap
		    // (approximated value from non-transformed distance values)
		    dist_order_count = (int)((double)dist_order_count * refine_clip_par);
		    // no need to flush, since we won't be touching them
		    // ever again!
		    //for (int i=(int)(dist_order_count); i<dist_order.length; i++) {
		    //    corresp[dist_order[i]] = -1; // no corresp point for you!
		    //}
		}
		
 		if (refine_clamp) {
		    
		    // kill all points farther away than a multiple (= refine_clamp_par) of
		    // the last reported mean distance and a multiple (= refine_clamp_par) of
		    // the current median distance

		    double med_distance = refine_clamp_par * sort_dist(dist_order_count/2);
		    double last_distance = refine_clamp_par * last_error;
		    dist_order_write = 0;
		    for (int i=0; i<dist_order_count; i++) {
			double di = sort_dist(i);
			if ((di <= last_distance) && (di <= med_distance)) {
			    // write back work point index into list
			    dist_order[dist_order_write++] = dist_order[i];
			}
		    }
		    System.out.println("refine_clamp");
		    System.out.println("last_stdDev: " + last_stdDev + "\n" 
                                     + "lastError: " + last_error + "\n" 
                                     + "workPoints.length: " + workPoints.length + "\n"
                                     + "dist_order_count: " + dist_order_count + "\n" 
                                     + "dist_order_write: " + dist_order_write);
		    dist_order_count = dist_order_write; // update
		
		}
                

                if (refine_sd) {        // sd = standardDeviation
		    // kill all points farther away than a multiple of the standard deviation from the current mean distance
                    double distance = refine_sd_par * last_stdDev;
		    dist_order_write = 0;
		    for (int i=0; i<dist_order_count; i++) {
			double di = sort_dist(i);
			if (di <= distance) {
			    // write back work point index into list
			    dist_order[dist_order_write++] = dist_order[i];
			}
		    }
                    System.out.println("refine_sd");
		    System.out.println("last_stdDev: " + last_stdDev + "\n"  
                                     + "lastError: " + last_error + "\n" 
                                     + "workPoints.length: " + workPoints.length + "\n"
                                     + "dist_order_count: " + dist_order_count+ "\n" 
                                     + "dist_order_write: " + dist_order_write);
                    System.out.println("");
		    dist_order_count = dist_order_write; // update
		
		}

		if (refine_unique) {
		    // allow only the closest work point per model point
		    // EXTREMELY SLOW!
		    for (int i=0; i<model_unique.length; i++) {
			model_unique[i] = 0;
		    }
		    
		    dist_order_write = 0;
		    for (int i=0; i<dist_order_count; i++) {
			int m = corresp[dist_order[i]]; // model point in question
			if (model_unique[m] == 0) {
			    model_unique[m] = 1;
			    dist_order[dist_order_write++] = dist_order[i];
			}
			//if (m != -1) {
			//    for (int j=i+1; j<dist_order_count; j++) {
			//	if (corresp[dist_order[j]] == m) {
			//	    corresp[dist_order[j]] = -1;
			//	}
			//    }
			//}
		    }
		    
		    dist_order_count = dist_order_write;
		}		

	    } // point refining on/off
            
            // KHK else is commented ... it is the same as refine clamp. This part could be deleted
	    else {          
// 		double last_distance = 2.5 * last_error;        // last_error = mean squared distance
// 		dist_order_write = 0;
// 		for (int i=0; i<dist_order_count; i++) {
// 		    double di = sort_dist(i); // = Double.POSITIVE_INFINITY;
// 		    if (di <= last_distance) {
// 			dist_order[dist_order_write++] = dist_order[i];
// 		    } 
// 		}
// 		dist_order_count = dist_order_write; // update
	    }

            // KHK: ab hier Berechung des Centroids 
            
            //*************************** KHK explanation

            /* this method implements formula (3) of Kanatani/Horn and others
             *  the weight is set to "1". It is a scaling factor only.
             * KHK exactly the same method is implemented for the corresponding landmark registration
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
             *  Point1'^T = (x'1, y'1, z'1) 
             * 
             *  The mark "'" means the corresponding point in the second image
             * 
             *  The multiplication of the two 3x1 vectors results 
             *  in a 3x3 matrix = correlationMatrixK.
             *  The matrices off all points are summarized and result in K
             * 
             */
             //*************************** KHK explanation
            
            
            
	    // - estimate rotation and translation
	    dataCenter.set(0.0, 0.0, 0.0);
	    modelCenter.set(0.0, 0.0, 0.0);

	    rMatrix.setZero();          // GMatrix
	    addMatrix.setZero();        // GMatrix

            // Berechnung des Centroids indem alle Punkte summiert werden und durch die Zahl der Punkte dividiert wird
            // jede einzelne Achse (x,y,z) wird separat addiert und dividiert (-> Division = .scale)
            
	    for (int i=0; i<dist_order_count; i++) {
		int j = dist_order[i]; // = i if unsorted
		//if (corresp[i] >= 0) {
		    dataCenter.add(dataPoints[j]);
		    modelCenter.add(modelPoints[corresp[j]]);
		    //c++;
		    //}
	    }

	    if (dist_order_count > 0) {
		dataCenter.scale(1.0 / ((double) dist_order_count));     // scale = hier Division durch Zahl der addierten Punkte
		modelCenter.scale(1.0 / ((double) dist_order_count));    // d.h. wir haben den Mittelwert aller Punkte = Centroid!
	    }
	    
	    if (run_verbose) System.out.println("Centers:");
	    if (run_verbose) System.out.println("  "+modelCenter);
	    if (run_verbose) System.out.println("  "+dataCenter);

	    
	    GVector dataCenterVector = new GVector(dataCenter); //  * Constructs a new GVector (javax.vecmath) and copies the initial values
                                                                //  * from the specified tuple - in this case a Point3d.  
                                                                //  * Point 3d = java.lang.Object +--javax.vecmath.Tuple3d  +--javax.vecmath.Point3d
	    GVector modelCenterVector = new GVector(modelCenter);
	    GVector dPoint = new GVector(3);
	    GVector mPoint = new GVector(3);
	    for (int i=0; i<dist_order_count; i++) {
		int j = dist_order[i]; // = i if unsorted
		//if (corresp[j] >= 0) {
		    dPoint.set(dataPoints[j]);                 
		    mPoint.set(modelPoints[corresp[j]]);        // ACHTUNG diese Zeile ist enorm wichtig!!!!
                                                                // hier wird bei jedem Durchgang eine bessere Zuordnung für die 
                                                                // Punktelisten verwendet. Dadurch wird die Rotation und Translation 
                                                                // auch immer bessere Ergebnisse erzielen.
                                                                // Es ist aber immer nur eine Rotation/Translation... es wird hier nichts akkumuliert!!

		    dPoint.sub(dataCenterVector);               // Punkte Liste - von den einzelnen Punkten wird der Centroid abgezogen
                                                                //  * Sets the value of this vector to the vector difference of itself
                                                                // * and vector (this = this - vector).
		    mPoint.sub(modelCenterVector);
		    addMatrix.mul(mPoint, dPoint);              // addMatrix ist als GMatrix 3x3 definiert
                                                                // addMatrix entspricht Kanatani's im Prinzip: vectorPoint1 * vectorPoint'1^T
                                                                // mPoint = row vector, dPoint = column vector
                                                                // Vektormultiplikation, mPoint, dPoint sind Vektoren
                                                                 /* Computes the outer product of the two vectors; multiplies the
                                                                 * the first vector by the transpose of the second vector and places
                                                                 * the matrix result into this matrix.  This matrix must be
                                                                 * be as big or bigger than getSize(v1)xgetSize(v2).
                                                                 */
                    
		    rMatrix.add(addMatrix);                      /**
                                                                 * Sets the value of this matrix to sum of itself and matrix m1.
                                                                 * @param m1 the other matrix (rMatrix wird bei ca. Z 386 mit Null initialisiert
                                                                 * 
                                                                 *  rMatrix entspricht: K = SUM(vectorPoint1 * vectorPoint'1^T) in der Kanatani Nomenklatur
                                                                 */  

	    }
	    

	    if (run_verbose) {
                System.out.println("Finished to calculate the correlation matrix K: Accumulated centers and rotation Matrix");
            }
            
            // if (dist_order_count >= 3) {     //... war Originalaufruf
            // Modifikation für minimum_valid_points
            System.out.println("Minimum valid points: "+ minimum_valid_points);
	    if (dist_order_count >= minimum_valid_points) {
		if (run_verbose) {
                    System.out.println("Minimizing...");
                    System.out.println("Kanatani's Correlation Matrix K: ");
                    System.out.println(rMatrix.toString());
                }

                //KHK: ab hier SVD singular value decomposition 
		// use Jama, not javax.vecmath (buggy)
		
		rMatrix.getRow(0, r_array[0]);              // getRow(int row, double[] array) Places the values of the specified row into the array parameter.
		rMatrix.getRow(1, r_array[1]);
		rMatrix.getRow(2, r_array[2]);
		
		Matrix rJamaMatrix = new Matrix(r_array, 3, 3);
		
		Matrix uJama;
		Matrix wJama;
		Matrix vJama;
		if (run_verbose) {
                    System.out.println("  computing SUV");
                }
		SingularValueDecomposition svdJama
		    = new SingularValueDecomposition(rJamaMatrix);
		uJama = svdJama.getU();
		wJama = svdJama.getS(); // diagonal matrix
		vJama = svdJama.getV();
		
		GMatrix U = new GMatrix(3,3);
		GMatrix W = new GMatrix(3,3);
		GMatrix V = new GMatrix(3,3);
		
		for (int i=0; i<3; i++) {
		    for (int j=0; j<3; j++) {
			U.setElement(i, j, uJama.get(i, j));
			W.setElement(i, j, wJama.get(i, j));
			V.setElement(i, j, vJama.get(i, j));
		    }
		}
		
                if (run_verbose){
                    System.out.println("Matrices U, W, V:");
                    System.out.println(U.toString());
                    System.out.println(W.toString());
                    System.out.println(V.toString());

                    System.out.println("Reconstructed rMatrix:");
                    
                    //W.mulTransposeRight(W, V);
                    //W.mul(W, V);
                    //rMatrix.mul(U, W);
                    //rMatrix.mulTransposeLeft(U, W);
                    System.out.println(rMatrix.toString());
                }
		
		if (run_verbose) {
                    System.out.println("  extracting rot and trans");
                }
		//rMatrix.mul(U, V);
		rMatrix.mulTransposeRight(U, V);  // mulTransposeRight(Matrix3d m1, Matrix3d m2) 
                                                  // Multiplies matrix m1 times the transpose of matrix m2, and places the result into this.
                                                  // also: U*V^T
		rMatrix.get(rotMatrix);           // rMatrix.get(Matrix3d m1) Places the values in the upper 3x3 of this GMatrix into the matrix m1.
                
		// rMatrix = GMatrix - die Korrelationsmatrix K
                // in rotMatrix ist jetzt die neue, aktuelle Rotationsmatrix für diesen einen Korrekturschritt enthalten
                
                /* Places the values in the upper 3x3 of this GMatrix into
                 * the matrix m1.
                 * @param m1  The matrix that will hold the new values
                 *  public final void get(Matrix3d m1)
                 */  
   
		// Test rotMatrix, since GMatrix has no determinant...
		if (rotMatrix.determinant() < 0.0) {
		    /*if (run_verbose) */System.out.println("Negative determinant!");
		    U.setElement(0, 2, -U.getElement(0, 2));
		    U.setElement(1, 2, -U.getElement(1, 2));
		    U.setElement(2, 2, -U.getElement(2, 2));
		    //rMatrix.mul(U, V);
		    rMatrix.mulTransposeRight(U, V);
		    rMatrix.get(rotMatrix);
		}
		
		//rotMatrix.transpose();    // das wäre dann die Umkehrfunktion der Rotationsmatrix R^T = R^-1

		
                
                // KH wichtig: wie bei meiner Vorwärtstransformation wird einer der beiden 
                // Centroidvektoren rotiert, bevor die Differenz als Translationsvektor berechnet wird.
		rotMatrix.transform(dataCenter);
		transVector.sub(modelCenter, dataCenter);           // werden die bisherigen Inhalte des Translationsvektor hier wirklich sauber überschrieben?
                                                                    // falls ja, dann wäre jetzt hier die Translation dieses einen Matching-Schrittes enthalten.
		
		// Apply motion (preliminary)
				
		for (int i=0; i<dataPoints.length; i++) {
		    // Apply Transformation
		    rotMatrix.transform(dataPoints[i], workPoints[i]); // rotMatrix ist 3x3 Matrix
                                                                         /**
                                                                         * Multiply this matrix by the tuple t and and place the result
                                                                         * into the tuple "result" (result = this*t).
                                                                         * @param t  the tuple to be multiplied by this matrix
                                                                         * @param result  the tuple into which the product is placed
                                                                         *  public final void transform(Tuple3d t, Tuple3d result)
                                                                         * 
                                                                         * Hier werden also die NEUEN WorkPoints (alte werden überschreiben) 
                                                                         * generiert, indem anhand der unveränderten dataPoints mit der neuen, besseren
                                                                         * Transformation die workPoints an neuer Stelle für die folgende Berechnung
                                                                         * der nearest neighbor Struktur verwendet werden!!
                                                                         * 
                                                                         */
   
		    workPoints[i].add(transVector);
		}
		if (run_verbose) {
                    System.out.println("Applied new and better transform");
                }
                /* for (int i=0; i<dataPoints.length;i++){
                    System.out.println("DataPoints/WorkPoints: "+dataPoints[i]+" "+workPoints[i]);
                }
                 */

                //************************************************************************
                // KHK: ab hier Abstandsmass berechnet , Fehler der minimiert werden soll
                
		//c = 0;
		error = 0.0;
		for (int i=0; i<dist_order_count; i++) {
		    int j = dist_order[i]; // = i if unsorted
		    //if (corresp[j] >= 0) {
                        // KHK error metric: error = mean of squared distance 
			error += workPoints[j].distanceSquared(modelPoints[corresp[j]]);     // KH: distanceSquared is a method of tuple3d from vecmath
			//c++;
			//}
		}
		if (dist_order_count>0) {
                    error /= (double) dist_order_count;
                }

                // variance
                double variance = 0;
                
                for (int i=0; i<dist_order_count; i++) {
		    int j = dist_order[i]; // = i if unsorted
		        variance += (error - (workPoints[j].distanceSquared(modelPoints[corresp[j]]))) *
                                (error - (workPoints[j].distanceSquared(modelPoints[corresp[j]]))) ;
                }
                variance = variance/dist_order_count;
                double stdDev = Math.sqrt(variance);
           
		// Now the current error is known
		
                if (run_verbose){
                    System.out.println("Result of the current iteration:");
                    System.out.println("Rotation Matrix of step:" + steps);
                    System.out.println(rotMatrix.toString());
                    System.out.println("Translation Vector:"+transVector);
                } 
		
                // KHK Interaktion mit Applet auskommentieren
                // evtl. als Abfrage:
                // if (!callFromMatching) ... then belassen, else ignorieren und unten weiterarbeiten

                 
		// if (run_verbose) System.out.println("Set new transform in scene");
		      

                // check abort conditions
				
		//error *= visualScale;
                

                System.out.println("");
		System.out.println("Error: "+error+" with "+dist_order_count+" points");
                System.out.println("stdDef: "+stdDev);
                System.out.println(""); 

		if (error < ERROR_BOUND) {
		    finished = true;
		}
		if (Math.abs(error-last_error) / error < ERROR_DIFF) {          // KHK: relative improvement less then ERROR_DIFF, at present: < 0.1% change
		    finished = true;
		}
		last_error = error;
                last_stdDev = stdDev;
		
	    } else { // if (c > 3)
		// not enough corresponding points
		finished = true;
	    }
	} // while(!finished)
        
        // Ende der Hauptschleife!!

	// finished, or aborted.
        
        System.out.println("");
        System.out.println("*******************************");
        System.out.println("*******************************");
        System.out.println();
	System.out.println("Finished after "+steps+" steps.");
        System.out.println();
        System.out.println("*******************************");
        System.out.println("*******************************");
        
        double[][] transformationMatrix = new double[4][4];
        
        transformationMatrix[0][0]=rotMatrix.getElement(0,0);
        transformationMatrix[1][0]=rotMatrix.getElement(1,0);
        transformationMatrix[2][0]=rotMatrix.getElement(2,0);
        transformationMatrix[3][0]=0.0;
        
        transformationMatrix[0][1]=rotMatrix.getElement(0,1);
        transformationMatrix[1][1]=rotMatrix.getElement(1,1);
        transformationMatrix[2][1]=rotMatrix.getElement(2,1);
        transformationMatrix[3][1]=0.0;
        
        transformationMatrix[0][2]=rotMatrix.getElement(0,2);
        transformationMatrix[1][2]=rotMatrix.getElement(1,2);
        transformationMatrix[2][2]=rotMatrix.getElement(2,2);
        transformationMatrix[3][2]=0.0;
        
        transformationMatrix[0][3]=transVector.x;
        transformationMatrix[1][3]=transVector.y;
        transformationMatrix[2][3]=transVector.z;
        transformationMatrix[3][3]=1.0;
               
        //if (run_verbose){
        if (true){
            System.out.println("Transformation Matrix formated for TransformJ: ");
            System.out.println(transformationMatrix[0][0]+"\t"+ transformationMatrix[0][1]+"\t"+ transformationMatrix[0][2]+"\t"+ transformationMatrix[0][3]);
            System.out.println(transformationMatrix[1][0]+"\t"+ transformationMatrix[1][1]+"\t"+ transformationMatrix[1][2]+"\t"+ transformationMatrix[1][3]);
            System.out.println(transformationMatrix[2][0]+"\t"+ transformationMatrix[2][1]+"\t"+ transformationMatrix[2][2]+"\t"+ transformationMatrix[2][3]);        
            System.out.println(transformationMatrix[3][0]+"\t"+ transformationMatrix[3][1]+"\t"+ transformationMatrix[3][2]+"\t"+ transformationMatrix[3][3]);
        }


        return transformationMatrix;
       
       

    } // run()
 // KHK: Ende der Run Methode
    


  
    private double sort_dist(int i) {
	int j = dist_order[i];
	if (corresp[j] < 0) {
            return Double.POSITIVE_INFINITY;
        }
	return workPoints[j].distanceSquared(modelPoints[corresp[j]]);
    }

    
    
    // Sort an int-array according to workPoint - modelPoint distance
    // this was the original quicksort algorithm
    
      // this is my interpretation of a quicksort algorithm from 
      // just as in KDNode.java I had to replace the original QS algorith 
      // with my interpretation of QS as I always got stack overflow errors
      // with the original version
    
      // from: http://algs4.cs.princeton.edu/23quicksort/Quick3way.java.html
      // quicksort the arraypoints] using 3-way partitioning
      
      private void sortKH(int left, int right) { 
        
        // System.out.println("     (sorting from "+left+" to "+right+")");
        if (right <= left) return;
        int lt = left, gt = right;
        int temp1 = dist_order[left];
        int i = left;
        int cmp;
        while (i <= gt) {
	    
		cmp = compareTo_KH(dist_order[i],temp1);
	    

            if      (cmp < 0) exch(dist_order, lt++, i++);
            else if (cmp > 0) exch(dist_order, i, gt--);
            else              i++;
        }

        //dist_order[left..lt-1] < temp = dist_order[lt..gt] < dist_order[gt+1..right]. 
        sortKH(left, lt-1);
        sortKH(gt+1, right);

	assert isSorted(dist_order, left, right);
	    
    }


      
   /***********************************************************************
    *  Helper sorting functions
    ***********************************************************************/
    // compareTo_KH(Point3d v, Point3d w)
    // replaces: 
    // ---> dist_order[i].compareTo(temp)
    // ---> v.compareTo(w) 
    
    private int compareTo_KH(int v, int w){
	int result;
	if (v < w){
		result = -1;
		return result;
	}
        else if (v == w) {
		result = 0;
		return result;		
        }
        else { 
		result = 1;
		return result;
	}
    }

        
    // exchange dist_order[i] and dist_order[j]
    private void exch(int[] dist_order, int i, int j) {
       int swap = dist_order[i];
       dist_order[i] = dist_order[j];
       dist_order[j] = swap;
    }


   /***********************************************************************
    *  Check if array is sorted - useful for debugging
    ***********************************************************************/

    private boolean isSorted(int[] dist_order, int left, int right) {
        for (int i = left + 1; i <= right; i++) {
            if (dist_order[i] < dist_order[i-1]) return false;
        }
        return true;
    }
    
    // KHK todo als eigenes Objekt mit Rückgabe von stdDev und mean ... anlegen a la ParameterICP
    // Methode funktioniert, wird im Moment nicht genutzt
    private double getStdDevSquaredEuklid(int dist_order_count, Point3d[] workPoints, Point3d[] modelPoints, int[] corresp, Matrix3d rotMatrix, Vector3d transVector){
     
                for (int i=0; i<dataPoints.length; i++) {
                    workPoints[i] = new Point3d(dataPoints[i]);
                    // hier: erst Rotation, dann Translation...
                    rotMatrix.transform(workPoints[i]);
                    workPoints[i].add(transVector);
                    // System.out.println("DataPoints und Workpoints:"+dataPoints[i]+" "+workPoints[i]);
                }
        
		double error = 0.0;
		for (int i=0; i<dist_order_count; i++) {
		    int j = dist_order[i]; // = i if unsorted
		    //if (corresp[j] >= 0) {
                        // KHK error metric: error = mean of squared distance 
			error += workPoints[j].distanceSquared(modelPoints[corresp[j]]);     // KH: distanceSquared is a method of tuple3d from vecmath
			//c++;
			//}
		}
		if (dist_order_count>0) {
                    error /= (double) dist_order_count;
                }

                // variance
                double variance = 0;
                
                for (int i=0; i<dist_order_count; i++) {
		    int j = dist_order[i]; // = i if unsorted
		        variance += (error - (workPoints[j].distanceSquared(modelPoints[corresp[j]]))) *
                                (error - (workPoints[j].distanceSquared(modelPoints[corresp[j]]))) ;
                }
                variance = variance/dist_order_count;
                
                double stdDev = Math.sqrt(variance);
                return stdDev;
    }
    
    // Methode funktioniert, wird im Moment nicht genutzt
    private double getMeanSquaredEuklid(int dist_order_count, Point3d[] workPoints, Point3d[] modelPoints, int[] corresp, Matrix3d rotMatrix, Vector3d transVector){
        
                for (int i=0; i<dataPoints.length; i++) {
                    workPoints[i] = new Point3d(dataPoints[i]);
                    // hier: erst Rotation, dann Translation...
                    rotMatrix.transform(workPoints[i]);
                    workPoints[i].add(transVector);
                    // System.out.println("DataPoints und Workpoints:"+dataPoints[i]+" "+workPoints[i]);
                }
        
        
		double error = 0.0;
		for (int i=0; i<dist_order_count; i++) {
		    int j = dist_order[i]; // = i if unsorted
		    //if (corresp[j] >= 0) {
                        // KHK error metric: error = mean of squared distance 
			error += workPoints[j].distanceSquared(modelPoints[corresp[j]]);     // KH: distanceSquared is a method of tuple3d from vecmath
			//c++;
			//}
		}
		if (dist_order_count>0) {
                    error /= (double) dist_order_count;
                }

                // System.out.println("inside getMeanSquaredEuklid: error =  " + error);
                return error;
    }
    
}
