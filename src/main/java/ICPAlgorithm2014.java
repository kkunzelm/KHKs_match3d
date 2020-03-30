import Jama.Matrix;
import Jama.SingularValueDecomposition;
import datastruct.KDNode;
import vecmath.*;

//todo Root Mean Square ... durchdenken.

// datastruct.KDNode.java muss im gleichen Verzeichnis sein!

// hier wurde die rechenzeitintensive Berechnung der nächsten Nachbarn durch einen Suchbaum optimiert
// siehe z. B. file:///KHKsData/usr2/Recherchen/3d-matching-und-microCT-news/ICP-Erlangen-Java/Kd-tree.htm 
// Quelle: http://www9.informatik.uni-erlangen.de:81/sfb603/Saeulen/Optimierung/Allgemein/applet2 - Zugriff 2008 
// (2014 nicht mehr zugängig)

class ICPAlgorithm2014 {

	private static final double ERROR_BOUND = 0.0001;
	private static final double ERROR_DIFF = 0.00001; // 0.1% change

	private Matrix3d landmarkRotation;
	private Vector3d landmarkTranslation;
	private Point3d[] targetModelPoints; // Target model
	private Point3d[] baseDataPoints; // Base 3D model: is matched to target
	private KDNode modelTree;
	private Point3d modelCorner0;
	private Point3d modelCorner1;
	private Point3d[] baseWorkPoints;
	private int[] correspondingPoints;
	private int[] distanceOrder;
	private int[] modelUnique;
	int validDistancesCtr; // used distances in distanceOrder
	private int minimumValidPoints; // Mindestanzahl gültiger Pixel
	private boolean refine;
	private boolean refineUnique; // klingt als ob jeder Punkt berücksichtigt wird
	private boolean refineAPrioriLandmark;

	// clipping ab bestimmter Position im Abstandshistogram, z.B. nur 75 % der niedrigeren Abstände berücksichtigt,
	// Rest der Punktepaare verworfen
	private boolean refineClip;
	private double refineClipParameter;

	// gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
	private boolean refineClamp;
	private double refineClampParameter;

	private boolean refineSd; // Punkte mit Abstand von sd_par mal stdDev + mean werden ausgewählt
	private double refineSdParameter;

	private boolean refineSparse; // nur jeder xte Punkte wird verwendet
	private double refineSparseParameter;

	private int iterationSteps;

	boolean runVerbose = false;

	Point3d landmarkCentroidTarget;
	Point3d landmarkCentroidBase;


	/** call this class like this:
	 * ICPAlgorithm2012 instanceOfICP;
	 * instanceOfICP = new ICPAlgorithm2012(modelPoints, dataPoints, rotationMatrix, translationVector);
	 * instanceOfICP.runICP();
	 *
	 * Constructor
 	**/
	public ICPAlgorithm2014() {
		boolean firstRun = true;
	}

	public void init(Point3d[] m, Point3d[] d, KDNode mT, Point3d mC0, Point3d mC1, Matrix3d rotationAfterPrealignment,
					 Vector3d translationAfterPrealignment, ParameterICP refineParameters,
					 double[] landmarkCentroidTarget, double[] landmarkCentroidBase) {

		// assign the important values
		// 3D data
		targetModelPoints = m;
		baseDataPoints = d;

		modelTree = mT;
		modelCorner0 = mC0;
		modelCorner1 = mC1;

		// results from landmark based prealignment
		landmarkRotation = rotationAfterPrealignment;
		landmarkTranslation = translationAfterPrealignment;

		this.landmarkCentroidTarget = new Point3d(landmarkCentroidTarget[0], landmarkCentroidTarget[1], landmarkCentroidTarget[2]);
		this.landmarkCentroidBase = new Point3d(landmarkCentroidBase[0], landmarkCentroidBase[1], landmarkCentroidBase[2]);

		// diese Parameter könnte man in der Maske, in der die Bilder fürs Matching aufgerufen werden
		// setzen/editierbar machen.
		refine = true;
		refineUnique = refineParameters.getRefineUnique(); // nearest points are used
		refineClip = refineParameters.getRefineClip(); // true;
		// the value is arbitrary here... 0.75 was used for one special example only
		refineClipParameter = refineParameters.getRefineClipPar(); // = 0.75;
		refineClamp = refineParameters.getRefineClamp(); // statistically determined outlier test
		refineClampParameter = refineParameters.getRefineClampPar(); // multiples of standard deviations which are
																	// clipped/clamped
		refineSd = refineParameters.getRefineSd();
		refineSdParameter = refineParameters.getRefineSdPar();
		refineSparse = refineParameters.getRefineSparse();
		refineSparseParameter = refineParameters.getRefineSparsePar(); // data reduction, here 50 %
		minimumValidPoints = refineParameters.getMinValidPoints();

		refineAPrioriLandmark = refineParameters.getRefineAPrioriLandmark();

		// some preparations
		distanceOrder = new int[d.length];
		if (refineUnique) {
			modelUnique = new int[m.length];
		} else {
			modelUnique = null;
		}

	}

	/*##################################################################################################################
	*
	* 										MAIN ALGORITHM
	*
	##################################################################################################################*/

	// wichtigste methode... der eigentliche Kern des ICP Matching
	public double[][] runICP() {
		//--------------------------------------------------------------------------------------------------------------
		//	Initialize variables
		//--------------------------------------------------------------------------------------------------------------

		// start calculating
		if (runVerbose)	System.out.println("Calculating...");

		// init rotation & translation matrices with results from landmark based prealignment
		Matrix3d rotationMatrix = landmarkRotation;
		Vector3d translationVector = landmarkTranslation;

		// KHK for debugging:
		// rotationMatrix.setIdentity();

		if (runVerbose)	{
			System.out.println("Rotation Matrix:" + rotationMatrix);
			System.out.println(rotationMatrix.toString());
			System.out.println("Translation Vector:" + translationVector);
		}

		Point3d baseDataCenter = new Point3d();
		Point3d targetModelCenter = new Point3d();

		GMatrix centroidRotationMatrix = new GMatrix(3, 3); // GMatrix = general matrix

		double error;
		double stdDev;
		double lastError = Double.POSITIVE_INFINITY;
		double lastStdDev = Double.POSITIVE_INFINITY;


		//--------------------------------------------------------------------------------------------------------------
		//	First transformation of Work Points
		//--------------------------------------------------------------------------------------------------------------
		/*
		 * as far as I understand the algorithm, the procedure is to find the corresponding point-pairs via modelPoints
		 * and workPoints (workPoints is a copy of dataPoints) based on the workPoints to which the rotation and
		 * translation is applied the next set of corresponding point-pairs is determined in the next iteration the
		 * rotation and translation, however, is always determined in total using the modelPoints and dataPoints.
		 */
		baseWorkPoints = transformBaseWorkPoints(rotationMatrix, translationVector);


		//--------------------------------------------------------------------------------------------------------------
		//	Iterations
		//--------------------------------------------------------------------------------------------------------------
		boolean iterationsFinished = false;

		while (!iterationsFinished) {

			iterationSteps++;

			//----------------------------------------------------------------------------------------------------------
			//	find corresponding points
			//----------------------------------------------------------------------------------------------------------
			correspondingPoints = new int[baseDataPoints.length];

			for (int i = 0; i < baseWorkPoints.length; i++) {
				correspondingPoints[i] = findNearestNeighborIndex(baseWorkPoints[i]);
			}

			if (runVerbose)
				System.out.println("Computed corresponding points");

			//----------------------------------------------------------------------------------------------------------
			//	refine corresponding point pairs
			//----------------------------------------------------------------------------------------------------------
			/*
			 * distanceOrder lists the indices of the workPoints (dataPoints) with increasing distance to their
			 * corresponding points.
			 */
			if ((distanceOrder == null) || (distanceOrder.length < baseWorkPoints.length)) {
				distanceOrder = new int[baseWorkPoints.length];
				System.out.println("Warning: Redoing distanceOrder");
			}

			/* Find valid point pairs:
			 * validDistancesCtr describes how many "used" data/work points are actually in the list.
			 * All beyond this count have no corresponding model point.
			 * distOrderWrite is a writing index
			 */
			int distOrderWrite = 0;
			for (int i = 0; i < distanceOrder.length; i++) {
				if (correspondingPoints[i] >= 0) {
					distanceOrder[distOrderWrite++] = i;
				}
			}
			validDistancesCtr = distOrderWrite;

			// KHK das war der Originalaufruf: sort(0, validDistancesCtr-1);
			sortDistanceOrderKH(0, validDistancesCtr - 1); // Aufruf meiner Variante

			/* Here, the set of corresponding points can be further trimmed down to improve convergence criteria, such
			* as rejecting pairs with too high a distance (multiple of average distance, or a certain percentage of
			* points
			*/
			if (refine) {
				// eliminate specific point pairs & adjust validDistancesCtr & distanceOrder
				refineCorrespondingPointPairs(lastError, lastStdDev);
			}


			//----------------------------------------------------------------------------------------------------------
			//	find centroids & calculate centroid rotation matrix
			//----------------------------------------------------------------------------------------------------------
			calcCentroids(baseDataCenter, targetModelCenter);

			centroidRotationMatrix = calcCentroidRotationMatrix(baseDataCenter, targetModelCenter);

			//----------------------------------------------------------------------------------------------------------
			//	Minimize distances
			//----------------------------------------------------------------------------------------------------------
			// if (validDistancesCtr >= 3) { //... war Originalaufruf: Modifikation für minimumValidPoints
			System.out.println("Minimum valid points: " + minimumValidPoints);
			if (validDistancesCtr >= minimumValidPoints) {
				if (runVerbose) {
					System.out.println("Minimizing...");
					System.out.println("Kanatani's Correlation Matrix K: ");
					System.out.println(centroidRotationMatrix.toString());
				}

				//------------------------------------------------------------------------------------------------------
				//	calculate SVD singular value decomposition
				//------------------------------------------------------------------------------------------------------
				calcSingularValueComposition(rotationMatrix, translationVector, baseDataCenter,
						targetModelCenter, centroidRotationMatrix);

				//------------------------------------------------------------------------------------------------------
				//	Apply motion (preliminary)
				//------------------------------------------------------------------------------------------------------
				applySVDTransformation(rotationMatrix, translationVector);

				//------------------------------------------------------------------------------------------------------
				//	Calculate Error
				//------------------------------------------------------------------------------------------------------
				// KHK: ab hier Abstandsmass berechnet: Fehler der minimiert werden soll
				error = getError();

				// get standard deviation
				stdDev = getStdDev(rotationMatrix, translationVector, error);

				System.out.println();
				System.out.println("Error: " + error + " with " + validDistancesCtr + " points");
				System.out.println("stdDef: " + stdDev);
				System.out.println();

				//------------------------------------------------------------------------------------------------------
				//	Evaluate abortion criteria
				//------------------------------------------------------------------------------------------------------
				if (error < ERROR_BOUND) {
					iterationsFinished = true;
				}
				if (Math.abs(error - lastError) / error < ERROR_DIFF) {
					// KHK: relative improvement less then ERROR_DIFF, at present: < 0.1% change
					iterationsFinished = true;
				}

				lastError = error;
				lastStdDev = stdDev;

			} else {
				//------------------------------------------------------------------------------------------------------
				//	Not enough corresponding points
				//------------------------------------------------------------------------------------------------------
				iterationsFinished = true;
			}
		} // while(!finished) // finished, or aborted.


		//------------------------------------------------------------------------------------------------------
		//	Apply landmark transformation again for a-priori refinements
		//------------------------------------------------------------------------------------------------------
		if (refineAPrioriLandmark) {
			//------------------------------------------------------------------------------------------------------
			//	calculate centroids from landmarks
			//------------------------------------------------------------------------------------------------------
			centroidRotationMatrix = calcCentroidRotationMatrix(landmarkCentroidBase, landmarkCentroidTarget);

			//------------------------------------------------------------------------------------------------------
			//	calculate SVD singular value decomposition
			//------------------------------------------------------------------------------------------------------
			calcSingularValueComposition(rotationMatrix, translationVector, landmarkCentroidBase,
					landmarkCentroidTarget, centroidRotationMatrix);

			//------------------------------------------------------------------------------------------------------
			//	Apply motion (preliminary)
			//------------------------------------------------------------------------------------------------------
			applySVDTransformation(rotationMatrix, translationVector);
		}


		//--------------------------------------------------------------------------------------------------------------
		//	Return Transformation Matrix
		//--------------------------------------------------------------------------------------------------------------

		System.out.println();
		System.out.println("*******************************");
		System.out.println("*******************************");
		System.out.println();
		System.out.println("Finished after " + iterationSteps + " steps.");
		System.out.println();
		System.out.println("*******************************");
		System.out.println("*******************************");

		return buildTransformationMatrix(rotationMatrix, translationVector);

	} // KHK: Ende der Run Methode


	/*##################################################################################################################
	*
	* 										MAIN CALCULATIONS
	*
	##################################################################################################################*/

	private Point3d[] transformBaseWorkPoints(Matrix3d rotationMatrix, Vector3d translationVector)
	{
		Point3d[] transformedWorkPoints = new Point3d[baseDataPoints.length];
		Point3d workPoint;

		// *************************************************
		// creating new work point & apply initial transform
		// *************************************************

		for (int i = 0; i < transformedWorkPoints.length; i++) {
			workPoint = new Point3d(baseDataPoints[i]);

			rotationMatrix.transform(workPoint);
			workPoint.add(translationVector);

			transformedWorkPoints[i] = workPoint;
		}

		if (runVerbose)
			System.out.println("Finished initial transform");

		return transformedWorkPoints;
	}

	private int findNearestNeighborIndex(Point3d baseWorkPoint) {

		// KHK den Teil mit x0, x1 verstehe ich nicht.
		// So wie ich das verstehe, müssen alle Werte gleich Null werden und werden sie auch.
		// Soll das sicherstellen, dass es keine Werte außerhalb der Bounding-Box gibt?

		// x0..z1 are the distances OUTSIDE the bounding box
		double x0 = modelCorner0.x - baseWorkPoint.x;
		if (x0 < 0.0)
			x0 = 0.0;

		double x1 = baseWorkPoint.x - modelCorner1.x;
		if (x1 < 0.0)
			x1 = 0.0;

		double y0 = modelCorner0.y - baseWorkPoint.y;
		if (y0 < 0.0)
			y0 = 0.0;

		double y1 = baseWorkPoint.y - modelCorner1.y;
		if (y1 < 0.0)
			y1 = 0.0;

		double z0 = modelCorner0.z - baseWorkPoint.z;
		if (z0 < 0.0)
			z0 = 0.0;

		double z1 = baseWorkPoint.z - modelCorner1.z;
		if (z1 < 0.0)
			z1 = 0.0;

		return modelTree.findNearest(
				baseWorkPoint,
				Double.POSITIVE_INFINITY,
				(x0 * x0) + (x1 * x1),
				(y0 * y0) + (y1 * y1),
				(z0 * z0) + (z1 * z1)
		);
	}

	private void refineCorrespondingPointPairs(double lastError, double lastStdDev) {
		// KHK todo: clip_gradient implementieren

		int distOrderWrite;// REFINE HERE

		if (refineSparse) {
			// kill a certain percentage of points. This works WITHOUT sorting!
			double sparseAccum = 0.0;
			distOrderWrite = 0;

			for (int i = 0; i < validDistancesCtr; i++) {
				sparseAccum += refineSparseParameter;
				if (sparseAccum < 1.0) {
					// write back work point index into list
					distanceOrder[distOrderWrite++] = distanceOrder[i];
				} else {
					sparseAccum -= 1.0;
				}
			}

			validDistancesCtr = distOrderWrite; // update
		}

		if (refineClip) {
			// kill percentage of points with highest distance
			// The two faces have about 75% overlap (approximated value from non-transformed distance values)
			validDistancesCtr = (int) ((double) validDistancesCtr * refineClipParameter);
			// no need to flush, since we won't be touching them ever again!
		}

		if (refineClamp) {
			// kill all points farther away than a multiple (= refine_clamp_par) of the last reported mean distance and
			// a multiple (= refine_clamp_par) of the current median distance

			double medianDistance = refineClampParameter * getDistanceToCorrespondingPoint(validDistancesCtr / 2);
			double lastDistance = refineClampParameter * lastError;
			distOrderWrite = 0;

			for (int i = 0; i < validDistancesCtr; i++) {
				double di = getDistanceToCorrespondingPoint(i);
				if ((di <= lastDistance) && (di <= medianDistance)) {
					// write back work point index into list
					distanceOrder[distOrderWrite++] = distanceOrder[i];
				}
			}

			System.out.println("refineClamp");
			System.out.println("lastStdDev: " + lastStdDev + "\n" + "lastError: " + lastError + "\n"
					+ "workPoints.length: " + baseWorkPoints.length + "\n" + "validDistancesCtr: "
					+ validDistancesCtr + "\n" + "distOrderWrite: " + distOrderWrite);

			validDistancesCtr = distOrderWrite; // update
		}

		if (refineSd) {
			// sd = standardDeviation
			// kill all points farther away than a multiple of the standard deviation from the current mean distance
			double distance = refineSdParameter * lastStdDev;
			distOrderWrite = 0;

			for (int i = 0; i < validDistancesCtr; i++) {
				double di = getDistanceToCorrespondingPoint(i);
				if (di <= distance) {
					// write back work point index into list
					distanceOrder[distOrderWrite++] = distanceOrder[i];
				}
			}

			System.out.println("refineSd");
			System.out.println("lastStdDev: " + lastStdDev + "\n" + "lastError: " + lastError + "\n"
					+ "baseWorkPoints.length: " + baseWorkPoints.length + "\n" + "validDistancesCtr: "
					+ validDistancesCtr + "\n" + "distOrderWrite: " + distOrderWrite);
			System.out.println();

			validDistancesCtr = distOrderWrite; // update
		}

		if (refineUnique) {
			// allow only the closest work point per model point. EXTREMELY SLOW!
			distOrderWrite = 0;

			for (int i = 0; i < validDistancesCtr; i++) {
				int m = correspondingPoints[distanceOrder[i]]; // model point in question
				if (modelUnique[m] == 0) {
					modelUnique[m] = 1;
					distanceOrder[distOrderWrite++] = distanceOrder[i];
				}
			}

			validDistancesCtr = distOrderWrite;
		}
	}

	private void calcCentroids(Point3d baseDataCenter, Point3d targetModelCenter) {
		// - estimate rotation and translation
		baseDataCenter.set(0.0, 0.0, 0.0);
		targetModelCenter.set(0.0, 0.0, 0.0);

		// Berechnung des Centroids indem alle Punkte summiert werden und durch die Zahl der Punkte dividiert wird
		// jede einzelne Achse (x,y,z) wird separat addiert und dividiert (-> Division = .scale)

		for (int i = 0; i < validDistancesCtr; i++) {
			int j = distanceOrder[i]; // = i if unsorted

			baseDataCenter.add(baseDataPoints[j]);
			targetModelCenter.add(targetModelPoints[correspondingPoints[j]]);
		}

		if (validDistancesCtr > 0) {
			// hier Division durch Zahl der addierten Punkte d.h. wir haben den Mittelwert aller Punkte = Centroid!
			baseDataCenter.scale(1.0 / ((double) validDistancesCtr));
			targetModelCenter.scale(1.0 / ((double) validDistancesCtr));
		}

		if (runVerbose)	{
			System.out.println("Centers:");
			System.out.println("  " + targetModelCenter);
			System.out.println("  " + baseDataCenter);
		}
	}

	private GMatrix calcCentroidRotationMatrix(Point3d baseDataCenter, Point3d targetModelCenter) {
		// *************************** KHK explanation
		/*
		 * this method implements formula (3) of Kanatani/Horn and others the weight is set to "1". It is a scaling
		 * factor only. KHK exactly the same method is implemented for the corresponding landmark registration
		 * Formula (3):
		 *
		 * K = SUM(vectorPoint1 * vectorPoint'1^T) --> ^T means transposed
		 *
		 * in detail:
		 * Each point consists of 3 coordinates x, y, z Each point is treated as a vector. The coordinates of these
		 * vectors are arranged vertically: Point1 = ( x1 ) ( y1 ) ( z1 )
		 *
		 * Point1'^T = (x'1, y'1, z'1)
		 * The mark "'" means the corresponding point in the second image
		 *
		 * The multiplication of the two 3x1 vectors results in a 3x3 matrix = correlationMatrixK.
		 * The matrices off all points are summarized and result in K
		 */
		// *************************** KHK explanation

		GMatrix rotMatrix = new GMatrix(3,3);
		GMatrix addMatrix = new GMatrix(3, 3);

		GVector dataCenterVector = new GVector(baseDataCenter);
		GVector modelCenterVector = new GVector(targetModelCenter);
		// * Constructs a new GVector (javax.vecmath) and copies the initial values from the specified tuple
		// - in this case a Point3d.

		GVector baseDataPoint = new GVector(3);
		GVector targetModelPoint = new GVector(3);

		for (int i = 0; i < validDistancesCtr; i++) {
			int j = distanceOrder[i]; // = i if unsorted
			// if (corresp[j] >= 0) {
			baseDataPoint.set(baseDataPoints[j]);
			targetModelPoint.set(targetModelPoints[correspondingPoints[j]]);
			// ACHTUNG diese Zeile ist enorm wichtig!!!!
			// hier wird bei jedem Durchgang eine bessere Zuordnung für die Punktelisten verwendet. Dadurch wird die
			// Rotation und Translation auch immer bessere Ergebnisse erzielen. Es ist aber immer nur eine Rotation/
			// Translation... es wird hier nichts akkumuliert!!

			// Sets the value of this vector to the vector difference of itself and vector (this = this - vector).
			baseDataPoint.sub(dataCenterVector);
			targetModelPoint.sub(modelCenterVector);
			// Punkte Liste - von den einzelnen Punkten wird der Centroid abgezogen

			/*
			 * Computes the outer product of the two vectors; multiplies the the first vector by the transpose of the
			 * second vector and places the matrix result into this matrix. This matrix must be be as big or bigger than
			 * getSize(v1) x getSize(v2).
			 */
			addMatrix.mul(targetModelPoint, baseDataPoint);
			// addMatrix ist als GMatrix 3x3 definiert
			// addMatrix entspricht Kanatani's im Prinzip: vectorPoint1 * vectorPoint'1^T
			// targetModelPoint = row vector, baseDataPoint = column vector
			// Vektormultiplikation, targetModelPoint, baseDataPoint sind Vektoren

			/*
			 * Sets the value of this matrix to sum of itself and matrix m1.
			 * the other matrix (rotMatrix wird bei ca. Z 386 mit Null initialisiert
			 * rotMatrix entspricht: K = SUM(vectorPoint1 * vectorPoint'1^T) in der Kanatani Nomenklatur
			 */
			rotMatrix.add(addMatrix);

		}

		if (runVerbose) {
			System.out.println(
					"Finished to calculate the correlation matrix K: Accumulated centers and rotation Matrix"
			);
		}

		return rotMatrix;
	}

	private void calcSingularValueComposition(Matrix3d rotationMatrix, Vector3d translationVector,
											  Point3d baseDataCenter, Point3d targetModelCenter,
											  GMatrix centroidRotationMatrix)
	{
		// use Jama, not javax.vecmath (buggy)
		if (runVerbose)
			System.out.println("  computing SUV");

		// rotation matrix, as two-dim. array
		double[][] rotationArray = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

		// getRow(int row, double[] array) Places the values of the specified row into the array parameter.
		centroidRotationMatrix.getRow(0, rotationArray[0]);
		centroidRotationMatrix.getRow(1, rotationArray[1]);
		centroidRotationMatrix.getRow(2, rotationArray[2]);

		Matrix rJamaMatrix = new Matrix(rotationArray, 3, 3);
		SingularValueDecomposition svdJama = new SingularValueDecomposition(rJamaMatrix);

		Matrix uJama = svdJama.getU();
		Matrix wJama = svdJama.getS(); // diagonal matrix
		Matrix vJama = svdJama.getV();

		GMatrix U = new GMatrix(3, 3);
		GMatrix W = new GMatrix(3, 3);
		GMatrix V = new GMatrix(3, 3);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				U.setElement(i, j, uJama.get(i, j));
				W.setElement(i, j, wJama.get(i, j));
				V.setElement(i, j, vJama.get(i, j));
			}
		}

		if (runVerbose) {
			System.out.println("Matrices U, W, V:");
			System.out.println(U.toString());
			System.out.println(W.toString());
			System.out.println(V.toString());

			System.out.println("Reconstructed centroidRotationMatrix:");
			System.out.println(centroidRotationMatrix.toString());
		}

		if (runVerbose) {
			System.out.println("  extracting rot and trans");
		}

		centroidRotationMatrix.mulTransposeRight(U, V);
		// mulTransposeRight(Matrix3d m1, Matrix3d m2)
		// Multiplies matrix m1 times the transpose of matrix m2, and places the result into this. also: U*V^T

		/*
		 * Places the values in the upper 3x3 of this GMatrix into the matrix m1.
		 * @param m1 The matrix that will hold the new values public final void get(Matrix3d m1)
		 */
		centroidRotationMatrix.get(rotationMatrix);
		// centroidRotationMatrix.get(Matrix3d m1) Places the values in the upper 3x3 of this GMatrix into the matrix m1
		// centroidRotationMatrix = GMatrix - die Korrelationsmatrix K in rotationMatrix ist jetzt die neue, aktuelle
		// Rotationsmatrix für diesen einen Korrekturschritt enthalten

		// Test rotationMatrix, since GMatrix has no determinant...
		if (rotationMatrix.determinant() < 0.0) {
			if (runVerbose)
				System.out.println("Negative determinant!");

			U.setElement(0, 2, -U.getElement(0, 2));
			U.setElement(1, 2, -U.getElement(1, 2));
			U.setElement(2, 2, -U.getElement(2, 2));

			centroidRotationMatrix.mulTransposeRight(U, V);
			centroidRotationMatrix.get(rotationMatrix);
		}

		// rotationMatrix.transpose(); // das wäre dann die Umkehrfunktion der Rotationsmatrix R^T = R^-1

		// KH wichtig: wie bei meiner Vorwärtstransformation wird einer der beiden Centroidvektoren rotiert, bevor die
		// Differenz als Translationsvektor berechnet wird.
		rotationMatrix.transform(baseDataCenter);
		translationVector.sub(targetModelCenter, baseDataCenter);
		// werden die bisherigen Inhalte des Translationsvektor hier wirklich sauber überschrieben?
		// falls ja, dann wäre jetzt hier die Translation dieses einen Matching-Schrittes enthalten.
	}

	private void applySVDTransformation(Matrix3d rotationMatrix, Vector3d translationVector) {
		for (int i = 0; i < baseDataPoints.length; i++) {
			// Apply Transformation
			rotationMatrix.transform(baseDataPoints[i], baseWorkPoints[i]);
			// rotationMatrix ist 3x3 Matrix

			/*
			 * Multiply this matrix by the tuple t and and place the result into the tuple: "result" (result = this * t)
			 *
			 * @param t
			 * 		the tuple to be multiplied by this matrix
			 * @param result
			 * 		the tuple into which the product is placed public final void transform(Tuple3d t, Tuple3d result)
			 *
			 * Hier werden also die NEUEN WorkPoints (alte werden überschreiben) generiert, indem anhand der
			 * unveränderten dataPoints mit der neuen, besseren Transformation die workPoints an neuer Stelle für die
			 * folgende Berechnung der nearest neighbor Struktur verwendet werden!!
			 */

			baseWorkPoints[i].add(translationVector);
		}
		if (runVerbose) {
			System.out.println("Applied new and better transform");
		}
	}

	private double getError() {
		double error;// c = 0;
		error = 0.0; // KHK error metric: error = mean of squared distance

		for (int i = 0; i < validDistancesCtr; i++) {
			int j = distanceOrder[i];
			error += baseWorkPoints[j].distanceSquared(targetModelPoints[correspondingPoints[j]]);
			// KH: distanceSquared is a method of tuple3d from vecmath
		}

		if (validDistancesCtr > 0) {
			error /= validDistancesCtr;
		}

		return error;
	}

	private double getStdDev(Matrix3d rotationMatrix, Vector3d translationVector, double error) {
		double variance = 0;

		for (int i = 0; i < validDistancesCtr; i++) {
			int j = distanceOrder[i];
			variance += (error - (baseWorkPoints[j].distanceSquared(targetModelPoints[correspondingPoints[j]])))
					* (error - (baseWorkPoints[j].distanceSquared(targetModelPoints[correspondingPoints[j]])));
		}

		variance = variance / validDistancesCtr;
		double stdDev = Math.sqrt(variance);

		// Now the current error is known
		if (runVerbose) {
			System.out.println("Result of the current iteration:");
			System.out.println("Rotation Matrix of step:" + iterationSteps);
			System.out.println(rotationMatrix.toString());
			System.out.println("Translation Vector:" + translationVector);
		}

		return stdDev;
	}

	private double[][] buildTransformationMatrix(Matrix3d rotationMatrix, Vector3d translationVector) {
		double[][] transformationMatrix = new double[4][4];

		transformationMatrix[0][0] = rotationMatrix.getElement(0, 0);
		transformationMatrix[1][0] = rotationMatrix.getElement(1, 0);
		transformationMatrix[2][0] = rotationMatrix.getElement(2, 0);
		transformationMatrix[3][0] = 0.0;

		transformationMatrix[0][1] = rotationMatrix.getElement(0, 1);
		transformationMatrix[1][1] = rotationMatrix.getElement(1, 1);
		transformationMatrix[2][1] = rotationMatrix.getElement(2, 1);
		transformationMatrix[3][1] = 0.0;

		transformationMatrix[0][2] = rotationMatrix.getElement(0, 2);
		transformationMatrix[1][2] = rotationMatrix.getElement(1, 2);
		transformationMatrix[2][2] = rotationMatrix.getElement(2, 2);
		transformationMatrix[3][2] = 0.0;

		transformationMatrix[0][3] = translationVector.x;
		transformationMatrix[1][3] = translationVector.y;
		transformationMatrix[2][3] = translationVector.z;
		transformationMatrix[3][3] = 1.0;

		// if (runVerbose){
		System.out.println("Transformation Matrix formated for TransformJ: ");
		System.out.println(transformationMatrix[0][0] + "\t" + transformationMatrix[0][1] + "\t"
				+ transformationMatrix[0][2] + "\t" + transformationMatrix[0][3]);
		System.out.println(transformationMatrix[1][0] + "\t" + transformationMatrix[1][1] + "\t"
				+ transformationMatrix[1][2] + "\t" + transformationMatrix[1][3]);
		System.out.println(transformationMatrix[2][0] + "\t" + transformationMatrix[2][1] + "\t"
				+ transformationMatrix[2][2] + "\t" + transformationMatrix[2][3]);
		System.out.println(transformationMatrix[3][0] + "\t" + transformationMatrix[3][1] + "\t"
				+ transformationMatrix[3][2] + "\t" + transformationMatrix[3][3]);
//		}

		return transformationMatrix;
	}

	/*##################################################################################################################
	*
	* 											Helper functions
	*
	##################################################################################################################*/

	private double getDistanceToCorrespondingPoint(int baseWorkPointIdx) {
		int j = distanceOrder[baseWorkPointIdx];

		if (correspondingPoints[j] < 0) {
			return Double.POSITIVE_INFINITY;
		}

		return baseWorkPoints[j].distanceSquared(targetModelPoints[correspondingPoints[j]]);
	}

	private void sortDistanceOrderKH(int left, int right) {
		// Sort an int-array according to workPoint - modelPoint distance
		// this was the original quicksort algorithm

		// this is my interpretation of a quicksort algorithm from just as in datastruct.KDNode.java
		// I had to replace the original QS algorith with my interpretation of QS as I always got stack overflow errors
		// with the original version

		// from: http://algs4.cs.princeton.edu/23quicksort/Quick3way.java.html
		// quicksort the arraypoints using 3-way partitioning

		if (right <= left)
			return;

		int lt = left, gt = right;
		int temp1 = distanceOrder[left];
		int i = left;
		int cmp;

		while (i <= gt) {

			cmp = compareToKH(distanceOrder[i], temp1);

			if (cmp < 0)
				exch(distanceOrder, lt++, i++);
			else if (cmp > 0)
				exch(distanceOrder, i, gt--);
			else
				i++;
		}

		// dist_order[left..lt-1] < temp = dist_order[lt..gt] < dist_order[gt+1..right].
		sortDistanceOrderKH(left, lt - 1);
		sortDistanceOrderKH(gt + 1, right);

		assert isSorted(distanceOrder, left, right);

	}

	private int compareToKH(int v, int w) {
		// compareTo_KH(Point3d v, Point3d w) replaces:
		// ---> dist_order[i].compareTo(temp)
		// ---> v.compareTo(w)
		int result;

		if (v < w) {
			result = -1;
			return result;
		} else if (v == w) {
			result = 0;
			return result;
		} else {
			result = 1;
			return result;
		}
	}

	private void exch(int[] distOrder, int i, int j) {
		// exchange distOrder[i] and distOrder[j]
		int swap = distOrder[i];

		distOrder[i] = distOrder[j];
		distOrder[j] = swap;
	}

	private boolean isSorted(int[] distOrder, int left, int right) {
		// Check if array is sorted - useful for debugging

		for (int i = left + 1; i <= right; i++) {
			if (distOrder[i] < distOrder[i - 1])
				return false;
		}

		return true;
	}

}
