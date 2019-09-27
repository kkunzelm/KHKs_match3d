/* package cz.fidentis.comparison.icp;

This class is part of Fidentis

FIDENTIS – Forensic Facial Identification Software
The project is based on multidisciplinary cooperation. Development is provided by two departments of Masaryk University: 

Human Computer Interaction Laboratory at Faculty of Informatics and Department of Anthropology at Faculty of Science.

Author: probably Katarína Furmanová
*/

// KHK singleton design pattern
// Aufruf z. B.:
/* import cz.fidentis.comparison.icp.Icp;
   Icp.setMaxIteration((Integer) xyz );
   Icp.instance().icp(modelData.getVerts(), compareData.getVerts(), (Float) xyz );
*/

    //*****************************************
    // Aufruf der Fidentis Implementierung.
    // wäre wegen der Quaternionen interessant,
    // hat aber ein Speicherproblem bei der nearestNeighbor-Bestimmung
    
    /*
    List<Vector3d> vectorListImg1;
    vectorListImg1 = getCorrespondingPointListFromImageData(imp1);

    List<Vector3d> vectorListImg2;
    vectorListImg2 = getCorrespondingPointListFromImageData(impTransformedImage);
    
    Icp.setMaxIteration((Integer) 1 );
    Icp.instance().icp(vectorListImg1, vectorListImg2, (Float) 0.05f );
    */


import Jama.Matrix;
import javax.vecmath.Vector3f;

import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA. User: ex3me Date: 2/22/13 Time: 6:56 PM To change
 * this template use File | Settings | File Templates.
 */
public class Icp {

    private static int maxIteration = 2;
    private static int currentIteration = 1;
    private static Icp unique;
   // private static Float prevMean;
    private static final int MAX_ERROR = 100;

    // KHK Singleton Pattern... takes care that only one instance of the class is active.
    // is this here the "lazzy instantiation" pattern ?

    
    public static Icp instance() {
        if (unique == null) {
            unique = new Icp();
        }
            
        System.out.println("KHK: arrived in ICP");
        return unique;
    }

    private Icp() {
    }

    public static void setMaxIteration(int maxIteration) {
        Icp.maxIteration = maxIteration;
            
        System.out.println("KHK: MaxIteration is set");
    }

    /**
     * Method to search for the nearest neighbour of point in mesh, within error
     * range. Uses method pointsInRange() to get all the points within the error
     * range. 
     *
     * Determines distance for all the points within the range and finds the one
     * with the lowest value, returning it as the nearest neighbour.
     *
     * @param point searching for nearest neighbour of the point
     * @param mesh searching in the mesh
     * @param error searching withing error range
     * @return point representing nearest neighbour
     */
    public Vector3f nearestNeighbour(Vector3f point, KdTree mesh, float error) {
        System.out.println("NearestNeighbour started");
        List<Vector3f> inRange;
        inRange = mesh.pointsInRange(point, error);

        double minDist;
        double curDist;
        Vector3f nearestPoint;

        nearestPoint = inRange.get(0);
        //computes distance from point to point
        minDist = distancePoints(point, nearestPoint);

        for (Vector3f p : inRange) {
            curDist = distancePoints(point, p);  // cur = current ?
            if (curDist < minDist) {
                minDist = curDist;
                nearestPoint = p;
            }
        }

        return nearestPoint;
    }

    /**
     * Computes distance between point1 and point2 in space.
     *
     * @param point1
     * @param point2
     * @return distance between point1 and point2
     */
    public double distancePoints(Vector3f point1, Vector3f point2) {
        return Math.sqrt(Math.pow(point1.getX() - point2.getX(), 2)
                + Math.pow(point1.getY() - point2.getY(), 2)
                + Math.pow(point1.getZ() - point2.getZ(), 2));
    }

    private float distanceCoordinates(float p1, float p2) {
        return p2 - p1;
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

    private Matrix sumMatrixComp(Vector3f p) {
        return new Matrix(new double[][]{{0, -p.getX(), -p.getY(), -p.getZ()},
            {p.getX(), 0, p.getZ(), -p.getY()},
            {p.getY(), -p.getZ(), 0, p.getX()},
            {p.getZ(), p.getY(), -p.getX(), 0}});
    }

    private Matrix sumMatrixMain(Vector3f p) {
        return new Matrix(new double[][]{{0, -p.getX(), -p.getY(), -p.getZ()},
            {p.getX(), 0, -p.getZ(), p.getY()},
            {p.getY(), p.getZ(), 0, -p.getX()},
            {p.getZ(), -p.getY(), p.getX(), 0}});
    }

    private Quaternion conjugateQ(Quaternion q) {
        return new Quaternion(-q.getX(), -q.getY(), -q.getZ(), q.getW());
    }

    private Quaternion multiply(Quaternion q1, Quaternion q2) {

        float x, y, z, w; // KHK w ist der reale Anteil, x,y,z sind die Komponenten von i,j,k

        w = q1.getW() * q2.getW() - q1.getX() * q2.getX() - q1.getY() * q2.getY() - q1.getZ() * q2.getZ();
        x = q1.getW() * q2.getX() + q2.getW() * q1.getX() + q1.getY() * q2.getZ() - q1.getZ() * q2.getY();
        y = q1.getW() * q2.getY() + q2.getW() * q1.getY() - q1.getX() * q2.getZ() + q1.getZ() * q2.getX();
        z = q1.getW() * q2.getZ() + q2.getW() * q1.getZ() + q1.getX() * q2.getY() - q1.getY() * q2.getX();
        return new Quaternion(x, y, z, w);
    }

    /**
     * ICP with rotations Transforms points in compF as close to mainF as
     * possible. Uses k-D Tree to search for nearest neighbours, and to find the
     * pairing between both meshes. Uses quaternions for rotations.
     *
     * Modifies the points in parameter Lists.
     *
     * @param mainF - Main Face
     * @param compF - Compared Face which will be aligned to Main Face
     * @param error - When difference of two meshes drops below the error rate,
     * the computation stops.
     */
    public void icp(List<Vector3f> mainF, List<Vector3f> compF, float error) {
        currentIteration = 1;
        // KHK warum ist prevMean ein Objekt?
        Float prevMean = null;
        
        //XYZ2DEM_ImporterHack hack0 = new XYZ2DEM_ImporterHack();
        //hack0.khkDisplayXYZ(compF); 
        
        System.out.println("Before KdTree");
        KdTree tree1 = new KdTree();
        tree1.putAll(mainF);
        System.out.println("After KdTree");
        float meanX , meanY , meanZ , meanD, x, y, z;
        // near und comp2 sind später Kopien von mainF und compF aber mit nearest point Zuordnung
        List<Vector3f> near = new ArrayList<Vector3f>();
        List<Vector3f> comp2 = new ArrayList<Vector3f>();
        
        int step = 1 + compF.size()/5000;
        //int step = 1;
        System.out.println("compF.size() und Step size: " + compF.size() + ", "+ step);
        System.out.println("compF.get(0)" + compF.get(0));
        //System.out.println("compF" + compF);
        do {
            meanX = 0;
            meanY = 0;
            meanZ = 0;
            meanD = 0;
            
            System.out.println("Inside do-loop");
            near.clear();
            comp2.clear();
            
            if(currentIteration > 1) {
                // nach der ersten Iteration ist meanD in der Schleife mit einem Wert gefüllt. Vorher nicht (= null).
                prevMean = meanD;
                System.out.println("First iteration");
            }
            
             for(int i = 0; i < compF.size(); i = i + step ) {
                System.out.println("Compute distance between points and list compF  " + compF.get(i));
                Vector3f nn = nearestNeighbour(compF.get(i), tree1, 3f);


                // KHK könnte man auch hier berechnen, ist nur Abstand zwischen den Punkten koordinatenweise
                x = distanceCoordinates(compF.get(i).getX(), nn.getX());
                y = distanceCoordinates(compF.get(i).getY(), nn.getY());
                z = distanceCoordinates(compF.get(i).getZ(), nn.getZ());

                meanX = meanX + x;
                meanY = meanY + y;
                meanZ = meanZ + z;

                // KHK distancePoints = euklidscher Abstand - error metric
                meanD = (float) (meanD + distancePoints(compF.get(i), nn));
                near.add(nn);
                comp2.add(compF.get(i));
            }

            // KHK mittlere Abstand zwischen den beiden meshes koordinatenweise
            meanX = meanX / compF.size();
            meanY = meanY / compF.size();
            meanZ = meanZ / compF.size();

            float[] mainCenter, compCenter;

            mainCenter = computeCentroid(near);
            compCenter = computeCentroid(comp2);
             
            System.out.println("Centroids done");
            // KHK centroid wird in den Ursprung verschoben - und zwar die Kopien der meshes
            List<Vector3f> relativeM = relativeCord(near, mainCenter);
            List<Vector3f> relativeC = relativeCord(comp2, compCenter);

            Matrix m = new Matrix(4, 4);
            Matrix tmp;

            for (int i = 0; i < relativeC.size(); i++) {
                // Beginn KHK:
                // sumMatrixComp und sumMatrixMain sind Klassenmethoden
                // siehe Seite 10 hier: http://www.cs.iastate.edu/~cs577/handouts/quaternion.pdf 
                // (oder Download in ./Recherchen/Matching/quaternion-anleitung_mega-genial.pdf)
                //
                // create a symmetric positive definite matrix
                // Matrix A = Matrix.random(N, N);
                // A = A.transpose().times(A);
                // Ende KHK
                
                tmp = sumMatrixComp(relativeC.get(i)).transpose().times(sumMatrixMain(relativeM.get(i)));
                // m is symmetric
                m = m.plus(tmp);
            }
            // Beginn KHK:
            // If A (or here m) is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal 
            // and the eigenvector matrix V is orthogonal. I.e. A = V.times(D.times(V.transpose())) and V.times(V.transpose()) 
            // equals the identity matrix.
            // von: http://introcs.cs.princeton.edu/java/95linear/Eigenvalues.java.html 
            // nur als Beispiel wegen der Analogie zu unten
            // compute the spectral decomposition
            // EigenvalueDecomposition e = A.eig();
            // Matrix V = e.getV();
            // Matrix D = e.getD();
            
            // EigenvalueDecomposition:            
            Matrix eigD = m.eig().getD();  // Return the block diagonal eigenvalue matrix
            Matrix eigM = m.eig().getV();  // Return the eigenvector matrix
            
            System.out.println("Decomposition done");

            // The unit quaternion q that maximizes (14) is the eigenvector that corresponds to the largest eigenvalue
            // of the matrix M. It describes the optimal rotation for (12), i.e, for data registration.
            int max = 0;
            for (int i = 0; i < 4; i++) {
                if (eigD.get(max, max) <= eigD.get(i, i)) {
                    max = i;
                }
            }

            Quaternion q = new Quaternion((float) eigM.get(1, max), (float) eigM.get(2, max), (float) eigM.get(3, max), (float) eigM.get(0, max));
            q.normalize();
            Quaternion qCopy;
            Vector3f p;
            System.out.println("before rotation!");
            // die Rotation wird auf die Daten des Originalmesh angewendet!
            for (int i = 0; i < compF.size(); i++) {
                p = compF.get(i);
                Quaternion point = new Quaternion(p.getX(), p.getY(), p.getZ(), 0);

                qCopy = multiply(point, conjugateQ(q));
                qCopy = multiply(q, qCopy);

                // KHK + meanX,Y,Z ist die Anwendung der Translation via mittlerer Abstand zwischen den beiden Punktewolken
                // ist das richtig? hier wird keine Rotation des Translationsvektors berücksichtigt.
                // important: to calculate the translation the centroid of the Targetimage has to be Rotated!!!
                // T = Sc - Tc*R
                // Sc = Source image center, Tc = center of target image
                // Könnte es sein, dass diese Korrektur vernachlässigt werden kann, da ja iterativ 
                // laufend rotiert und translatiert wird. Dadurch ist die Änderung am Schluss nur noch minimal.
                // KHK Ende Kommentar.
                
                p.setX(qCopy.getX() + meanX);
                p.setY(qCopy.getY() + meanY);
                p.setZ(qCopy.getZ() + meanZ);
            }
            
            System.out.println("Size compF before hack: " + compF.size());
            XYZ2DEM_ImporterHack hack = new XYZ2DEM_ImporterHack();
            hack.khkDisplayXYZ(compF); 
            
            currentIteration++;
          } while (currentIteration <= maxIteration
                && (prevMean == null || Math.abs(prevMean - meanD) > error));
    }
}
