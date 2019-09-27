/* package cz.fidentis.comparison.icp;

This class is part of Fidentis

FIDENTIS – Forensic Facial Identification Software
The project is based on multidisciplinary cooperation. Development is provided by two departments of Masaryk University: 

Human Computer Interaction Laboratory at Faculty of Informatics and Department of Anthropology at Faculty of Science.

Author: Bachelor thesis of Zuzana Ferková (http://is.muni.cz/th/374016/fi_b/bakalarka.pdf)
Advisor: Mgr. Igor Chalás
*/

import javax.vecmath.Vector3f;
import java.util.ArrayList;
import java.util.Set;
import java.util.List;
import java.util.HashSet;
import java.util.TreeMap;

/**
 * Created with IntelliJ IDEA.
 * User: ex3me
 * Date: 2/9/13
 * Time: 9:49 PM
 * Class representing KdTree data structure.
 */
public class KdTree {

    //Orders points based on their x-values, then orders y-values based on x-values, etc.
    private TreeMap<Float,TreeMap<Float, Set<Float>>> kdTree;
    
    // Eine TreeMap implementiert NavigableMap und speichert eine sortierte Abfolge von Element-Paaren, 
    // bei denen jeweils ein Wert (Value) einem Schlüssel (Key) zugewiesen ist. 
    // Die Sortierung erfolgt anhand des Keys.
    
    /* Befehle:
    TreeMap.put(K key, V value) Fügt ein Wertepaar in eine TreeMap ein. Wenn der Key bereits existiert, wird statt dessen der zugehörige Wert überschrieben, sodass ein Key nie doppelt vorhanden sein kann.
    TreeMap.size()              Gibt die Anzahl der Einträge in einer TreeMap zurück.
    TreeMap.keySet()            Returns a Set view of the keys contained in this map.
    TreeMap.subMap(K fromKey, boolean fromInclusive, K toKey, boolean toInclusive) 
                                Returns a view of the portion of this map whose keys range from fromKey to toKey.
    */
    
    private static final int MAX_ERROR = 100;

    public KdTree() {
        kdTree = new TreeMap<Float,TreeMap<Float, Set<Float>>>();
    }

    /**
     * Adds a single point into a data structure.
     * If there is a point with same x,y values, adds z value to the appropriate set.
     * If there is a point with same x value only, creates new Tree for y value and set for z value
     * If there is not a point with same x value, creates new Tree for x,y values and set for z value
     *      *
     * @param point - point to be added
     */
    public void put(Vector3f point){
        if(kdTree.containsKey(point.getX())){
            if(kdTree.get(point.getX()).containsKey(point.getY())){
                kdTree.get(point.getX()).get(point.getY()).add(point.getZ());
            }else{
                Set<Float> list = new HashSet<Float>();
                list.add(point.getZ());
                kdTree.get(point.getX()).put(point.getY(), list);
            }
        }else {
            TreeMap<Float,Set<Float>> tree = new TreeMap<Float, Set<Float>>();
            Set<Float> list = new HashSet<Float>();
            list.add(point.getZ());
            tree.put(point.getY(),list);
            kdTree.put(point.getX(),tree);
        }
    }

    /**
     * Uses method put to add all points from the list into a data structure.
     *
     * @param points - list of points to be added
     */
    public void putAll(List<Vector3f> points){
        for(Vector3f point: points){
         this.put(point);
        }
    }


    /**
     * Method to determine whether there is point in a data structure or not.
     * Uses methods defined in MapTree class.
     *
     * @param point
     * @return true - if data structure contains point, false if data structure doesn't contain point
     */
    public boolean containsPoint(Vector3f point){
       return kdTree.containsKey(point.getX()) &&
              kdTree.get(point.getX()).containsKey(point.getY()) &&
              kdTree.get(point.getX()).get(point.getY()).contains(point.getZ());

    }

    /**
     * Method returns list of points that are within an error range of point.
     * Searches a data structure to find all the points that are in desired radius.
     *
     * @param point
     * @param error - defines the range
     * @return list of points within an error range of point
     */
    public List<Vector3f> pointsInRange(Vector3f point, float error){
        List<Vector3f> pointsInRange = new ArrayList<Vector3f>();

    // while finished when both conditions are false
    while(pointsInRange.size() == 0 || error > MAX_ERROR){   // MAX_ERROR is currently set to 100
        for(Float x : this.kdTree.subMap(point.getX() - error, true, point.getX() + error, true).keySet()){
            //System.out.println("point.getX: " + point.getX() + ", error: " + error);
            for(Float y : this.kdTree.get(x).subMap(point.getY() - error, true, point.getY() + error, true).keySet()){
                //System.out.println("   point.getY: " + point.getY() + ", error: " + error);
              for(Float z : this.kdTree.get(x).get(y)){
                  if(z >= point.getZ() - error && z <= point.getZ() + error){
                      //System.out.println("      point.getZ: " + point.getZ() + ", error: " + error);
                      pointsInRange.add(new Vector3f(x,y,z));
                     // System.out.println("pointsInRange: "+pointsInRange);
                  }
              }
            }
         }
        error = error + 3f;
       }
        return pointsInRange;
        }
    }

