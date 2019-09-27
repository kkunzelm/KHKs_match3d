// KHK: status 30.4.2014
// KHK: I downloaded this file
// KHK: good ... the limited search depth, bad: the objects cannot be retrieved on coordinates, just based on neighborhood
// KHK: should be modified, but I do not need quadtrees anylonger at the moment


import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.

    In a typical game, you’ll start by creating the quadtree and passing the bounds of the screen.

    Quadtree quad = new Quadtree(0, new Rectangle(0,0,600,600));

    At every frame, you’ll insert all objects into the quadtree by first clearing the quadtree then using the insert method for every object.


    quad.clear();
    for (int i = 0; i < allObjects.size(); i++) {
      quad.insert(allObjects.get(i));
    }

    Once all objects have been inserted, you’ll go through each object and retrieve a list of objects it could possibly collide with. You'll then check for collisions between each object in the list and the initial object, using a collision detection algorithm.


    List returnObjects = new ArrayList();
    for (int i = 0; i < allObjects.size(); i++) {
      returnObjects.clear();
      quad.retrieve(returnObjects, objects.get(i));

      for (int x = 0; x < returnObjects.size(); x++) {
        // Run collision detection algorithm between objects
      }
    }




 */

/**
 * from http://gamedevelopment.tutsplus.com/tutorials/quick-tip-use-quadtrees-to-detect-likely-collisions-in-2d-space--gamedev-374
 * @author kkunzelm
 */


public class Quadtree {
 
  private final int MAX_OBJECTS = 18;
  private int max_levels;            // todo change this later to a meaningful value (z.B. Anzahl ist Breite/2 oder Höhe/2 je nachdem, was größer ist).
   
  private int level;
  private List objects;
  private Rectangle bounds;
  private Quadtree[] nodes;
 
 /*
  * Constructor
  */
  public Quadtree(int pLevel, Rectangle pBounds) {
   level = pLevel;
   objects = new ArrayList();
   bounds = pBounds;
   nodes = new Quadtree[4];
   max_levels = (int)pBounds.getHeight()/2;
   if(pBounds.getHeight() < pBounds.getWidth()){
       max_levels = (int)pBounds.getWidth()/2;
   }
  }
  
 /*
 * Clears the quadtree
 */
 public void clear() {
   objects.clear();
 
   for (int i = 0; i < nodes.length; i++) {
     if (nodes[i] != null) {
       nodes[i].clear();
       nodes[i] = null;
     }
   }
 }
 
  /*
 * Splits the node into 4 subnodes
 */
 private void split() {
   int subWidth = (int)(bounds.getWidth() / 2);
   int subHeight = (int)(bounds.getHeight() / 2);
   int x = (int)bounds.getX();
   int y = (int)bounds.getY();
 
   nodes[0] = new Quadtree(level+1, new Rectangle(x + subWidth, y, subWidth, subHeight));
   nodes[1] = new Quadtree(level+1, new Rectangle(x, y, subWidth, subHeight));
   nodes[2] = new Quadtree(level+1, new Rectangle(x, y + subHeight, subWidth, subHeight));
   nodes[3] = new Quadtree(level+1, new Rectangle(x + subWidth, y + subHeight, subWidth, subHeight));
 }
 
 /*
 * Determine which node the object belongs to. -1 means
 * object cannot completely fit within a child node and is part
 * of the parent node
 */
 private int getIndex(Triangle3D pTri) {
   int index = -1;
   double verticalMidpoint = bounds.getX() + (bounds.getWidth() / 2);
   double horizontalMidpoint = bounds.getY() + (bounds.getHeight() / 2);
 
   double[] boundingBox = pTri.getBoundingBox();
   
   // Object can completely fit within the top quadrants
   boolean topQuadrant = (boundingBox[2] < horizontalMidpoint && boundingBox[3] < horizontalMidpoint);
   // Object can completely fit within the bottom quadrants
   boolean bottomQuadrant = (boundingBox[2] > horizontalMidpoint);
 
   // Object can completely fit within the left quadrants
   if (boundingBox[0] < verticalMidpoint && boundingBox[1] < verticalMidpoint) {
      if (topQuadrant) {
        index = 1;
      }
      else if (bottomQuadrant) {
        index = 2;
      }
    }
    // Object can completely fit within the right quadrants
    else if (boundingBox[0] > verticalMidpoint) {
     if (topQuadrant) {
       index = 0;
     }
     else if (bottomQuadrant) {
       index = 3;
     }
   }
   System.out.println("Index von QuadTree (via Triangle): " + index);
   return index;
 }
 
 private int getIndex(double[] gridPosition) {            // gridPosition[0] = x, gridPosition[1] = y
   int index = -1;
   double verticalMidpoint = bounds.getX() + (bounds.getWidth() / 2);
   double horizontalMidpoint = bounds.getY() + (bounds.getHeight() / 2);
 
   // Object can completely fit within the top quadrants
   boolean topQuadrant = (gridPosition[1] < horizontalMidpoint);
   // Object can completely fit within the bottom quadrants
   boolean bottomQuadrant = (gridPosition[1] > horizontalMidpoint);
 
   // Object can completely fit within the left quadrants
   if (gridPosition[0] < verticalMidpoint) {
      if (topQuadrant) {
        index = 1;
      }
      else if (bottomQuadrant) {
        index = 2;
      }
    }
    // Object can completely fit within the right quadrants
    else if (gridPosition[0] > verticalMidpoint) {
     if (topQuadrant) {
       index = 0;
     }
     else if (bottomQuadrant) {
       index = 3;
     }
   }
   //System.out.println("Index von QuadTree (via x,y Pos): " + index);
   return index;
 }
 
 
 
 /*
 * Insert the object into the quadtree. If the node
 * exceeds the capacity, it will split and add all
 * objects to their corresponding nodes.
 */
 public void insert(Triangle3D pTri) {
   if (nodes[0] != null) {
     int index = getIndex(pTri);
 
     if (index != -1) {
       nodes[index].insert(pTri);
 
       return;
     }
   }
 
   objects.add(pTri);
 
   if (objects.size() > MAX_OBJECTS && level < max_levels) {
      if (nodes[0] == null) { 
         split(); 
      }
 
     int i = 0;
     while (i < objects.size()) {
       int index = getIndex((Triangle3D)objects.get(i));
       if (index != -1) {
         nodes[index].insert((Triangle3D)objects.remove(i));
       }
       else {
         i++;
       }
     }
   }
 }
 
 /*
 * Return all objects that could collide with the given object
 */
 public List retrieve(List returnObjects, Triangle3D pTri) {
   int index = getIndex(pTri);
   if (index != -1 && nodes[0] != null) {
     nodes[index].retrieve(returnObjects, pTri);
   }
 
   returnObjects.addAll(objects);
 
   return returnObjects;
 }
 
  public List retrieve(List returnObjects, double[] gridPosition) {
   int index = getIndex(gridPosition);
   if (index != -1 && nodes[0] != null) {
     nodes[index].retrieve(returnObjects, gridPosition);
   }
 
   returnObjects.addAll(objects);
 
   return returnObjects;
 }
 
}
