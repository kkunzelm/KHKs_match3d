// KHK: status 30.4.2014
// KHK: I downloaded this file with the accompaning files 
// KHK: I used it at an intermediate devolopment state for my target to source transformation
//      it worked fine, but was very slow
//      I learned a lot from these files. Especially the math behind raytracing!! Many thanks to the author.
//      I also like the data structures he implemented: vec3, Triangle3d...
//      This implementation of Triangle3d was modified by me to add several methods to return values (like boundingBox ...)


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 * siehe u.a. hier: http://www.blackpawn.com/texts/pointinpoly/default.html und 
 * hier:            http://www.uninformativ.de

        http://www.uninformativ.de/?section=news&ndo=single&newsid=117
        http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf
 *
 */

/**
 *
 * @author 
 */

import static java.lang.Math.*;

public class Triangle3D {
    
    double alpha1;
    private Vec3 n, uBeta , uGamma ;
    private double d, kBeta , kGamma ;
    private Vec3 vertex1, vertex2, vertex0;
    /** Vorberechnen statischer Werte . */
    public void makeTriangle3D(Vec3[] vertices){

        vertex0 = vertices[0];
        vertex1 = vertices[1];
        vertex2 = vertices[2];
        
        // Kantenvektoren .
        Vec3 b = vertices[1].minus(vertices[0]) ;
        Vec3 c = vertices[2].minus(vertices[0]) ;

        // Ebenenparameter der durch das Dreieck aufgespannten Ebene.
        n = b.cross(c);
        n.normalize ();
        d = n.dot(vertices[0]) ;

        // Berechne die grossen Klammern , sodass beim Schnitttest weniger
        // Arbeit getan werden muss .
        double bb = b. dot (b);
        double bc = b. dot (c);
        double cc = c. dot (c);
        double D = 1.0 / (cc * bb - bc * bc);
        double bbD = bb * D;
        double bcD = bc * D;
        double ccD = cc * D;
        uBeta = b. times ( ccD ). minus (c. times ( bcD ));
        uGamma = c. times ( bbD ). minus (b. times ( bcD ));
        kBeta = -vertices [0]. dot ( uBeta );
        kGamma = -vertices [0]. dot ( uGamma );
        
        /*
        another resource:
        //https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
        Transcribed from Christer Ericson's Real-Time Collision Detection (which, incidentally, is an excellent book):

        // Compute barycentric coordinates (u, v, w) for
        // point p with respect to triangle (a, b, c)
        void Barycentric(Point a, Point b, Point c, float &u, float &v, float &w)
        {
            Vector v0 = b - a, v1 = c - a, v2 = p - a;
            float d00 = Dot(v0, v0);
            float d01 = Dot(v0, v1);
            float d11 = Dot(v1, v1);
            float d20 = Dot(v2, v0);
            float d21 = Dot(v2, v1);
            float denom = d00 * d11 - d01 * d01;
            v = (d11 * d20 - d01 * d21) / denom;
            w = (d00 * d21 - d01 * d20) / denom;
            u = 1.0f - v - w;
        }
        
        
        */
        
        
        
    }

    public double getAlpha(){
        return alpha1;
    }
    
    public Vec3 getNormal(){
        return n;
    }
    
    public Vec3 getVertex0(){
        return vertex0;
    }
    
    public Vec3 getVertex1(){
        return vertex1;
    }
    public Vec3 getVertex2(){
        return vertex2;
    }
    
    public Vec3 getBaryCenter(){
        double[] bb = getBoundingBox(); // bb = boundingBox
        Vec3 center = new Vec3();
        center.x = (bb[0]+(bb[1]-bb[0])/2);
        center.y = (bb[2]+(bb[3]-bb[2])/2);
        center.z = (bb[4]+(bb[5]-bb[4])/2);
        
        return center;
        
    }
    
    public double[] getBoundingBox(){
        // bounding box
        double x0, x1, y0, y1, z0, z1;
        x0 = y0 = z0 = Double.POSITIVE_INFINITY;
        x1 = y1 = z1 = Double.NEGATIVE_INFINITY;
        

            if (vertex0.x < x0){
                x0 = vertex0.x;
            } 
            if (vertex0.x > x1){
                x1 = vertex0.x;
            }
            if (vertex0.y < y0){
                y0 = vertex0.y;
            } 
            if (vertex0.y > y1){
                y1 = vertex0.y;
            }
            if (vertex0.z < z0){
                z0 = vertex0.z;
            } 
            if (vertex0.z > z1){
                z1 = vertex0.z;
            }
        
            if (vertex1.x < x0){
                x0 = vertex1.x;
            } 
            if (vertex1.x > x1){
                x1 = vertex1.x;
            }
            if (vertex1.y < y0){
                y0 = vertex1.y;
            } 
            if (vertex1.y > y1){
                y1 = vertex1.y;
            }
            if (vertex1.z < z0){
                z0 = vertex1.z;
            } 
            if (vertex1.z > z1){
                z1 = vertex1.z;
            }
 
            if (vertex2.x < x0){
                x0 = vertex2.x;
            } 
            if (vertex2.x > x1){
                x1 = vertex2.x;
            }
            if (vertex2.y < y0){
                y0 = vertex2.y;
            } 
            if (vertex2.y > y1){
                y1 = vertex2.y;
            }
            if (vertex2.z < z0){
                z0 = vertex2.z;
            } 
            if (vertex2.z > z1){
                z1 = vertex2.z;
            }
            
        double[] boundingBox = new double[6];
        
        boundingBox[0]=x0 ;
        boundingBox[1]=x1 ;
        boundingBox[2]=y0 ;
        boundingBox[3]=y1 ;
        boundingBox[4]=z0 ;
        boundingBox[5]=z1 ;
        return boundingBox;
    }
    
    /** Fuehrt einen Schnitttest mit einem Strahl durch . */
    public Vec3 intersectionTest( Ray r){
        // Schnitttest Ray -> Ebene 
        // KHK: "Strahl parallel zur Ebene" oder "Strahl in Ebene" Test
        double rn = r.direction.dot (n);
        if ( abs(rn) < 1e-15){              // KHK Bedingung: vec_u dot vec_n darf nicht null sein (ist im Nenner nicht definiert) 
            return null;
        }

        // Wie weit hat es der Ray von seinem Ursprung zum Schnittpunkt ?
        // KHK: positive Richtung/negative Richtung - alpha > 0 => Punkt über grid, alpha < 0 Punkt unter grid
        alpha1 = (d - r.origin.dot(n)) / rn;         // alpha ist Skalierungsfaktor für Strahl vec_q = vec_p + alpha * vec_u (vec_p = Aufpunkt, vec_u = Richtung genormt, vec_q = Geradengleichung in Normalenform

        // KHK bei Raytracing nur positive Werte erwünscht, ich möchte alpha in beide Richtungen
        /*if ( alpha1 <= 0.0){        
             return null;
        }
 */
        // Schnittpunkt q mit Ebene gefunden , liegt er im Dreieck ?
        Vec3 q = r.evaluate(alpha1);                // q = Abstand des Punktes zur Ebene
        
        double beta = uBeta.dot(q) + kBeta ;
        if ( beta < 0.0){
            return null;
        }

        double gamma = uGamma.dot(q) + kGamma ;
        if ( gamma < 0.0){
            return null;
        }

        double alpha = 1 - beta - gamma ;
        if ( alpha < 0.0){
            return null;
        }
                

        // Es gibt also einen Schnittpunkt und dieser ist q. alpha , beta
         // und gamma koeonnten fuer Texturierung oder aehnliches weiter -
         // verwendet werden .
         return q;
     }
 }