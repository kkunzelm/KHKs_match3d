/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import datastruct.Ray;
import datastruct.Triangle3D;
import datastruct.Vec3;

/**
 *
 * @author kkunzelm
 */
public class KHKs_RaytracerTest {
    
    public static void main(String[] args){

      // vertices of triangle
        Vec3 A = new Vec3();
        Vec3 B = new Vec3();    
        Vec3 C = new Vec3();

        // sample vertices
        A.x = 2;
        A.y = 2;
        A.z = 4;

        B.x = 4;
        B.y = 2;
        B.z = 6;

        C.x = 2;
        C.y = 4;
        C.z = 4;

        // accumulate vertex array
        Vec3[] vertices = new Vec3[3];

        vertices[0] = A;
        vertices[1] = B;
        vertices[2] = C;

        // make triangle
        Triangle3D sampleTriangle = new Triangle3D();
        sampleTriangle.makeTriangle3D(vertices);

        // datastruct.Ray
        Ray sampleRay;
        Vec3 origin = new Vec3();
        Vec3 direction = new Vec3();

        // sample origin
        origin.x = 3;
        origin.y = 3;
        origin.z = 0;

        // sample direction
        direction.x = 0;
        direction.y = 0;
        direction.z = 1;

        sampleRay = new Ray(origin, direction);

        // do intersection test

        Vec3 result = sampleTriangle.intersectionTest(sampleRay);


        if (result != null){
            System.out.println("datastruct.Ray hits triangle at:" + result.toString());
        }
        else {
            System.out.println("Triangle missed");
        }
    }
    
}
