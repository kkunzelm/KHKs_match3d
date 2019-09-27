/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */



import javax.vecmath.*;
/**
 *
 * @author kkunzelm
 */
public class khsTransform3d {
    public Matrix3d rotationMatrix;
    public Vector3d translationVector;

    
    public khsTransform3d(Matrix3d rotMatrix, Vector3d transVector){
        this.rotationMatrix = rotMatrix;
        this.translationVector = transVector;

    }
    
}
