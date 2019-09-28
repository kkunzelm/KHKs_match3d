/*
Format Specifications

An ASCII STL file begins with the line:

 solid name

where name is an optional string. The file continues with any number of triangles, each represented as follows:

solid name
 facet normal n1 n2 n3
   outer loop
     vertex v11 v12 v13
     vertex v21 v22 v23
     vertex v31 v32 v33
   endloop
 endfacet
 endsolid name

where n1-n3 and v11-v33 are floating point numbers in sign-mantissa'e'-sign-exponent format and concludes with:

 endsolid name

The structure of the format suggests that other possibilities exist (eg Facets with more than one 'loop' or loops with other than three vertices) but in practice, all facets are simple triangles.

White space (spaces, tabs, newlines) may be used anywhere in the file except within numbers or words. The spaces between 'facet' and 'normal' and between 'outer' and 'loop' are required.

The Facet Normal

In both ASCII and binary versions of STL, the facet normal should be a unit vector pointing outwards from the solid object. In most software this may be set to (0,0,0) and the software will automatically calculate a normal based on the order of the triangle vertices using the 'right hand rule'. 
 
25.4.2014: Modification of SimpleSTL_Writer

In this version I used the datastruct.Vec3 datatype to generate the triangles.
The file is used for debugging purposes to analyse the triangles - in context with the T2sTransformationRangeData debugging.

*/


import java.io.*;

import datastruct.Triangle3D;
import datastruct.Vec3;
import ij.*;
import ij.process.*;
import ij.measure.*;
import ij.plugin.filter.PlugInFilter;
import ij.io.*;
import java.util.ArrayList;
import java.util.List;



// Saves an image described by an ImageProcessor object as a tab-delimited text file. */
// STL File Format
// generated very simple as triangle strips 
// for use with our match3d images
// input image has heigth data encoded as grey value z = I(x,y)
// input data should be on regular grid already 

public class SimpleSTL_WriterVec3Based implements PlugInFilter {
    
	private ImageProcessor ip;
	private Calibration cal;
        double scaleX;
        double scaleY;
	
        int width;
	int height;
        
        static ImagePlus imp;

	private String name;
	private String directory;
        

        
        public int setup(String arg, ImagePlus imp) {

            return DOES_32;
	}
        
        public void run(ImageProcessor ip) {
            imp = WindowManager.getCurrentImage();
            
            FileInfo fi = imp.getFileInfo();
            scaleX = fi.pixelWidth; 
            scaleY = fi.pixelHeight;

            this.ip = ip;


		String path = getPath("STL-File", ".stl");
                
                try {

			//TextEncoder file = new TextEncoder(imp.getProcessor(), cal, precision);
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
			//file.
                        write(out);
			out.close();
		}
		catch (IOException e) {
			showErrorMessage(e);
		}
	}
        
        
        
                String getPath(String type, String extension) {
		name = imp.getTitle();
		SaveDialog sd = new SaveDialog("Save as "+type, name, extension);
		name = sd.getFileName();
		if (name==null)
			return null;
		directory = sd.getDirectory();
		imp.startTiming();
		String path = directory+name;
		return path;
	}
		
        
	/** Saves the image as a text file. */
	public void write(DataOutputStream out) throws IOException {
           
		PrintWriter pw = new PrintWriter(out);
		boolean calibrated = cal!=null && cal.calibrated();    // calibrated wird im MOment nicht verwendet. Ich gehe den Umweg über FileInfo und scaleX/Y
		if (calibrated)
			ip.setCalibrationTable(cal.getCTable());
		else
			ip.setCalibrationTable(null);
		width = ip.getWidth();
		height = ip.getHeight();

		//IJ.showStatus("Exporting as text...");
		

                // KHK old: validPoints = readVerticesForObjectFile();
                
                
                List triangleList = makeTriangleList(imp);
                
                pw.println("solid khk");
                for (Object get : triangleList) {
                    Triangle3D currentTriangle = (Triangle3D)get;
                    double temp = currentTriangle.getVertex0().x;
                    pw.println("facet normal "+currentTriangle.getNormal().x+" "+currentTriangle.getNormal().y+" "+currentTriangle.getNormal().z);
                           pw.println("outer loop");
                           // todo check: if scaleX and Y are different, then it could be that I have to swap scaleX and scaleY here --- see comment line 260ff
                             pw.println("vertex "+ currentTriangle.getVertex0().x+" "+currentTriangle.getVertex0().y+" "+currentTriangle.getVertex0().z);
                             pw.println("vertex "+ currentTriangle.getVertex1().x+" "+currentTriangle.getVertex1().y+" "+currentTriangle.getVertex1().z);
                             pw.println("vertex "+ currentTriangle.getVertex2().x+" "+currentTriangle.getVertex2().y+" "+currentTriangle.getVertex2().z);
                           pw.println("endloop");
                         pw.println("endfacet");
		
		}
	        pw.println("endsolid khk");
                pw.println();
              
		pw.close();
		IJ.showProgress(1.0);
		//IJ.showStatus("");
    
                    
                }
                
	     
		void showErrorMessage(IOException e) {
		String msg = e.getMessage();
		if (msg.length()>100)
			msg = msg.substring(0, 100);
		IJ.error("FileSaver", "An error occured writing the file.\n \n" + msg);
	}
                
    public static List makeTriangleList(ImagePlus imp){
        
        ImageProcessor ip = imp.getProcessor();
        
        int width = ip.getWidth();
	int height = ip.getHeight();
        
        FileInfo fi = imp.getFileInfo();
        double scaleX = fi.pixelWidth; 
        double scaleY = fi.pixelHeight;
                
        // List with triangles
        List triangleList = new ArrayList();      
   
       // get range data and triangulate them
       //
       // mentally arrange the lines like:
       
       // first line:  0 2 4 6 8 ...
       // second line: 1 3 5 7 9 ...
       // etc.
       //
       // Then make Triangles always in the same rotation direction - clockwise or alternatively counterclockwise (here clockwise)
       // Triangle 1: 021  = ----
       //                    | / 
       //                    |/ 
       // Triangle 2: 231  =  /| 
       // ...                / |
       //                   ----
       //
       // translate 0,2,4..., 1,3,5... in k, k+1, l, k+1
       // this should give a triangulated surface

        for (int k=0; k<width-1; k++) {
            for (int l=0; l<height-1; l++) {
                // first triangle 
                if(ip.getPixelValue(k,l)>0 && ip.getPixelValue(k+1,l)>0 && ip.getPixelValue(k,l+1)>0){          // pixels with value zero are not used
                    
                    // vertices of triangle
                    Vec3 A = new Vec3();
                    Vec3 B = new Vec3();
                    Vec3 C = new Vec3();
                    
                    // first vertex
                    A.x = k*scaleX;
                    A.y = l*scaleY;
                    A.z = ip.getPixelValue(k,l);

                    // second vertex
                    B.x = (k+1)*scaleX;
                    B.y = l*scaleY;
                    B.z = ip.getPixelValue(k+1,l);

                    // third vertex
                    C.x = k*scaleX;
                    C.y = (l+1)*scaleY;
                    C.z = ip.getPixelValue(k,l+1);
                    
                    // accumulate temporary vertex array
                    Vec3[] vertices = new Vec3[3];
                    
                    // vertex array
                    vertices[0] = A;
                    vertices[1] = B;
                    vertices[2] = C;
                    
                    // make triangle
                    Triangle3D triangle = new Triangle3D();
                    triangle.makeTriangle3D(vertices);
                    
                    // add current triangle to triangleList
                    triangleList.add(triangle);
                }
                
                // second triangle with same k,l,k+1,l+1
                if(ip.getPixelValue(k+1,l)>0 && ip.getPixelValue(k+1,l+1)>0 && ip.getPixelValue(k,l+1)>0){
                    
                    // vertices of triangle
                    Vec3 A = new Vec3();
                    Vec3 B = new Vec3();    
                    Vec3 C = new Vec3();
                    
                    // first vertex
                    A.x = (k+1)*scaleX;
                    A.y = l*scaleY;
                    A.z = ip.getPixelValue(k+1,l);

                    // second vertex
                    B.x = (k+1)*scaleX;
                    B.y = (l+1)*scaleY;
                    B.z = ip.getPixelValue(k+1,l+1);

                    // third vertex
                    C.x = k*scaleX;
                    C.y = (l+1)*scaleY;
                    C.z = ip.getPixelValue(k,l+1);
                    
                    // accumulate temporary vertex array
                    Vec3[] vertices = new Vec3[3];
                    
                    // vertex array
                    vertices[0] = A;
                    vertices[1] = B;
                    vertices[2] = C;
                    
                    // make triangle
                    Triangle3D triangle = new Triangle3D();
                    triangle.makeTriangle3D(vertices);
                    
                    // add current triangle to triangleList
                    triangleList.add(triangle);
                }
            }
        }

        
        return triangleList;
        
    }
   
}