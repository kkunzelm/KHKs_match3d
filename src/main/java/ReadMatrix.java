import ij.IJ;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;

// probably from package imagescience.transform; -> Transform.java
class ReadMatrix {

    private double[][] transformationMatrix = new double[4][4];

    double[][] run(final String file) {

        try {
            if (file == null || file.equals(""))
                throw new IllegalArgumentException("Empty matrix file name");

            transformationMatrix = load(file);

        } catch (OutOfMemoryError e) {
            IJ.error("Not enough memory for this operation");

        } catch (UnknownError e) {
            IJ.error("Could not create output image for some reason.\nPossibly there is not enough free memory");

        } catch (IllegalArgumentException e) {
            IJ.error(e.getMessage());

        }
        return transformationMatrix;
    }

    // from TransformJ TJ_Matrix.java
    private double[][] load(final String file) {

        // Read lines:
        final Vector<String> lines = new Vector<>(10, 10); // 10 lines first, will be incremented for 10 lines, if full
        String line;
        try {
            final BufferedReader br = new BufferedReader(new FileReader(file));
            line = br.readLine();
            while (line != null) {
                line = line.trim();
                if (!line.equals(""))
                    lines.add(line);
                line = br.readLine();
            }
            br.close();
        } catch (FileNotFoundException e) {
            throw new IllegalArgumentException("Unable to find " + file);
        } catch (IOException e) {
            throw new IllegalArgumentException("Error reading from " + file);
        }

        // Convert lines:
        if (lines.size() != 4)
            throw new IllegalArgumentException("File " + file + " does not contain a 4 x 4 matrix");
        String delim = "\t";
        line = lines.get(0);
        if (line.contains(","))
            delim = ",";
        else if (line.contains(" "))
            delim = " ";
        final double[][] matrix = new double[4][4];
        for (int r = 0; r < 4; ++r) {
            line = lines.get(r);
            final StringTokenizer st = new StringTokenizer(line, delim);
            if (st.countTokens() != 4)
                throw new IllegalArgumentException("File " + file + " does not contain a 4 x 4 matrix");
            for (int c = 0; c < 4; ++c) {
                try {
                    matrix[r][c] = Double.parseDouble(st.nextToken());
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Error reading element (" + r + "," + c + ") in " + file);
                }
            }
        }

        // Store matrix:
        return matrix;
    }
}
