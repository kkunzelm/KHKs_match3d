import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

class WriteMatrix {

    void run(final double[][] transformationMatrix) {

        // basic idea from TransformJ and imagescience but modified for my purposes

        String filename = "";

        final String prefix = "";
        final String delim = "\t";
        final String postfix = "\n";

        StringBuilder sb = new StringBuilder();
        for (int r = 0; r < 4; ++r) {
            sb.append(prefix);
            for (int c = 0; c < 4; ++c) {
                // formatiert noch nicht richtig
                sb.append(d2s(transformationMatrix[r][c]));
                // sb.append(fmt.d2s(transformationMatrix[r][c]));
                if (c < 3)
                    sb.append(delim);
            }
            sb.append(postfix);
        }

        // Set System Look&Feel
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (UnsupportedLookAndFeelException | ClassNotFoundException | InstantiationException
                | IllegalAccessException e) {
            // handle exception
        }
        JFileChooser saveDialog = new JFileChooser();

        FileNameExtensionFilter zipExtensionFilter = new FileNameExtensionFilter("ZIP File(*.zip)", "zip");
        saveDialog.addChoosableFileFilter(zipExtensionFilter);
        int saveDialogReturn = saveDialog.showSaveDialog(null);

        if (saveDialogReturn == JFileChooser.APPROVE_OPTION) {
            filename = saveDialog.getSelectedFile().getAbsolutePath();
        }

        try {
            final BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
            bw.write(sb.toString());
            bw.close();
        } catch (IOException e) {
            throw new IllegalArgumentException("Error writing to " + filename);
        }
    }

    // double to formated string (from imagescience/utility/Formatter.java)

    /**
     * Returns a {@code String} representation of a {@code double} value.
     *
     * @param d the {@code double} value to be represented.
     * @return a new {@code String} object containing a string representation of
     * {@code d}. The maximum number of decimals used in representing
     * {@code d} can be specified with method . The value of {@code d} is
     * rounded to the specified maximum number of decimals. The returned
     * string will contain less than the maximum number of decimals if
     * {@code d} can be represented exactly that way. In particular, if
     * {@code d} is equal to an integer value, the returned string
     * represents that integer value, without decimals and preceding decimal separator symbol. The string returned when
	 * {@code Double.isNaN(d)}
     * yields {@code true} can be specified with method . Similarly, the string returned when
	 * {@code Double.isInfinite(d)} yields {@code true} can be specified with method. The returned string is "0" if the
     * absolute value of {@code d} is less than the limit set with method .
     */
    private String d2s(final double d) {

        final DecimalFormat edf = new DecimalFormat("0.#E0", new DecimalFormatSymbols(Locale.US));
        final DecimalFormat ndf = new DecimalFormat("0.#", new DecimalFormatSymbols(Locale.US));
        ndf.setMaximumFractionDigits(8);

        double limit = 1.0E-12;
        double dflobo = 1.0E-10; // war original 0.1;

        String nan = "NaN";
        String inf = "Inf";

        if (Double.isNaN(d))
            return nan;
        if (Double.isInfinite(d))
            return (d < 0) ? ("-" + inf) : ("+" + inf); // Abfrage pos oder neg infinity

        final long dr = Math.round(d);
        if (dr == d)
            return String.valueOf(dr);

        // Conditional operator is also known as the ternary operator.
        // This operator consists of three operands and is used to evaluate Boolean expressions.
        // The goal of the operator is to decide which value should be assigned to the variable.
        // The operator is written as:

        // variable x = (expression) ? value if true : value if false

        final double da = (d < 0) ? -d : d;
        if (da < limit)
            return "0";

        String ds;
        if (da < dflobo || da > 10000000)
            ds = edf.format(d);
        else
            ds = ndf.format(d);

        return ds;
    }
}
