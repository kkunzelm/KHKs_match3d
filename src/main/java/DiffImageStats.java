import java.util.Arrays;

class DiffImageStats {
    double mean;
    int ctr;
    double[] differences;

    private double minDiff;
    private double maxDiff;
    private double variance;
    private double stddev;
    private double hi;
    private double lo;

    DiffImageStats(double mean, int ctr, double[] differences){
        this.mean = mean;
        this.ctr = ctr;
        this.differences = differences;

        mean = mean / ctr;
        double var = 0.0;

        for (int i = 0; i < ctr; i++) {
            if (!Double.isNaN(differences[i]))
                var += (differences[i] - mean) * (differences[i] - mean);
        }

        variance = var / (ctr - 1);
        stddev = Math.sqrt(variance);

        lo = mean - 1.96 * stddev;
        hi = mean + 1.96 * stddev;

        System.out.println("######### IMAGE STATS ############");
        System.out.println("average          = " + mean);
        System.out.println("sample variance  = " + variance);
        System.out.println("sample stddev    = " + stddev);
        System.out.println("95% approximate confidence interval");
        System.out.println("[ " + lo + ", " + hi + " ]");
        System.out.println("##################################");

        Arrays.sort(differences);
        minDiff = differences[0];
        maxDiff = differences[differences.length - 1];
    }

    public double getMean(){
        return mean;
    }

    public double getCtr() {
        return ctr;
    }

    public double getMinDiff() {
        return minDiff;
    }

    public double getMaxDiff() {
        return maxDiff;
    }

    public double getVariance() {
        return variance;
    }

    public double getStddev() {
        return stddev;
    }

    public double getHi() {
        return hi;
    }

    public double getLo() {
        return lo;
    }
}
