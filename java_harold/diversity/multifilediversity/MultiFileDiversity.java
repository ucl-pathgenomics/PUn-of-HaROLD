// 
// Decompiled by Procyon v0.5.36
// 

package multifilediversity;

import java.io.IOException;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.Iterator;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.exception.MathIllegalStateException;
import java.util.Arrays;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Random;

public class MultiFileDiversity
{
    static Random random;
    static int nSites;
    ArrayList<Sample> samples;
    static double fracSampled;
    HashMap<String, double[]> storedValuesHash;
    
    public static void main(final String[] args) {
        final MultiFileDiversity mfd = new MultiFileDiversity();
        mfd.run(args);
    }
    
    MultiFileDiversity() {
        this.samples = new ArrayList<Sample>();
        this.storedValuesHash = new HashMap<String, double[]>();
    }
    
    void run(final String[] args) {
        final CommandLineValues commandLineValues = new CommandLineValues(args);
        if (commandLineValues.getStoredValues()) {
            this.storedValuesHash = commandLineValues.getStoredValuesHash();
        }
        MultiFileDiversity.fracSampled = commandLineValues.getFracSampled();
        this.readData(commandLineValues.getDataFileList());
        final MultivariateOptimizer optimize = (MultivariateOptimizer)new BOBYQAOptimizer(10, 1.0E-4, 1.0E-6);
        final double[] params = { 2.3, 0.4, 0.016, 0.001, 0.04 };
        final double[] lb = { 0.001, 0.001, 1.0E-10, 1.0E-10, 0.001 };
        final double[] ub = { 100.0, 100.0, 100.0, 100.0, 0.999 };
        for (final Sample sample : this.samples) {
            double[] best = new double[5];
            final OptimiseAlpha oa = new OptimiseAlpha(sample);
            if (this.storedValuesHash.containsKey(sample.name)) {
                best = this.storedValuesHash.get(sample.name);
            }
            else {
                final OptimizationData[] parm = { (OptimizationData)new InitialGuess(params), (OptimizationData)new MaxEval(100000), (OptimizationData)GoalType.MINIMIZE, (OptimizationData)new ObjectiveFunction((MultivariateFunction)oa), (OptimizationData)new SimpleBounds(lb, ub) };
                try {
                    best = optimize.optimize(parm).getPoint();
                }
                catch (MathIllegalStateException mis) {
                    System.out.println("Convergence problem, taking best so far");
                    best = Arrays.copyOf(oa.bestSoFar, oa.bestSoFar.length);
                }
            }
            final double[] alphas = Arrays.copyOf(best, 4);
            final double F = best[4];
            final GammaValuesPoly gammaValuesPoly = new GammaValuesPoly(alphas);
            sample.setAlpha(F, gammaValuesPoly);
            System.out.format("kkk\t%s\t%.8f\t%.8f\t%.8g\t%.8g\t%.8g\n", sample.name, best[0], best[1], best[2], best[3], F);
        }
        double errorRate = commandLineValues.getErrorRate();
        if (errorRate < 0.0) {
            System.out.println("Optimising error");
            final OptimiseError oe = new OptimiseError(this.samples);
            errorRate = oe.fmin(1.0E-5, 0.1, 1.0E-6);
        }
        for (final Sample sample2 : this.samples) {
            sample2.getDivergence(errorRate);
        }
    }
    
    public void readData(final File dataFileList) {
        try {
            final FileReader file = new FileReader(dataFileList);
            final BufferedReader buff = new BufferedReader(file);
            boolean eof = false;
            while (!eof) {
                final String line = buff.readLine();
                if (line == null) {
                    eof = true;
                }
                else {
                    final Sample newSample = new Sample(line);
                    this.samples.add(newSample);
                }
            }
            buff.close();
            file.close();
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
    }
    
    static {
        MultiFileDiversity.random = new Random();
        MultiFileDiversity.nSites = 235650;
        MultiFileDiversity.fracSampled = 1.0;
    }
    
    private class OptimiseError
    {
        ArrayList<Sample> samples;
        int iIter;
        
        OptimiseError(final ArrayList<Sample> samples) {
            this.iIter = 0;
            this.samples = samples;
        }
        
        public double value(final double errorRate) {
            double llh = 0.0;
            for (final Sample sample : this.samples) {
                llh += sample.calculateLogLikelihood(errorRate);
            }
            ++this.iIter;
            if (this.iIter % 1 == 0) {
                System.out.println(this.iIter + "\t" + errorRate + "\t" + llh);
            }
            return -llh;
        }
        
        private double fmin(double a, double b, final double tol) {
            final double c = 0.5 * (3.0 - Math.sqrt(5.0));
            double d = 0.0;
            double eps = 1.2E-16;
            double tol2 = eps + 1.0;
            eps = Math.sqrt(eps);
            double w;
            double x;
            double v = x = (w = a + c * (b - a));
            double e = 0.0;
            double fv;
            double fw;
            double fx = fw = (fv = this.value(x));
            final double tol3 = tol / 3.0;
            double xm = 0.5 * (a + b);
            tol2 = eps * Math.abs(x) + tol3;
            for (double t2 = 2.0 * tol2; Math.abs(x - xm) > t2 - 0.5 * (b - a); xm = 0.5 * (a + b), tol2 = eps * Math.abs(x) + tol3, t2 = 2.0 * tol2) {
                double r;
                double p;
                double q = p = (r = 0.0);
                if (Math.abs(e) > tol2) {
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0) {
                        p = -p;
                    }
                    else {
                        q = -q;
                    }
                    r = e;
                    e = d;
                }
                if (Math.abs(p) < Math.abs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
                    d = p / q;
                    final double u = x + d;
                    if (u - a < t2 || b - u < t2) {
                        d = tol2;
                        if (x >= xm) {
                            d = -d;
                        }
                    }
                }
                else {
                    if (x < xm) {
                        e = b - x;
                    }
                    else {
                        e = a - x;
                    }
                    d = c * e;
                }
                double u;
                if (Math.abs(d) >= tol2) {
                    u = x + d;
                }
                else if (d > 0.0) {
                    u = x + tol2;
                }
                else {
                    u = x - tol2;
                }
                final double fu = this.value(u);
                if (fx <= fu) {
                    if (u < x) {
                        a = u;
                    }
                    else {
                        b = u;
                    }
                }
                if (fu <= fx) {
                    if (u < x) {
                        b = x;
                    }
                    else {
                        a = x;
                    }
                    v = w;
                    fv = fw;
                    w = x;
                    fw = fx;
                    x = u;
                    fx = fu;
                }
                else if (fu <= fw || w == x) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
            return x;
        }
    }
    
    private class OptimiseAlpha implements MultivariateFunction
    {
        Sample sample;
        int iIter;
        double[] bestSoFar;
        double topLlh;
        
        OptimiseAlpha(final Sample sample) {
            this.iIter = 0;
            this.bestSoFar = new double[5];
            this.topLlh = -1.0E7;
            this.sample = sample;
        }
        
        public double value(final double[] point) {
            final double[] alphas = Arrays.copyOf(point, 4);
            final double F = point[4];
            final GammaValuesPoly gammaValuesPoly = new GammaValuesPoly(alphas);
            final double llh = this.sample.calculateLogLikelihood(F, gammaValuesPoly);
            if (llh > this.topLlh) {
                this.topLlh = llh;
                this.bestSoFar = Arrays.copyOf(point, point.length);
            }
            ++this.iIter;
            if (this.iIter % 100 == 0) {
                final double sumAlpha = point[0] + point[1] + point[2] + point[3];
                final double probMatch = (point[0] * (1.0 + point[0]) + point[1] * (1.0 + point[1]) + point[2] * (1.0 + point[2]) + point[3] * (1.0 + point[3])) / ((sumAlpha + 1.0) * (sumAlpha + 2.0));
                final double var = F * (1.0 - probMatch);
                System.out.println(this.iIter + "\t" + Arrays.toString(point) + "\t" + llh + "\t" + var);
            }
            return -llh;
        }
    }
}
