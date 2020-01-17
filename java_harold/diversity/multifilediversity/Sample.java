// 
// Decompiled by Procyon v0.5.36
// 

package multifilediversity;

import java.util.Iterator;
import java.io.IOException;
import java.util.Collection;
import java.util.Arrays;
import java.util.Collections;
import org.apache.commons.math3.distribution.BinomialDistribution;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.ArrayList;

public class Sample
{
    ArrayList<CountSet> countsArray;
    String name;
    GammaValuesPoly gammaValuesPoly;
    double F;
    double[] avgReadDepth;
    
    Sample(final String fileName) {
        this.countsArray = new ArrayList<CountSet>();
        this.avgReadDepth = new double[2];
        this.name = fileName;
        final HashMap<String, CountSet> countHash = new HashMap<String, CountSet>();
        final int[] biasResults = new int[2];
        try {
            final FileReader file = new FileReader(fileName);
            final BufferedReader buff = new BufferedReader(file);
            String line = buff.readLine();
            boolean eof = false;
            while (!eof) {
                line = buff.readLine();
                if (line == null) {
                    eof = true;
                }
                else {
                    final int[] iCounts = new int[5];
                    final String[] words = line.split(",");
                    final int totCount = 0;
                    if (MultiFileDiversity.fracSampled < 0.0) {
                        for (int iBase = 0; iBase < 4; ++iBase) {
                            iCounts[iBase] = Integer.parseInt(words[3 + iBase]);
                            final int[] array = iCounts;
                            final int n = 4;
                            array[n] += iCounts[iBase];
                        }
                    }
                    else {
                        for (int iBase = 0; iBase < 4; ++iBase) {
                            final BinomialDistribution bd = new BinomialDistribution(Integer.parseInt(words[3 + iBase]), MultiFileDiversity.fracSampled);
                            iCounts[iBase] = bd.sample();
                            final int[] array2 = iCounts;
                            final int n2 = 4;
                            array2[n2] += iCounts[iBase];
                        }
                    }
                    final double[] avgReadDepth = this.avgReadDepth;
                    final int n3 = 0;
                    avgReadDepth[n3] += iCounts[4];
                    final double[] avgReadDepth2 = this.avgReadDepth;
                    final int n4 = 1;
                    ++avgReadDepth2[n4];
                    final Integer[] counts = { iCounts[0], iCounts[1], iCounts[2], iCounts[3], iCounts[4] };
                    if (counts[4] <= 3) {
                        continue;
                    }
                    Arrays.sort(counts, 0, 4, Collections.reverseOrder());
                    final String summary = Arrays.toString(counts);
                    if (!countHash.containsKey(summary)) {
                        final CountSet newCountSet = new CountSet(counts);
                        countHash.put(summary, newCountSet);
                    }
                    else {
                        countHash.get(summary).incrementCopies();
                    }
                }
            }
            System.out.println("*** Finished file " + this.name);
            this.countsArray = new ArrayList<CountSet>(countHash.values());
            System.out.println(this.countsArray.size());
            buff.close();
            file.close();
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
    }
    
    void setAlpha(final double F, final GammaValuesPoly gammaValuesPoly) {
        this.gammaValuesPoly = gammaValuesPoly;
        this.F = F;
    }
    
    double calculateLogLikelihood(final double errorRate) {
        double totalLogLikelihood = 0.0;
        for (final CountSet countSet : this.countsArray) {
            final int[] count = countSet.count;
            final double term1 = Math.log(1.0 - this.F) + count[0] * Math.log(1.0 - errorRate) + (count[4] - count[0]) * Math.log(errorRate);
            final double betaTopPoly = this.gammaValuesPoly.getLGammaAlphaPlusN(0, count[0]) + this.gammaValuesPoly.getLGammaAlphaPlusN(1, count[1]) + this.gammaValuesPoly.getLGammaAlphaPlusN(2, count[2]) + this.gammaValuesPoly.getLGammaAlphaPlusN(3, count[3]) - this.gammaValuesPoly.getLGammaSumAlphaPlusN(count[4]);
            final double betaBottomPoly = this.gammaValuesPoly.lGammaAlpha[0] + this.gammaValuesPoly.lGammaAlpha[1] + this.gammaValuesPoly.lGammaAlpha[2] + this.gammaValuesPoly.lGammaAlpha[3] - this.gammaValuesPoly.lGammaSumAlpha;
            final double term2 = Math.log(this.F) + betaTopPoly - betaBottomPoly;
            final double bigTerm = Math.max(term1, term2);
            final double smallTerm = Math.min(term1, term2);
            final double diff = -Math.abs(term1 - term2);
            final double logLikelihood = bigTerm + Math.log(1.0 + Math.exp(diff));
            totalLogLikelihood += countSet.nCopies * logLikelihood;
        }
        return totalLogLikelihood;
    }
    
    void getDivergence(final double errorRate) {
        final double[] avgD = new double[4];
        for (final CountSet countSet : this.countsArray) {
            final int[] count = countSet.count;
            final double term1 = Math.log(1.0 - this.F) + count[0] * Math.log(1.0 - errorRate) + (count[4] - count[0]) * Math.log(errorRate);
            final double betaTopPoly = this.gammaValuesPoly.getLGammaAlphaPlusN(0, count[0]) + this.gammaValuesPoly.getLGammaAlphaPlusN(1, count[1]) + this.gammaValuesPoly.getLGammaAlphaPlusN(2, count[2]) + this.gammaValuesPoly.getLGammaAlphaPlusN(3, count[3]) - this.gammaValuesPoly.getLGammaSumAlphaPlusN(count[4]);
            final double betaBottomPoly = this.gammaValuesPoly.lGammaAlpha[0] + this.gammaValuesPoly.lGammaAlpha[1] + this.gammaValuesPoly.lGammaAlpha[2] + this.gammaValuesPoly.lGammaAlpha[3] - this.gammaValuesPoly.lGammaSumAlpha;
            final double term2 = Math.log(this.F) + betaTopPoly - betaBottomPoly;
            final double bigTerm = Math.max(term1, term2);
            final double probVar = Math.exp(term2 - bigTerm) / (Math.exp(term1 - bigTerm) + Math.exp(term2 - bigTerm));
            final double[] alpha = this.gammaValuesPoly.alpha;
            final double sumAlpha = this.gammaValuesPoly.sumAlpha;
            final double diverg = ((alpha[0] + count[0] + 1.0) * (alpha[0] + count[0]) + (alpha[1] + count[1] + 1.0) * (alpha[1] + count[1]) + (alpha[2] + count[2] + 1.0) * (alpha[2] + count[2]) + (alpha[3] + count[3] + 1.0) * (alpha[3] + count[3])) / ((sumAlpha + count[4] + 1.0) * (sumAlpha + count[4]));
            final double[] frac = { 1.0 * count[0] / count[4], 1.0 * count[1] / count[4], 1.0 * count[2] / count[4], 1.0 * count[3] / count[4] };
            final double[] array = avgD;
            final int n = 0;
            array[n] += countSet.nCopies * probVar * (1.0 - diverg);
            if (frac[0] < 0.99) {
                final double[] array2 = avgD;
                final int n2 = 1;
                array2[n2] += probVar * countSet.nCopies * (count[4] / (count[4] - 1.0)) * (1.0 - frac[0] * frac[0] - frac[1] * frac[1] - frac[2] * frac[2] - frac[3] * frac[3]);
                final double[] array3 = avgD;
                final int n3 = 2;
                array3[n3] += countSet.nCopies * (count[4] / (count[4] - 1.0)) * (1.0 - frac[0] * frac[0] - frac[1] * frac[1] - frac[2] * frac[2] - frac[3] * frac[3]);
            }
            final double[] array4 = avgD;
            final int n4 = 3;
            array4[n4] += countSet.nCopies;
        }
        double localFrac = 1.0;
        if (MultiFileDiversity.fracSampled > 0.0) {
            localFrac = MultiFileDiversity.fracSampled;
        }
        System.out.println(this.name + "\t" + localFrac + "\t" + this.avgReadDepth[0] / this.avgReadDepth[1] + "\t" + errorRate + "\t" + avgD[0] / avgD[3] + "\t" + avgD[1] / avgD[3] + "\t" + avgD[2] / avgD[3]);
    }
    
    double calculateLogLikelihood(final double F, final GammaValuesPoly gammaValuesPoly) {
        double totalLogLikelihood = 0.0;
        for (final CountSet countSet : this.countsArray) {
            final int[] count = countSet.count;
            final int lowerLimit = 1 + count[4] / 500;
            double logLikelihood = 0.0;
            if (count[4] - count[0] > lowerLimit) {
                logLikelihood = Math.log(F) + this.evaluateDirichletPoly(count[0], count[1], count[2], count[3], count[4], gammaValuesPoly);
            }
            else {
                final ArrayList<Double> termArray = new ArrayList<Double>();
                termArray.add(Math.log(1.0 - F));
                double maxVal = Math.log(1.0 - F);
                for (int count2 = 0; count2 <= lowerLimit; ++count2) {
                    for (int count3 = 0; count3 <= Math.min(count2, lowerLimit - count2); ++count3) {
                        for (int count4 = 0; count4 <= Math.min(count3, lowerLimit - count2 - count3); ++count4) {
                            final int count5 = count[4] - count2 - count3 - count4;
                            final double val = Math.log(F) + this.evaluateDirichletPoly(count5, count2, count3, count4, count[4], gammaValuesPoly);
                            termArray.add(val);
                            maxVal = Math.max(maxVal, val);
                        }
                    }
                }
                double likelihood = 0.0;
                for (final double val2 : termArray) {
                    likelihood += Math.exp(val2 - maxVal);
                }
                logLikelihood = maxVal + Math.log(likelihood);
            }
            totalLogLikelihood += countSet.nCopies * logLikelihood;
        }
        return totalLogLikelihood;
    }
    
    double evaluateDirichletPoly(final int count0, final int count1, final int count2, final int count3, final int count4, final GammaValuesPoly gammaValuesPoly) {
        final double logFact = gammaValuesPoly.getLGammaN(count4 + 1) - gammaValuesPoly.getLGammaN(count0 + 1) - gammaValuesPoly.getLGammaN(count1 + 1) - gammaValuesPoly.getLGammaN(count2 + 1) - gammaValuesPoly.getLGammaN(count3 + 1);
        final double betaTopPoly = gammaValuesPoly.getLGammaAlphaPlusN(0, count0) + gammaValuesPoly.getLGammaAlphaPlusN(1, count1) + gammaValuesPoly.getLGammaAlphaPlusN(2, count2) + gammaValuesPoly.getLGammaAlphaPlusN(3, count3) - gammaValuesPoly.getLGammaSumAlphaPlusN(count4);
        final double betaBottomPoly = gammaValuesPoly.lGammaAlpha[0] + gammaValuesPoly.lGammaAlpha[1] + gammaValuesPoly.lGammaAlpha[2] + gammaValuesPoly.lGammaAlpha[3] - gammaValuesPoly.lGammaSumAlpha;
        return logFact + betaTopPoly - betaBottomPoly;
    }
    
    private class CountSet
    {
        int[] count;
        int nCopies;
        
        CountSet(final Integer[] inputCount) {
            this.count = new int[5];
            this.nCopies = 1;
            for (int i = 0; i < 5; ++i) {
                this.count[i] = inputCount[i];
            }
            this.nCopies = 1;
        }
        
        void incrementCopies() {
            ++this.nCopies;
        }
    }
}
