package refineHaplotypes;

import java.util.Iterator;
import java.util.Arrays;
import java.util.ArrayList;

public class Parameters
{
    double totalLogLikelihood;
    double totPenalty;
    double[][][] probBase;
    double[] hapFreq;
    int[][] iMaxBase;
    int nHaplo;
    double adjLogLikelihood;
    ArrayList<Integer>[] polySites;
    int nActive;
    
    Parameters(final double totalLogLikelihood, final double totPenalty, final double[] hapFreq, final double[][][] probBase, final int[][] iMaxBase, final ArrayList<Integer>[] polySites) {
        this.totalLogLikelihood = 0.0;
        this.totPenalty = 0.0;
        this.probBase = null;
        this.hapFreq = null;
        this.iMaxBase = null;
        this.nHaplo = 0;
        this.adjLogLikelihood = 0.0;
        this.nActive = 0;
        this.totalLogLikelihood = totalLogLikelihood;
        this.totPenalty = totPenalty;
        this.adjLogLikelihood = totalLogLikelihood - totPenalty;
        this.hapFreq = hapFreq;
        this.probBase = probBase;
        this.iMaxBase = iMaxBase;
        this.nHaplo = hapFreq.length;
        this.polySites = polySites;
        this.nActive = 0;
        for (int iHaplo = 0; iHaplo < this.nHaplo; ++iHaplo) {
            if (hapFreq[iHaplo] > 1.0E-10) {
                ++this.nActive;
            }
        }
    }
    
    Parameters recombineHaplo(final int nSites, final ArrayList<Integer> activeSites) {
        final double[] newHapFreq = Arrays.copyOf(this.hapFreq, this.nHaplo);
        final double[][][] newProbBase = new double[this.nHaplo][nSites][4];
        for (final int iSite : activeSites) {
            for (int iHaplo = 0; iHaplo < this.nHaplo; ++iHaplo) {
                newProbBase[iHaplo][iSite] = Arrays.copyOf(this.probBase[iHaplo][iSite], 4);
            }
        }
        boolean ok = false;
        int iTry = 0;
        int iHaplo = 0;
        int jHaplo = 0;
        while (!ok) {
            iHaplo = RefineHaplotypes.random.nextInt(this.nHaplo);
            jHaplo = (iHaplo + RefineHaplotypes.random.nextInt(this.nHaplo - 1) + 1) % this.nHaplo;
            final double ranForm = RefineHaplotypes.random.nextInt(4);
            final int iBreak = RefineHaplotypes.random.nextInt(nSites);
            final int jBreak = iBreak + Math.round((float)Math.round(RefineHaplotypes.nd.sample()));
            final int lBreak = Math.min(iBreak, jBreak);
            final int rBreak = Math.max(iBreak, jBreak);
            int nSwitch = 0;
            if (this.hapFreq[iHaplo] > 0.01 && this.hapFreq[jHaplo] > 0.01) {
                for (final int iSite2 : activeSites) {
                    if (iSite2 > lBreak && iSite2 < rBreak) {
                        if (ranForm == 2.0) {
                            newProbBase[iHaplo][iSite2] = Arrays.copyOf(this.probBase[jHaplo][iSite2], 4);
                            newProbBase[jHaplo][iSite2] = Arrays.copyOf(this.probBase[jHaplo][iSite2], 4);
                        }
                        else if (ranForm == 3.0) {
                            newProbBase[iHaplo][iSite2] = Arrays.copyOf(this.probBase[iHaplo][iSite2], 4);
                            newProbBase[jHaplo][iSite2] = Arrays.copyOf(this.probBase[iHaplo][iSite2], 4);
                        }
                        else {
                            newProbBase[iHaplo][iSite2] = Arrays.copyOf(this.probBase[jHaplo][iSite2], 4);
                            newProbBase[jHaplo][iSite2] = Arrays.copyOf(this.probBase[iHaplo][iSite2], 4);
                        }
                        ++nSwitch;
                    }
                }
            }
            ++iTry;
            ok = (nSwitch > 0 || iTry > 10);
        }
        System.out.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, iHaplo, jHaplo));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, iHaplo, jHaplo));
        final Parameters newParameters = new Parameters(0.0, 0.0, newHapFreq, newProbBase, null, null);
        return newParameters;
    }
    
    Parameters compressHaplo(final int iHaplo, final int jHaplo, final int nSites) {
        System.out.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, iHaplo, jHaplo));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, iHaplo, jHaplo));
        final double[] contractedHapFreq = new double[this.nHaplo - 1];
        final double[][][] contractedProbBase = new double[this.nHaplo - 1][nSites][4];
        int iPoint = 0;
        for (int kHaplo = 0; kHaplo < this.nHaplo; ++kHaplo) {
            if (kHaplo != iHaplo && kHaplo != jHaplo) {
                contractedHapFreq[iPoint] = this.hapFreq[kHaplo];
                for (int iSite = 0; iSite < nSites; ++iSite) {
                    contractedProbBase[iPoint][iSite] = Arrays.copyOf(this.probBase[kHaplo][iSite], 4);
                }
                ++iPoint;
            }
        }
        contractedHapFreq[this.nHaplo - 2] = this.hapFreq[iHaplo] + this.hapFreq[jHaplo];
        for (int iSite2 = 0; iSite2 < nSites; ++iSite2) {
            for (int iBase = 0; iBase < 4; ++iBase) {
                contractedProbBase[this.nHaplo - 2][iSite2][iBase] = 0.5 * (this.probBase[iHaplo][iSite2][iBase] + this.probBase[jHaplo][iSite2][iBase]);
            }
        }
        final int contractedNHaplo = this.nHaplo - 1;
        final Parameters newParameters = new Parameters(0.0, 0.0, contractedHapFreq, contractedProbBase, null, null);
        return newParameters;
    }
    
    Parameters expandHaplo(final int expandHaplo, final int nSites) {
        System.out.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, expandHaplo));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, expandHaplo));
        final double[] expandedHapFreq = new double[this.nHaplo + 1];
        final double[][][] expandedProbBase = new double[this.nHaplo + 1][nSites][4];
        for (int iHaplo = 0; iHaplo < this.nHaplo; ++iHaplo) {
            expandedHapFreq[iHaplo] = this.hapFreq[iHaplo];
            for (int iSite = 0; iSite < nSites; ++iSite) {
                expandedProbBase[iHaplo][iSite] = Arrays.copyOf(this.probBase[iHaplo][iSite], 4);
            }
        }
        for (int iSite2 = 0; iSite2 < nSites; ++iSite2) {
            expandedProbBase[this.nHaplo][iSite2] = Arrays.copyOf(this.probBase[expandHaplo][iSite2], 4);
        }
        double avgFrac = 0.0;
        for (final int iSite3 : this.polySites[expandHaplo]) {
            int maxBase = -1;
            double maxBaseProb = -1.0;
            for (int iBase = 0; iBase < 4; ++iBase) {
                if (this.probBase[expandHaplo][iSite3][iBase] > maxBaseProb) {
                    maxBaseProb = this.probBase[expandHaplo][iSite3][iBase];
                    maxBase = iBase;
                }
                expandedProbBase[expandHaplo][iSite3][iBase] = 0.2 * this.probBase[expandHaplo][iSite3][iBase];
                expandedProbBase[this.nHaplo][iSite3][iBase] = 0.8 * this.probBase[expandHaplo][iSite3][iBase];
            }
            expandedProbBase[expandHaplo][iSite3][maxBase] = 0.8 * this.probBase[expandHaplo][iSite3][maxBase];
            expandedProbBase[this.nHaplo][iSite3][maxBase] = 0.2 * this.probBase[expandHaplo][iSite3][maxBase];
            double summ = expandedProbBase[expandHaplo][iSite3][0] + expandedProbBase[expandHaplo][iSite3][1] + expandedProbBase[expandHaplo][iSite3][2] + expandedProbBase[expandHaplo][iSite3][3];
            final double[] array = expandedProbBase[expandHaplo][iSite3];
            final int n = 0;
            array[n] /= summ;
            final double[] array2 = expandedProbBase[expandHaplo][iSite3];
            final int n2 = 1;
            array2[n2] /= summ;
            final double[] array3 = expandedProbBase[expandHaplo][iSite3];
            final int n3 = 2;
            array3[n3] /= summ;
            final double[] array4 = expandedProbBase[expandHaplo][iSite3];
            final int n4 = 3;
            array4[n4] /= summ;
            summ = expandedProbBase[this.nHaplo][iSite3][0] + expandedProbBase[this.nHaplo][iSite3][1] + expandedProbBase[this.nHaplo][iSite3][2] + expandedProbBase[this.nHaplo][iSite3][3];
            final double[] array5 = expandedProbBase[this.nHaplo][iSite3];
            final int n5 = 0;
            array5[n5] /= summ;
            final double[] array6 = expandedProbBase[this.nHaplo][iSite3];
            final int n6 = 1;
            array6[n6] /= summ;
            final double[] array7 = expandedProbBase[this.nHaplo][iSite3];
            final int n7 = 2;
            array7[n7] /= summ;
            final double[] array8 = expandedProbBase[this.nHaplo][iSite3];
            final int n8 = 3;
            array8[n8] /= summ;
            avgFrac += maxBaseProb;
        }
        avgFrac /= this.polySites[expandHaplo].size();
        expandedHapFreq[expandHaplo] = avgFrac * this.hapFreq[expandHaplo];
        expandedHapFreq[this.nHaplo] = (1.0 - avgFrac) * this.hapFreq[expandHaplo];
        final int expandedNHaplo = this.nHaplo + 1;
        final Parameters newParameters = new Parameters(0.0, 0.0, expandedHapFreq, expandedProbBase, null, null);
        return newParameters;
    }
}