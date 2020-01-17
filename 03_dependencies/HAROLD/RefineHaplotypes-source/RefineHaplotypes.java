package refineHaplotypes;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import java.io.File;
import java.util.Collection;
import java.io.IOException;
import java.io.BufferedReader;
import pal.datatype.DataType;
import java.io.Reader;
import pal.alignment.AlignmentReaders;
import java.io.FileReader;
import pal.datatype.Nucleotides;
import java.util.Arrays;
import java.util.Iterator;
import java.util.HashMap;
import org.apache.commons.math3.distribution.NormalDistribution;
import java.util.Random;
import java.io.PrintWriter;
import pal.alignment.Alignment;
import java.util.ArrayList;

public class RefineHaplotypes
{
    ReadArguments readArguments;
    int nSites;
    int nSamples;
    int nRegions;
    private boolean[][] basePresent;
    private int[] nBasePresent;
    private int[][][] siteCount;
    private int[] consensusSeq;
    private ArrayList<Integer> activeSites;
    Alignment alignment;
    Alignment refAlignment;
    ArrayList<Read> readList;
    static double errorRate;
    static byte minBaseQual;
    static int minMappingQual;
    double minHapFreq;
    double minEValue;
    double minReadDepth;
    double minReads;
    int nOptBaseMaxIter;
    double nParamsPerHaplo;
    boolean iterate;
    int maxIter;
    int maxHaplo;
    int maxRecombine;
    boolean expand;
    String outputFileNameTag;
    String sampleTag;
    static PrintWriter alignmentWriter;
    static PrintWriter logWriter;
    double smallValue;
    double bigValue;
    char[] dna;
    static Random random;
    static NormalDistribution nd;
    private static final HashMap<Character, Integer> DNA_CHAR_TO_INT_DNA_HASH;
    
    public static void main(final String[] args) {
        final RefineHaplotypes m = new RefineHaplotypes();
        m.run(args);
    }
    
    private RefineHaplotypes() {
        this.activeSites = new ArrayList<Integer>();
        this.readList = new ArrayList<Read>();
        this.minHapFreq = 0.0;
        this.minEValue = 0.001;
        this.minReadDepth = -1.0;
        this.minReads = 20.0;
        this.nOptBaseMaxIter = 1000;
        this.nParamsPerHaplo = 1.0;
        this.iterate = false;
        this.maxIter = -1;
        this.maxHaplo = 10;
        this.maxRecombine = 1000;
        this.expand = false;
        this.smallValue = 0.0;
        this.bigValue = 1.0;
        this.dna = new char[] { 'A', 'C', 'G', 'T', '-', 'X', 'a', 'c', 'g', 't', 'x' };
        this.fillHash();
    }
    
    private void run(final String[] args) {
        final Parameters startParameters = this.readData(args);
        this.process(startParameters);
        this.finish();
    }
    
    public void process(final Parameters startParameters) {
        Parameters bestParameters = this.optimiseBases(startParameters, 0, 0.0);
        this.printStuff(bestParameters, 0);
        System.out.println(invokedynamic(makeConcatWithConstants:(DDD)Ljava/lang/String;, bestParameters.totalLogLikelihood, bestParameters.totPenalty, bestParameters.adjLogLikelihood));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(DDD)Ljava/lang/String;, bestParameters.totalLogLikelihood, bestParameters.totPenalty, bestParameters.adjLogLikelihood));
        Parameters nextBest = bestParameters;
        if (this.maxRecombine > 0) {
            for (int iTry = 0; iTry < this.maxRecombine; ++iTry) {
                if (nextBest.nActive > 1) {
                    Parameters recombineParameters = nextBest.recombineHaplo(this.nSites, this.activeSites);
                    recombineParameters = this.optimiseBases(recombineParameters, 0, nextBest.adjLogLikelihood);
                    if (recombineParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                        nextBest = recombineParameters;
                    }
                    System.out.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, iTry, recombineParameters.adjLogLikelihood));
                    RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, iTry, recombineParameters.adjLogLikelihood));
                }
            }
            if (nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                bestParameters = nextBest;
            }
            System.out.println(invokedynamic(makeConcatWithConstants:(DDD)Ljava/lang/String;, bestParameters.totalLogLikelihood, bestParameters.totPenalty, bestParameters.adjLogLikelihood));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(DDD)Ljava/lang/String;, bestParameters.totalLogLikelihood, bestParameters.totPenalty, bestParameters.adjLogLikelihood));
        }
        if (this.iterate) {
            nextBest = bestParameters;
            int bigloop = 1;
            while (bigloop < this.maxIter) {
                if (bestParameters.nActive > 1) {
                    for (int iHaplo = 0; iHaplo < bestParameters.nHaplo - 1; ++iHaplo) {
                        for (int jHaplo = iHaplo + 1; jHaplo < bestParameters.nHaplo; ++jHaplo) {
                            if (bestParameters.hapFreq[iHaplo] > 1.0E-6 && bestParameters.hapFreq[jHaplo] > 1.0E-6) {
                                Parameters condenseParameters = bestParameters.compressHaplo(iHaplo, jHaplo, this.nSites);
                                condenseParameters = this.optimiseBases(condenseParameters, bigloop, nextBest.adjLogLikelihood);
                                if (condenseParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                    nextBest = condenseParameters;
                                }
                                System.out.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, bigloop, condenseParameters.adjLogLikelihood));
                                RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, bigloop, condenseParameters.adjLogLikelihood));
                            }
                        }
                    }
                }
                if (this.expand && bestParameters.nHaplo < this.maxHaplo) {
                    for (int iHaplo = 0; iHaplo < bestParameters.nHaplo; ++iHaplo) {
                        if (bestParameters.hapFreq[iHaplo] > 1.0E-6) {
                            Parameters expandParameters = bestParameters.expandHaplo(iHaplo, this.nSites);
                            expandParameters = this.optimiseBases(expandParameters, bigloop, nextBest.adjLogLikelihood);
                            System.out.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, iHaplo, expandParameters.adjLogLikelihood));
                            if (expandParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                nextBest = expandParameters;
                            }
                        }
                    }
                }
                if (bigloop == 1 || nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                    for (int iTry2 = 0; iTry2 < this.maxRecombine; ++iTry2) {
                        if (nextBest.nActive > 1) {
                            Parameters recombineParameters2 = nextBest.recombineHaplo(this.nSites, this.activeSites);
                            recombineParameters2 = this.optimiseBases(recombineParameters2, bigloop, nextBest.adjLogLikelihood);
                            if (recombineParameters2.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                nextBest = recombineParameters2;
                            }
                            System.out.println(invokedynamic(makeConcatWithConstants:(IID)Ljava/lang/String;, bigloop, iTry2, recombineParameters2.adjLogLikelihood));
                            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(IID)Ljava/lang/String;, bigloop, iTry2, recombineParameters2.adjLogLikelihood));
                        }
                    }
                }
                System.out.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, bigloop, nextBest.adjLogLikelihood));
                this.printStuff(nextBest, bigloop);
                if (nextBest.adjLogLikelihood - bestParameters.adjLogLikelihood < 1.0 || nextBest.nActive == 10) {
                    if (nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                        bestParameters = nextBest;
                        break;
                    }
                    break;
                }
                else {
                    bestParameters = nextBest;
                    ++bigloop;
                }
            }
        }
        final int[][] iMaxBase = bestParameters.iMaxBase;
        for (int iHaplo = 0; iHaplo < bestParameters.nHaplo - 1; ++iHaplo) {
            for (int jHaplo = iHaplo + 1; jHaplo < bestParameters.nHaplo; ++jHaplo) {
                if (bestParameters.hapFreq[iHaplo] > 1.0E-6 && bestParameters.hapFreq[jHaplo] > 1.0E-6) {
                    final double[] hamming = new double[2];
                    for (final int iSite : this.activeSites) {
                        if (iMaxBase[iHaplo][iSite] < 4 && iMaxBase[jHaplo][iSite] < 4) {
                            final double[] array = hamming;
                            final int n = 1;
                            ++array[n];
                            if (iMaxBase[iHaplo][iSite] == iMaxBase[jHaplo][iSite]) {
                                continue;
                            }
                            final double[] array2 = hamming;
                            final int n2 = 0;
                            ++array2[n2];
                        }
                    }
                    System.out.println(invokedynamic(makeConcatWithConstants:(IID)Ljava/lang/String;, iHaplo, jHaplo, hamming[0] / (1.0E-10 + hamming[1])));
                }
            }
        }
        RefineHaplotypes.logWriter.println("FINISHED");
        RefineHaplotypes.logWriter.println("\nHaplotype frequencies");
        for (int iHaplo = 0; iHaplo < bestParameters.nHaplo; ++iHaplo) {
            if (bestParameters.hapFreq[iHaplo] > 1.0E-10) {
                RefineHaplotypes.logWriter.format("%d\t%.6f\n", iHaplo, bestParameters.hapFreq[iHaplo]);
            }
        }
        RefineHaplotypes.logWriter.println("Active sites");
        for (final int iSite2 : this.activeSites) {
            RefineHaplotypes.logWriter.print(iSite2);
            for (int iHaplo2 = 0; iHaplo2 < bestParameters.nHaplo; ++iHaplo2) {
                if (bestParameters.hapFreq[iHaplo2] > 1.0E-10) {
                    RefineHaplotypes.logWriter.format("\t%.6f,%.6f,%.6f,%.6f", bestParameters.probBase[iHaplo2][iSite2][0], bestParameters.probBase[iHaplo2][iSite2][1], bestParameters.probBase[iHaplo2][iSite2][2], bestParameters.probBase[iHaplo2][iSite2][3]);
                }
            }
            RefineHaplotypes.logWriter.println();
        }
        RefineHaplotypes.alignmentWriter.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.refAlignment.getIdentifier(0).toString()));
        for (int iSite3 = 1; iSite3 < this.nSites; ++iSite3) {
            RefineHaplotypes.alignmentWriter.print(this.refAlignment.getData(0, iSite3 - 1));
        }
        RefineHaplotypes.alignmentWriter.println();
    }
    
    Parameters optimiseBases(final Parameters parameters, final int bigLoop, final double valToBeat) {
        Parameters bestParameters;
        Parameters nextParameters = bestParameters = this.assignToHaplotypes(parameters, 0, bigLoop);
        final double[] values = new double[this.nOptBaseMaxIter];
        double bestValue = nextParameters.adjLogLikelihood;
        values[0] = bestValue;
        for (int iIter = 1; iIter < this.nOptBaseMaxIter; ++iIter) {
            nextParameters = this.assignToHaplotypes(nextParameters, iIter, bigLoop);
            values[iIter] = Math.max(nextParameters.adjLogLikelihood, bestValue);
            if (iIter > 5 && Math.abs(nextParameters.adjLogLikelihood - bestValue) < 1.0) {
                bestParameters = nextParameters;
                break;
            }
            if (nextParameters.adjLogLikelihood > bestValue) {
                bestParameters = nextParameters;
                bestValue = bestParameters.adjLogLikelihood;
            }
            if (Math.abs(valToBeat) > 1.0E-9 && iIter > 5) {
                final double timeTillBest = (valToBeat - bestValue) / (0.2 * (values[iIter] - values[iIter - 5]));
                System.out.println(invokedynamic(makeConcatWithConstants:(DDD)Ljava/lang/String;, timeTillBest, valToBeat, bestValue));
                if (timeTillBest > 1000.0) {
                    break;
                }
            }
        }
        return bestParameters;
    }
    
    public Parameters assignToHaplotypes(final Parameters parameters, final int iTry, final int bigLoop) {
        final double[][][] probBase = parameters.probBase;
        double[] hapFreq = parameters.hapFreq;
        final int nHaplo = parameters.nHaplo;
        final HashMap<Read, double[]> logLikelihoodHash = new HashMap<Read, double[]>();
        for (final Read read : this.readList) {
            final double[] logLikelihood = new double[nHaplo];
            final HashMap<Integer, Integer> sigSiteHash = read.getSigSiteHash();
            for (int iHaplo = 0; iHaplo < nHaplo; ++iHaplo) {
                if (hapFreq[iHaplo] > 1.0E-10) {
                    for (final int iSite : sigSiteHash.keySet()) {
                        final double[] array = logLikelihood;
                        final int n = iHaplo;
                        array[n] += Math.log(probBase[iHaplo][iSite][sigSiteHash.get(iSite)]);
                    }
                }
            }
            logLikelihoodHash.put(read, logLikelihood);
        }
        double[] totalReads = new double[nHaplo];
        double previousLogLikelihood = -1.0E20;
        boolean finished = false;
        int iFreqIter = 0;
        final double[] newHapFreq = Arrays.copyOf(hapFreq, nHaplo);
        while (!finished) {
            double newTotalLogLikelihood = 0.0;
            totalReads = new double[nHaplo];
            for (final Read read2 : this.readList) {
                final double[] logLikelihood2 = logLikelihoodHash.get(read2);
                final HashMap<Integer, Integer> sigSiteHash2 = read2.getSigSiteHash();
                final double[] logProbHaplo = new double[nHaplo];
                double maxTerm = 0.0;
                double iMax = 0.0;
                maxTerm = -1.0E100;
                iMax = -1.0;
                for (int iHaplo2 = 0; iHaplo2 < nHaplo; ++iHaplo2) {
                    logProbHaplo[iHaplo2] = -1.0E10;
                    if (newHapFreq[iHaplo2] > 1.0E-10) {
                        logProbHaplo[iHaplo2] = Math.log(newHapFreq[iHaplo2]) + logLikelihood2[iHaplo2];
                        if (logProbHaplo[iHaplo2] > maxTerm) {
                            maxTerm = logProbHaplo[iHaplo2];
                            iMax = iHaplo2;
                        }
                    }
                }
                double summ = 0.0;
                final double[] probHaplo = new double[nHaplo];
                final double[] deltaLog = new double[nHaplo];
                final double sumProb = 0.0;
                for (int iHaplo3 = 0; iHaplo3 < nHaplo; ++iHaplo3) {
                    deltaLog[iHaplo3] = logProbHaplo[iHaplo3] - maxTerm;
                    if (deltaLog[iHaplo3] > -10.0) {
                        probHaplo[iHaplo3] = Math.exp(deltaLog[iHaplo3]);
                    }
                    summ += probHaplo[iHaplo3];
                }
                newTotalLogLikelihood += read2.getNCopies() * (maxTerm + Math.log(summ));
                for (int iHaplo3 = 0; iHaplo3 < nHaplo; ++iHaplo3) {
                    final double[] array2 = probHaplo;
                    final int n2 = iHaplo3;
                    array2[n2] /= summ;
                    final double[] array3 = totalReads;
                    final int n3 = iHaplo3;
                    array3[n3] += read2.getNCopies() * probHaplo[iHaplo3];
                }
            }
            double tot = 0.0;
            for (int iHaplo4 = 0; iHaplo4 < nHaplo; ++iHaplo4) {
                if ((this.minReadDepth < 0.0 || totalReads[iHaplo4] > this.minReadDepth * this.nSites / 250.0) && (this.minReads < 0.0 || totalReads[iHaplo4] > this.minReads)) {
                    tot += totalReads[iHaplo4];
                }
            }
            for (int iHaplo4 = 0; iHaplo4 < nHaplo; ++iHaplo4) {
                newHapFreq[iHaplo4] = 0.0;
                if ((this.minReadDepth < 0.0 || totalReads[iHaplo4] > this.minReadDepth * this.nSites / 250.0) && (this.minReads < 0.0 || totalReads[iHaplo4] > this.minReads)) {
                    newHapFreq[iHaplo4] = totalReads[iHaplo4] / tot;
                }
            }
            if (++iFreqIter > 100 || Math.abs(newTotalLogLikelihood - previousLogLikelihood) < 1.0) {
                finished = true;
            }
            else {
                previousLogLikelihood = newTotalLogLikelihood;
            }
        }
        hapFreq = newHapFreq;
        double totalLogLikelihood = 0.0;
        int nTotalAssigned = 0;
        final int[] nAssigned = new int[nHaplo];
        final int[] nPoly = new int[nHaplo];
        final double[][][] newProbBase = new double[nHaplo][this.nSites][4];
        final double[][][] reconstructSequence = new double[nHaplo][this.nSites][4];
        final ArrayList<Integer>[] polySites = (ArrayList<Integer>[])new ArrayList[nHaplo];
        final int[][] iMaxBase = new int[nHaplo][this.nSites];
        totalReads = new double[nHaplo];
        for (final Read read3 : this.readList) {
            final double[] logLikelihood3 = logLikelihoodHash.get(read3);
            final HashMap<Integer, Integer> sigSiteHash3 = read3.getSigSiteHash();
            final double[] logProbHaplo2 = new double[nHaplo];
            double maxTerm2 = 0.0;
            double iMax2 = 0.0;
            maxTerm2 = -1.0E100;
            iMax2 = -1.0;
            for (int iHaplo5 = 0; iHaplo5 < nHaplo; ++iHaplo5) {
                logProbHaplo2[iHaplo5] = -1.0E10;
                if (hapFreq[iHaplo5] > 1.0E-10) {
                    logProbHaplo2[iHaplo5] = Math.log(hapFreq[iHaplo5]) + logLikelihood3[iHaplo5];
                    if (logProbHaplo2[iHaplo5] > maxTerm2) {
                        maxTerm2 = logProbHaplo2[iHaplo5];
                        iMax2 = iHaplo5;
                    }
                }
            }
            double summ2 = 0.0;
            final double[] probHaplo2 = new double[nHaplo];
            final double[] deltaLog2 = new double[nHaplo];
            final double sumProb2 = 0.0;
            for (int iHaplo6 = 0; iHaplo6 < nHaplo; ++iHaplo6) {
                deltaLog2[iHaplo6] = logProbHaplo2[iHaplo6] - maxTerm2;
                if (deltaLog2[iHaplo6] > -10.0) {
                    probHaplo2[iHaplo6] = Math.exp(deltaLog2[iHaplo6]);
                }
                summ2 += probHaplo2[iHaplo6];
            }
            totalLogLikelihood += read3.getNCopies() * (maxTerm2 + Math.log(summ2));
            for (int iHaplo6 = 0; iHaplo6 < nHaplo; ++iHaplo6) {
                final double[] array4 = probHaplo2;
                final int n4 = iHaplo6;
                array4[n4] /= summ2;
                final double[] array5 = totalReads;
                final int n5 = iHaplo6;
                array5[n5] += read3.getNCopies() * probHaplo2[iHaplo6];
            }
            for (final int iSite2 : read3.getSigSiteHash().keySet()) {
                final int iSeq = read3.getSigSiteHash().get(iSite2);
                for (int iHaplo7 = 0; iHaplo7 < nHaplo; ++iHaplo7) {
                    final double[] array6 = reconstructSequence[iHaplo7][iSite2];
                    final int n6 = iSeq;
                    array6[n6] += probHaplo2[iHaplo7] * read3.getNCopies();
                }
            }
        }
        double tot2 = 0.0;
        for (int iHaplo2 = 0; iHaplo2 < nHaplo; ++iHaplo2) {
            if ((this.minReadDepth < 0.0 || totalReads[iHaplo2] > this.minReadDepth * this.nSites / 250.0) && (this.minReads < 0.0 || totalReads[iHaplo2] > this.minReads)) {
                tot2 += totalReads[iHaplo2];
            }
        }
        for (int iHaplo2 = 0; iHaplo2 < nHaplo; ++iHaplo2) {
            polySites[iHaplo2] = new ArrayList<Integer>();
            newHapFreq[iHaplo2] = 0.0;
            Arrays.fill(iMaxBase[iHaplo2], 4);
            for (final int iSite3 : this.activeSites) {
                Arrays.fill(newProbBase[iHaplo2][iSite3], 0.25);
            }
            if ((this.minReadDepth < 0.0 || totalReads[iHaplo2] > this.minReadDepth * this.nSites / 250.0) && (this.minReads < 0.0 || totalReads[iHaplo2] > this.minReads)) {
                newHapFreq[iHaplo2] = totalReads[iHaplo2] / tot2;
                for (final int iSite3 : this.activeSites) {
                    final double summ3 = reconstructSequence[iHaplo2][iSite3][0] + reconstructSequence[iHaplo2][iSite3][1] + reconstructSequence[iHaplo2][iSite3][2] + reconstructSequence[iHaplo2][iSite3][3];
                    if (summ3 > 0.99999) {
                        final double[] fracBase = new double[5];
                        for (int iBase = 0; iBase < 4; ++iBase) {
                            final double fracObsBase = reconstructSequence[iHaplo2][iSite3][iBase] / summ3;
                            if (fracObsBase > RefineHaplotypes.errorRate) {
                                fracBase[iBase] = (fracObsBase - RefineHaplotypes.errorRate) / (1.0 - 4.0 * RefineHaplotypes.errorRate);
                                final double[] array7 = fracBase;
                                final int n7 = 4;
                                array7[n7] += fracBase[iBase];
                            }
                        }
                        for (int iBase = 0; iBase < 4; ++iBase) {
                            newProbBase[iHaplo2][iSite3][iBase] = (1.0 - 3.0 * RefineHaplotypes.errorRate) * (fracBase[iBase] / fracBase[4]) + RefineHaplotypes.errorRate * (1.0 - fracBase[iBase] / fracBase[4]);
                            if (reconstructSequence[iHaplo2][iSite3][iBase] / summ3 > 0.6) {
                                ++nTotalAssigned;
                                final int[] array8 = nAssigned;
                                final int n8 = iHaplo2;
                                ++array8[n8];
                                iMaxBase[iHaplo2][iSite3] = iBase;
                            }
                        }
                        if (iMaxBase[iHaplo2][iSite3] != 4) {
                            continue;
                        }
                        final int[] array9 = nPoly;
                        final int n9 = iHaplo2;
                        ++array9[n9];
                        polySites[iHaplo2].add(iSite3);
                    }
                }
            }
        }
        final double totPenalty = nHaplo * this.nParamsPerHaplo;
        System.out.println(invokedynamic(makeConcatWithConstants:(IIDDD)Ljava/lang/String;, iTry, nTotalAssigned, totalLogLikelihood, totPenalty, totalLogLikelihood - totPenalty));
        for (int iHaplo8 = 0; iHaplo8 < nHaplo; ++iHaplo8) {
            if (hapFreq[iHaplo8] > 0.001) {
                System.out.println(invokedynamic(makeConcatWithConstants:(IDIIID)Ljava/lang/String;, iHaplo8, hapFreq[iHaplo8], nAssigned[iHaplo8], nPoly[iHaplo8], nAssigned[iHaplo8] + nPoly[iHaplo8], totalReads[iHaplo8]));
            }
        }
        final Parameters newParameters = new Parameters(totalLogLikelihood, totPenalty, newHapFreq, newProbBase, iMaxBase, polySites);
        return newParameters;
    }
    
    void printStuff(final Parameters inputParameters, final int iTry) {
        final int[][] iMaxBase = inputParameters.iMaxBase;
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(IDDD)Ljava/lang/String;, iTry, inputParameters.totalLogLikelihood, inputParameters.totPenalty, inputParameters.adjLogLikelihood));
        System.out.println(invokedynamic(makeConcatWithConstants:(IDDD)Ljava/lang/String;, iTry, inputParameters.totalLogLikelihood, inputParameters.totPenalty, inputParameters.adjLogLikelihood));
        final double[] hapFreq = inputParameters.hapFreq;
        for (int iHaplo = 0; iHaplo < inputParameters.nHaplo; ++iHaplo) {
            if (hapFreq[iHaplo] > 1.0E-6) {
                System.out.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, iHaplo, hapFreq[iHaplo]));
                RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(ID)Ljava/lang/String;, iHaplo, hapFreq[iHaplo]));
                RefineHaplotypes.alignmentWriter.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;II)Ljava/lang/String;, this.sampleTag, iHaplo, iTry));
                for (int iSite = 1; iSite < this.nSites; ++iSite) {
                    if (this.activeSites.contains(iSite)) {
                        RefineHaplotypes.alignmentWriter.print(this.dna[iMaxBase[iHaplo][iSite]]);
                    }
                    else {
                        RefineHaplotypes.alignmentWriter.print(this.dna[this.consensusSeq[iSite]]);
                    }
                }
                RefineHaplotypes.alignmentWriter.println();
            }
        }
        RefineHaplotypes.alignmentWriter.flush();
        System.out.println();
        RefineHaplotypes.logWriter.println();
    }
    
    public Parameters readData(final String[] args) {
        this.readArguments = new ReadArguments(args);
        int nHaplo = 0;
        double[] hapFreq = null;
        double[][][] probBase = null;
        this.outputFileNameTag = this.readArguments.getTag();
        this.minReadDepth = this.readArguments.getMinReadDepth();
        this.minReads = this.readArguments.getMinReads();
        RefineHaplotypes.errorRate = this.readArguments.getErrorRate();
        this.iterate = this.readArguments.getIterate();
        this.maxIter = this.readArguments.getMaxIter();
        this.maxHaplo = this.readArguments.getMaxHaplo();
        this.maxRecombine = this.readArguments.getMaxRecombine();
        this.expand = this.readArguments.getExpand();
        this.openFiles();
        this.readArguments.printStuff();
        try {
            final Nucleotides dataType = new Nucleotides();
            FileReader in = new FileReader(this.readArguments.getHapAlignmentFile());
            this.alignment = AlignmentReaders.readFastaSequences((Reader)in, (DataType)dataType);
            in.close();
            nHaplo = this.alignment.getSequenceCount();
            in = new FileReader(this.readArguments.getRefSeqFile());
            this.refAlignment = AlignmentReaders.readFastaSequences((Reader)in, (DataType)dataType);
            in.close();
            this.nSites = this.refAlignment.getSiteCount() + 1;
            this.nRegions = this.nSites / 1000;
            this.siteCount = new int[this.nSites][3][5];
            hapFreq = new double[nHaplo];
            probBase = new double[nHaplo][this.nSites][4];
            Arrays.fill(this.consensusSeq = new int[this.nSites], 4);
            FileReader file = null;
            BufferedReader buff = null;
            String line = "";
            String[] words = null;
            if (this.readArguments.getHapFreqFile().exists()) {
                file = new FileReader(this.readArguments.getHapFreqFile());
                buff = new BufferedReader(file);
                line = buff.readLine();
                words = line.split("\\t");
                for (int iHaplo = 0; iHaplo < nHaplo; ++iHaplo) {
                    hapFreq[iHaplo] = Double.parseDouble(words[iHaplo + 1]);
                }
                buff.close();
                file.close();
            }
            else {
                System.out.println("Initialising hapfreq with even values");
                for (int iHaplo = 0; iHaplo < nHaplo; ++iHaplo) {
                    hapFreq[iHaplo] = 1.0 / nHaplo;
                }
            }
            file = new FileReader(this.readArguments.getBaseFreqFile());
            buff = new BufferedReader(file);
            for (boolean eof = false; !eof; eof = true) {
                line = buff.readLine();
                if (line != null) {
                    words = line.split("\\t");
                    final int iSite = Integer.parseInt(words[0]);
                    for (int iHaplo2 = 0; iHaplo2 < nHaplo; ++iHaplo2) {
                        final String[] freqs = words[iHaplo2 + 4].split(",");
                        for (int iBase = 0; iBase < 4; ++iBase) {
                            final double rawProb = Double.parseDouble(freqs[iBase]);
                            probBase[iHaplo2][iSite][iBase] = rawProb * (1.0 - 3.0 * RefineHaplotypes.errorRate) + RefineHaplotypes.errorRate * (1.0 - rawProb);
                        }
                    }
                }
            }
            final File sampleFile = this.readArguments.getReadsFile();
            this.sampleTag = sampleFile.getName().replaceAll(".*\\/", "").replaceAll(".bam", "").replaceAll(".sam", "").replaceAll(".BAM", "").replaceAll(".SAM", "");
            this.readFromFile(sampleFile);
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
        this.basePresent = new boolean[this.nSites][4];
        this.nBasePresent = new int[this.nSites];
        for (int iSite2 = 0; iSite2 < this.nSites; ++iSite2) {
            final double[] pValues = this.checkStrandCount(this.siteCount[iSite2]);
            if (Double.isNaN(pValues[1]) || pValues[1] * this.nSites > this.minEValue) {
                int iMax = -1;
                int nMax = 0;
                int nTot = 0;
                for (int iBase2 = 0; iBase2 < 4; ++iBase2) {
                    nTot += this.siteCount[iSite2][2][iBase2];
                    if (this.siteCount[iSite2][2][iBase2] > nMax) {
                        iMax = iBase2;
                        nMax = this.siteCount[iSite2][2][iBase2];
                    }
                }
                if (nMax >= 1.999) {
                    this.basePresent[iSite2][iMax] = true;
                    this.consensusSeq[iSite2] = iMax;
                    final double noise = Math.max(1.999, RefineHaplotypes.errorRate * nMax);
                    for (int iBase3 = 0; iBase3 < 4; ++iBase3) {
                        if (iBase3 != iMax && this.siteCount[iSite2][2][iBase3] > noise) {
                            this.basePresent[iSite2][iBase3] = true;
                        }
                    }
                }
            }
        }
        for (int iSite2 = 0; iSite2 < this.nSites; ++iSite2) {
            for (int iBase4 = 0; iBase4 < 4; ++iBase4) {
                if (this.basePresent[iSite2][iBase4]) {
                    final int[] nBasePresent = this.nBasePresent;
                    final int n = iSite2;
                    ++nBasePresent[n];
                }
            }
            if (this.nBasePresent[iSite2] > 1) {
                this.nParamsPerHaplo += this.nBasePresent[iSite2] - 1.0;
                this.activeSites.add(iSite2);
            }
        }
        System.out.print(this.readList.size());
        final HashMap<String, Read> compressReadListHash = new HashMap<String, Read>();
        for (final Read read : this.readList) {
            read.setSignificantSites(this.activeSites);
            final String tag = read.getSigSiteTag();
            if (compressReadListHash.containsKey(tag)) {
                compressReadListHash.get(tag).addCopies();
            }
            else {
                if (tag.length() <= 0) {
                    continue;
                }
                compressReadListHash.put(tag, read);
            }
        }
        this.readList = new ArrayList<Read>(compressReadListHash.values());
        System.out.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.readList.size()));
        final Parameters parameters = new Parameters(0.0, 0.0, hapFreq, probBase, null, null);
        return parameters;
    }
    
    private void readFromFile(final File inputSamOrBamFile) throws IOException {
        final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
        for (final SAMRecord samRecord : reader) {
            if (samRecord.getMappingQuality() >= RefineHaplotypes.minMappingQual) {
                final Read newRead = new Read(samRecord);
                this.readList.add(newRead);
                final int[] limits = newRead.getLimits();
                final int[] sequence = newRead.getSequence();
                final boolean[] siteExists = newRead.getSiteExists();
                for (int iSite = limits[0]; iSite < limits[1]; ++iSite) {
                    final int lSite = iSite - limits[0];
                    if (siteExists[lSite]) {
                        final int[] array = this.siteCount[iSite][2];
                        final int n = sequence[lSite];
                        ++array[n];
                        if (newRead.getNegativeStrand()) {
                            final int[] array2 = this.siteCount[iSite][1];
                            final int n2 = sequence[lSite];
                            ++array2[n2];
                        }
                        else {
                            final int[] array3 = this.siteCount[iSite][0];
                            final int n3 = sequence[lSite];
                            ++array3[n3];
                        }
                    }
                }
            }
        }
    }
    
    double[] checkStrandCount(final int[][] sCount) {
        final double[] pValues = new double[2];
        final int[] iCount = { sCount[0][0] + sCount[0][1] + sCount[0][2] + sCount[0][3], sCount[1][0] + sCount[1][1] + sCount[1][2] + sCount[1][3] };
        final int total = iCount[0] + iCount[1];
        int iMin = 0;
        if (iCount[1] < iCount[0]) {
            iMin = 1;
        }
        final BinomialDistribution bd = new BinomialDistribution(total, 0.5);
        pValues[0] = bd.cumulativeProbability(iCount[iMin]);
        pValues[1] = 100.0;
        if (iCount[0] > 2 && iCount[1] > 2) {
            for (int iBase = 0; iBase < 4; ++iBase) {
                final int sumBase = sCount[0][iBase] + sCount[1][iBase];
                if (sumBase > 0 && sumBase < (total + 1) / 2) {
                    final double realTerm = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], sCount[0][iBase]) * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], sCount[1][iBase]);
                    final double[] p = new double[2];
                    for (int iTerm = 0; iTerm <= Math.min(sumBase, iCount[0]); ++iTerm) {
                        if (sumBase - iTerm <= iCount[1]) {
                            final double term = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], iTerm) * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], sumBase - iTerm);
                            if (term <= realTerm + 1.0E-5) {
                                final double[] array = p;
                                final int n = 0;
                                array[n] += term;
                            }
                            else {
                                final double[] array2 = p;
                                final int n2 = 1;
                                array2[n2] += term;
                            }
                        }
                    }
                    pValues[1] = Math.min(p[0] / (p[0] + p[1]), pValues[1]);
                }
            }
        }
        return pValues;
    }
    
    private void fillHash() {
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('A', 0);
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('C', 1);
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('G', 2);
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('T', 3);
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('-', 4);
        RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.put('N', 4);
    }
    
    public static int getDNA(final char DNA) {
        return RefineHaplotypes.DNA_CHAR_TO_INT_DNA_HASH.get(DNA);
    }
    
    public void openFiles() {
        try {
            final File logFile = new File(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.outputFileNameTag));
            if (logFile.exists()) {
                logFile.delete();
            }
            RefineHaplotypes.logWriter = new PrintWriter(logFile);
            final File alignmentOutputFile = new File(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.outputFileNameTag));
            if (alignmentOutputFile.exists()) {
                alignmentOutputFile.delete();
            }
            RefineHaplotypes.alignmentWriter = new PrintWriter(alignmentOutputFile);
        }
        catch (IOException ioe) {
            System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ioe.toString()));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ioe.toString()));
            System.exit(1);
        }
    }
    
    public void finish() {
        try {
            RefineHaplotypes.logWriter.close();
            RefineHaplotypes.alignmentWriter.close();
        }
        catch (Exception ex) {
            System.exit(1);
        }
    }
    
    static {
        RefineHaplotypes.errorRate = 0.002;
        RefineHaplotypes.minBaseQual = 30;
        RefineHaplotypes.minMappingQual = 10;
        RefineHaplotypes.alignmentWriter = null;
        RefineHaplotypes.logWriter = null;
        RefineHaplotypes.random = new Random();
        RefineHaplotypes.nd = new NormalDistribution(0.0, 10.0);
        DNA_CHAR_TO_INT_DNA_HASH = new HashMap<Character, Integer>();
    }
}