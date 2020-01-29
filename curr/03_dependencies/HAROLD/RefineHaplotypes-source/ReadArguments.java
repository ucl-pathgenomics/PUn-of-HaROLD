package refineHaplotypes;

import java.util.Arrays;
import java.io.File;

public class ReadArguments
{
    private String baseFreqFileName;
    private String hapAlignmentFileName;
    private String hapFreqFileName;
    private String readsFileName;
    private String refSeqFileName;
    private String tag;
    private File baseFreqFile;
    private File hapAlignmentFile;
    private File hapFreqFile;
    private File readsFile;
    private File refSeqFile;
    private String[] arguments;
    private double minReadDepth;
    private double minReads;
    private double errorRate;
    private boolean iterate;
    private int maxIter;
    private int maxHaplo;
    private int maxRecombine;
    private boolean expand;
    
    ReadArguments(final String[] args) {
        this.baseFreqFileName = "";
        this.hapAlignmentFileName = "";
        this.hapFreqFileName = "";
        this.readsFileName = "";
        this.refSeqFileName = "";
        this.tag = "";
        this.minReadDepth = -1.0;
        this.minReads = 20.0;
        this.errorRate = 0.002;
        this.iterate = false;
        this.maxIter = 10;
        this.maxHaplo = 10;
        this.maxRecombine = 20;
        this.expand = false;
        this.arguments = Arrays.copyOf(args, args.length);
        int iArg = 0;
        while (iArg < args.length) {
            if (args[iArg].startsWith("-")) {
                if (args[iArg].equals("-bam") || args[iArg].equals("-sam")) {
                    this.readsFileName = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-baseFreq")) {
                    this.baseFreqFileName = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-hapFreq")) {
                    this.hapFreqFileName = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-hapSeq")) {
                    this.hapAlignmentFileName = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-refSeq")) {
                    this.refSeqFileName = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-tag")) {
                    this.tag = args[iArg + 1];
                    iArg += 2;
                }
                else if (args[iArg].equals("-minReadDepth")) {
                    this.minReadDepth = Double.parseDouble(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-minReads")) {
                    this.minReads = Double.parseDouble(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-errorRate")) {
                    this.errorRate = Double.parseDouble(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-iterate")) {
                    this.iterate = true;
                    ++iArg;
                }
                else if (args[iArg].equals("-maxIterate")) {
                    this.maxIter = Integer.parseInt(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-maxHaplo")) {
                    this.maxHaplo = Integer.parseInt(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-Expand")) {
                    this.expand = true;
                    ++iArg;
                }
                else {
                    if (!args[iArg].equals("-maxRecombine")) {
                        continue;
                    }
                    this.maxRecombine = Integer.parseInt(args[iArg + 1]);
                    iArg += 2;
                }
            }
        }
        if (this.maxIter < 1) {
            this.iterate = false;
        }
        boolean notOK = false;
        if (this.readsFileName.length() == 0) {
            System.out.println("Must supply BAM or SAM file name, \"-bam\" or \"-sam\"");
            notOK = true;
        }
        if (this.baseFreqFileName.length() == 0) {
            System.out.println("Must supply base frequency file name, \"-baseFreq\"");
            notOK = true;
        }
        if (this.hapAlignmentFileName.length() == 0) {
            System.out.println("Must supply hap alignment file name, \"-hapSeq\"");
            notOK = true;
        }
        if (this.refSeqFileName.length() == 0) {
            System.out.println("Must supply reference sequence file name, \"-refSeq\"");
            notOK = true;
        }
        if (this.tag.length() == 0) {
            System.out.println("Must supply output file tag, \"-tag\"");
            notOK = true;
        }
        if (!notOK) {
            this.baseFreqFile = new File(this.baseFreqFileName);
            this.hapAlignmentFile = new File(this.hapAlignmentFileName);
            this.hapFreqFile = new File(this.hapFreqFileName);
            this.readsFile = new File(this.readsFileName);
            this.refSeqFile = new File(this.refSeqFileName);
            if (!this.baseFreqFile.exists()) {
                System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.baseFreqFileName));
                notOK = true;
            }
            if (!this.hapAlignmentFile.exists()) {
                System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.hapAlignmentFileName));
                notOK = true;
            }
            if (!this.readsFile.exists()) {
                System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.readsFileName));
                notOK = true;
            }
            if (!this.refSeqFile.exists()) {
                System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.refSeqFileName));
                notOK = true;
            }
        }
        if (notOK) {
            System.exit(1);
        }
    }
    
    public File getBaseFreqFile() {
        return this.baseFreqFile;
    }
    
    public File getHapAlignmentFile() {
        return this.hapAlignmentFile;
    }
    
    public File getHapFreqFile() {
        return this.hapFreqFile;
    }
    
    public File getReadsFile() {
        return this.readsFile;
    }
    
    public File getRefSeqFile() {
        return this.refSeqFile;
    }
    
    public String getTag() {
        return this.tag;
    }
    
    public double getMinReadDepth() {
        return this.minReadDepth;
    }
    
    public double getMinReads() {
        return this.minReads;
    }
    
    public double getErrorRate() {
        return this.errorRate;
    }
    
    public boolean getIterate() {
        return this.iterate;
    }
    
    public int getMaxIter() {
        return this.maxIter;
    }
    
    public int getMaxHaplo() {
        return this.maxHaplo;
    }
    
    public int getMaxRecombine() {
        return this.maxRecombine;
    }
    
    public boolean getExpand() {
        return this.expand;
    }
    
    public void printStuff() {
        System.out.print("Arguments:");
        RefineHaplotypes.logWriter.print("Arguments:");
        for (int iArg = 0; iArg < this.arguments.length; ++iArg) {
            System.out.print(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.arguments[iArg]));
            RefineHaplotypes.logWriter.print(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, this.arguments[iArg]));
        }
        System.out.println();
        RefineHaplotypes.logWriter.println();
        System.out.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.minReadDepth));
        System.out.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.minReads));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.minReadDepth));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.minReads));
        System.out.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.errorRate));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(D)Ljava/lang/String;, this.errorRate));
        System.out.println(invokedynamic(makeConcatWithConstants:(Z)Ljava/lang/String;, this.iterate));
        RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(Z)Ljava/lang/String;, this.iterate));
        if (this.iterate) {
            System.out.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxHaplo));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxHaplo));
            System.out.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxIter));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxIter));
            System.out.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxRecombine));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, this.maxRecombine));
            System.out.println(invokedynamic(makeConcatWithConstants:(Z)Ljava/lang/String;, this.expand));
            RefineHaplotypes.logWriter.println(invokedynamic(makeConcatWithConstants:(Z)Ljava/lang/String;, this.expand));
        }
        System.out.println();
        RefineHaplotypes.logWriter.println();
    }
}