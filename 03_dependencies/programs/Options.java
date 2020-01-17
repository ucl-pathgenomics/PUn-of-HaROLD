package cluster_rg;

//oc options for multiple steps in the process

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;

@Command(name = "richards-haplotype-model", footer = "Copyright (c) 2018 Richard A Goldstein", description = "", version = "1.0")
# the 
public class Options {
    
    @Option(names = {"-a", "--initial-alpha"}, arity = "2", description = "") // initial paramaters for the alpha distribution
    double[] initialAlphaParams = new double[]{Constants.DEFAULT_ALPHA_0, Constants.DEFAULT_ALPHA_1}; 
    
    @Option(names = {"-A", "--fix-alpha"}, required = false, description = "Fix alpha parameters")
    boolean fixAlpha = false;    

    @Option(names = {"--alpha-frac"}, required = false, description = "Fraction of sites to use to optimise error parameters") //between 1 (100%) or 0 0% for all bases to read for calculations, i.e. cmv use 10% , noro use 100% as small
    double alpha_frac = 1.0;
    
    // TODO: -c and -n should be n-arity inputs - make arrays
    @Option(names = {"-c", "--count-file"}, arity = "1..*", required = true, description = "file containing list of count files") //file containinig the names of the count files, .standcount.csv
    File[] countFile;
    
    /// OC - opportunity to optimise, at multiple points this is anoother.

    // TODO: -c and -n should be n-arity inputs - make arrays
    @Option(names = {"-f", "--initial-freq-file"}, arity = "1..*", description = "optional file containing hap frequency values") //frequuency of the different haplotypes for the different samples
    File[] initialFreqFile = null;  

    @Option(names = {"-g", "--gamma-cache"}, description = "", hidden = true) // a lot of time is spent otpimisng gamma funcitons, cache of gamma functionsthat have been evaluated prviisuly. larger the size the more efficient, as can store more gamma parameters
    int gammaCache = 1000; //10000

    @Option(names = {"-h", "-?", "--help"}, usageHelp = true, description = "give this help list")
    protected boolean helpRequested;

    @Option(names = {"-H", "--printHaplotypes"}, required = false, description = "Print haplotypes") // used for haplotype refinement
    boolean printHaplotypes = false;    

    @Option(names = {"-L", "--printLikelihoods"}, required = false, description = "Print likelihoods")
    boolean printLikelihoods = false;        

    @Option(names = {"-n", "--haplotypes"}, arity = "1..*", required = true, description = "number of haplotypes") //num of haplotypes to give
    int[] haplotypes;

    @Option(names = {"-N", "--noOpt"}, required = false, description = "Process without optimising") // don't optimise any of the features, all features. Faster more accurate to 
    boolean process = false;   

    @Option(names = {"-o", "--optimiser"}, required = false, description = "Optimiser for haplotype frequencies") //choose an optimier from a list of opportunities, 
    String optimiser = "BOBYQ";      

    @Option(names = {"-p", "--prefix"}, required = false, description = "Results file prefix") //just a name
    String prefix = "Results";

    @Option(names = {"-s", "--seed"}, description = "") //random seed
    long randomSeed = System.currentTimeMillis();

    @Option(names = {"--threads"}) //as many threads as you have samples i.e per patient. error parameters ar the same for all smaples.
    int threads = 1;

    @Option(names = {"--tol"})
    double tol = Constants.DEFAULT_TOL;
               
    @Option(names = {"-v", "--verbose"}, description = "")
    boolean verbose = false;

    @Option(names = {"-V", "--version"}, versionHelp = true, description = "")
    boolean versionRequested;

}

