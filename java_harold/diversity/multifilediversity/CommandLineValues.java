// 
// Decompiled by Procyon v0.5.36
// 

package multifilediversity;

import java.io.IOException;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.io.File;

public class CommandLineValues
{
    private double fracSampled;
    private double inputErrorRate;
    private File dataFileList;
    private boolean storedValues;
    HashMap<String, double[]> storedValuesHash;
    
    CommandLineValues(final String[] args) {
        this.fracSampled = -1.0;
        this.inputErrorRate = -1.0;
        this.dataFileList = null;
        this.storedValues = false;
        this.storedValuesHash = new HashMap<String, double[]>();
        int iArg = 0;
        while (iArg < args.length) {
            if (args[iArg].startsWith("-")) {
                if (args[iArg].equals("-sampled")) {
                    this.fracSampled = Double.parseDouble(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-errorRate")) {
                    this.inputErrorRate = Double.parseDouble(args[iArg + 1]);
                    iArg += 2;
                }
                else if (args[iArg].equals("-dataFileList")) {
                    this.dataFileList = new File(args[iArg + 1]);
                    iArg += 2;
                }
                else {
                    if (!args[iArg].equals("-storedParameters")) {
                        continue;
                    }
                    this.storedValues = true;
                    final File storedValuesList = new File(args[iArg + 1]);
                    if (storedValuesList.exists()) {
                        try {
                            final FileReader file = new FileReader(storedValuesList);
                            final BufferedReader buff = new BufferedReader(file);
                            boolean eof = false;
                            while (!eof) {
                                final String line = buff.readLine();
                                if (line == null) {
                                    eof = true;
                                }
                                else {
                                    final String[] words = line.split("\\t");
                                    final double[] values = new double[5];
                                    for (int iVal = 0; iVal < 5; ++iVal) {
                                        values[iVal] = Double.parseDouble(words[iVal + 1]);
                                    }
                                    this.storedValuesHash.put(words[0], values);
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
                    else {
                        System.out.println(args[iArg + 1] + " not valid");
                        System.exit(1);
                    }
                    iArg += 2;
                }
            }
        }
        boolean notOK = false;
        if (this.dataFileList == null) {
            System.out.println("Must supply output data file list, \"-dataFileList\"");
            notOK = true;
        }
        else if (!this.dataFileList.exists()) {
            System.out.println("File " + this.dataFileList.getName() + " does not exist");
            notOK = true;
        }
        if (notOK) {
            System.exit(1);
        }
    }
    
    boolean getStoredValues() {
        return this.storedValues;
    }
    
    HashMap<String, double[]> getStoredValuesHash() {
        return this.storedValuesHash;
    }
    
    double getFracSampled() {
        return this.fracSampled;
    }
    
    double getErrorRate() {
        return this.inputErrorRate;
    }
    
    File getDataFileList() {
        return this.dataFileList;
    }
}
