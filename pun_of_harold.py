# PUN of HaROLD - Pretty Useful ruNtime of HaROLD
# I AM HaROLD - Intuitive AutoMation for HaROLD
# 
#####
# author:   Oscar Charles - github:ojcharles
# Nov 2019
# descr: automates the optimal way of reconstructing haplotypes using the (R Goldstein 2019, biorxiv) method.
# runs harold, produces formatted outputs, applies logic automatically
#####
## notes
#1  takes 01_input .bam files into strandcount
#2  runs haplotype reconstruction on all files, in a process chain from 2-n number of haplotypes. Each time piping old inputs into the new process
#3  detects when a local maximum number of haplotypes has been hit, and stops this process.
#4  Runs haplotype refinement on each file in parallel, using the results from the best number of haplotypes run in 2.

#instructions
runvars = {}
runvars['max-haplo'] = 7        # default:7 CMV
runvars['alpha-frac'] = 0.01    # default: 0.01 CMV     fraction of variable sites to use in the error calculation. higher = slower
runvars['command'] = "qsub"     # command requied to run a shell script on your computer, or submit a job on a cluster. e.g. $sh script.sh $qsub script.sh
runvars['conda_activate_loc'] = '/share/apps/anaconda/bin/activate'
runvars['conda_env'] = '/home/ocharles/.conda/envs/harold-oct19'
runvars['java'] = "java" # how to execute java on target machine

# functions
runvars['max-haplo'] = runvars['max-haplo'] + 1 # fixes inclusive
import numpy as np
import pandas as pd
import os
import re
import time
import csv
import matplotlib.pyplot as plt 
from glob import glob
from os import listdir
from os.path import isfile, join
from datetime import datetime
from shutil import copyfile
from distutils.dir_util import copy_tree
import warnings
warnings.filterwarnings("ignore")

def has_haplo_finished(haplo_num, alphafrac, path_curr):
    fin_file = glob(path_curr+"/ResultsHaplo.fasta") #list files called haplo sequences as last step 
    if len(fin_file) == 1:
        return(True)
    else:
        return(False)

def checkfilelist(path_proj):
    filelist = glob(path_proj+"/02_strandcount/*.txt")
    if len(filelist) ==1:
        return(True)
    else:
        print("Requires a single .txt file in input folder, defining the input files. has found 0 or more than 1!")
        return(False)

def formatHaploFreqOutput(path_curr):
    #hello
    '''takes harold output and returns formatted
    tab separated frequency file ready for harold'''
    infile = "/Results.log"
    outfile = "/freq_table.tsv"
    lines = []
    with open(path_curr + infile, "r") as f:
        #extract lines with frequency table
        for num, line in enumerate(f, 1):
            lines.append(line)
            if 'Main: Final total likelihood' in line:
                end = num
                LL = line
            if 'Haplotype frequencies' in line:
                start = num
    text = lines[start:end-2]
    with open(path_curr + outfile, "w") as w:
        for line in text:
            w.write(line)
    #dont write a file but return the i row table of LL
    LL_num = re.search("-[0-9]{1,20}.[0-9]{1,20}", LL)
    LL_out = float(LL_num.group())
    return(LL_out)

def create_strandcount_filelist():
    infile = glob(path_proj+"/01_bam/*.txt")[0]
    outfile = os.path.join(path_proj, "02_strandcount", "filelist.txt")
    with open(infile, "r") as f:
        #print(f.read())
        newText=f.read().replace('.bam', 'strandcount.csv')
        print(newText)
    with open(outfile, "w") as f:
        f.write(newText) 
        print('filelist script written')
        pass

def create_harold_script(numhaplo, alphafrac, path_curr, path_proj, runvars, refinebam, refinebestharold):
    #create shell script for execution  
    infile = path_proj + "/harold.sh" # from path up 1
    outfile = path_curr + "/init.sh"
    ref_file = glob(os.path.join(path_proj, "01_bam", "*.fasta"))[0]
    with open(infile, "r") as f:
        newText=f.read().replace('${dir_proj}', path_proj)
        newText = newText.replace('${v_haplonum}', str(numhaplo))
        newText = newText.replace('${v_alphafrac}', str(alphafrac))
        newText = newText.replace('${dir_harold}', os.path.join(path_proj, '03_dependencies'))
        newText = newText.replace('${dir_curr}', path_curr)
        newText = newText.replace('${infile}', ''.join(glob(path_curr +"/*.txt")))
        newText = newText.replace('${command}', runvars['command'])
        newText = newText.replace('${conda_activate_loc}', runvars['conda_activate_loc'])
        newText = newText.replace('${conda_env}', runvars['conda_env'])
        newText = newText.replace('${v_java}', runvars['java'])
        newText = newText.replace('${v_threadnum}', str(runvars['threads']))
        newText = newText.replace('${refine.bam}', str(refinebam))
        newText = newText.replace('${refine.best_harold_dir}', str(refinebestharold))
        newText = newText.replace('${ref.fasta}', str(ref_file))
        if(numhaplo <= 2):
            newText = newText.replace('${init_freq}', "")
        else:
            #create frequency file
            #formatHaploFreqOutput("Results.log", dir) # needs to be at end
            #read previous run frequency table
            newText = newText.replace('${init_freq}', "-f " +  path_proj+"/h"+str(numhaplo-1)+"_a"+str(alphafrac)+"/freq_table.tsv")
        pass
    with open(outfile, "w") as f:
        f.write(newText) 
        print('Shell script written')
        pass

def estimate_coef(x, y): 
    # number of observations/points 
    n = np.size(x) 
    # mean of x and y vector 
    m_x, m_y = np.mean(x), np.mean(y) 
    # calculating cross-deviation and deviation about x 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
    # calculating regression coefficients 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
    return(b_0, b_1) 

def plot_regression_line(x, y, b, filename): 
    # plotting the actual points as scatter plot 
    plt.scatter(x, y, color = "m", marker = "o", s = 30) 
    # predicted response vector 
    y_pred = b[0] + b[1]*x 
    # plotting the regression line 
    plt.plot(x, y_pred, color = "g") 
    # putting labels 
    plt.xlabel('x') 
    plt.ylabel('y') 
    # function to show plot 
    #plt.show() 
    plt.savefig('LogLikelihood.png')

def bam2strand(files_bam, path_proj):
    # this could / should just be one shell script 
    # as 1 bam file -> bam refined -> strandcount -> diversity - probably have a script already
    for i in range(len(files_bam)):
        file_bam = files_bam[i]
        file_bam_base = os.path.basename(os.path.splitext(file_bam)[0]) #base filename
        ## make bam file clc like
        file_clean = str(os.path.splitext(file_bam)[0])+"_ready"+str(os.path.splitext(file_bam)[1])
        command = "samtools view -h -b -G69 "+str(file_bam)+" | samtools view -h -b -G133 > "+str(file_clean)
        os.system(command)
        ## clc like bam file to strandcount
        # java program will outpuput filename.log and filename.strandcount.csv to path_project need to clean up
        path_curr = os.path.join(path_proj,"02_strandcount")
        command = "java -Xms4096m -Xmx8192m -cp 03_dependencies/HAROLD/Cluster_RG/dist/lib/htsjdk-unspecified-SNAPSHOT.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/picocli-3.6.0.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/pal-1.5.1.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/cache2k-all-1.0.2.Final.jar:03_dependencies/HAROLD/Cluster_RG/dist/lib/commons-math3-3.5.jar:03_dependencies/MakeReadCount.jar makereadcount.MakeReadCount " + str(file_clean)
        os.system(command)
        os.remove(os.path.join(path_proj , os.path.basename(os.path.splitext(file_clean)[0])+".log")) #dont need this
        os.rename(os.path.join(path_proj , os.path.basename(os.path.splitext(file_clean)[0])+".strandcount.csv"), os.path.join(path_curr , os.path.basename(os.path.splitext(file_clean)[0])+".strandcount.csv")) #strandcount.csv
        os.remove(file_clean) # remove very large bam file


# *****************************actual process ******************************************
print("*************** Punofharold Started ***************")
path_proj = os.getcwd()
r_alphafrac = str(runvars['alpha-frac'])

## bam 2 strand ##
files_bam = glob(os.path.join(path_proj, "01_bam", '*.bam'))
files_strand = glob(os.path.join(path_proj, "02_strandcount", "*.csv"))
if len(files_strand) == 0:
    print("Strandcount files not found, producing now.")
    bam2strand(files_bam, path_proj)

create_strandcount_filelist()

## haplotype reconstruction ##
print("** Haplotype Reconstruction **")
for r_num_haplo in range(2,7):#, runvars['max-haplo']):
    print("haplotypes - "+str(r_num_haplo)+" // alphafrac - "+r_alphafrac)
    # setup
    path_curr = path_proj+"/h"+str(r_num_haplo)+"_a"+r_alphafrac
    # create folder for run and copy files
    if os.path.exists(path_curr) == False:
        os.mkdir(path_curr)
    # copy files required for exection
    fromDirectory = path_proj+"/02_strandcount"
    toDirectory = path_curr
    copy_tree(fromDirectory, toDirectory)
    runvars['threads'] = len(glob(fromDirectory+"\*"))
    create_harold_script(r_num_haplo, r_alphafrac, path_curr, path_proj, runvars, " ", " ")

    # run shell script to start job
    command =  str(runvars['command']) + " " + path_curr + "/init.sh" # calls the shell file
    #os.system(command)
    
    r_time = 0
    while has_haplo_finished(r_num_haplo, r_alphafrac, path_curr) == False:
        #print("waiting for Harold to finish")
        r_wait = 5 # seconds
        time.sleep(5) # add to top as instruction?
        r_time = r_time + r_wait
        pass
    print("process took "+str(r_time/60)+" minutes to complete")
    #time.sleep(20) # allow any connections to close
    
    # clean up frequency crate dataframe - to query for best haplo
    LL = formatHaploFreqOutput(path_curr) #creates file & returns LL for run
    if r_num_haplo == 2:
        df_LL = pd.DataFrame([[r_num_haplo,LL]], columns=['num_haplo', 'LL'])
        pass
    else:
        df2 = pd.DataFrame([[r_num_haplo,LL]], columns=['num_haplo', 'LL'])
        df_LL = df_LL.append(df2)
        pass
    #print("LL and frequency data cleaned")

    # have we plateaued? if so go to refinement, else run with num_haplo ++
    if r_num_haplo >= 4: #oc hardcoded minimmum agreeable number of haplotypes to test
        #get last 3 points
        X = np.array(df_LL[["num_haplo"]])
        Y = np.array(df_LL[["LL"]])
        X = X[-3:]
        Y = Y[-3:]

        # estimating coefficients - get slope
        b = estimate_coef(X, Y) 
        slope = b[1]

        # have we found the best haplotype?
        if slope < 0: #then we have likely found a good maximum LL num haplo, pipe this into refinement
            refinebestharold = path_proj+"/h"+str(r_num_haplo)+"_a"+r_alphafrac
            print("---Maximum Likelihood Maxima has been reached ---")
            print("Linear Regression:\nintercept = {}  \nslope      = {}".format(b[0], b[1]))
            path_curr = path_curr = path_proj+"/refinement"
            if os.path.exists(path_curr) == False:
                os.mkdir(path_curr)
            out_best_haplo = int(df_LL.max()[0])
            # plotting regression line - saves png
            plot_regression_line(np.array(df_LL[["num_haplo"]]), np.array(df_LL[["LL"]]), b, 'LogLikelihood_all.png') #why two?
            plot_regression_line(X, Y, b, "LogLikelihood_best-decision.png")
            # save LL table
            df_LL.to_csv("LogLikelihood.csv", index = False)


            ## refinement ##
            # n files = n jobs submitted, punofharold will sent them in parallell then cease.
            print("** Haplotype Refinement started **")
            for i in range(len(files_bam)):
                file_bam = files_bam[i]
                print("fired off refinement script - i = "+str(i)+" && filename = "+str(file_bam))
                
                path_curr = os.path.join(path_proj, "refinement", "file-"+str(i))
                if os.path.exists(path_curr) == False:
                    os.mkdir(path_curr)

                create_harold_script(r_num_haplo, r_alphafrac, path_curr, path_proj, runvars, file_bam, refinebestharold)
                command =  str(runvars['command']) + " " + path_curr + "/init.sh" # calls the shell file

                #refinement could take 48hours + ,so at this point kill the python script 
                next
            print("all refinement files submitted in paralell")
            print("script completed at") # get time
            print("*************** Analysis completed ***************")





#def handle_refinement(path_curr):
#   
#    return(False)
    

    
# if __name__ == "__main__": main()# if being called as a python scipt, not debugging:



######################## en code wasteland #############################
# example raise expection
# raise Exception('x should not exceed 5. The value of x was: {}'.format(x))
