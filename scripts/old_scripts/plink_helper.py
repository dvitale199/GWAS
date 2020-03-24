import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys
from contextlib import redirect_stderr, redirect_stdout

# parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
# parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

# args = parser.parse_args()

class Logger:
    def __init__(self, geno_path, out_path):
        self.geno_path = geno_path
        self.out_path = out_path
        print("PROCESSING THE FOLLOWING GENOTYPES:", geno_path)
        
    def logging(self):

        log_out_name = self.geno_path + "_GWAS_QC.log"
        log = open(log_out_name, "a", newline='\n')
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("RUNNING PLINK COMMANDS FOR " + self.geno_path)
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
            
        return log
    
    def rm_log(self):
        log_out_name = self.geno_path + "_GWAS_QC.log"
        if os.path.exists(log_out_name):
            os.remove(log_out_name)
            
            
class Driver(Logger):
    def __init__(self, geno_path, out_path):
        super().__init__(geno_path, out_path)
        
    def run_cmds(self, cmds_list, step):
        log = self.logging()
        log.write(step + " WITH THE FOLLOWING COMMANDS:")
        log.write("\n")
        log.write("\n")
        
        
        for cmd in cmds_list:
            log.write(cmd)
            log.write("\n")
            log.write("\n")
            
            # check length of list of files ending in ".log"
            logs_count = len(sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime))
            subprocess.run(cmd, shell=True)
            
            # check length of new list of files ending in ".log". if longer than original, append newest .log file to running logfile
            # this indicates new log created
            new_logs_count = len(sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime))
            
            if new_logs_count > logs_count:
                cmd_log = sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime)[-1]
                
                new_log = open(cmd_log, "r")
                new_log_read = new_log.read()
                new_log.close()
                
                log.write(new_log_read)
                log.write("\n")
                log.write("***********************************************")
                log.write("\n")
                log.write("\n")

        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        
        
        
#QC and data cleaning
class QC(Driver):
    def __init__(self, geno_path, out_path):
        super().__init__(geno_path, out_path)
        # create new names for each step
        self.geno_het = geno_path + "_het"
        self.geno_call_rate = self.geno_het + "_call_rate"
        self.geno_sex = self.geno_call_rate + "_sex"
        self.geno_relatedness = self.geno_sex + "_relatedness"
        self.geno_variant = self.geno_relatedness + "_variant"
        self.geno_final = self.geno_variant + "_final"
    
    def het_pruning(self):
        geno_path = self.geno_path
        out_path = self.out_path
        
        step = "NOW PRUNING FOR HETEROZYGOSITY"

        bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
        bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
        bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"
        bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
        bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
        bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
        bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_het"

        cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
        
        self.run_cmds(cmds, step)
        
    def call_rate_pruning(self):
        geno_path = self.geno_het
        out_path = self.out_path
        step = "PRUNING FOR CALL RATE"
        print(step)
        bash1 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
        bash2 = "mv " + geno_path + "_after_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"

        cmds = [bash1, bash2]

        self.run_cmds(cmds, step)