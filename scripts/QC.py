#### THINGS TO ADD #####
# -stdout/stderr output to log after each step
# -CLEANUP DURING EACH STEP
# -NEED TO CREATE MORE ORGANIZED DIRECTORY LAYOUT (OR JUST DELETE ALL EXTRANEOUS FILES)
# -FIX LOGGING TO SINGLE FILE
# -METHOD TO CHECK FILE PATH AT INPUT
# -COMMENTS!!!!!!!

#### SEPARATE SCRIPTS #####
# -AUTOMATE IMPUTATION VIA API
# -GET PLINK AND HAIL GWAS RUNNING


import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys
from contextlib import redirect_stderr, redirect_stdout

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')
parser.add_argument('--rare', default=False, action="store_true", help='Pruning toggle for rare variants. If --rare is used, final MAF pruning (0.01) will not be conducted, otherwise, rare variants will be pruned')


args = parser.parse_args()

geno = args.geno
out = args.out
rare = args.rare


# This will eventually need to be a class which checks for new log outputs in each step and appends contents between plink runs
def logger(geno_name):
    out = geno_name + "_GWAS_QC.log"
    if os.path.exists(out):
        # if log already exists, delete and then open new one in "append" mode
        bash_clear_log = "rm " + out
        subprocess.run(bash_clear_log, shell=True)
        log = open(out, "a", newline='\n')
    else:
        #if log does not exist, open new one in "append" mode
        log = open(out, "a", newline='\n')
        
    return log

# open log
log = logger(geno)

log.write("RUNNING QC FOR " + geno)
log.write("\n")
log.write("***********************************************")
log.write("\n")
log.write("***********************************************")
log.write("\n")
log.write("***********************************************")
log.write("\n")


def run_cmds(cmds_list, pruning_step, outpath, log=log):
        log.write(pruning_step + " WITH THE FOLLOWING COMMANDS:")
        log.write("\n")
        log.write("\n")
        
        
        for cmd in cmds_list:
            log.write(cmd)
            log.write("\n")
            
            # check length of list of files ending in ".log"
            logs_count = len(sorted(glob.glob(outpath + '*.log'), key=os.path.getmtime))
            subprocess.run(cmd, shell=True)
            
            # check length of new list of files ending in ".log". if longer than original, append newest .log file to running logfile
            # this indicates new log created
            new_logs_count = len(sorted(glob.glob(outpath + '*.log'), key=os.path.getmtime))
            
            if new_logs_count > logs_count:
                cmd_log = sorted(glob.glob(outpath + '*.log'), key=os.path.getmtime)[-1]
                
                new_log = open(cmd_log, "r")
                new_log_read = new_log.read()
                new_log.close()
                
                log.write(new_log_read)
                log.write("\n")
                log.write("***********************************************")
                
                
                
                    

        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")  

print("PROCESSING THE FOLLOWING GENOTYPES:", geno)

#QC and data cleaning
def het_pruning(geno_path, out_path):
#     print("NOW PRUNING FOR HETEROZYGOSITY")
    step = "NOW PRUNING FOR HETEROZYGOSITY"
    print(step)
    bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
    bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
    bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"
    bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
    bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
    bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_het"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    run_cmds(cmds, step, out_path)

    
def call_rate_pruning(geno_path, out_path):
    step = "PRUNING FOR CALL RATE"
    print(step)
    bash1 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
    bash2 = "mv " + geno_path + "_after_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"
    
    cmds = [bash1, bash2]
    
    run_cmds(cmds, step, out_path)


def sex_check(geno_path, out_path):
    step = "CHECKING SEXES"
    print(step)
    bash1 = "plink --bfile " + geno_path + " --check-sex 0.25 0.75 --maf 0.05 --out " + out_path + "gender_check1"
    bash2 = "plink --bfile "+ geno_path + " --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out " + out_path + "gender_check2"
    bash3 = "grep 'PROBLEM' " + out_path + "gender_check1.sexcheck > " + out_path + "problems1.txt"
    bash4 = "grep 'PROBLEM' " + out_path + "gender_check2.sexcheck > " + out_path + "problems2.txt"
    bash5 = "cat " + out_path + "problems1.txt " + out_path + "problems2.txt > " + out_path + "GENDER_FAILURES.txt"
    bash6 = "cut -f 1,2 " + out_path + "GENDER_FAILURES.txt > " + out_path + "samples_to_remove.txt"
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "samples_to_remove.txt --make-bed --out " + geno_path + "_sex"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    run_cmds(cmds, step, out_path)
    
    
# MAY NEED TO ADD FUNCTION FOR ANCESTRY OUTLIERS (PCR FOR RELATEDNESS)     
         
def relatedness_pruning(geno_path, out_path):
    step = "RELATEDNESS PRUNING"
    print(step)
    bash1 = "gcta --bfile " + geno_path + " --make-grm --out " + out_path + "GRM_matrix --autosome --maf 0.05" 
    bash2 = "gcta --grm-cutoff 0.125 --grm " + out_path + "GRM_matrix --out " + out_path + "GRM_matrix_0125 --make-grm"
    bash3 = "plink --bfile " + geno_path + " --keep " + out_path + "GRM_matrix_0125.grm.id --make-bed --out " + geno_path + "_relatedness"
    
    cmds = [bash1, bash2, bash3]
    
    run_cmds(cmds, step, out_path)
       
        
##variant checks
def variant_pruning(geno_path, out_path):
    step = "VARIANT-LEVEL PRUNING"
    print(step)
    # variant missingness
    bash1 = "plink --bfile " + geno_path + " --make-bed --out " + geno_path + "_geno --geno 0.05"
    
    #missingness by case control (--test-missing), using P > 1E-4
    bash2 = "plink --bfile " + geno_path + "_geno --test-missing --out " + out_path + "missing_snps" 
    bash3 = "awk '{if ($5 <= 0.0001) print $2 }' " + out_path + "missing_snps.missing > " + out_path + "missing_snps_1E4.txt"
    bash4 = "plink --bfile " + geno_path + "_geno --exclude " + out_path + "missing_snps_1E4.txt --make-bed --out " + geno_path + "_missing1"
    
    #missingness by haplotype (--test-mishap), using P > 1E-4
    bash5 = "plink --bfile " + geno_path + "_missing1 --test-mishap --out " + geno_path + "_missing_hap" 
    bash6 = "awk '{if ($8 <= 0.0001) print $9 }' " + out_path + "missing_hap.missing.hap > " + out_path + "missing_haps_1E4.txt"
    bash7 = "sed 's/|/\/g' " + out_path + "missing_haps_1E4.txt > " + out_path + "missing_haps_1E4_final.txt"
    bash8 = "plink --bfile " + geno_path + " --exclude " + out_path + "missing_haps_1E4_final.txt --make-bed --out " +  geno_path + "_missing2"
    
    
    ###### THIS DOES NOT WORK WITHOUT PHENOTYPES!!!!!!!!
    #HWE from controls only using P > 1E-4
    bash9 = "plink --bfile " + geno_path + "_missing2 --filter-controls --hwe 1E-4 --write-snplist --out " + out_path + "plink"
    bash10 = "plink --bfile " + geno_path + "_missing2 --extract " + out_path + "plink.snplist --make-bed --out " + geno_path + "_HWE"

    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7, bash8, bash9, bash10]
    
    run_cmds(cmds, step, out_path)
    
    exts = [".bed",".bim",".fam",".log",".hh"]
    for ext in exts:
        shutil.move(geno_path + "_HWE" + ext, geno_path + "_variant" + ext)

    
def rare_prune(geno_path, out_path, log=log, rare=rare):
    # OPTIONAL STEP: if --rare flag included when running, rare variants will be left alone, otherwise they will be pruned with --maf 0.01
    bash = "plink --bfile " + geno_path + " --maf 0.01 --make-bed --out " + geno_path + "_MAF"
    
    if rare:
        print("SKIPPING FINAL MAF PRUNING (0.01)... RARE VARIANTS LEFT ALONE")
        log.write("SKIPPING FINAL MAF PRUNING (0.01)... RARE VARIANTS LEFT ALONE")
        log.write("\n")
        log.write("\n")

        exts = [".bed",".bim",".fam",".log",".hh"]
        for ext in exts:
            shutil.move(geno_path + ext, geno_path + "_final" + ext)

        log.write("MOVED " + geno_path + " to " + geno_path + "_final")
        log.write("\n")    

    else:
        print("RARE VARIANTS (MAF <= 0.O1) PRUNING")
        log.write("RARE VARIANTS (MAF <= 0.O1) PRUNING WITH THE FOLLOWING COMMANDS:")
        log.write("\n")
        log.write("\n")
        subprocess.run(bash, shell=True)
        log.write(bash)
        log.write("\n")

        exts = [".bed",".bim",".fam",".log",".hh"]
        for ext in exts:
            shutil.move(geno_path + "_MAF" + ext, geno_path + "_final" + ext)
            
        log.write("MOVED " + geno_path + "_MAF to " + geno_path + "_final")
        log.write("\n")

    
# create new names for each step
geno_het = geno + "_het"
geno_call_rate = geno_het + "_call_rate"
geno_sex = geno_call_rate + "_sex"
geno_relatedness = geno_sex + "_relatedness"
geno_variant = geno_relatedness + "_variant"
geno_final = geno_variant + "_final"

# run pruning steps
het_pruning(geno, out)
call_rate_pruning(geno_het, out)
sex_check(geno_call_rate, out)
relatedness_pruning(geno_sex, out)
variant_pruning(geno_relatedness, out)
rare_prune(geno_variant, out)


log.close()




###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#more code to be added
"""

# post-imputation
plink2 --double-id --vcf chr1.dose.vcf.gz --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg38.txt --make-pgen --out chr1_dose_imputed

######## RUN plink_cleanup.py SCRIPT HERE ########

#extract phenos from .fam file
cat /data/CARD/PD/genotype_data/DUTCH/Dutch_gwas.fam | awk '{ print $1,$2,$6 }'

plink2 --pfile chr1_dose_imputed_fixedIDs_sex --indep-pairwise 1000 10 0.02 --autosome --pheno dutch_phenos.txt --out pruned_data

"""