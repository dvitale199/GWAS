#### THINGS TO ADD #####
# -stdout/stderr output to log after each step
# -CLEANUP DURING EACH STEP
# -NEED TO CREATE MORE ORGANIZED DIRECTORY LAYOUT
# -FIX LOGGING TO SINGLE FILE
# -Make toggle for rare variants
# -Need method to check file path

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

args = parser.parse_args()

geno = args.geno
out = args.out



def logging(geno_name, out_name):
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


#### MAY WANT TO THINK ABOUT CHANGING INITIAL LOGGING! THIS MAY BE CONFUSING USING "log" BOTH WITHIN AND OUTSIDE OF FUNCTIONS
log = logging(geno, out)
log.write("RUNNING QC FOR " + geno)
log.write("\n")
log.write("***********************************************")
log.write("\n")
log.write("***********************************************")
log.write("\n")
log.write("***********************************************")
log.write("\n")


print("PROCESSING THE FOLLOWING GENOTYPES:", geno)

#QC and data cleaning
def het_pruning(geno_path, out_path, log=log):
    print("NOW PRUNING FOR HETEROZYGOSITY")
    
    bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
    bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
    bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"

    # bash for now. convert these to python
    bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
    bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
    bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
    
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_het"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    log.write("PRUNING FOR HETEROZYGOSITY WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")  
        subprocess.run(cmd, shell=True)
        
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")    
        
def call_rate_pruning(geno_path, out_path, log=log):
    print("PRUNING FOR CALL RATE")
    
    bash1 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
    bash2 = "mv " + geno_path + "_after_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"
    
    cmds = [bash1, bash2]
    
    log.write("\n")
    log.write("PRUNING FOR CALL RATE WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n") 
        

def sex_check(geno_path, out_path, log=log):
    print("CHECKING SEXES")
    bash1 = "plink --bfile " + geno_path + " --check-sex 0.25 0.75 --maf 0.05 --out " + out_path + "gender_check1"
    bash2 = "plink --bfile "+ geno_path + " --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out " + out_path + "gender_check2"
    bash3 = "grep 'PROBLEM' " + out_path + "gender_check1.sexcheck > " + out_path + "problems1.txt"
    bash4 = "grep 'PROBLEM' " + out_path + "gender_check2.sexcheck > " + out_path + "problems2.txt"
    bash5 = "cat " + out_path + "problems1.txt " + out_path + "problems2.txt > " + out_path + "GENDER_FAILURES.txt"
    bash6 = "cut -f 1,2 " + out_path + "GENDER_FAILURES.txt > " + out_path + "samples_to_remove.txt"
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "samples_to_remove.txt --make-bed --out " + geno_path + "_sex"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    log.write("\n")
    log.write("PRUNING FOR SEX WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")

###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# MAY NEED TO ADD FUNCTION FOR ANCESTRY OUTLIERS (PCR FOR RELATEDNESS)     
         
def relatedness_pruning(geno_path, out_path, log=log):
    print("RELATEDNESS PRUNING")
    bash1 = "gcta --bfile " + geno_path + " --make-grm --out " + out_path + "GRM_matrix --autosome --maf 0.05" 
    bash2 = "gcta --grm-cutoff 0.125 --grm " + out_path + "GRM_matrix --out " + out_path + "GRM_matrix_0125 --make-grm"
    bash3 = "plink --bfile " + geno_path + " --keep " + out_path + "GRM_matrix_0125.grm.id --make-bed --out " + geno_path + "_relatedness"
    
    cmds = [bash1, bash2, bash3]
    
    
    log.write("\n")
    log.write("PRUNING FOR RELATEDNESS WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
        
        
###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################

##variant checks
def variant_pruning(geno_path, out_path, rare=False, log=log):
    print("VARIANT-LEVEL PRUNING")
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
    
    
    # OPTIONAL STEP: the following may not be used if you want to use specific rare variants, otherwise, rare variants will be removed here
    bash11 = "plink --bfile " + geno_path + "_HWE --maf 0.01 --make-bed --out " + geno_path + "_HWE_MAF"

    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7, bash8, bash9, bash10]
    
    log.write("\n")
    log.write("VARIANT-LEVEL PRUNING WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)

    if rare:
        print("rare")
        log.write("SKIPPING FINAL MAF PRUNING (0.01)... RARE VARIANTS LEFT ALONE")
#         bash_mv = "for f in "  + geno_path + "_HWE*; do mv $f " + geno_path + "_variant; done"
        
        exts = [".bed",".bim",".fam"]
        for ext in exts:
            shutil.move(geno_path + "_HWE" + ext, geno_path + "_variant" + ext)
        
        
#         subprocess.run(bash_mv, shell=True)
        log.write("MOVED " + geno_path + "_HWE to " + geno_path + "_variant")
    
    else:
        log.write("PRUNING RARE VARIANTS (MAF 0.O1):")
        subprocess.run(bash11, shell=True)
        log.write(bash11)
        
        exts = [".bed",".bim",".fam",".log",".hh"]
        for ext in exts:
            shutil.move(geno_path + "_HWE_MAF" + ext, geno_path + "_variant" + ext)
#         bash_mv = "mv " + geno_path + "_MAF " + geno_path + "_variant"
#         subprocess.run(bash_mv, shell=True)
        log.write("MOVED " + geno_path + "_HWE_MAF to " + geno_path + "_variant")
        
        
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")
    log.write("***********************************************")
    log.write("\n")

geno_het = geno + "_het"
geno_call_rate = geno_het + "_call_rate"
geno_sex = geno_call_rate + "_sex"
geno_relatedness = geno_sex + "_relatedness"
geno_variant = geno_relatedness + "_variant"

# run pruning steps
het_pruning(geno, out)
call_rate_pruning(geno_het, out)
sex_check(geno_call_rate, out)
relatedness_pruning(geno_sex, out)
variant_pruning(geno_relatedness, out)




log.close()





# check if geno_het output exists, if not, run on original geno data (there may not have been anything pruned in previous step and thus, no new file created)
# if os.path.exists(geno_het + ".bim"):
#     call_rate_pruning(geno_het, out)
# else:
#     call_rate_pruning(geno, out)

# # check if geno_call_rate exists, if not, check geno_het, if not, use original
# if os.path.exists(geno_call_rate + ".bim"):
#     sex_check(geno_call_rate, out)
# elif os.path.exists(geno_het + ".bim"):
#     sex_check(geno_het, out)
# else:
#     sex_check(geno, out)

# if os.path.exists(geno_sex + ".bim"):
#     relatedness_pruning(geno_sex, out)  
# elif os.path.exists(geno_call_rate + ".bim"):
#     relatedness_pruning(geno_call_rate)
# elif os.path.exists(geno_het + ".bim"):
#     relatedness_pruning(geno_het, out)
# else:
#     relatedness_pruning(geno, out)

# if os.path.exists(geno_relatedness + ".bim"):
#     variant_pruning(geno_relatedness, out)
# elif os.path.exists(geno_sex + ".bim"):
#     variant_pruning(geno_sex, out)  
# elif os.path.exists(geno_call_rate + ".bim"):
#     variant_pruning(geno_call_rate)
# elif os.path.exists(geno_het + ".bim"):
#     variant_pruning(geno_het, out)
# else:
#     variant_pruning(geno, out)




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