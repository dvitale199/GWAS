#### THINGS TO ADD #####
# -CLEANUP DURING EACH STEP
# -FIX LOGGING TO SINGLE FILE
# -FIX AWFUL CONDITIONAL METHOD
# -FIGURE OUT WHERE PHENOS ARE FOR ADNI DATA AND MERGE
# -FIGURE OUT WHAT IS FAILING ON VARIANT-LEVEL
# -AUTOMATE IMPUTATION VIA API
# -GET PLINK AND HAIL GWAS RUNNING



import subprocess
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

args = parser.parse_args()

print("PROCESSING THE FOLLOWING GENOTYPES:", args.geno)

def logging(geno_path, out_path):
    out_name = geno_path + "_GWAS_QC.log"
    log = open(out_name, "a", newline='\n')
    
    return log
    
#QC and data cleaning
def het_pruning(geno_path, out_path):
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
    
    log = logging(geno_path, out_path)
    log.write("RUNNING QC FOR " + geno_path)
    log.write("\n")
    log.write("PRUNING FOR HETEROZYGOSITY WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("***********************************************")
    log.write("***********************************************")
    log.write("***********************************************")
    log.close()
        
        
def call_rate_pruning(geno_path, out_path):
    print("PRUNING FOR CALL RATE")
    
    bash1 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
    bash2 = "mv " + geno_path + "_after_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"
    
    cmds = [bash1, bash2]
    
    log = logging(geno_path, out_path)
    log.write("\n")
    log.write("PRUNING FOR CALL RATE WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("***********************************************")
    log.write("***********************************************")
    log.write("***********************************************")
    log.close()    
        

def sex_check(geno_path, out_path):
    print("CHECKING SEXES")
    bash1 = "plink --bfile " + geno_path + " --check-sex 0.25 0.75 --maf 0.05 --out " + out_path + "gender_check1"
    bash2 = "plink --bfile "+ geno_path + " --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out " + out_path + "gender_check2"
    bash3 = "grep 'PROBLEM' " + out_path + "gender_check1.sexcheck > " + out_path + "problems1.txt"
    bash4 = "grep 'PROBLEM' " + out_path + "gender_check2.sexcheck > " + out_path + "problems2.txt"
    bash5 = "cat " + out_path + "problems1.txt " + out_path + "problems2.txt > " + out_path + "GENDER_FAILURES.txt"
    bash6 = "cut -f 1,2 " + out_path + "GENDER_FAILURES.txt > " + out_path + "samples_to_remove.txt"
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "samples_to_remove.txt --make-bed --out " + geno_path + "_sex"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    log = logging(geno_path, out_path)
    log.write("\n")
    log.write("PRUNING FOR SEX WITH THE FOLLOWING COMMANDS:")
    log.write("\n")
    log.write("\n")
    
    for cmd in cmds:
        log.write(cmd)
        log.write("\n")
        subprocess.run(cmd, shell=True)
        
    log.write("***********************************************")
    log.write("***********************************************")
    log.write("***********************************************")
    log.close()
 






###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# MAY NEED TO ADD FUNCTION FOR ANCESTRY OUTLIERS (PCR FOR RELATEDNESS)     
        
        
def relatedness_pruning(geno_path, out_path):
    print("RELATEDNESS PRUNING")
    bash1 = "gcta --bfile " + geno_path + " --make-grm --out " + out_path + "GRM_matrix --autosome --maf 0.05" 
    bash2 = "gcta --grm-cutoff 0.125 --grm " + out_path + "GRM_matrix --out " + out_path + "GRM_matrix_0125 --make-grm"
    bash3 = "plink --bfile " + geno_path + " --keep " + out_path + "GRM_matrix_0125.grm.id --make-bed --out " + geno_path + "_relatedness"
    
    cmds = [bash1, bash2, bash3]
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

        
        
###########################################################################################################################################################################################################################################################################################################################################################################################################################################################################

##variant checks
def variant_pruning(geno_path, out_path):
    print("VARIANT-LEVEL PRUNING")
    # variant missingness
    "plink --bfile " + geno_path + " --make-bed --out " + geno_path + "_geno --geno 0.05"
    
    #missingness by case control (--test-missing), using P > 1E-4
    "plink --bfile " + geno_path + " --test-missing --out " + out_path + "missing_snps" 
    "awk '{if ($5 <= 0.0001) print $2 }' " + out_path + "missing_snps.missing > " + out_path + "missing_snps_1E4.txt"
    "plink --bfile " + geno_path + " --exclude " + out_path + "missing_snps_1E4.txt --make-bed --out " + geno_path + "_missing1"
    
    #missingness by haplotype (--test-mishap), using P > 1E-4
    "plink --bfile " + geno_path + "_missing1 --test-mishap --out " + geno_path + "missing_hap" 
    "awk '{if ($8 <= 0.0001) print $9 }' " + out_path + "missing_hap.missing.hap > " + out_path + "missing_haps_1E4.txt"
    "sed 's/|/\/g' " + out_path + "missing_haps_1E4.txt > " + out_path + "missing_haps_1E4_final.txt"
    "plink --bfile " + geno_path + " --exclude " + out_path + "missing_haps_1E4_final.txt --make-bed --out " +  geno_path + "_missing2"
    
    #HWE from controls only using P > 1E-4
    "plink --bfile " + geno_path + "_missing2 --filter-controls --hwe 1E-4 --write-snplist"
    "plink --bfile " + geno_path + "_missing2 --extract " + out_path + "plink.snplist --make-bed --out " + geno_path + "_variant"
    
    # OPTIONAL STEP: the following may not be used if you want to use specific rare variants, otherwise, rare variants will be removed here
    "plink --bfile " + geno_path + " --maf 0.01 --make-bed --out " + geno_path + "_MAF"



        

geno = args.geno
out = args.out
geno_het = geno + "_het"
geno_call_rate = geno_het + "_call_rate"
geno_sex = geno_call_rate + "_sex"
geno_relatedness = geno_sex + "_relatedness"
geno_variant = geno_relatedness + "_variant"




# NEED TO FIX THIS MONSTROSITY BELOW... SOMETHING LIKE:
# [LIST OF STEPS]
# for each step after het_pruning:
    #new elif statement

# het pruning
het_pruning(geno, out)

# check if geno_het output exists, if not, run on original geno data (there may not have been anything pruned in previous step and thus, no new file created)
if os.path.exists(geno_het + ".bim"):
    call_rate_pruning(geno_het, out)
else:
    call_rate_pruning(geno, out)

# check if geno_call_rate exists, if not, check geno_het, if not, use original
if os.path.exists(geno_call_rate + ".bim"):
    sex_check(geno_call_rate, out)
elif os.path.exists(geno_het + ".bim"):
    sex_check(geno_het, out)
else:
    sex_check(geno, out)

if os.path.exists(geno_sex + ".bim"):
    relatedness_pruning(geno_sex, out)  
elif os.path.exists(geno_call_rate + ".bim"):
    relatedness_pruning(geno_call_rate)
elif os.path.exists(geno_het + ".bim"):
    relatedness_pruning(geno_het, out)
else:
    relatedness_pruning(geno, out)

if os.path.exists(geno_relatedness + ".bim"):
    variant_pruning(geno_relatedness, out)
elif os.path.exists(geno_sex + ".bim"):
    variant_pruning(geno_sex, out)  
elif os.path.exists(geno_call_rate + ".bim"):
    variant_pruning(geno_call_rate)
elif os.path.exists(geno_het + ".bim"):
    variant_pruning(geno_het, out)
else:
    variant_pruning(geno, out)




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