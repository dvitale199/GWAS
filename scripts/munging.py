import subprocess
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

args = parser.parse_args()

print("GENOTYPE FILES", args.geno)

#QC and data cleaning
def het_pruning(geno_path, out_path):

    bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
    bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
    bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"

    # bash for now. convert these to python
    bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
    bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
    bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
    
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_het"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    for cmd in cmds:
        subprocess.run(cmd, shell=True)
        
def call_rate_pruning(geno_path, out_path):
    
    bash1 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
    bash2 = "mv " + geno_path + "_after_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"
    
    cmds = [bash1, bash2]
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

def sex_check(geno_path, out_path):
    bash1 = "plink --bfile " + geno_path + " --check-sex 0.25 0.75 --maf 0.05 --out " + out_path + "gender_check1"
    bash2 = "plink --bfile "+ geno_path + "--chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out " + out_path + "gender_check2"
    bash3 = "grep 'PROBLEM' " + out_path + "gender_check1.sexcheck > " + out_path + "problems1.txt"
    bash4 = "grep 'PROBLEM' " + out_path + "gender_check2.sexcheck > " + out_path + "problems2.txt"
    bash5 = "cat " + out_path + "problems1.txt " + out_path + "problems2.txt > " + out_path + "GENDER_FAILURES.txt"
    bash6 = "cut -f 1,2 " + out_path + "GENDER_FAILURES.txt > " + out_path + "samples_to_remove.txt"
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "samples_to_remove.txt --make-bed --out " + geno_path + "_sex"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

 ################################################################################################################################################################################################################################################################################################################################################
# MAY NEED TO ADD FUNCTION FOR ANCESTRY OUTLIERS (PCR FOR RELATEDNESS)     
        
        
def relatedness_pruning(geno_path, out_path):
    bash1 = "gcta --bfile " + geno_path + " --make-grm --out " + out_path + "GRM_matrix --autosome --maf 0.05" 
    bash2 = "gcta --grm-cutoff 0.125 --grm " + out_path + "GRM_matrix --out " + out_path + "GRM_matrix_0125 --make-grm"
    bash3 = "plink --bfile " + geno_path + " --keep " + out_path + "GRM_matrix_0125.grm.id --make-bed --out " + geno_path + "_relatedness"
    
    cmds = [bash1, bash2, bash3]
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

        
        
################################################################################################################################################################################################################################################################################################################################################

##variant checks
def variant_missingness(geno_path, out_path):
    # variant missingness
    "plink --bfile " + geno_path + " --make-bed --out " + geno_path + "_geno --geno 0.05"
    
    #missingness by case control (--test-missing), using P > 1E-4
    "plink --bfile " + geno_path + " --test-missing --out " + out_path + "missing_snps" 
    "awk '{if ($5 <= 0.0001) print $2 }'" + out_path + "missing_snps.missing > " + out_path + "missing_snps_1E4.txt"
    "plink --bfile " + geno_path + " --exclude " + out_path + "missing_snps_1E4.txt --make-bed --out " + geno_path + "_missing1"
    
    #missingness by haplotype (--test-mishap), using P > 1E-4
    "plink --bfile " + geno_path + "_missing1 --test-mishap --out " + geno_path + "missing_hap" 
    "awk '{if ($8 <= 0.0001) print $9 }' " + out_path + "missing_hap.missing.hap > " + out_path + "missing_haps_1E4.txt"
    "sed 's/|/\/g' " + out_path + "missing_haps_1E4.txt > " + out_path + "missing_haps_1E4_final.txt"
    "plink --bfile " + geno_path + " --exclude " + out_path + "missing_haps_1E4_final.txt --make-bed --out " +  geno_path + "_missing2"
    
    #HWE from controls only using P > 1E-4
    "plink --bfile " + geno_path + "_missing2 --filter-controls --hwe 1E-4 --write-snplist"
    "plink --bfile " + geno_path + "_missing2 --extract " + out_path + "plink.snplist --make-bed --out " + geno_path + "_HWE"
    



        

geno = args.geno
out = args.out
geno_het = geno + "_het"
geno_call_rate = geno_het + "_call_rate"
geno_sex = geno_call_rate + "_sex"
        
het_pruning(geno, out)
call_rate_pruning(geno_het, out)
sex_check(geno_call_rate, out)
relatedness_pruning(geno_sex, out)

# THIS IS A TEST FOR GIT