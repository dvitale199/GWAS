import subprocess
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Arguments for Genotyping Imputation (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

args = parser.parse_args()

geno = args.geno
out = args.out




################################################
#DATA PREP

def prep_data(geno_path, out_path):
    os.chdir(out_path)
    # download file to check
    bash1 = "wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip -P " + out_path
     
    bash2 = "unzip " + out_path + "HRC-1000G-check-bim.v4.2.5.zip -d " + out_path
    
    bash3 = "wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz -P " + out_path
    bash4 = "gunzip " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
    # make your .frq file
    bash5 = "plink --bfile " + geno_path + " --freq --out " + geno_path

    bash6 = "perl " + out_path + "HRC-1000G-check-bim.pl -b " + geno_path + ".bim -f " + geno_path + ".frq -r " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"

    # then run to fix your data
    bash7 = "sh " + out_path + "Run-plink.sh"
    
    cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

def vcf_prep(geno_path, out_path):
    # then make vcf files

    for i in range(1,24):
        
        bash1 = "plink --bfile " + geno_path + "-updated-chr" + str(i) + " --recode vcf --chr " + str(i) + " --out " + geno_path + "_chr" + str(i)
        print(bash1)
        subprocess.run(bash1, shell=True)

    ## then sort and zip
    for i in range(1,24):
        
        bash2 = "vcf-sort " + geno_path + "_chr" + str(i) + ".vcf | bgzip -c > pre_impute_" + geno_path + "_chr" + str(i) + ".vcf.gz"
        print(bash2)
        subprocess.run(bash2, shell=True)
    # and then you are ready to submit to the imputation server


# prep_data(geno, out)
vcf_prep(geno, out)