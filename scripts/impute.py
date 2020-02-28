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
    # download file to check
    bash1 = "wget " + "http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip -P " + out_path
     
    bash2 = "unzip " + out_path + "HRC-1000G-check-bim.v4.2.5.zip -d " + out_path
    ######### ADD LINE BELOW FOR REF
    # wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
    
    # make your .frq file
    bash3 = "plink --bfile " + geno_path + " --freq --out " + geno_path

    bash4 = "perl " + out_path + "HRC-1000G-check-bim.pl -b " + geno_path + ".bim -f " + geno_path + ".frq -r " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"

    # then run to fix your data
    bash5 = "sh " + out_path + "Run-plink.sh"
    
    cmds = [bash1, bash2, bash3, bash4, bash5]
    
    for cmd in cmds:
        subprocess.run(cmd, shell=True)

        
        
prep_data(geno, out)        
#     # then make vcf files

#     for i in range(1,24):
#       plink --bfile YOURFILE-updated-chr$chnum --recode vcf --chr $chnum --out YOURFILE$chnum 


#     ## then sort and zip

#     for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
#       do
#         vcf-sort YOURFILE$chnum.vcf | bgzip -c >  pre_impute_YOURFILE_$chnum.vcf.gz
#     done

#     # and then you are ready to submit to the imputation server