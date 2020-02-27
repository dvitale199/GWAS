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
     
    bash2 = "unzip " + out_path + "HRC-1000G-check-bim.v4.2.5.zip"
    
    # make your .frq file
    bash3 = "plink --bfile " + geno_path " --freq --out " geno_path

    bash4 = "perl " + out_path + HRC-1000G-check-bim.pl -b $FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# then run to fix your data
sh Run-plink.sh

# then make vcf files

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile YOURFILE-updated-chr$chnum --recode vcf --chr $chnum --out YOURFILE$chnum 
done

## then sort and zip

module load vcftools

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	vcf-sort YOURFILE$chnum.vcf | bgzip -c >  pre_impute_YOURFILE_$chnum.vcf.gz
done

# and then you are ready to submit to the imputation server