import subprocess
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

args = parser.parse_args()

print("GENOTYPE FILES", args.geno)

geno = args.geno
out = args.out

#QC and data cleaning
def het_pruning(geno_path, out_path):

    bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
    bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
    bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"

    # bash for now. convert these to python
    bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
    bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
    bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
    
    bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_after_heterozyg"
    
    cmd_1 = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
    
    for cmd in cmd_1:
        subprocess.run(cmd, shell=True)
        
het_pruning(geno, out)