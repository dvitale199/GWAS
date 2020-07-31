import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

#local imports
from plink_helper.plink_driver import Driver

#QC and data cleaning
class gwas(Driver):
    def __init__(self, geno_path):
        super().__init__(geno_path)
        
        
    def pca(geno_path, out_path, exclusion_region):

       step= "GENERATING PCs"
       print(step)

        #Make sure to use high-quality SNPs 
        bash1 = 'plink --bfile ' + geno_path + ' --maf 0.01 --geno 0.05 --hwe 1E-6 --exclude range ' + exclusion_region + ' --make-bed --out filter'
        # Prune out unnecessary SNPs (only need to do this to generate PCs)
        bash2 = 'plink --bfile filter --indep-pairwise 1000 10 0.02 --out prune'
        # Keep only pruned SNPs (only need to do this to generate PCs)
        bash3 = 'plink --bfile filter --extract prune.prune.in --make-bed --out prune' 
        # Generate PCs 
        bash4 = 'plink --bfile prune --pca 20 --out ' + out_path
        bash5 = 'rm filter* prune*'

        cmds = [bash1, bash2, bash3, bash4, bash5]
        
        self.run_cmds(cmds, step)