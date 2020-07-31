import subprocess
import argparse
import pandas as pd
import os
import glob
import shutil
import sys

#local imports
from plink_helper.plink_driver import Driver


class Ancestry(Driver):
    def __init__(self, geno_path):
        super().__init__(geno_path)
        
    def pca(self, geno_path, ref_path):
        out_path = self.out_path
        step = "PCA for Ancestry Comparison"
        
        bash1 = "plink --bfile " + geno_path + " --bmerge " + ref_path + " --out " + out_path + "bin_snplis --make-bed"
        bash2 = "plink --bfile " + geno_path + " --flip bin_snplis-merge.missnp --make-bed --out " + geno_path + "_flip"
        bash3 = "plink --bfile " + geno_path + "_flip --bmerge " + ref_path + " --out " + out_path + "bin_snplis --make-bed"
        bash4 = "plink --bfile " + geno_path + "_flip --exclude " + out_path + "bin_snplis-merge.missnp --out " + geno_path + "_flip_pruned --make-bed"
        bash5 = "plink --bfile " + geno_path + "_flip_pruned --bmerge " + ref_path + " --out " + out_path + "bin_snplis --make-bed"
        bash6 = "plink --bfile " + out_path + "bin_snplis --geno 0.01 --out " + out_path +  "pca --make-bed --pca 4"
        
        cmds = [bash1, bash2, bash3, bash4, bash5, bash6]
        
        self.run_cmds(cmds, step)
        
    def 