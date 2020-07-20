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
class QC(Driver):
    def __init__(self, geno_path, rare=False):
        super().__init__(geno_path)
        # create new names for each step
        self.geno_call_rate = geno_path + "_call_rate"
        self.geno_het =  self.geno_call_rate + "_het"
        self.geno_sex = self.geno_het + "_sex"
        self.geno_relatedness = self.geno_sex + "_relatedness"
        self.geno_variant = self.geno_relatedness + "_variant"
        self.geno_final = self.geno_variant + "_final"
        self.tmp_file_list = [self.geno_call_rate, self.geno_het, self.geno_sex, self.geno_relatedness, self.geno_variant]
        self.rare = rare

    def call_rate_pruning(self):
            geno_path = self.geno_path
            out_path = self.out_path

            step = "PRUNING FOR CALL RATE"
            print(step)
            bash1 = "awk '{print $1,$2,$6}' " + geno_path + '.fam > ' + geno_path + '.phenos'
            bash2 = "plink --bfile " + geno_path + " --mind 0.05 --make-bed --out " + geno_path + "_call_rate"
            bash3 = "mv " + geno_path + "_call_rate.irem " + out_path + "CALL_RATE_OUTLIERS.txt"

            cmds = [bash1, bash2, bash3]

            self.run_cmds(cmds, step)
            
            
    def het_pruning(self):
        geno_path = self.geno_call_rate
        out_path = self.out_path
        
        step = "PRUNING FOR HETEROZYGOSITY"
        print(step)
        
        
        bash1 = "plink --bfile " + geno_path + " --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out " + out_path + "pruning"
        bash2 = "plink --bfile " + geno_path + " --extract " + out_path + "pruning.prune.in --make-bed --out " + out_path + "pruned_data"
        bash3 = "plink --bfile " + out_path + "pruned_data --het --out " + out_path + "prunedHet"
        bash4 = "awk '{if ($6 <= -0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers1.txt" 
        bash5 = "awk '{if ($6 >= 0.15) print $0 }' " + out_path + "prunedHet.het > " + out_path + "outliers2.txt" 
        bash6 = "cat " + out_path + "outliers2.txt " + out_path + "outliers1.txt > " + out_path + "HETEROZYGOSITY_OUTLIERS.txt"
        bash7 = "plink --bfile " + geno_path + " --remove " + out_path + "HETEROZYGOSITY_OUTLIERS.txt --make-bed --out " + geno_path + "_het"

        cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]
        
        self.run_cmds(cmds, step)

        
    def sex_check(self):
        geno_path = self.geno_het
        out_path = self.out_path
        
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

        self.run_cmds(cmds, step)


    # MAY NEED TO ADD FUNCTION FOR ANCESTRY OUTLIERS (PCR FOR RELATEDNESS)     

    def relatedness_pruning(self):
        geno_path = self.geno_sex
        out_path = self.out_path
        
        step = "RELATEDNESS PRUNING"
        print(step)

        bash1 = "gcta --bfile " + geno_path + " --make-grm --out " + out_path + "GRM_matrix --autosome --maf 0.05" 
        bash2 = "gcta --grm-cutoff 0.125 --grm " + out_path + "GRM_matrix --out " + out_path + "GRM_matrix_0125 --make-grm"
        bash3 = "plink --bfile " + geno_path + " --keep " + out_path + "GRM_matrix_0125.grm.id --make-bed --out " + geno_path + "_relatedness"

        cmds = [bash1, bash2, bash3]

        self.run_cmds(cmds, step)


    ##variant checks
    def variant_pruning(self):
        geno_path = self.geno_relatedness
        out_path = self.out_path
        
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
        bash11 = "#### moved " + geno_path + "_HWE to " + geno_path + "_variant #####"
        cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7, bash8, bash9, bash10, bash11]

        self.run_cmds(cmds, step)

        exts = [".bed",".bim",".fam",".log"]
        for ext in exts:
            
            shutil.move(geno_path + "_HWE" + ext, geno_path + "_variant" + ext)


    def rare_prune(self):

        geno_path = self.geno_variant
        out_path = self.out_path
        rare = self.rare
        # OPTIONAL STEP: if --rare flag included when running, rare variants will be left alone, otherwise they will be pruned with --maf 0.01
        bash = "plink --bfile " + geno_path + " --maf 0.01 --make-bed --out " + geno_path + "_MAF"
        
        if rare:
            print("SKIPPING FINAL MAF PRUNING (0.01)... RARE VARIANTS LEFT ALONE")
             # call logging to do some custom logging
            log = self.logging()
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
            # call logging to do some custom logging
            log = self.logging()
            log.write("RARE VARIANTS (MAF <= 0.O1) PRUNING WITH THE FOLLOWING COMMANDS:")
            log.write("\n")
            log.write("\n")
            subprocess.run(bash, shell=True)
            log.write(bash)
            log.write("\n")

            new_log = open(geno_path + ".log", "r")
            new_log_read = new_log.read()
            new_log.close()

            log.write(new_log_read)
            log.write("\n")
            log.write("***********************************************")
            log.write("\n")
            log.write("\n")

            exts = [".bed",".bim",".fam",".log"]
            for ext in exts:
                shutil.move(geno_path + "_MAF" + ext, geno_path + "_final" + ext)

            log.write("MOVED " + geno_path + "_MAF to " + geno_path + "_final")
            log.write("\n")
            log.close()
            
            
    def cleanup(self):
        print("CLEANING UP THE DIRECTORY OF INTERMEDIATE FILES")
        print("***********************************************")
        print()
            
        save_files = [self.geno_path + ext for ext in ['.bim','.bed','.fam','.log','.hh', '.phenos', '.PLINK_STEPS.log']] + [self.geno_final + ext for ext in ['.bim','.bed','.fam','.hh']] + [self.out_path + 'imputed']
        all_files = glob.glob(self.out_path + '*')
        rm_files = [x for x in all_files if x not in save_files]
        for file in rm_files:
            os.remove(file)
        
