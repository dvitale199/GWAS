import subprocess
import argparse
import os
import requests
import json

#local imports
from plink_helper.plink_driver import Driver


# parser = argparse.ArgumentParser(description='Arguments for Genotyping Imputation (data in Plink .bim/.bam/.fam format)')
# parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

# args = parser.parse_args()

# geno = args.geno
# out = args.out
geno = '/data/vitaled2/test_data/PDBP/PDBP_het_call_rate_sex_relatedness_variant_final'
out = '/data/vitaled2/test_data/PDBP/'



################################################

class Impute(Driver):
    def __init__(self, geno_path):
        super().__init__(geno_path)
        self.vcf_list_for_impute = [geno_path + "_pre_impute" + "_chr" + str(i) + ".vcf.gz" for i in range(1,24)]
        
        
    def impute_prep_data(self):
        geno_path = self.geno_path
        out_path = self.out_path
        step = "PREP PLINK FILES FOR IMPUTATION"
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

#         for cmd in cmds:
#             subprocess.run(cmd, shell=True)
        self.run_cmds(cmds, step)
            

    def impute_make_vcf(self):
        geno_path = self.geno_path
        out_path = self.out_path
        # then make vcf files
        step1 = "RECODE PLINK FILES TO VCF"
        
        cmds1 = ["plink --bfile " + geno_path + "-updated-chr" + str(i) + " --recode vcf --chr " + str(i) + " --out " + geno_path + "_chr" + str(i) for i in range(1,24)]   
#         for i in range(1,24):

#             bash1 = "plink --bfile " + geno_path + "-updated-chr" + str(i) + " --recode vcf --chr " + str(i) + " --out " + geno_path + "_chr" + str(i)
#             print(bash1)
#             subprocess.run(bash1, shell=True)
        self.run_cmds(cmds1, step1)
        ## then sort and zip
        step2 =  "vcf-sort AND bgzip VCFS"
        cmds2 = ["vcf-sort " + geno_path + "_chr" + str(i) + ".vcf | bgzip -c > " + geno_path + "_pre_impute" + "_chr" + str(i) + ".vcf.gz" for i in range(1,24)]
#         for i in range(1,24):

#             bash2 = "vcf-sort " + geno_path + "_chr" + str(i) + ".vcf | bgzip -c > " + geno_path + "_chr" + str(i) + "_pre_impute.vcf.gz"
#             print(bash2)
#             subprocess.run(bash2, shell=True)
        self.run_cmds(cmds2, step2)

        # and then you are ready to submit to the imputation server
        
    def impute(self, key, input_population='eur', vcf_list=None):
        geno_path = self.geno_path
        vcf_list = self.vcf_list_for_impute
        
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see Authentication)
        headers = {'X-Auth-Token' : key}

        open_vcfs = [open(vcf,'rb') for vcf in vcf_list]
        
        # submit new job
#         vcf1 = '/data/vitaled2/test_data/PDBP/PDBP_het_call_rate_sex_relatedness_variant_final_pre_impute_chr1.vcf.gz'
#         vcf2 = '/data/vitaled2/test_data/PDBP/PDBP_het_call_rate_sex_relatedness_variant_final_pre_impute_chr2.vcf.gz';
        files = {'input-files': open_vcfs[0],
                'input-files': open_vcfs[1],
                'input-files': open_vcfs[2],
                'input-files': open_vcfs[3],
                'input-files': open_vcfs[4],
                'input-files': open_vcfs[5],
                'input-files': open_vcfs[6],
                'input-files': open_vcfs[7],
                'input-files': open_vcfs[8],
                'input-files': open_vcfs[9],
                'input-files': open_vcfs[10],
                'input-files': open_vcfs[11],
                'input-files': open_vcfs[12],
                'input-files': open_vcfs[13],
                'input-files': open_vcfs[14],
                'input-files': open_vcfs[15],
                'input-files': open_vcfs[16],
                'input-files': open_vcfs[17],
                'input-files': open_vcfs[18],
                'input-files': open_vcfs[19],
                'input-files': open_vcfs[20],
                'input-files': open_vcfs[21],
                'input-files': open_vcfs[22],
                }
    
        data = {'input-mode' : 'imputation',
                 'input-files-source': 'file-upload',
                'input-refpanel': 'apps@hrc-r1.1',
                 'input-phasing': 'eagle',
                 'input-population': input_population}

        r = requests.post(url + "/jobs/submit/minimac4", files=files, headers=headers, data=data)
        if r.status_code != 200:
            raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))

        # print message
        print(r.json()['message'])
        print(r.json()['id'])

        

imputer = Impute(geno)
# imputer.impute_prep_data()
# imputer.impute_make_vcf()
imputer.impute(key='hehehehe')
