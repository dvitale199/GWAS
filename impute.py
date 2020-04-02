import subprocess
import argparse
import os
import requests
import json
import time

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
        
    def check_impute_status(self, _key, _id):
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see authentication)
        headers = {'X-Auth-Token' : _key }

        # get all jobs
        r = requests.get(url + "/jobs", headers=headers)
        if r.status_code != 200:
            raise Exception('GET /jobs/ {}'.format(r.status_code))
        
        status = r.json()
        for stat in status['data']:
            if stat['id'] == _id:
                if stat['state'] == 1:
                    print("Launching Job:", stat['id'])
                elif stat['state'] == 2:
                    print("Running Job:", stat['id'])
                elif stat['state'] == 3:
                    print(stat['id'], "returned state '3', have a look at jobs on the web front for more information")
                elif stat['state'] == 5:
                    print(stat['id'], "has failed. consult docs on data input to ensure your vcfs are correct")
                elif stat['state'] == 4:
                    print(stat['id'], "COMPLETED!")
                
                return stat['state']
            
            else:
                pass
        
        
        
#         for stat in status['data']:
#             if stat['id'] == _id:
#                 if stat['complete']:
#                     print(stat['id'],'is complete')
#                     #insert script to pull results!!!!
#                 else:
#                     print(stat['id'],"still running!")
        
#                 return stat['complete']

#             else:
#                 pass
                
                
        # print all jobs
#         for job in r.json():
#             print('{} [{}]'.format(job['id'], job['state']))

    def pull_imputed_data(self):
        pass
   

    def impute(self, key, input_population='eur', pw='imputer', vcf_list=None):
        geno_path = self.geno_path
        vcf_list = self.vcf_list_for_impute
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see Authentication)
        headers = {'X-Auth-Token' : key}

        open_vcfs = [open(vcf, 'rb') for vcf in vcf_list]
        
        files = set([('input-files-upload', vcf) for vcf in open_vcfs])

        data = {'input-mode' : 'imputation',
                'input-files-source': 'file-upload',
                'input-password': pw,
                'input-refpanel': 'apps@hrc-r1.1',
                'input-phasing': 'eagle',
                'input-population': input_population}

        r = requests.post(url + "/jobs/submit/minimac4", files=files, headers=headers, data=data)
        if r.status_code != 200:
            raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))
        
        impute_id = r.json()['id']
        message = r.json()['message']

        print(message)
        print(impute_id)
        print('***************************')
        print('* * * * * * * * * * * * * *') 
        
        imp_state = 0
        while imp_state < 3:
            time.sleep(5)
            print('***************************')
            time.sleep(5)
            print('***************************')
            time.sleep(5)
            print('***************************')
            time.sleep(5)
            print('***************************')
            time.sleep(5)
            print('***************************')
            time.sleep(5)
            os.system('cls')
            imp_state = check_impute_status(key, impute_id)
        
    
    

imputer = Impute(geno)
# imputer.impute_prep_data()
# imputer.impute_make_vcf()
imputer.impute(key='***REMOVED***')
