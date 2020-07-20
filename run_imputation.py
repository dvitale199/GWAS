import argparse

#local imports
from gwas.imputation import Impute

parser = argparse.ArgumentParser(description='Arguments for Genotyping Imputation (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--key', type=str, default='nope', help='Please input imputation server key!')
#for now, getting the "out" location from --geno
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

args = parser.parse_args()

geno = args.geno
key = args.key

imputer = Impute(geno)
imputer.impute_prep_data()
imputer.impute_make_vcf()
imputer.impute(key=key)

#not quite ready to use!!!! will delete necessary files
# imputer.cleanup()