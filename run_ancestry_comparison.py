import argparse

#local imports
from gwas.ancestry import Ancestry

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')
parser.add_argument('--ref', type=str, default='ref_data/HAPMAP_hg19_new', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam [default: ref_data/HAPMAP_hg19_new].')

args = parser.parse_args()

geno_name = args.geno
ref_name = args.ref


ancestry = Ancestry(geno_name)
ancestry.pca(geno_name, ref_name)