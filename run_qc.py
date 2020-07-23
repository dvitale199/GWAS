import argparse

#local imports
from gwas.qc import QC

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')
parser.add_argument('--rare', default=False, action="store_true", help='Pruning toggle for rare variants. If --rare is used, final MAF pruning (0.01) will not be conducted, otherwise, rare variants will be pruned')

args = parser.parse_args()

geno_name = args.geno
rare_flag = args.rare

# now make filenames:
geno_call_rate = geno_name + "_call_rate"
geno_het =  geno_call_rate + "_het"
geno_sex = geno_het + "_sex"
geno_relatedness = geno_sex + "_relatedness"
geno_variant = geno_relatedness + "_variant"
geno_final = geno_variant + "_final"

# INSTANTIATE QC WITH INPUT NAME AND OUTPUT NAME     
qc = QC(geno_name, rare=rare_flag)

# NOW RUN COMMANDS
# FIRST, CLEAR EXISTING LOGFILE
qc.rm_log()

# run het pruning
# qc.call_rate_pruning(geno_name)
# qc.het_pruning(geno_call_rate)
# qc.sex_check(geno_het)
# qc.relatedness_pruning(geno_sex)
qc.variant_pruning(geno_relatedness)
qc.rare_prune(geno_variant)

# need to fix cleanup
# qc.cleanup()

print("DONE!!!!!!!")