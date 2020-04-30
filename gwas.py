
# local imports
from plink_helper.plink_driver import Driver


pheno_path = '/data/vitaled2/test_data/PDBP/PDBP.phenos'
imputed_geno_path = '/data/vitaled2/test_data/PDBP/imputed/'

class GWAS(Driver):
    def __init__(self, geno_path):
        super().__init__(geno_path)

    def vcf_to_plink(self, phenos):
#         vcf_to_plink_cmds = ['plink --vcf chr' + i + '.dose.vcf.gz --make-bed --out ' for i in ]
        pass
    
    def pca(self, exclusion_regions_file):
        geno_path = self.geno_path
        out_path = self.out_path
        geno_filtered = geno_path + '_filtered'
        geno_filtered_pruned = geno_filtered + '_pruned'
        bash1 = 'plink --bfile ' + geno_path + ' --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude ' + exclusion_regions_files + ' --make-bed --out ' + geno_path + '_filtered'
        # Prune snps 
        bash2 = 'plink --bfile ' + geno_path + ' --indep-pairwise 1000 10 0.02 --autosome --out ' + out_path + 'pruned_data'
        # Extract pruned SNPs and only these variants will be used for PC calculation
        bash3 = 'plink --bfile ' + geno_filtered + ' --extract pruned_data.prune.in --make-bed --out ' + geno_filtered_pruned
        # Calculate PCs
        bash4 = 'plink --bfile ' + geno_filtered_pruned + ' --pca --out PCA'
        
        
    