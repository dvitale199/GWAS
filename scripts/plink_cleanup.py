import pandas as pd
#this is a short script to fix FID and IIDs of post-imputation .psam file
sam = pd.read_csv("/data/vitaled2/PD/DUTCH/chr1_dose_imputed.psam" , sep='\t')
fam = pd.read_csv("/data/CARD/PD/genotype_data/DUTCH/Dutch_GWAS.fam", sep=' ', header=None)
# fix FID and IID
sam_new = pd.DataFrame()
sam_new[["#FID","IID"]] = sam.IID.str.split('_',1,expand=True)

#merge and add "SEX" column to sam_new
merged =  sam_new.merge(fam, how='left', left_on=["#FID","IID"], right_on=[0,1])
sam_new["SEX"] = merged[4]

#make "pheno" file
pheno = merged[["#FID","IID",5]]

sam_new.to_csv("/data/vitaled2/PD/DUTCH/chr1_dose_imputed_fixedIDs_sex.psam", sep=' ', index=False)
pheno.to_csv("/data/vitaled2/PD/DUTCH/dutch_phenos.txt", sep=' ', index=False, header=False)
# now need to copy .pgen and .pvar to files with matching name!