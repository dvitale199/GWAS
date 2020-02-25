#create merge list
ls /data/CARD/AD/ADNI/sequence_data_hg38/genotypes/*.bim | sed 's/[.].*$//' > /data/vitaled2/test_data/ADNI/merge.list

# this will fail if there are multiallelic positions
plink --allow-extra-chr --biallelic-only --keep-allele-order --autosome --make-bed --merge-list /data/vitaled2/test_data/ADNI/merge.list --out /data/vitaled2/test_data/ADNI/adni

#when fail occurs, run again with --exclude adni-merge.missnp
for i in {1..22} X; do plink --bfile /data/CARD/AD/ADNI/sequence_data_hg38/genotypes/adni.nov2018.chr$i --allow-extra-chr --keep-allele-order --autosome --exclude /data/vitaled2/test_data/ADNI/adni-merge.missnp --make-bed --out /data/vitaled2/test_data/ADNI/chr${i}_multi_prune; done

#now create merge list again with new pruned files
ls /data/vitaled2/test_data/ADNI/*.bim | sed 's/[.].*$//' > /data/vitaled2/test_data/ADNI/merge.list

# now try merge again
plink --allow-extra-chr --biallelic-only --keep-allele-order --autosome --make-bed --bmerge-list /data/vitaled2/test_data/ADNI/merge.list --out /data/vitaled2/test_data/ADNI/clean/adni