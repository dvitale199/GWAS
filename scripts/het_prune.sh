#!/bin/bash
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet
awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
plink --bfile $FILENAME --remove all_outliers.txt --make-bed --out $FILENAME_after_heterozyg