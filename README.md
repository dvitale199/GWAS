# This repo is no longer being developed and has been replaced by: https://github.com/dvitale199/GenoTools




# GWAS Pipeline- Quality Control, Imputation, GWAS, and more!
This repo provides easy execution for automated GWAS Pipelines which are used by the Laboratory of Neurogenetics of the National Institute on Aging, National Institutes of Health. The basic pipeline was written by Cornelis Blauwendraat (Laboratory of Neurogenetics) and was automated by Dan Vitale (Data Tecnica International/LNG) This is a work in progress but will currently run automated Quality Control and and Imputation (Michigan Imputation Server) given genotypes in Plink (.bed/.bim/.fam) format.

#### Go have a look at GWAS_QC.ipynb to see how to run. You must have Plink (v1.9), GCTA, vcf-sort, and Perl loaded up in path to run this!

### Coming Soon
- Fully automated GWAS execution (PLINK) along with diagnostics and visualizations (Scree plots, QQ plots, Manhattan plots) which were automated as an effort in the International Parkinson's Disease Genetics Consortium (IPDGC) 2020 Virtual Hackathon.
- TOPMed Imputation Server API access
- Dynamic Loading Plink, GCTA, and other tools used in the repo
- Automated Meta-analysis of GWAS summary stats
- Ability to direct output to directory of choice (currently outputs to dir containing the genotypes
- Post-imputation directory cleanup (function exists but is currently not working)
- Full testing and tuning with new, uncleaned genotypes (note: at this point, the pipeline has only been tested on pre-cleaned datasets)

# Some other important information
The functions in this repository utilize functions in plink_logger.py and plink_driver.py (in plink_helper directory). These are just methods for running plink (and other commands through bash) and compiling the steps and logs into a single log file
