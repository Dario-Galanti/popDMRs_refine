# DMRs_VCA
Scripts to perform Variance Component Analysis of sequence context-specific DMRs (CpG, CHG or CHH), to quantify cis, trans and environment variance explained.

This repository provides tools to quantify DMR variance explained by three predictors: i) cis genetic variance, ii) trans genetic variance and iii) environment of origin.
For each DMR, it runs three mixed models, each one with one of three random factor matrices capturing cis, trans and environment. <br/>
i) The cis-IBS matrix is calculated in the script based on all vairants within 50kb from the DMR center, also using [PLINK](https://zzz.bwh.harvard.edu/plink/).<br/>
ii) The trans-IBS matrix can be obtained with [PLINK](https://zzz.bwh.harvard.edu/plink/). <br/>
iii) The environmental matrix can be obtained as euclidean distance between locations of origin using environmental variables and should be normalized between 0 and 1 and inverted to a similarity matrix.

SCRIPTS DESCRIPTION: <br/>
[1a_DMRs_VCA.sh](https://github.com/Dario-Galanti/popDMRs_refine_VCA/blob/main/DMRs_VCA/1a_DMRs_VCA.sh)<br/>
Iterate through all DMRs, calculate cis-IBS matrix and run 1b_DMRs_VCA_heritabiliy.R for variance decomposition.


[1b_DMRs_VCA_heritability.R](https://github.com/Dario-Galanti/popDMRs_refine_VCA/blob/main/DMRs_VCA/1b_DMRs_VCA_heritability.R)<br/>
Run variance decomposition of an individual DMR.


[2_DMRs_h2xLocxCont_withunexplained.sh](https://github.com/Dario-Galanti/popDMRs_refine_VCA/blob/main/DMRs_VCA/2_DMRs_h2xLocxCont_withunexplained.sh)<br/>
Summarize the results: 1) calculate context-specific and genomic feature specific variance explained. 2) DMR binning based on the major predictor.
