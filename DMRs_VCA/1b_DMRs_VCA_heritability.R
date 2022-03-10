## Author: Dario Galanti Apr 2021
## Aim: Calculate heritability of individual DMRs based on trans_ibs, cis_ibs and environmental distance using "heritability" package
## Input is described below
## Run: Rscript 1b_DMRs_VCA_heritability.R cis_ibs.kinship trans_ibs.kinship climate_mx.csv DMR_pheno.txt Chr1:245-481
## Note: If possible, run the script recursively without closing the R environment. This will keep the trans and env matrixes open and save some time.


## MODELS USED: This version of the script runs the following models for each DMR:
# marker_h2(my_data$DMR_meth, my_data$env)    using "heritability"
# marker_h2(my_data$DMR_meth, my_data$cis)    using "heritability"      -> rarely failing
# marker_h2(my_data$DMR_meth, my_data$trans)  using "heritability"

### SETUP
library("dplyr")
library("data.table")
library(heritability)
library(tidyr)       # for function "unite()"
library(stringr)

args <- commandArgs(trailingOnly = TRUE)


### DEFINE INPUT FILES
cis_mx <- args[1]   # cis IBS matrix in csv format with column names but no rownames.
trans_mx <- args[2] # trans IBS matrix in csv format with column names but no rownames.
climate_mx <- args[3] # similarity matrix based on environmental variables, normalized between 0 and 1. Can be obtained from a normal distance matrix with "clim_norm <- 1-(clim_mx-min(clim_mx))/(max(clim_mx)-min(clim_mx))"
DMR_pheno <- args[4] # 2 columns txt file with sample names and methylation values. Filename should contain context
cont <- str_extract(DMR_pheno, "C[p,H][G,H]")
DMR <- args[5] # String with DMR coordinates eg. Chr1:13823-13935
DMRvec <- strsplit(DMR,":|-")[[1]]
### DEFINE OUTPUT FILE. IF FIRST DMR, PRINT HEADERS AND CLEAR ENVIRONMENT
fout <- args[6] # Output file name
headers <- paste("chr","st","end","cis_h2","trans_h2","env_h2",sep="\t")
if(file.exists(fout)==FALSE){write(headers, file=fout, sep="\t");rm(trans_ibs,clim_norm)}

### 1) READING DMR PHENOTYPE AND ADD COLUMNS
my_data <- read.table(DMR_pheno, header = TRUE)
#my_data$id <- my_data$sample           # Add one column with samples for one of the two ibd matrixes
my_data$Pop <- substr(my_data$sample,4,8)
### 2) CLEAN DMR PHENOTYPE
my_data <- my_data[!is.na(my_data$DMR_meth),]

### 2) READ TRANS IBS MATRIX
## Only reload matrix if not already loaded in the environment
if(exists("clim_norm")==FALSE){
  trans_ibs <- as.matrix(fread(trans_mx))
  rownames(trans_ibs) <- colnames(trans_ibs) #Because the kinshim matrix is printed out without rownames
  vec <- colnames(trans_ibs) %in% my_data$sample
  trans_ibs <- trans_ibs[vec,vec]
}

### 3) READ CIS IBS MATRIX
# NB: When the cis_ibs is made of few markers only, you may end up with a non positive definite matrix (a singular matrix with determinant = 0)!!!
# That means that at least one of your variables can be expressed as a linear combination of the others.
# In our case some individuals (rows and cols) are simply identical
# This will trigger an error. See https://stats.stackexchange.com/questions/30465/what-does-a-non-positive-definite-covariance-matrix-tell-me-about-my-data#:~:text=The%20covariance%20matrix%20is%20not,a%20subset%20of%20the%20others.
cis_ibs <- as.matrix(fread(cis_mx))
rownames(cis_ibs) <- colnames(cis_ibs) #Because the kinshim matrix is printed out without rownames
vec <- colnames(cis_ibs) %in% my_data$sample
cis_ibs <- cis_ibs[vec,vec]


### 4) READ ENVIRONMENTAL MATRIX (reverse and normalise only if not already done)
## NB: climate_mx is a distance matrix, not similarity like IBS. So we normalize between 0 and 1 and reverse it!!!
## Only reload matrix if not already loaded in the environment
if(exists("clim_norm")==FALSE){
  clim_norm <- as.matrix(fread(climate_mx, sep = ","))
  rownames(clim_norm) <- colnames(clim_norm)
  #clim_norm <- 1-(clim_mx-min(clim_mx))/(max(clim_mx)-min(clim_mx))
}

### 5) CALCULATE HERITABILITIES
try(cis_h2 <- marker_h2(my_data$DMR_meth, my_data$sample, covariates = NULL, cis_ibs, max.iter = 100))
if( exists("cis_h2")==FALSE ){cis_h2 <- "NA"}else{cis_h2 <- round(cis_h2$h2,digits=4)} # Sometimes the cis_ibs is made of too few markers and the model doesn't converge
trans_h2 <- marker_h2(my_data$DMR_meth, my_data$sample, covariates = NULL, trans_ibs, max.iter = 100)
trans_h2 <- round(trans_h2$h2,digits=4)
env_h2 <- marker_h2(my_data$DMR_meth, my_data$Pop, covariates = NULL, clim_norm, max.iter = 100)
env_h2 <- round(env_h2$h2,digits=4)

### 6) APPEND RESULTS TO FOUT
newline <- c(DMRvec[1], DMRvec[2], DMRvec[3], cis_h2, trans_h2, env_h2)
newline <- paste(newline, collapse="\t")
write(newline, file=fout, append=TRUE, sep="\t")




