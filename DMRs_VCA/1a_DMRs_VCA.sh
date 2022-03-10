#!/bin/bash

### Aim: Perform variant decomposition of DMRs using R package heritability or lme4qtl
### Author: Dario Galanti Apr 2021
### Input arg 1): position-sorted unionbed file (sort -k1,1V -k2,2n), with average methylation of individual DMRs for all individuals. Use extr_mergedDMRs_meth.sh to obtain it.
### Input 2): plink_base -> genomic variants file in plink format (map and ped)
### Input 3): trans_ibs -> trans IBS matrix, tab separated, with column names but no rownames
### Input 4): climate_mx -> similarity matrix based on environmental variables, normalized between 0 and 1, in csv format, with column names but no rownames.
### Input 4) NB: Right now climate_mx headers have to be population codes that can be extracted from sample names (r script line 39), but with a little change they can also be sample names (r script env model should be changed to "marker_h2(my_data$DMR_meth, my_data$sample, covariates = NULL, clim_norm, max.iter = 100)"
### Run: bash 1a_DMRs_VCA.sh allspls_5MEF_20ED_joined_146bpDMRs_CpG_merged.bed
### Run on Epi: sbatch --partition test --cpus-per-task 4 --mem 10G --time 60:00:00 --wrap "bash 1a_DMRs_VCA.sh ../merged_DMRs/CpG/allspls_5MEF_20ED_joined_146bpDMRs_CpG_merged.bed"

## IMPORTANT FOR RUNNING THE SCRIPT!!!!
# RUNNING TIME: Using DMRs_VCA_heritability.R --> 5:50min for 20 DMRs. -> ~ 5h for 1000 DMRs -> 50h for 10000
# The script can be run for the whole unionbed file (for all DMRs) or for a subset (can be defined on lines 46 and 47). Decision should be made based on the running time required and available
# If running for subsets, remember to change output file name. And combine output files in the end

## Define tools Epi
plink=~/conda/vcftools/bin/plink
Rscript=~/conda/R4/bin/Rscript		# Use conda R4 to access libraries

## Define input and outputs
unionbed=$1			# /scr/episan/RP07/DMRs/pop_DMRs/merged_DMRs/CpG/allspls_5MEF_20ED_joined_146bpDMRs_CpG_merged.bed
plink_base=/scr/episan/RP07/Ta_SNPs/GATK_v3_imputed/Ta_v3_vrts_1MAF_imputed_GWAspls_withref		# We use mild MAF filtered variants for making cis-IBS matrix
trans_ibs=/scr/episan/RP07/GWAS/input_v3/Ta_v3_vrts_1MAF_imp_8prun_GWAspls_withref.kinship		# We use mild MAF but pruned trans-IBS matrix
climate_mx=/scr/episan/RP07/DMRs/pop_DMRs/VCA_mergedDMRs/Clim_dist_mx_allbiovars_93-18_norm.csv
cont=$(basename $unionbed | grep -o 'C[pH][HG]')	# Extract context
rscript_vers=/scr/episan/RP07/DMRs/pop_DMRs/VCA_mergedDMRs/DMRs_VCA_heritability.R
fout=${cont}_DMRs_h2s.txt
cis_reg=50000		# Define size of cis region (Dubin et al.2015: all variants within 50 kb from the DMR were defined as cis-acting)

# PROCESS:
# Iterate through unionbed file. In each line:
# 1) Calculate DMR midpoint and cis-IBS matrix from variants located 50kb up and downstream
# 2) make txt file with "samples" "DMRmeth" columns for the DMR in the current line
# 3) Run R script that runs models to explin DMR methylation variance and appends results of each DMR to fout

### CLEAN START
rm $fout	# Otherwise new lines will be appended to old file

### Prepare sample list (NB: We do this outside the loop, so the loop can start from further positions)
echo sample > ${cont}_spls.txt
head -1 $unionbed | cut -f4- | tr "\t" "\n" | tr "_" "-" >> ${cont}_spls.txt	# plink crushes with "_"

### Iterate through unionbed file
tail -n+2 $unionbed | while read line							# NB: SKIP FIRST LINE WITH HEADERS!!!!!! # Run for the whole file
#tail -n+2 $unionbed | head -10000 | tail -10000 | while read line	# NB: SKIP FIRST LINE WITH HEADERS!!!!!! # Run for a subset of the file
do
 chr=$(echo $line | cut -d" " -f1)
 st=$(echo $line | cut -d" " -f2)
 end=$(echo $line | cut -d" " -f3)
 cis_st=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);st=(mid-cis_r);cis_st=(st > 0 ? st : 0);print cis_st}')
 cis_end=$(echo $line | awk -v cis_r=$cis_reg '{mid=int(($2+$3)/2);cis_end=(mid+cis_r);print cis_end}')
 DMR=${chr}:${st}-${end}
 #echo $chr $st $end $middle $cis_st $cis_end
 $plink --file $plink_base --allow-extra-chr --chr ${chr} --from-bp ${cis_st} --to-bp ${cis_end}  --distance square ibs flat-missing --out ${cont}_cis_plink_ibs
 echo $(cat ${cont}_cis_plink_ibs.mibs.id | cut -f1 | tr "\n" "\t") | tr " " "\t" > ${cont}_cis_ibs.kinship
 cat ${cont}_cis_plink_ibs.mibs >> ${cont}_cis_ibs.kinship
 rm ${cont}_cis_plink_ibs.*
 echo DMR_meth > ${cont}_tmp_line.txt
 echo $line | cut -d" " -f4- | tr " " "\n" >> ${cont}_tmp_line.txt
 paste ${cont}_spls.txt ${cont}_tmp_line.txt > ${cont}_DMR_pheno.txt
 rm ${cont}_tmp_line.txt
 $Rscript $rscript_vers ${cont}_cis_ibs.kinship $trans_ibs $climate_mx ${cont}_DMR_pheno.txt $DMR $fout
done

### Cleanup
rm ${cont}_spls.txt ${cont}_cis_ibs.kinship ${cont}_DMR_pheno.txt

### OPTIONAL: DMR binning based on major predictor
## NB: NAs in the cis column are converted to 0 for the DMR binning or NA results higher than any other number.
## NAs are generated when the model does not converge, so probably very low amount of variation is explained
fout2=${cont}_DMRs_h2s_kind.txt
awk '{OFS="\t";if(NR>1){if($4=="NA"){cis=0}else{cis=$4};if(cis>0.1 || $5>0.1 || $6>0.1){if(cis>$5 && cis>$6){kind="cis"}else{if($5>$6){kind="trans"}else{kind="env"}}} \
else{kind="unexplained"}}else{kind="major_predictor"};print $0,kind}' ${fout} > ${fout2}
rm ${fout}
