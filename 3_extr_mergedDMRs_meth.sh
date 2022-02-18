#!/bin/bash

### Aim: Extract methylation of set of non-overlapping DMRs using Adam's average_over_bed.py 
### Author: Dario Galanti Mar 2021
### Input 1): position-sorted unionbedg file (properly sorted, so "sort -k1,1V -k2,2n")
### Input 2): merged_DMRs file or any bed file with sorted, non-overlapping regions (columns after the third will be ignored)
### Run: bash 03_extr_mergedDMRs_meth.sh unionbedg.bed region_file
### Run: sbatch --partition test --cpus-per-task 2 --mem 20G --time 08:00:00 --wrap "bash 3_extr_mergedDMRs_meth.sh METHunion_filtered_CpG.bed merged_DMRs/CpG/joined_146bpDMRs_CpG_merged.bed"

## Dependencies: This script is using Adam Nunn's average_over_bed.py script (https://github.com/EpiDiverse/ewas/blob/master/bin/average_over_bed.py)

### Steps:
### 1) Intersect region_bed and unionbedg and discard DMRs with less than x=5 Cs covered
### 2) Extract methylation of individual DMRs -> Adam's average_over_bed.py
### 3) Filter DMRs by MEF and ED

# DEFINE INPUT, OUTPUT, CONTEXT AND INDEX FILE
unionbed=$1										# /scr/episan/RP07/bam_aligned/WGBS/unionbed_v3_3cov_0.25NAs/CpG_unionbed_v3_0.25NAs.bed
cont=$(echo $unionbed | grep -o 'C[pH][HG]')	# Extract context
regions=$2										# /scr/episan/RP07/DMRs/pop_DMRs/merged_DMRs/CpG/joined_146bpDMRs_CpG_merged.bed
x=5				# NB: Minimum number of Cs per region that should be present in the unionbed, regions with less Cs are discarded (NB: Context specific, so don't be too restrictive!!)
wDir=$(dirname $regions)
covered_regions=${wDir}/covered_$(basename $regions)
#subset_unionbed=${wDir}/${cont}_DMRs_unionbed_v3_0.25NAs.bed
index=${wDir}/${cont}_index.txt		# Will be created! Single column txt file containing all scaffold names.
average_over_bed=/scr/episan/RP07/region_meth/gene_meth/gene_methIMP/average_over_bed.py
fout1=${wDir}/allspls_$(basename $regions)

# IMPORTANT: Define MEF and ED to further filter DMRs
MEF=0.05	# Define Minor Epiallele Frequency (proportion of samples which need to have differential methylation from the others)
# Define Epiallele Difference (minimum methylation difference to define different epialleles) NB: Use 20 for CpG and 15 for CHG and CHH
if [ ${cont} == "CpG" ];then ED=20;else ED=15;fi

fout2=${wDir}/allspls_5MEF_${ED}ED_$(basename $regions)

## PRE-STEP: Make an index file from the unionbed file, containing a single column with all scaffolds
#tail -n+2 METHunion_filtered_CpG.bed | cut -f1 | uniq > $index

## NB: All DMRs should be covered in the unionbed file, but we should still make sure of that
## 1) Intersect region_bed and unionbedg to discard uncovered DMRs (less than "x" Cs) and report coverage of covered ones for further filtering
#echo Intersecting inputs and discarding uncovered regions
awk 'NR!=1{OFS="\t";print $1,$2,$3}' $unionbed | bedtools intersect -a $regions -b stdin -c | awk -v x=$x '{if($NF>=x){OFS="\t";print $1,$2,$3,$NF}}' | sort -k1,1V -k2,2n > $covered_regions

## 2) Extract regions methylation with average_over_bed.py
tail -n+2 $unionbed | cut -f1 | uniq > $index
python3 $average_over_bed $covered_regions $unionbed $index $fout1
rm $index

## 3) Filter DMRs by MEF and ED
awk -v MEF=$MEF -v ED=$ED 'OFS="\t"{if(NR==1){print}else{for(i=4;i<=NF;i++){if($i!="NA"){list[i]=$i}};asort(list);
spls=length(list); MEC=int(MEF*spls+0.999); MECdiff=(list[(spls-MEC+1)]-list[MEC]);
if(MECdiff>ED){print}; delete list}}' $fout1 > $fout2

## Cleanup
rm $covered_regions

