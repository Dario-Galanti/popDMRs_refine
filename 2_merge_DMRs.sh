#!/bin/bash

## Author: Dario Galanti 2020 (modified from Cristian Peña https://github.com/EpiDiverse/scripts/blob/master/merge_DMRs.sh)
## Aim: This script merges all bed files (included in an input directory) in one unique bed file. In this way each context can be summarized.
## Input directory: run script in DMR pipe output directory (outdir/${context}/${comparison}/ )
## Run: bash 2_merge_DMRs.sh context outputDirectory
## Run: bash 2_merge_DMRs.sh CpG merged_DMRs
## Dependencies: bedtools (https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)


cont=$1
output=$2/$1
mkdir -p ${output}
FILES=$(ls ${cont}/*.bed | tr "\n" " ")
NAMES=$(ls ${cont}/*.bed | cut -f2 -d/ | rev | cut -f2- -d . | rev | tr "\n" " ")

cat ${cont}/*.bed > ${output}/DMRs.bed
bedtools sort -i ${output}/DMRs.bed > ${output}/DMRs_sorted.bed
bedtools merge -i ${output}/DMRs_sorted.bed -c 4 -o count > ${output}/DMRs_${cont}_merged.bed
echo $(wc -l ${output}/DMRs_${cont}_merged.bed) merged DMRs created

rm ${output}/DMRs.bed
rm ${output}/DMRs_sorted.bed
