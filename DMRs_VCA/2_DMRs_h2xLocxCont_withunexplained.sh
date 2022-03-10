#!/bin/bash

### Aim: Calculate average heritability and number of mergedDMRs overlapping with genes, TEs, promoters and intergenic regions
### Author: Dario Galanti May 2021
### Input 1): DMRs VCA files for all contexts from the DMRs_VCA.sh script. Only the CpG file is given as input argument
### Input 2): genomic features bed files. To obtain them check https://github.com/Dario-Galanti/region_meth
### Run: bash 2_DMRs_h2xLocxCont_withunexplained.sh CpG_DMRs_h2s.txt

# DEFINE INPUT FILES
gene_regions="/scr/episan/RP07/region_meth/gene_meth/gene_meth/Ta_gene_regions.bed"
TE_regions="/scr/episan/RP07/region_meth/TE_meth/TE_meth/Ta_TE_regions.bed"
prom_regions="/scr/episan/RP07/region_meth/promoter_meth/prom_meth/Ta_promoter_regions_merged.bed"
intergenic_regions="/scr/episan/RP07/region_meth/intergenic_meth/intergenic_meth/Ta_intergenic.bed"

# DEFINE DMR INPUT FILES
CpG_DMRs=$1
CHG_DMRs=$(echo $CpG_DMRs | sed 's/CpG/CHG/g')
CHH_DMRs=$(echo $CpG_DMRs | sed 's/CpG/CHH/g')

# OUTPUT DIR AND FILE
outDir=DMRs_h2_summary
fout=${outDir}/summary_withunexplained_$(basename $CpG_DMRs | sed 's/CpG_//g')
mkdir -p ${outDir}

# MAKE FILE ARRAYS
reg_arr=( $gene_regions $TE_regions $prom_regions $intergenic_regions )
cont_arr=( $CpG_DMRs $CHG_DMRs $CHH_DMRs )

# INTERSECT ALL CONTEXT AND REGION FILES
echo -e location"\t"context"\t"tot_dmrs"\t"cis_h2"\t"trans_h2"\t"env_h2"\t"cis_DMRs"\t"trans_DMRs"\t"env_DMRs"\t"unexplained_DMRs > $fout
for cont in ${cont_arr[@]};
do
	for reg in ${reg_arr[@]};
	do
		context=$(echo $cont | grep -o 'C[pH][HG]' | head -1)
		region=$(basename $reg .bed | cut -d"_" -f2)
		
		info_DMRs=$(tail -n+2 $cont | bedtools intersect -a stdin -b $reg -f 0.05 -wa \
		 | awk 'BEGIN{cis=0;tr=0;env=0;unexp=0} {cis_h2+=$4;tr_h2+=$5;env_h2+=$6;if($4=="NA"){c=0}else{c=$4};
		 if(c>0.1 || $5>0.1 || $6>0.1){if(c>$5 && c>$6){cis++}else{if($5>$6){tr++}else{env++}}}else{unexp++}} \
		 END{OFS="\t";print NR,(cis_h2/NR),(tr_h2/NR),(env_h2/NR),cis,tr,env,unexp}')
		echo -e $region"\t"$context"\t"$info_DMRs | tr " " "\t" >> $fout
		rm ${context}_${region}_cis.tr.env_h2.bed
	done
done


