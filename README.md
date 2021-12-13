# popDMRs_refine
Workflow for downstream analysis of EpiDiverse [DMR pipeline](https://github.com/EpiDiverse/dmr) results to merge and refine comparison-specific DMRs from different pairwise-comparisons.

The EpiDiverse [DMR pipeline](https://github.com/EpiDiverse/dmr) is a great tool for calling context-specific (CpG, CHG or CHH) Differentially Methylated Regions between groups of samples. This is a downstream workflow that refines DMRs for large datasets where the pipeline needs to perform a large amount of pairwise comparisons between groups. It includes 1) joining of supposedly "fragmented" DMRs (pairwise-specific short and close DMRs), 2) merging of DMRs from different pairwise comparisons, 3) methylation extraction of the newly obtained merged-DMRs using a unionbed file and final DMR filtering.

WORKFLOW DESCRIPTION: <br/>

[1_join_fragmented_DMRs.py](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/1_join_fragmented_DMRs.py) <br/>
[DMR pipeline](https://github.com/EpiDiverse/dmr) fragmentation acts differently depending on context (CpG, CHG or CHH), sometimes leading to CHH-DMRs being particularly "fragmented". This python script joins all comparison-specific DMRs that are closer than a user-defined distance (default 146bp as fragmentation default in the DMR pipeline) and have the same directionality (higher methylation in the same comparison group). 
DMR length distribution (left figure below) should be checked to see if very small DMRs are more abundant and if there are differences between contexts. But especially the distribution of distance between subsequent DMRs (right figure below) should be observed to see if very close DMRs are over-represented. This can be done with [DMR_distance_distr.py](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/DMR_distance_distr.py).
An over-representation of short and close DMRs as in the figures below suggests excessive "fragmentation".

![image](https://user-images.githubusercontent.com/58292612/145868384-42509bda-5227-4e93-8f83-56f7cdf6fabf.png)


[2_merge_DMRs.sh](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/2_merge_DMRs.sh) <br/>
Merging DMRs from all pairwise comparisons in the [DMR pipeline](https://github.com/EpiDiverse/dmr) output directory to optain DMRs for the whole dataset (merged-DMRs).

[3_extr_mergedDMRs_meth.sh](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/3_extr_mergedDMRs_meth.sh) <br/>
Extract methylation of merged-DMRs from a unionbed file and filter merged-DMRs. The filtering is based on two user-defined parameters: 1) a minimum Methylation Difference (MD) and 2) a minimum proportion of samples which need to show the aforementioned difference, which we call Minor Epiallele Frequency (MEF), similarly to Minor Allele Frequency filtering often used in population genetics.
This script needs [average_over_bed.py](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/average_over_bed.py).

