# popDMRs_refine
Workflow for downstream analysis of EpiDiverse [DMR pipeline](https://github.com/EpiDiverse/dmr) results to merge and refine comparison-specific DMRs from different pairwise-comparisons.

The EpiDiverse [DMR pipeline](https://github.com/EpiDiverse/dmr) is a great tool for calling context-specific (CpG, CHG or CHH) Differentially Methylated Regions between groups of samples. This is a downstream workflow that refines DMRs for large datasets where the pipeline needs to perform a large amount of pairwise comparisons between groups. In includes 1) joining of supposedly "fragmented" DMRs (pairwise-specific short and close DMRs), 2) merging of DMRs from different pairwise comparisons, 3) methylation extraction of the newly obtained DMRs using a unionbed file and final DMR filtering.

WORKFLOW DESCRIPTION: <br/>

[01_join_fragmented_DMRs.py](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/01_join_fragmented_DMRs.py) <br/>
[DMR pipeline](https://github.com/EpiDiverse/dmr) fragmentation acts differently depending on context (CpG, CHG or CHH), sometimes leading to CHH-DMRs being particularly "fragmented". This python script joins all comparison-specific DMRs that are closer than a user-defined distance (default 146bp as fragmentation default in the DMR pipeline) and have the same directionality (higher methylation in the same comparison group). 
DMR length distribution can be compared between the 3 contexts with reg_length_distr_bycont.py and distance between subsequent DMRs with DMR_distance_distr.py.
An enrichment of short and close DMRs as in the figure below suggests excessive "fragmentation".

![image](https://user-images.githubusercontent.com/58292612/121940472-6e80a580-cd4e-11eb-964f-25de4ee85b5e.png)


[02_merge_DMRs.sh](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/02_merge_DMRs.sh) <br/>
Merging of DMRs from different pairwise comparisons in the DMR pipeline directory structure to optain DMRs for the whole dataset.

[03_extr_mergedDMRs_meth.sh](https://github.com/Dario-Galanti/popDMRs_refine/blob/main/03_extr_mergedDMRs_meth.sh) <br/>
Extract methylation from the newly merged DMRs from a unionbed file and filter merged DMRs. The filtering uses a minimum methylation difference and a minimum proportion of samples which need to show the aforementioned difference, similarly to Minor Allele Frequency filtering often used in population genetics.
This script needs [average_over_bed.py](https://github.com/Dario-Galanti/EpiWGBS_downstream/blob/main/average_over_bed.py).

