#!/usr/bin/python3
## Author: Dario Galanti, Feb 2021
## Input: 3 bed files with DMRs (1 per context). Either merged_DMRs or DMRs from a specific pw comparison
## Hypothesis: If "fragmented_DMRs" were called instead of whole DMRs, there should be a strong over-representation of really close DMRs
## Aim: Plot distribution of distance between DMRs for all context.
## Run: python3 DMR_distance_distr.py CpG_bedfile (NB: there needs to be CHG and CHH corresponding files with same dir structure) plotname.png
## Run: sbatch --partition test --cpus-per-task 2 --mem 10G --time 02:00:00 --wrap "python3 DMR_distance_distr.py merged_DMRs/CpG/DMRs_CpG_merged.bed pop_merged_DMRs_distance.png"

## NB: There is no need to differentiate between different scaffolds since we will filter out really high or negative distances!

## Matplotlib code
## https://www.machinelearningplus.com/plots/matplotlib-histogram-python-examples/
## https://matplotlib.org/3.1.1/gallery/statistics/hist.html

## Dependencies: conda activate Python (matplotlib and seaborn)

## Import modules
import sys
import matplotlib.pyplot as plt				# To draw SD histogram
from matplotlib import colors
import seaborn as sns

## Input & output files
fin_CpG = str(sys.argv[1])		# NB: there needs to be CHG and CHH corresponding files with same dir structure
fin_CHG = fin_CpG.replace("CpG","CHG")
fin_CHH = fin_CpG.replace("CpG","CHH")
fout = str(sys.argv[2])

## Define function to extract length distribution list.
def reg_bed_dist(fin):
	## Store lengths into list
	global dist_list
	dist_list = []
	end = 1000000 #Whatever high value will discard the first line
	with open(fin, "r") as fin:
		for line in fin:
			line = line.strip()
			line = line.split("\t")
			dist = (int(line[1])-end)
			if 0 < dist < 5000:			# We only look at DMRs which are less than 5kb apart
				dist_list.append(dist)
				#if dist < 500: print(dist,line[0],line[1], sep=" ")
			end = int(line[2])

## Make plot
sns.set_style("white")
kwargs = dict(hist_kws={'alpha':.4}, kde_kws={'linewidth':2})	#The lower alpha the more transparent and blury
plt.figure(figsize=(10,7), dpi= 100)

binsize=10
#plot CHH
reg_bed_dist(fin_CHH)
bins=int(max(dist_list)/binsize)
sns.distplot(dist_list, color="#1E90FF", label="CHH", **kwargs, bins=bins, kde=False)
# plot CpG
reg_bed_dist(fin_CpG)
bins=int(max(dist_list)/binsize)
sns.distplot(dist_list, color="#CD5B45", label="CpG", **kwargs, bins=bins, kde=False)
#plot CHG
reg_bed_dist(fin_CHG)
bins=int(max(dist_list)/binsize)
sns.distplot(dist_list, color="#EEC900", label="CHG", **kwargs, bins=bins, kde=False)

plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('Distance between successive DMRs', fontsize=18)
plt.ylabel('Density', fontsize=18)
plt.xlim(0,1000)				# Zooming in	
#plt.ylim(0,0.08)	
plt.legend(fontsize=18)
plt.savefig(fout)

