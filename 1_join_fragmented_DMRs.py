#!/usr/bin/python3
## Author: Dario Galanti, Feb 2021
## Input: DMR pipe output file
## Aim: Merge DMRs closer then "x" bp that have the same directionality (both hyper in the same group)
## Output: bed file with same format as input but some DMRs merged
## Run: python3 1_join_fragmented_DMRs.py input.bed output.bed
## Run: for f in pop_DMRs/CpG/R*.bed;do fout=$(echo $f | rev | cut -c5- | rev)_joined_146bp.bed; python3 1_join_fragmented_DMRs.py $f $fout;done

## Purpose: Especially in CHH (higher density) DMRs are sometimes detected in smaller chuncks instead of detecting full DMRs
## IMPORTANT: Bare in mind that the DMR stats from this script are not very accurate because they don't make use of the original unionbed files!!!
## Cs is only obtained summing up Cs from merged DMRs but does not count any Cs in the region between the merged DMRs
## The average difference suffers from the same imprecision.
## The merged_DMR q-value is taken from the original DMR with the minimum q-value
## For these reasons, these stats are not used for further analysis.

## Import modules
import sys

## Define input, output and min dist to merge DMRs
fin = str(sys.argv[1])
fout = str(sys.argv[2])
index = "/scr/epi/genomes/thlaspi_arvense/thlaspi.fa.fai"
min_dist = 146												## IMPORTANT!!! Defines the minimum distance to merge DMRs with same directionality


line_num = 0
with open(fin, "r") as fin, open(fout, "w") as fout:
	for line in fin:
		line = line.strip()
		line = line.split("\t")
		line_num += 1
	## Skip first line
		if (line_num > 1):
		## Calculate distance from previous DMR and directionality of new DMR
			dist = (int(line[1]) - end)
			new_direction = ("-" if float(line[4]) < 0 else "+")
		## NEW LINE SHOULD BE MERGED TO PREVIOS BLOCK!!! Update few info
			if (line[0] == scaff) and (dist < min_dist) and (direction == new_direction):
				Cs += int(line[3])
				diff += float(line[4]) * int(line[3])
				sig = min(sig, float(line[5]))
		## NEW LINE SHOULD BE SEPARATED FROM PREVIOUS BLOCK!!! Print previous block and then update info
			else:
			## First print previous block
				diff = (diff/Cs)
				length = (end - st)
				print(scaff,st,end,Cs,diff,sig,length, sep="\t", file=fout)
			## Then update start and other info
				scaff = line[0]
				st = int(line[1])
				Cs = int(line[3])
				diff = float(line[4]) * int(line[3])
				direction = ("-" if float(line[4]) < 0 else "+")
				sig = float(line[5])
				
	## First line: only store DMR info
		else:
			scaff = line[0]
			st = int(line[1])
			Cs = int(line[3])
			diff = float(line[4]) * int(line[3])
			direction = ("-" if float(line[4]) < 0 else "+")
			sig = float(line[5])
		## The original DMR end is the only info that has to be updated every line!!!
		end = int(line[2])
	


