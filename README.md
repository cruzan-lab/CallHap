# CallHap
A bioinformatics pipeline for processing of population-level pooled chloroplast DNA

Last Update Date: 03/13/2018
By: Elizabeth Hendrickson

CHANGE LOG: 

V0.02.01
	add --randHighFreq flag
		Generate multiple output files by randomly generating a specified number of highest frequency SNPs before building haplotypes
		Intended use of confirming deterministic output
	Misc. language/grammar updates; code cleaning and removed unnecessary functions
	Add allele count output file to VCF_Filt.py
	Update new haplotype names
	
V0.01.24
	Fixed date/time reporting

V0.01.13
 	Raised realignment limit for indel_realigner
 	Cleaned up unnecessary lines

V0.01.12
	Removed samtools reheader as it seemed to be crashing errors and may no longer be necessary

V0.01.11
	Moved reference identification to being defined per-sample
	Matches ConfigCreator_0.1.2

V0.01.10
	Removed de-duplication to preserve depth in pools 
	Also changed program to be able to run se and pe files at the same time for future convenience
	This feature may not be necessary, but for the purposes of running current data sets, it will remain implemented
