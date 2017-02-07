#!/bin/bash

#This is a template for all steps prior to the CallHap programs

# Version 0.01.22
# 
# Change Log:
#     V0.01.13
#          Raised realignment limit for indel_realigner
#          Cleaned up unnecessary lines
#     V0.01.12
#          Removed samtools reheader as it seemed to be caushing errors and may no longer be necessary.
#     V0.01.11
#          Moved reference identification to being defined per-sample
#          Matches ConfigCreator_0.1.2
#     V0.01.10
#          Removed deduplication to preserve depth in pools.  
#          Also changed program to be able to run se and pe files at the same time for future convenience.
#             This feature may not be necessary, but for the purposes of running current data sets, it will remain implemented.  

# Locations of various programs:
# Config file with paths to the programs
source $1
# Run-specific config file (created using config creator)
source $2

# Good sorted header for later header replacement
# This will eventually be in the config file.
auxileryHdr=LasburSCircHeadr.sam
startDate=date
echo "Running CallHap Preprocessing 1(Trimming) using program config $1 and run config $2 at $startDate.  "
for sampleIter in $(seq 0 $(($NumberSamples-1)) ); do
    if [ "${Mode[$sampleIter]}" = "pe" ]; then
        # PE Processing
        echo "Preprocessing"
        # Adapter trimming:
        nice -n 5 $cutadaptPath -a $inAdapter1 -A $inAdapter2 -o ${RGSMs[$sampleIter]}_R1_at.fastq.gz -p ${RGSMs[$sampleIter]}_R2_at.fastq.gz ${inputRead1[$sampleIter]} ${inputRead2[$sampleIter]}
        # Quality trimming:
        nice -n 5 $sicklePath pe -f ${RGSMs[$sampleIter]}_R1_at.fastq.gz -r ${RGSMs[$sampleIter]}_R2_at.fastq.gz -o ${RGSMs[$sampleIter]}_R1_at_qt.fastq.gz -p ${RGSMs[$sampleIter]}_R2_at_qt.fastq.gz -t sanger -s ${RGSMs[$sampleIter]}_extras.fastq.gz -q $minBaseQuality -g
        echo "Aligning"
    elif [ "${Mode[$sampleIter]}" = "se" ]; then
    #SE Processing
        echo "Preprocessing"
        # Adapter trimming with cutadapt:
        nice -n 5 $cutadaptPath -a $inAdapter1 -o ${RGSMs[$sampleIter]}_at.fastq.gz ${inputRead1[$sampleIter]}
        # Quality trimming with sickle:
        nice -n 5 $sicklePath se -f ${RGSMs[$sampleIter]}_at.fastq.gz -o ${RGSMs[$sampleIter]}_at_qt.fastq.gz -t sanger -q $minBaseQuality -g
    fi
done
$FastqcPah *_at.fastq.gz
$FastqcPath *_at_qt.fastq.gz
echo "Run ran from $startDate to $endDate"