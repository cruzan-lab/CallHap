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
echo "Running CallHap Preprocessing using program config $1 and run config $2 at $startDate.  "
for sampleIter in $(seq 0 $(($NumberSamples-1)) ); do
    if [ "${Mode[$sampleIter]}" = "pe" ]; then
        # PE Processing
        echo "Aligning"
        # Alignment with BWA:
        nice -n 5 $bwaPath mem -M ${Refs[$sampleIter]} ${inputRead1[$sampleIter]} ${inputRead2[$sampleIter]} > ${RGSMs[$sampleIter]}_pe.sam
        # Sort the alignment:
        # Also filter out unmapped reads
        echo "Sorting"
        nice -n 5 $samtoolsPath view -Sbu -F 4 ${RGSMs[$sampleIter]}_pe.sam | $samtoolsPath sort - ${RGSMs[$sampleIter]}_pe.sort
        nice -n 5 $samtoolsPath index ${RGSMs[$sampleIter]}_pe.sort.bam
        echo "Getting Statistics"
        # Get alignment statistics:
        nice -n 5 $samtoolsPath flagstat ${RGSMs[$sampleIter]}_pe.sort.bam > ${RGSMs[$sampleIter]}_pe.sort.flagstat.txt
        echo "Starting GATK Steps"
        # Now for the GATK steps
        # Add readgroups to prepare files for GATK:
        echo "Adding Readgroups"
        nice -n 5 java -jar $picardPath AddOrReplaceReadGroups INPUT=${RGSMs[$sampleIter]}_pe.sort.bam OUTPUT=${RGSMs[$sampleIter]}_pe.sort.rg.bam RGID=${RGSMs[$sampleIter]} RGLB=${RGLBs[$sampleIter]} RGPL=$inRGPL RGPU=${RGPUs[$sampleIter]} RGSM=${RGSMs[$sampleIter]}
        $samtoolsPath index ${RGSMs[$sampleIter]}_pe.sort.rg.bam
        # Perform local realingment around indels:
        nice -n 5 java -jar $GATKpath -T RealignerTargetCreator -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_pe.sort.rg.bam -o ${RGSMs[$sampleIter]}_pe.sort.rg.intervals
        nice -n 5 java -jar $GATKpath -T IndelRealigner -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_pe.sort.rg.bam -targetIntervals ${RGSMs[$sampleIter]}_pe.sort.rg.intervals -o ${RGSMs[$sampleIter]}_pe.sort.rg.ra.bam -dt NONE --maxReadsForRealignment 200000
        rm *.intervals
        rm *.rg.bam
        rm *.sort.bam
        rm *.sam
    elif [ "${Mode[$sampleIter]}" = "se" ]; then
    #SE Processing
        echo "Aligning"
        # Alignment with BWA:
        nice -n 5 $bwaPath mem -M ${Refs[$sampleIter]} ${inputRead1[$sampleIter]} > ${RGSMs[$sampleIter]}_se.sam
        echo "Sorting"
        # Sort the alignment:
        # Also filter out unmapped reads
        nice -n 5 $samtoolsPath view -Sbu -F 4 ${RGSMs[$sampleIter]}_se.sam | $samtoolsPath sort - ${RGSMs[sampleIter]}_se.sort
        nice -n 5 $samtoolsPath index ${RGSMs[$sampleIter]}_se.sort.bam
        echo "Getting Statistics"
        # Get alignment statistics:
        nice -n 5 $samtoolsPath flagstat ${RGSMs[$sampleIter]}_se.sort.bam > ${RGSMs[$sampleIter]}_se.sort.flagstat.txt
        # Now for the GATK steps
        # Add readgroups to prepare files for GATK:
        echo "Adding Readgroups"
        nice -n 5 java -jar $picardPath AddOrReplaceReadGroups INPUT=${RGSMs[$sampleIter]}_se.sort.bam OUTPUT=${RGSMs[$sampleIter]}_se.sort.rg.bam RGID=${RGSMs[$sampleIter]} RGLB=${RGLBs[$sampleIter]} RGPL=$inRGPL RGPU=${RGPUs[$sampleIter]} RGSM=${RGSMs[$sampleIter]} CREATE_INDEX=true
        $samtoolsPath index ${RGSMs[$sampleIter]}_se.sort.rg.bam
        # Perform local realingment around indels:
        nice -n 5 java -jar $GATKpath -T RealignerTargetCreator -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_se.sort.rg.bam -o ${RGSMs[$sampleIter]}_se.sort.rg.intervals 
        nice -n 5 java -jar $GATKpath -T IndelRealigner -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_se.sort.rg.bam -targetIntervals ${RGSMs[$sampleIter]}_se.sort.rg.intervals -o ${RGSMs[$sampleIter]}_se.sort.rg.ra.bam -dt NONE --maxReadsForRealignment 200000
        rm *.intervals
        rm *.rg.bam
        rm *.sort.bam
        rm *.sam
    fi
done

endDate=date
echo "Run ran from $startDate to $endDate"
