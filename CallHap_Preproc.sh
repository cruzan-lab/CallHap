#!/bin/bash

#This is a template for all steps prior to the CallHap programs

# Locations of various programs:
# Config file with paths to the programs
source $1
# Run-specific config file (created using config creator)
source $2

# Good sorted header for later header replacement
auxileryHdr=LasburSCircHeadr.sam
startDate=$(date)
echo "Running CallHap Preprocessing using program config $1 and run config $2 at $startDate"
for sampleIter in $(seq 0 $(($NumberSamples-1)) ); do
    if [ "${Mode[$sampleIter]}" = "pe" ]; then
        # PE Processing
        echo "Preprocessing"
        # Adapter trimming:
        nice -n 5 $cutadaptPath -a $inAdapter1 -A $inAdapter2 -o ${RGSMs[$sampleIter]}_R1_at.fastq.gz -p ${RGSMs[$sampleIter]}_R2_at.fastq.gz ${inputRead1[$sampleIter]} ${inputRead2[$sampleIter]}
        # Quality trimming:
        nice -n 5 $sicklePath pe -f ${RGSMs[$sampleIter]}_R1_at.fastq.gz -r ${RGSMs[$sampleIter]}_R2_at.fastq.gz -o ${RGSMs[$sampleIter]}_R1_at_qt.fastq.gz -p ${RGSMs[$sampleIter]}_R2_at_qt.fastq.gz -t sanger -s ${RGSMs[$sampleIter]}_extras.fastq.gz -q $minBaseQuality -g
        echo "Aligning"
        # Alignment with BWA:
        nice -n 5 $bwaPath mem -M ${Refs[$sampleIter]} ${RGSMs[$sampleIter]}_R1_at_qt.fastq.gz ${RGSMs[$sampleIter]}_R2_at_qt.fastq.gz > ${RGSMs[$sampleIter]}_pe.sam
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
        # Perform local realignment around indels:
        nice -n 5 java -jar $GATKpath -T RealignerTargetCreator -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_pe.sort.rg.bam -o ${RGSMs[$sampleIter]}_pe.sort.rg.intervals
        nice -n 5 java -jar $GATKpath -T IndelRealigner -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_pe.sort.rg.bam -targetIntervals ${RGSMs[$sampleIter]}_pe.sort.rg.intervals -o ${RGSMs[$sampleIter]}_pe.sort.rg.ra.bam -dt NONE --maxReadsForRealignment 200000
        rm *.fastq.gz
        rm *.intervals
        rm *.rg.bam
        rm *.sort.bam
        rm *.sam
    elif [ "${Mode[$sampleIter]}" = "se" ]; then
    #SE Processing
        echo "Preprocessing"
        # Adapter trimming with cutadapt:
        nice -n 5 $cutadaptPath -a $inAdapter1 -o ${RGSMs[$sampleIter]}_at.fastq.gz ${inputRead1[$sampleIter]}
        # Quality trimming with sickle:
        nice -n 5 $sicklePath se -f ${RGSMs[$sampleIter]}_at.fastq.gz -o ${RGSMs[$sampleIter]}_at_qt.fastq.gz -t sanger -q $minBaseQuality -g
        echo "Aligning"
        # Alignment with BWA:
        nice -n 5 $bwaPath mem -M ${Refs[$sampleIter]} ${RGSMs[$sampleIter]}_R1_at_qt.fastq.gz > ${RGSMs[$sampleIter]}_se.sam
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
        # Perform local realignment around indels:
        nice -n 5 java -jar $GATKpath -T RealignerTargetCreator -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_se.sort.rg.bam -o ${RGSMs[$sampleIter]}_se.sort.rg.intervals 
        nice -n 5 java -jar $GATKpath -T IndelRealigner -R ${Refs[$sampleIter]} -I ${RGSMs[$sampleIter]}_se.sort.rg.bam -targetIntervals ${RGSMs[$sampleIter]}_se.sort.rg.intervals -o ${RGSMs[$sampleIter]}_se.sort.rg.ra.bam -dt NONE --maxReadsForRealignment 200000
        rm *.fastq.gz
        rm *.intervals
        rm *.rg.bam
        rm *.sort.bam
        rm *.sam
    fi
done

endDate=$(date)
echo "Start time: $startDate to End time: $endDate"
