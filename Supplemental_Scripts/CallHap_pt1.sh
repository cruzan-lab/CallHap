#!/bin/bash 

#assign unique name for all output files relating to this run
samplename=

#reference genome
#indexed with dictionary
reference=/path/to/reference/{reference_genome}.fasta

#assign where you want to be working out of
workingDirectory=

#Adapter Sequences (check with FastQC if you're not sure)
a1=AGATCGGAAGAG
A2=AGATCGGAAGAG

#(1 / total number of PLs)
minAltFrac=0.03

#Program locations
program_config=/path/to/CallHap/program_config.sh
ConfigCreator=/path/to/CallHap/CallHap_ConfigCreator.py
Preprocessing=/path/to/CallHap/CallHap_Preproc.sh
FreeBayes=/path/to/freeBayes/freebayes


#Preprocessing
python "$ConfigCreator" --input "$samplename".csv --adapt1 "$a1" --adapt2 "$A2" --sequencer IlluminaHighSeq --minBaseQual 30 --minReadQual 30 --runID "$samplename"

#run specific samplename_config.sh is created
bash "$Preprocessing" "$program_config" "$samplename"_config.sh

#SNP Calling
ls -1 "$workingDirectory"/*SSL*.rg.ra.bam > "$samplename".txt
ls -1 "$workingDirectory"/*Pool*.rg.ra.bam >> "$samplename".txt

#FreeBayes
"$FreeBayes" -L "$samplename".txt -p 1 -f "$reference" -v "$samplename".vcf --use-best-n-alleles 2 --min-repeat-entropy 1 --no-partial-observations --min-alternate-fraction "$minAltFrac"