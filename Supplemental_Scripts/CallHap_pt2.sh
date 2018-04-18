#!/bin/bash 

#continue naming convention from the CallHap_pt1.sh for your output files.
samplename=
#Pools text file
#list all pools in one column separated by a "tab" and enter the number of indvls in the pool (usually 20)
pools=/path/to/pools/file/{pools}.txt

#reference genome
reference=/path/to/reference/genome/

#number of SSLs
SSL=
#number of Pools
PLs=
#quality score used
q=


#Program locations
VCF_filt=/vol/share/cruzan_lab/bioinformatics/programs/CallHap/CallHap-master_110617v2/CallHap_VCF_Filt.py
##CallHap=/vol/share/cruzan_lab/bioinformatics/programs/CallHap/CallHap-master/CallHap_HapCallr.py

#SNP Filtering at 100 depth incriments.
#it is recommended to check the depth of library and pool reads before increasing depth parameter to avoid dropping samples 
#quality score is currently set to 0.
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d300q"$q"_Haps.vcf -O "$samplename"_d300q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 300 -q "$q" -p "$pools" --dropLowDepth 100
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d400q"$q"_Haps.vcf -O "$samplename"_d400q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 400 -q "$q" -p "$pools" --dropLowDepth 100
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d500q"$q"_Haps.vcf -O "$samplename"_d500q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 500 -q "$q" -p "$pools" --dropLowDepth 100
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d600q"$q"_Haps.vcf -O "$samplename"_d600q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 600 -q "$q" -p "$pools" --dropLowDepth 100
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d700q"$q"_Haps.vcf -O "$samplename"_d700q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 700 -q "$q" -p "$pools" --dropLowDepth 100
python "$VCF_filt" -i "$samplename".vcf -o "$samplename"_d800q"$q"_Haps.vcf -O "$samplename"_d800q"$q"_Pools.vcf -n "$SSL" -N "$PLs" -d 800 -q "$q" -p "$pools" --dropLowDepth 100

#CallHap
#decide where you want to run this line, as well with which depth.
##python "$CallHap" --inputHaps "$samplename"_d600q20_Haps.vcf --inputFreqs "$samplename"_d400q20_Pools.vcf -o "$samplename" -p 20 -t 5 -l 2 --numRandom 100 –numTopRSS 3