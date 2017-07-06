# -*- coding: utf-8 -*-
"""
Created on Fri May 19 10:34:18 2017

@author: Brendan Kohrn
"""

#AugmentHaps.py
# Import necessary modules
import numpy as np
from argparse import ArgumentParser
from Modules.VCF_parser import *
from Modules.General import *
from Modules.IO import *

# Inputs: 

parser = ArgumentParser()
parser.add_argument(
    '-i','--inputHaps', 
    action="store", 
    dest="knownHaps", 
    help = "A VCF-formatted file containing the known haplotypes encoded \
            in the GT field.  GT must be present in the FORMAT field, and \
            ploidy must be 1.  ", 
    required=True
    )
parser.add_argument(
    "-n", "--newHaps", 
    action="store", 
    dest="inNewHaps", 
    help="A file containting the decimal identifiers for any new haplotypes \
          to be added, one per line.  ", 
    required=True
    )
parser.add_argument(
    '-o','--outFile', 
    action="store", 
    dest="outFile", 
    required=True, 
    help="A name for the output file.  "
    )
o = parser.parse_args()
#class o:
#    pass
## the initial haplotypes file
#o.knownHaps = "Lascal_DeNovoAlign_d600q20_Haps_Extended.vcf"
## A file containing those haplotypes to be added, one per line
#o.inNewHaps = "NewHapsToAdd.txt"
## A name for the output file
#o.outFile = "Rerun2/WithAddedHaps.vcf"

# Load haplotypes
KnownHaps, KnownNames = toNP_array(o.knownHaps, "GT")
# Invert haplotypes so that ref allele is 1
#KnownHaps = invertArray(KnownHaps)
# Find unique haplotypes
inHapArray = ExtendHaps(KnownHaps)

# Convert to np array
# Read in new haplotypes, new haplotype names
inNewHaps = open(o.inNewHaps, "rb")
newHapNames = []
NewDecHaps = []
NewHaps = []
for iterLine in inNewHaps:
    newHapNames.append("NH_%s" % iterLine.strip())
    NewDecHaps.append(int(iterLine.strip()))
NewHaps = [DecHapToNPHap(NewDecHaps[x]) 
            for x in xrange(len(NewDecHaps))]
# Append new haplotypes to old haplotypes
fullHaps = [np.copy(inHapArray[:,x]) for x in xrange(len(KnownNames))]
fullHaps.extend([invertArray(NewHaps[x].transpose()) for x in xrange(len(NewHaps))])
# Combine names fields
KnownNames.extend(newHapNames)

# Write out haplotypes
output3 = vcfWriter(
    o.outFile, 
    source="")
output3.writeHeader(KnownNames)
output3.setFormat("GT")

tmpVCF = vcfReader(o.knownHaps)
output3.importLinesInfo(
    tmpVCF.getData("chrom", lineTarget="a"),
    tmpVCF.getData("pos", lineTarget="a"), 
    tmpVCF.getData("ref", lineTarget="a"), 
    tmpVCF.getData("alt", lineTarget="a"), 
    tmpVCF.getData("qual", lineTarget="a")
    )
for hapIter in xrange(len(fullHaps)):
    # Add predicted SNP frequencies to VCF output
    output3.importSampleValues(list(fullHaps[hapIter]), KnownNames[hapIter])
output3.writeSamples()
# Close output files
output3.close()