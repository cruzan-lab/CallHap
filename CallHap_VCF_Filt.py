#!/usr/bin/env python
# CallHap_VCF_Filt.py
# By Brendan F. Kohrn
#
# This is the VCF filter used by the CallHap pipeline
#
# To ensure that they can be used by the main haplotype caller, the following variants are removed:
# A non-SNP variants 
# B variants with low depth or quality
# C variants that do not have a defined identity 
# D SNPs in close proximity to indels
#
# Filters have a limit: depth filter, indel proximity, quality filters that can be modified

import numpy as np
from argparse import ArgumentParser
import time
from Modules.VCF_parser import *
from Modules.IO import *
import sys

parser = ArgumentParser()
parser.add_argument(
    "-i","--inVCF", 
    action="store", 
    dest="inFile", 
    help="The input VCF file to be filtered.  All SSLs should be grouped \
          together in the first columns of the VCF, and all pools grouped \
          together afterwards, as in 'SSL1, SSL2, SSL3, ..., SSLN, Pool1, \
          Pool2, Pool3, ..., PoolM'", 
    required=True)
parser.add_argument(
    "-o", "--outHaps", 
    action="store", 
    dest="outHaps", 
    help="Output haplotypes VCF file name" , 
    required=True)
parser.add_argument(
    "-O", "--outPools", 
    action="store", 
    dest="outPools", 
    help="Output pools VCF file name", 
    required=True)
parser.add_argument(
    "-n", "--numSamps", 
    action="store", 
    dest="numSamps", 
    type=int, 
    help="The number of SSLs in your input VCF file", 
    required=True)
parser.add_argument(
    "-N", "--numPools", 
    action="store", 
    dest="numPools", 
    type=int, 
    help="The number of pools in your input VCF file", 
    required=True)
parser.add_argument(
    "-d", "--minDepth", 
    action="store", 
    dest="minDepth", 
    type=int, 
    help="The minimum depth to process a line, and the minimum average depth \
          to process a column, default: 500", 
    default = 500)
parser.add_argument(
    "--minCallPrev", 
    action="store", 
    dest="minCallPrev", 
    type=float, 
    help="The percentage of reads that must be of a given identity in a SSL, default: 0.9", 
    default=0.9)
parser.add_argument(
    "--minSnpPrev", 
    action="store", 
    dest="minSnpPrev", 
    type=float, 
    help="The percent of a single individual's worth of reads that must be of a \
          given identity to call a position as polymorphic based on pool \
          samples, default: 0.75", 
    default = 0.75)
parser.add_argument(
    "-p", "--poolSizes", 
    action="store", 
    dest="poolSizesFile", 
    help="A file detailing the number of individuals in each pooled library", 
    required=True)
parser.add_argument(
    "-q", "--minQual", 
    action="store", 
    dest="minQual", 
    type=int, 
    help="The minimum quality a given variant call must have to be processed, default: 100",
    default=100)
parser.add_argument(
    "--reportInterval", 
    action="store", 
    dest="rptInt", 
    type=int, 
    help="Report progress at this number of lines, default: 1000", 
    default=1000)
parser.add_argument(
    "--dropLowDepth", 
    action="store", 
    type=int,
    dest="dropLow", 
    help="This should be no greater than the requested minimum depth. Failure \
          to set this flag to something other than the default may result in \
          loss of variants based on incomplete coverage. \
          Default is to keep all columns, default: 0",
    default=0)
parser.add_argument(
    "--indelDist", 
    action="store", 
    dest="indelDist", 
    default=100, 
    help="How far away from indels to make SNPs, default: 100")

o = parser.parse_args()
print("Running CallHap VCF filter on %s at %s" % (time.strftime("%d/%m/%Y"),
                                                  time.strftime("%H:%M:%S")))
pyCommand = "python CallHap_VCF_Filt.py %s" % " ".join(sys.argv[1:])
print("Command = python CallHap_VCF_Filt.py %s" % " ".join(sys.argv[1:]))

print("\nOpening files...")
# Open input VCF file
inVCF = vcfReader(o.inFile)
# Open output files
outHaps = open(o.outHaps, 'wb')
outPools = open(o.outPools, 'wb')
# Write VCF version lines
outHaps.write("\n".join(inVCF.headInfo["headBlock"]))
outHaps.write("##fileDate=%s\n" % time.strftime("%Y%m%d"))
outHaps.write("##source=%s\n" % ("CallHap_VCF_parser"))
outHaps.write("##commandline=\"%s\"" % pyCommand)
outHaps.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
outHaps.write("\n".join(inVCF.headInfo["contig"]))
if len(inVCF.headInfo["contig"]) > 0: 
    outHaps.write("\n")

# Write command into header lines of both output files
outPools.write("\n".join(inVCF.headInfo["headBlock"]))
outPools.write("##fileDate2=%s\n" % time.strftime("%Y%m%d"))
outPools.write("##source2=%s\n" % ("CallHap_VCF_parser"))
outPools.write("##commandline2=\"%s\"" % pyCommand)
#outPools.write("".join(inVCF.headInfo["INFO"]))
outPools.write('##FORMAT=<ID=RF,Number=1,Type=Float,Description="Reference Frequency">\n')
outPools.write("\n".join(inVCF.headInfo["contig"]))
if len(inVCF.headInfo["contig"]) > 0: 
    outPools.write("\n")

# Generate poolSize related numbers:
poolSizes = []
inPoolSizes = open(o.poolSizesFile,"rb")
for line in inPoolSizes:
    poolSizes.append(int(line.strip().split()[1]))
inPoolSizes.close()
if len(poolSizes) != o.numPools:
    print("Error: Number of pool sizes does not match reported number of pools!")
    exit()
snpPrevCutoffs = [o.minSnpPrev/x for x in poolSizes]

# Check average depth of each column in input
print("Checking depth of input columns...")
vcfNames = inVCF.getNames()
# set depths float value from the input file
depths = [0. for x in xrange(o.numSamps + o.numPools)] 
# iterator for lines in input file
lines = 0 
# good depth boolean
goodDepth = [True for x in xrange(o.numSamps + o.numPools)] 
# generates a log of errors 
lineChecker = [] 
# good variant counter 
goodVarCtr = 0 
# counter for the maximum depth in input file
maxDepth = [0. for x in xrange(o.numSamps + o.numPools)] 
# Determine which columns pass the depth filter
for line in inVCF.lines:
    lineDPs = line.getData("DP","a")
    for iter1 in xrange(o.numSamps + o.numPools):
        if np.isnan(float(lineDPs[iter1])) == True:
            depths[iter1] += 0. 
        else:
            depths[iter1] += float(lineDPs[iter1]) 
            if float(lineDPs[iter1]) >= maxDepth[iter1]: 
                maxDepth[iter1] = float(lineDPs[iter1])
    lines += 1
for iter1 in xrange(o.numSamps + o.numPools):
    # compare the average depth of the column to the specified low depth value
    if depths[iter1]/lines >= o.dropLow: 
        goodDepth[iter1] = True
    else:
        goodDepth[iter1] = False 
    if depths[iter1]/lines <= o.minDepth:
        print("\tWarning: Sample %s is has too low of a depth (%s)%s" % 
              (vcfNames[iter1], depths[iter1]/lines, "; skipping" 
              if not goodDepth[iter1] else ""))
# Write column labels to both output files
outHaps.write(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % (
        "\t".join(
            [vcfNames[x] for x in  xrange(0,o.numSamps) 
                if goodDepth[x] == True])))
outPools.write(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % (
        "\t".join([vcfNames[x] for x in  xrange(o.numSamps, o.numSamps+o.numPools) 
                  if goodDepth[x] == True ])))

print("\nStarting analysis of lines...")
# setting counters
lineCounter = 0
indelLocs = [] 
outHapsLines = []
outPoolsLines = []
outLinesPoss = []
for line in inVCF.lines:
    # Print periodic progress reports
    if lineCounter % o.rptInt == 0 and lineCounter > 0:
        print("%s lines processed..." % lineCounter)
        print(line.getData("pos"))
        errors = [x[1] for x in lineChecker]
        print("\nReport:")
        print("%s lines processed" % lineCounter)
        print("%s lines passing initial filters" % 
              [x[0] for x in lineChecker].count(True))
        print("\t%s good variants" % goodVarCtr)
        print("%s lines failing initial filters" % 
              [x[0] for x in lineChecker].count(False))
        print("\t%s incomplete coverage" % 
            errors.count("incomplete coverage"))
        print("\t%s low depth" % errors.count("low depth"))
        print("\t%s too many differences between ref and one alt" % 
            errors.count("too many differences between ref and one alt"))
        print("\t%s different differences between two alts and the ref" % 
            errors.count("Different differences between two alts and the ref"))
        print("\t%s unequal ref and alt lengths" % 
            errors.count("unequal ref and alt lengths"))
        print("\t%s alt longer than 1 with ref length of 1" % 
            errors.count("alt longer than 1 with ref length of 1"))
        print("\t%s Low quality variant call" % 
            errors.count("low quality SNP call"))

    lineCounter += 1
    # Retrieve basic information about this line
    pos = line.getData("pos")
    lineRefCounts = line.getData("RO","a")
    lineDPs = line.getData("DP","a")
    checkDPs = [lineDPs[x] for x in xrange(len(lineDPs)) if goodDepth[x] == True]
    useLine = True
    # Check that there is depth in all samples for this line
    if np.nan in checkDPs:
        lineChecker.append((False, "incomplete coverage", pos))
        useLine = False
    # Check that this line passes the quality filter
    elif float(line.getData("qual")) < o.minQual:
        lineChecker.append((False, "low quality SNP call", pos))
        useLine = False
    # Check that all used columns in this line have adequate depth
    elif False in [
            True if int(checkDPs[x]) >= o.minDepth
            else False for x in xrange(len(checkDPs))
            ]:
        lineChecker.append((False, "low depth", pos))
        useLine = False
    # Check that the length of the reference is 1
    elif len(line.getData("ref")) > 1:
        # If the length of the reference is greater than 1, check if length of alt allele matches length of reference
        # Check if all alt alleles are the same length 
        if (
         max([len(x) for x in line.getData("alt")]) == len(line.getData("ref"))
         and len(set([len(y) for y in line.getData("alt")])) <= 1
         ):
            altValues = line.getData("alt")
            refValue = line.getData("ref")
            newAlt = []
            diffIdxs = []
            
            # If ref and alt alleles are the same length, check how many differences between them
            # If more than one difference, discard line
            for altValue in altValues:
                numDiffs = 0
                diffIdx = []
                for diffCounter in xrange(len(altValue)):
                    if refValue[diffCounter] != altValue[diffCounter]:
                        numDiffs += 1
                        diffIdx.append(diffCounter)
                if numDiffs > 1:
                    if useLine == True:
                        lineChecker.append(
                            (False, "too many differences between ref and one alt", pos))
                    useLine = False
                
                else:
                    # Calculate base pair differences between alt and reference
                    diffIdxs.append(diffIdx[0])
                    newAlt.append(altValue[diffIdx[0]])
            
            # Check that all alts have the same base pair different from the reference
            # If true, adjust the position accordingly
            # Otherwise, discard line
            if len(set(diffIdxs)) == 1:
                line.setElmt("pos", line.getData("pos") + diffIdx[0])
                line.setElmt("ref", refValue[diffIdxs[0]])
                line.setElmt("alt", newAlt)
            else:
                useLine = False
                lineChecker.append(
                    (False, "Different variations between two alternate alleles and the references", pos))
                # Keep track of this location as an indel
                indelLocs.append(int(line.getData("pos")))
        else:
            # If reference and alt are different lengths, discard the line
            useLine = False
            lineChecker.append((False, "unequal ref and alt lengths", pos))
            indelLocs.append(int(line.getData("pos")))
    
    # Check that the length of the alt allele is 1
    # If not, discard line
    elif max([len(x) for x in line.getData("alt")]) > 1:
        lineChecker.append(
        (False, "alt longer than 1 with ref length of 1", pos)
        )
        useLine = False
        # Keep track of this location as an indel
        indelLocs.append(int(line.getData("pos")))
    if useLine == True:
        # If the line passes filtering, count alternate alleles 
        lineChecker.append([True, "good line", pos])
        lineAltCounts = [[x] for x in line.getData("AO", "a")]
        sampIDs = []
        conflicted = False
        monomorphic = False
        
        # Determine the identity (Ref/Alt) of each sample (SSL) in this line
        for sampIter in xrange(o.numSamps):
            if goodDepth[sampIter] == True:
                sampTest = [int(lineRefCounts[sampIter])]
                sampTest.extend([int(x) for x in lineAltCounts[sampIter]])
                maxIter = None
                for iter1 in xrange(len(sampTest)):
                    if maxIter == None:
                        maxIter = iter1
                    elif sampTest[iter1] > sampTest[maxIter]:
                        maxIter = iter1
                sampIDs.append(maxIter)
                # Test if there are no reads in any of the identities for this 
                # sample/line
                if sum(sampTest) == 0:
                    if conflicted == False:
                        lineChecker[-1].append(
                            "conflicted because sum sampTest = 0 (line 167)"
                            )
                    conflicted = True
                # Test if the proportion of the most common identity in this 
                # sample/line is high enough to reliably call
                elif float(sampTest[maxIter])/sum(sampTest) < o.minCallPrev:
                    if conflicted == False:

                        lineChecker[-1].append(
                            "conflicted because of minCallPrev (line 170)"
                            )
                    conflicted = True
        
        # If no conflicts exist, determine reference allele 
        if conflicted == False:
            refAllele = None
            altAllele = None
            for iter1 in [0, 1, 2]:
                if refAllele == None:
                    refAllele = iter1
                elif sampIDs.count(iter1) > sampIDs.count(refAllele):
                    altAllele = refAllele
                    refAllele = iter1
                elif altAllele == None:
                    altAllele = iter1
                elif sampIDs.count(iter1) > sampIDs.count(altAllele):
                    altAllele = iter1
            
            # Calculate reference allele frequency in each pool
            poolFreqs = []
            for iter1 in xrange(o.numPools):
                if goodDepth[o.numSamps + iter1] == True:
                    if refAllele == 0:
                        poolFreqs.append(
                            float(lineRefCounts[o.numSamps + iter1])/
                            int(lineDPs[o.numSamps + iter1])
                            )
                    else:
                        poolFreqs.append(float(
                            lineAltCounts[o.numSamps + iter1][refAllele - 1])/
                            int(lineDPs[o.numSamps + iter1]))
            monomorphicSamps = False
            polymorphicPools = False
            
            # Check if locus is monomorphic in single samples
            if sampIDs.count(refAllele) == len(sampIDs):
                monomorphicSamps = True
                lineChecker[-1].append("monomorphicSamps")
            # Check if sample is polymorphic in pools
            if o.numPools > 0:
                testArr = []
                for freqIter in xrange(len(poolFreqs)):
                    if poolFreqs[freqIter] <= 1. - snpPrevCutoffs[freqIter]:
                        testArr.append(True)
                    else:
                        testArr.append(False)
                if True in testArr:
                    polymorphicPools = True
                    lineChecker[-1].append("polymorphicPools")
            
            # If either SSLs are polymorphic or SSLs are monomorphic and Pools are polymorphic, keep line
            if monomorphicSamps==False or (polymorphicPools == True and 
                                           monomorphicSamps == True):
                goodVarCtr += 1
                lineRef = line.getData(
                    "ref" if refAllele == 0 else "alt")[0 if refAllele == 0 
                                                        else refAllele - 1]
                lineAlt = line.getData(
                    "ref" if altAllele == 0 else "alt")[0 if altAllele == 0 
                                                        else altAllele - 1]
                linePos = line.getData("pos")
                lineChrom = line.getData("chrom")
                lineQual = line.getData("qual")

                outHapsLines.append(
                    "\n%s\t%s\t.\t%s\t%s\t%s\t.\t.\tGT\t%s" % (lineChrom, 
                        linePos, 
                        lineRef, 
                        lineAlt, 
                        lineQual, 
                        "\t".join(["0" if x == refAllele else "1" 
                                  for x in sampIDs]))
                    )
                outPoolsLines.append(
                    "\n%s\t%s\t.\t%s\t%s\t%s\t.\t.\tRF\t%s" % (lineChrom, 
                        linePos, 
                        lineRef, 
                        lineAlt, 
                        lineQual, 
                        "\t".join([str(x) for x in poolFreqs]))
                    )
                outLinesPoss.append(int(line.getData("pos")))

# Check all located variants for proximity to indels
finGoodVars = 0
for outIter in xrange(len(outLinesPoss)):
    currPos = outLinesPoss[outIter]
    useVar = True
    indelIter = 0
    indelPos = 0
    while indelPos <= currPos + o.indelDist and indelIter < len(indelLocs):
        indelPos = indelLocs[indelIter]
        if indelPos > currPos and indelPos - 10 < currPos:
            useVar = False
        elif indelPos < currPos and indelPos + 10 > currPos:
            useVar = False
        elif indelPos == currPos:
            useVar = False
        indelIter += 1
    
    # Write output files
    if useVar == True:
        outHaps.write(outHapsLines[outIter])
        outPools.write(outPoolsLines[outIter])
        finGoodVars += 1

outHaps.close()
outPools.close()

# Create Nexus output
finOutToNex,finOutNames = toNP_array(o.outHaps,"GT")
if finOutToNex.shape[0] != 0:
    NexusWriter(
        [vcfNames[x] for x in  xrange(0,o.numSamps) if goodDepth[x] == True], 
        finOutToNex, finOutToNex.shape[0], o.outHaps[:-4], "", o.outHaps)
    
#output report
errors = [x[1] for x in lineChecker]
print("\nFinal report:")
print("End time = %s %s" % 
    (time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S")))
print("%s lines processed" % lineCounter)
print("%s lines passing initial filters" % 
    [x[0] for x in lineChecker].count(True))
print("\t%s good variants" % goodVarCtr)
print("\t%s not within %s bp of an indel" % (finGoodVars,o.indelDist))
print("%s lines failing initial filters" % 
    [x[0] for x in lineChecker].count(False))
print("\t%s incomplete coverage" % errors.count("incomplete coverage"))
print("\t%s low depth" % errors.count("low depth"))
print("\t%s too many differences between ref and one alt" % 
    errors.count("too many differences between ref and one alt"))
print("\t%s different differences between two alts and the ref" % 
    errors.count("Different differences between two alts and the ref"))
print("\t%s unequal ref and alt lengths" % 
    errors.count("unequal ref and alt lengths"))
print("\t%s alt longer than 1 with ref length of 1" % 
    errors.count("alt longer than 1 with ref length of 1"))
print("\t%s Low quality variant call" % 
    errors.count("low quality SNP call"))
