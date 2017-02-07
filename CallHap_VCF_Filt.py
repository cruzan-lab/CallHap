import numpy as np
from argparse import ArgumentParser
import time
from Modules.VCF_parser import *
from Modules.IO import *

parser = ArgumentParser()
parser.add_argument("-i","--inVCF", action="store", dest="inFile", help="The input VCF file to be filtered.  All SSLs should be grouped together in the first columns of the VCF, and all pools grouped together afterwards, as in 'SSL1, SSL2, SSL3, ..., SSLN, Pool1, Pool2, Pool3, ..., PoolM'.  ", required=True)
parser.add_argument("-o", "--outHaps", action="store", dest="outHaps", help="A name for the output haplotypes VCF file." , required=True)
parser.add_argument("-O", "--outPools", action="store", dest="outPools", help="A name for the output pools VCF file.  ", required=True)
parser.add_argument("-n", "--numSamps", action="store", dest="numSamps", type=int, help="The number of SSLs in your input VCF file", required=True)
parser.add_argument("-N", "--numPools", action="store", dest="numPools", type=int, help="The number of pools in your input VCF file", required=True)
parser.add_argument("-d", "--minDepth", action="store", dest="minDepth", type=int, help="The minimum depth to process a line, and the minimum average depth to process a column.  ", default = 500)
parser.add_argument("--minCallPrev", action="store", dest="minCallPrev", type=float, help="The percentage of reads that must be of a given identity in a SSL to have that position be good.  ", default=0.9)
parser.add_argument("--minSnpPrev", action="store", dest="minSnpPrev", type=float, help="The percent of a single individuals worth of reads that must be of a given idetity to call a position as polymorphic based on pool samples", default = 0.75)
parser.add_argument("-p", "--poolSize", action="store", dest="poolSize", type=int, help="the number of idividuals in each pooled library.  ", required=True)
parser.add_argument("-q", "--minQual", action="store", dest="minQual", type=int, help="The minimum quality a given varient call must have to be processed.  ", default=100)
parser.add_argument("--reportInterval", action="store", dest="rptInt", type=int, help="Report progress at this number of lines", default=1000)
parser.add_argument("--dropLowDepth", action="store_true", dest="dropLow", help= "Automatically drop any samples with an average depth under the minimum depth.  ")
parser.add_argument("--indelDist", action="store", dest="indelDist", default=100, help="How far away from indels to make SNPs.  Defaults to 100")
o = parser.parse_args()

print("Running CallHap VCF filter on %s at %s" % (time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S")))
pyCommand = "python CallHap_VCF_Filt.py --inVCF %s --outHaps %s --outPools %s --numSamps %s --numPools %s --minDepth %s --minCallPrev %s --minSnpPrev %s --poolSize %s --minQual %s" % (o.inFile, o.outHaps, o.outPools, o.numSamps, o.numPools, o.minDepth, o.minCallPrev, o.minSnpPrev, o.poolSize, o.minQual)

print("Command = %s" % pyCommand)
print("\nOpening files...")
inVCF = vcfReader(o.inFile)
outHaps = open(o.outHaps, 'wb')
outPools = open(o.outPools, 'wb')
outHaps.write("##Command=\"%s\"" % pyCommand)
outPools.write("##Command=\"%s\"" % pyCommand)

# Check average depth of each column in input
print("\nChecking depth of input columns...")
depths = [0. for x in xrange(o.numSamps + o.numPools)]
lines = 0
goodDepth = [True for x in xrange(o.numSamps + o.numPools)]
lineChekcer = []
goodVarCtr = 0
for line in inVCF.lines:
    lineDPs = line.getData("DP","a")
    for iter1 in xrange(o.numSamps + o.numPools):
        if np.isnan(float(lineDPs[iter1])) == True:
            depths[iter1] += 0.
        else:
            depths[iter1] += float(lineDPs[iter1])
    lines += 1
for iter1 in xrange(o.numSamps + o.numPools):
    if depths[iter1]/lines >= o.minDepth: 
        goodDepth[iter1] = True
    else:
        goodDepth[iter1] = False

vcfNames = inVCF.getNames()
if False in goodDepth:
    badColumns = [x for x in range(len(goodDepth)) if goodDepth[x] == False]
    for badIter in badColumns:
        print("\tWarning: Sample %s is has too low of a depth (%s)" % (vcfNames[badIter], depths[badIter]/lines))
        if o.dropLow:
            print("; skipping\n")
        else:
            print("\n")
else:
    print("\tAll columns have greater than minimum average depth.  ")
if not o.dropLow:
    goodDepth = [True for x in xrange(o.numSamps + o.numPools)]

outHaps.write("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % ("\t".join([vcfNames[x] for x in  xrange(0,o.numSamps) if goodDepth[x] == True])))
outPools.write("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % ("\t".join([vcfNames[x] for x in  xrange(o.numSamps, o.numSamps+o.numPools) if goodDepth[x] == True ])))

print("\nStarting analysis of lines...")
lineCounter = 0
indelLocs = []
outHapsLines = []
outPoolsLines = []
outLinesPoss = []
for line in inVCF.lines:
    if lineCounter % o.rptInt == 0 and lineCounter > 0:
        print("%s lines processed..." % lineCounter)
        print(line.getData("pos"))
        badReasons = [x[1] for x in lineChekcer]
        print("\nReport:")
        print("%s lines processed" % lineCounter)
        print("%s lines passing initial filters" % [x[0] for x in lineChekcer].count(True))
        print("\t%s good varients" % goodVarCtr)
        print("%s lines failing initial filters" % [x[0] for x in lineChekcer].count(False))
        print("\t%s incomplete coverage" % badReasons.count("incomplete coverage"))
        print("\t%s low depth" % badReasons.count("low depth"))
        print("\t%s too many differences between ref and one alt" % badReasons.count("too many differences between ref and one alt"))
        print("\t%s different differences between two alts and the ref" % badReasons.count("Different differences between two alts and the ref"))
        print("\t%s unequal ref and alt lengths" % badReasons.count("unequal ref and alt lengths"))
        print("\t%s alt longer than 1 with ref length of 1" % badReasons.count("alt longer than 1 with ref length of 1"))
        print("\t%s Low quality varient call" % badReasons.count("low quality SNP call"))

    lineCounter += 1
    pos = line.getData("pos")
    lineRefCounts = line.getData("RO","a")
    lineDPs = line.getData("DP","a")
    useLine = True
    if np.nan in lineDPs:
        lineChekcer.append((False, "incomplete coverage", pos))
        useLine = False
    elif float(line.getData("qual")) < o.minQual:
        lineChekcer.append((False, "low quality SNP call", pos))
        useLine = False
    elif False in [True if int(lineDPs[x]) >= o.minDepth or goodDepth[x] == False else False for x in xrange(len(lineDPs))]:
        lineChekcer.append((False, "low depth", pos))
        useLine = False
    elif len(line.getData("ref")) > 1:
        if max([len(x) for x in line.getData("alt")]) == len(line.getData("ref")) and len(set([len(y) for y in line.getData("alt")])) <= 1:
            altValues = line.getData("alt")
            refValue = line.getData("ref")
            newAlt = []
            diffIdxs = []
            for altValue in altValues:
                numDiffs = 0
                diffIdx = []
                for diffCounter in xrange(len(altValue)):
                    if refValue[diffCounter] != altValue[diffCounter]:
                        numDiffs += 1
                        diffIdx.append(diffCounter)
                if numDiffs > 1:
                    if useLine == True:
                        lineChekcer.append((False, "too many differences between ref and one alt", pos))
                    useLine = False
                else:
                    diffIdxs.append(diffIdx[0])
                    newAlt.append(altValue[diffIdx[0]])
            if len(set(diffIdxs)) == 1:
                line.setElmt("pos", line.getData("pos") + diffIdx[0])
                line.setElmt("ref", refValue[diffIdxs[0]])
                line.setElmt("alt", newAlt)
            else:
                useLine = False
                lineChekcer.append((False, "Different differences between two alts and the ref", pos))
                indelLocs.append(int(line.getData("pos")))
        else:
            useLine = False
            lineChekcer.append((False, "unequal ref and alt lengths", pos))
            indelLocs.append(int(line.getData("pos")))
    elif max([len(x) for x in line.getData("alt")]) > 1:
        lineChekcer.append((False, "alt longer than 1 with ref length of 1", pos))
        useLine = False
        indelLocs.append(int(line.getData("pos")))
    if useLine == True:
        lineChekcer.append([True, "good line", pos])

        lineAltCounts = [[x] for x in line.getData("AO", "a")]
        sampIDs = []
        conflicted = False
        monomorphic = False
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
                if sum(sampTest) == 0:
                    if conflicted == False:
                        lineChekcer[-1].append("conflicted because sum sampTest = 0 (line 167)")
                    conflicted = True
                elif float(sampTest[maxIter])/sum(sampTest) < o.minCallPrev:
                    if conflicted == False:

                        lineChekcer[-1].append("conflicted because of minCallPrev (line 170)")
                    conflicted = True
        if conflicted == False:
            # Figure out the reference allele
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
            # Calculate ref allele frequency in each pool
            poolFreqs = []
            for iter1 in xrange(o.numPools):
                if goodDepth[o.numSamps + iter1] == True:
                    if refAllele == 0:
                        poolFreqs.append(float(lineRefCounts[o.numSamps + iter1])/int(lineDPs[o.numSamps + iter1]))
                    else:
                        poolFreqs.append(float(lineAltCounts[o.numSamps + iter1][refAllele - 1])/int(lineDPs[o.numSamps + iter1]))
            monomorphicSamps = False
            polymorphicPools = False
            # Check if locus is monomorphic in single samples
            if sampIDs.count(refAllele) == len(sampIDs):
                monomorphicSamps = True
                lineChekcer[-1].append("monomorphicSamps")
            if o.numPools > 0:
                if min(poolFreqs) <= 1. - (o.minSnpPrev/o.poolSize):
                    polymorphicPools = True
                    lineChekcer[-1].append("polymorphicPools")
            if monomorphicSamps==False or (polymorphicPools == True and monomorphicSamps == True):
                goodVarCtr += 1
                lineRef = line.getData("ref" if refAllele == 0 else "alt")[0 if refAllele == 0 else refAllele - 1]
                lineAlt = line.getData("ref" if altAllele == 0 else "alt")[0 if altAllele == 0 else altAllele - 1]
                linePos = line.getData("pos")
                lineChrom = line.getData("chrom")
                lineQual = line.getData("qual")
                
                outHapsLines.append("\n%s\t%s\t.\t%s\t%s\t%s\t.\t.\tGT\t%s" % (lineChrom, linePos, lineRef, lineAlt, lineQual, "\t".join(["0" if x == refAllele else "1" for x in sampIDs])))
                outPoolsLines.append("\n%s\t%s\t.\t%s\t%s\t%s\t.\t.\tRF\t%s" % (lineChrom, linePos, lineRef, lineAlt, lineQual, "\t".join([str(x) for x in poolFreqs])))
                outLinesPoss.append(int(line.getData("pos")))
finGoodVars = 0
for outIter in xrange(len(outLinesPoss)):
    currPos = outLinesPoss[outIter]
    useVar = True
    indelIter = 0
    indelPos = 0
    while indelPos <= currPos + o.indelDist and indelIter < len(indelLocs): # Replace with user defined variable o.readLength
        indelPos = indelLocs[indelIter]
        if indelPos > currPos and indelPos - 10 < currPos:
            useVar = False
        elif indelPos < currPos and indelPos + 10 > currPos:
            useVar = False
        elif indelPos == currPos:
            useVar = False
        indelIter += 1
    if useVar == True:
        outHaps.write(outHapsLines[outIter])
        outPools.write(outPoolsLines[outIter])
        finGoodVars += 1

outHaps.close()
outPools.close()

# Create Nexus output
finOutToNex,finOutNames = toNP_array(o.outHaps,"GT")
if finOutToNex.shape[0] != 0:
    NexusWriter([vcfNames[x] for x in  xrange(0,o.numSamps) if goodDepth[x] == True], finOutToNex,finOutToNex.shape[0], o.outHaps[:-4],"",o.outHaps)
#output report
badReasons = [x[1] for x in lineChekcer]
print("\nFinal report:")
print("End time = %s %s" % (time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S")))
print("%s lines processed" % lineCounter)
print("%s lines passing initial filters" % [x[0] for x in lineChekcer].count(True))
print("\t%s good varients" % goodVarCtr)
print("\t%s not within %s bp of an indel" % (finGoodVars,o.indelDist))
print("%s lines failing initial filters" % [x[0] for x in lineChekcer].count(False))
print("\t%s incomplete coverage" % badReasons.count("incomplete coverage"))
print("\t%s low depth" % badReasons.count("low depth"))
print("\t%s too many differences between ref and one alt" % badReasons.count("too many differences between ref and one alt"))
print("\t%s different differences between two alts and the ref" % badReasons.count("Different differences between two alts and the ref"))
print("\t%s unequal ref and alt lengths" % badReasons.count("unequal ref and alt lengths"))
print("\t%s alt longer than 1 with ref length of 1" % badReasons.count("alt longer than 1 with ref length of 1"))
print("\t%s Low quality varient call" % badReasons.count("low quality SNP call"))
