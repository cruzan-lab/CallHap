import numpy as np
from VCF_parser import *

def ExtendHaps(origHaps):
    allHapsToArray = [origHaps]
    allHapsToArray.extend([np.array([[1 for x in xrange(int(origHaps.shape[1]))]])])
    return(np.concatenate(allHapsToArray, axis=0))

def UniqueHaps(inHaps, inNames):
    '''Find unique haplotypes and reduce the haplotypes and their names; 
    still needs some work to fix merged names'''
    remove = [False for n in range(int(inHaps.shape[1]))]
    for iterx in range(int(inHaps.shape[1])-1):
        for y in range(iterx+1, inHaps.shape[1]):
            if not remove[y]:
                if np.ma.all(inHaps[:,iterx] == inHaps[:,y]):
                    remove[y] = True
    toRemove = [iterx for iterx in range(len(remove)) if remove[iterx]]
    names = [inNames[iterx] for iterx in range(len(inNames)) if not remove[iterx]]
    return(np.delete(inHaps, toRemove, 1), names)

def outputProt(UniqueNames, bestFreqs, bestArray, poolSize, poolNames, population, outFile):
    decHaps = []
    hapNames = UniqueNames[:]
    newHapNumber = 1
    for haplotypeIter in xrange(len(bestFreqs[0])):
        if haplotypeIter >= len(UniqueNames):
            hapNames.append("NewHap_%s" % str(newHapNumber).zfill(2))
            newHapNumber += 1
        decHaps.append(int("1"+"".join([str(int(x)) for x in bestArray[:, haplotypeIter]]),2))
    indivHaps = []
    for haplotypeIter in xrange(len(bestFreqs[0])):
        for indivIter in xrange(int(bestFreqs[0][haplotypeIter])):
            
            indivHaps.append(decHaps[haplotypeIter])
    for individual in xrange(poolSize):
        outLine = " ".join(["%s_%s" % (poolNames[population], str(individual).zfill(len(str(poolSize)))), str(population),str(indivHaps[individual])])
        outFile.write("%s\n" % outLine)

def NexusWriter(myHapNames, finSolution, numSNPs, outPrefix, outIdx, knownHaps, snpsToRemove=[]):
    outFile3 = open("%s_%s_haps.nex" % (outPrefix, outIdx), 'wb')
    outFile3.write("##NEXUS\n")
    outFile3.write("Begin Data;\n")
    outFile3.write("\tDimensions ntax=%s nchar=%s;\n" % (finSolution.shape[1], numSNPs))
    outFile3.write("\tFormat datatype=DNA missing=N gap=-;\n")
    outFile3.write("\tMatrix\n")
    
    finVCF = vcfReader(knownHaps)
    refAlleles = []
    altAlleles = []
    for line in finVCF.lines:
        if line.getData("pos") not in snpsToRemove:
            refAlleles.append(line.getData("ref"))
            altAlleles.append(line.getData("alt")[0])
    
    for hap in xrange(len(myHapNames)):
        outFile3.write("%s\t%s\n" % (myHapNames[hap], "".join([refAlleles[x] if finSolution[x, hap] == 1 else altAlleles[x] for x in xrange(numSNPs)])))
    outFile3.write(";\n")
    outFile3.write("End;\n")
    outFile3.close()