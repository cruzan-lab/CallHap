#!/usr/bin/env python
# CallHap General.py
# By Brendan Kohrn
# 3/20/2017
#
# This script contains general functions for CallHap

import numpy as np
import math
import decimal as dec

def comparePotHaps(potHapSetA, potHapSetB, numInitialHaps):
    '''Check if all haplotypes in two haplotype sets are the same'''
    # If two haplotype sets are different lengths, they are different
    if len(potHapSetA) != len(potHapSetB):
        return(False)
    else:
        return(all(np.all(x==y) for x,y in zip(potHapSetA[numInitialHaps:], 
                                               potHapSetB[numInitialHaps:])))

def average(inList): 
    ''' Take the average value of a list'''
    # Make sure the list has length
    if len(inList) == 0:
        raise("Error in Average: %s" % inList)
    return(float(sum(inList))/len(inList))

def npDecZeros(rows, cols=0):
    '''Create a numpy array of Decimal(0) values'''
    if cols == 0:
        outArray = np.zeros(rows,dtype=dec.Decimal)
        for rowIter in xrange(rows):
            outArray[rowIter] = dec.Decimal(outArray[rowIter])
    else:
        outArray = np.zeros((rows,cols), dtype=dec.Decimal)
        for rowIter in xrange(rows):
            for colIter in xrange(cols):
                outArray[rowIter,colIter] = dec.Decimal(outArray[rowIter, 
                                                                 colIter])
    
    return(outArray)

def npToDecNp(inArray):
    '''Convert a numpy array of floats to a numpy array of Decimal numbers to 
    avoid rounding errors'''
    outArray = np.array(inArray, dtype=dec.Decimal)
    for elmnt, value in np.ndenumerate(outArray):
        outArray[elmnt] = dec.Decimal(outArray[elmnt])
    return(outArray)

def copy(inArr, elmntType = "int"):
    '''Copy a list (particularly of numpy arrays).'''
    if elmntType == "nparray":
        return([np.copy(x) for x in inArr])
    else:
        return([x for x in inArr])

def AIC_from_RSS(RSS, numHaps, numSNPs):
    '''Calculate AIC from RSS values'''
    AIC = 2 * numHaps + (numSNPs * math.log10(RSS/numSNPs))
    return(AIC)

def AICc_from_RSS(RSS, numHaps, numSNPs):
    '''Calculate AICc from RSS values'''
    AIC = 2 * numHaps + (numSNPs * math.log10(RSS/numSNPs)) + (2*numHaps * 
        (numHaps + 1))/(numSNPs - numHaps - 1)
    return(AIC)

def invertArray(inArray):
    '''Invert an array of 0s and 1s (such as the Haplotypes array) or an array 
    between 0 and 1 (such as the SNP Freqs array).'''
    OutArray = 1 - inArray
    return(OutArray)

def residuals(inSol, inData, inFreqs, poolSize):
    '''Calculate residuals for one particular least-squares solution of Ax=b'''
    calculated = np.sum((inSol * inData)/poolSize, 1)
    resid = np.subtract(inFreqs, calculated)
    return(resid)

def ArrayHaps(origHaps, newHaps):
    allHapsToArray = [origHaps]
    allHapsToArray.extend(newHaps)
    return(np.concatenate(allHapsToArray, axis=1))

def numDiffs(inHap1, inHap2):
    if inHap1.shape != inHap2.shape:
        raise
    else:
        in1 = inHap1.ravel()
        in2 = inHap2.ravel()
        diffCounter = sum([0 if in1[x] == in2[x] else 1 
                          for x in xrange(len(in1)) ])
        return(diffCounter)

def areEqual(inHap1, inHap2):
    if inHap1.shape != inHap2.shape:
        return(False)
    else:
        return(np.all(inHap1 == inHap2))

def FindLegalSnpsByNetwork(inHaps, testHapIdx):
    closestHaps = []
    closestDiffs = []
    notClosest = []
    numSnps = len(inHaps[0])
    distances=[numSnps - np.sum(a==inHaps[testHapIdx]) for a in inHaps]
    # Determine the distance between this haplotype and every other haplotype 
    # in number of SNPs different
    # Sort by closeness
    distIters = sorted(range(len(distances)), key=lambda x: distances[x])
    # For each haplotype, from closest to furthest away, check if it shares a 
    # difference in the target SNP 
    # with another haplotype in closestHaps
    for hapIter in distIters:
        if hapIter != testHapIdx:
            if closestHaps == []:
                # If no haplotype is closest yet, this one is the closest
                closestHaps.append(hapIter)
                closestDiffs.append([])
                notClosest.append([])
                for x in xrange(numSnps):
                    if inHaps[hapIter][x] == inHaps[testHapIdx][x]:
                        pass
                    else:
                        closestDiffs[-1].append(x)
            else:
                # Otherwise, test to see if this haplotype shares a different 
                # SNP with any closer haplotype
                diffBranch = True
                for closeHap in xrange(len(closestHaps)):
                    for difSnp in xrange(len(closestDiffs[closeHap])):
                        tmpPointer = closestDiffs[closeHap][difSnp]
                        if (inHaps[hapIter][tmpPointer] == 
                                inHaps[testHapIdx][tmpPointer]):
                            notClosest[closeHap].append(difSnp)
                        else:
                            diffBranch = False
                if diffBranch:
                    closestHaps.append(hapIter)
                    notClosest.append([])
                closestDiffs.append([])
                for x in xrange(numSnps):
                    if inHaps[hapIter][x] == inHaps[testHapIdx][x]:
                        pass
                    else:
                        closestDiffs[-1].append(x)
        CanChange = []
        for hap in xrange(len(closestHaps)):
            CanChange.extend(closestDiffs[hap])
    return(closestHaps[0],CanChange)
    
def ValidSnpsFromPhylogeny(inHaps):
    countDiffs = [[(a!=b) for a in inHaps] for b in inHaps]
    diffSnps = [[[b for b in xrange(len(countDiffs[x][a])) 
                  if countDiffs[x][a][b] == True] 
                  for a in xrange(len(countDiffs[x]))] 
                  for x in xrange(len(countDiffs))]
    # Find adgacend haplotypes for each haplotype
    validSnps = []
    nextHaps = []
    for hap in xrange(len(inHaps)):
        nextHaps.append([])
        validSnps.append([])
        minDistOrder = sorted(range(len(inHaps)), 
                              key=lambda x: len(diffSnps[hap][x]))
#        print("DEBUG")
        for hap2 in minDistOrder:
            if hap != hap2:

                if len(nextHaps[-1]) > 0:
                    isAdj = True
                    for closeHap in nextHaps[-1]:
                        if len(np.intersect1d(diffSnps[hap][closeHap], 
                                              diffSnps[hap][hap2])) != 0:
                            isAdj = False
                            validSnps[-1] = list(np.setdiff1d(validSnps[-1], 
                                diffSnps[closeHap][hap2]))
                    if isAdj == True:
                        nextHaps[-1].append(hap2)
                        validSnps[-1].extend(diffSnps[hap][hap2])
                else:
                    nextHaps[-1].append(hap2)
                    validSnps[-1].extend(diffSnps[hap][hap2])
    return(validSnps)
    
def DecHapToNPHap(decHap):
    '''Convert a decimal haplotype back into a numpy array'''
    binHap = bin(decHap)[2:]
    binHap = np.array([dec.Decimal(x) for x in binHap[1:]])
    return(binHap)