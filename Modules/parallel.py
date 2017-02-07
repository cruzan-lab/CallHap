import numpy as np
from functools import partial
from MakeHaplotypes import *
from CallHap_LeastSquares import *
from General import *
import sys

def easyConcat(listHaps):
    return(np.concatenate([x[np.newaxis].transpose() for x in listHaps], axis=1))

def massFindFreqs(inHaps, inSnpFreqs, p):
    mySLSqs = []
    myFreqs = []
    for poolIter in xrange(inSnpFreqs.shape[1]):
        tmpSol = Find_Freqs(inHaps, inSnpFreqs[:, poolIter], p)
        mySLSqs.append(tmpSol[1])
        myFreqs.append(tmpSol[0])
    myAIC = sum(mySLSqs)/len(mySLSqs)
    return(mySLSqs, myFreqs, myAIC)

def massFindFreqs2(inHaps, inSnpFreqs, p):
    mySLSqs = []
    myFreqs = []
    tmpSol = Find_Freqs(inHaps, inSnpFreqs, p)
    mySLSqs.append(tmpSol[1])
    myFreqs.append(tmpSol[0])
    myAIC = sum(mySLSqs)/len(mySLSqs)
    return(mySLSqs, myFreqs, myAIC)

def analyzePotHap(potSolIter, SnpIdx, localPotHapFreqs, AdjPoolCorrSnps, localPotHaplotypes, numHaplotypes, usedCorrSnps):
    newPotSols = []
    newPotHaps = []
    newNumHaplotypes = []
        
    # Figure out which haplotypes this CorrSnp frequency could come from based on frequency
    potSources = []
    for SrcIter in xrange(len(localPotHapFreqs[potSolIter])):
        if localPotHapFreqs[potSolIter][SrcIter] >= AdjPoolCorrSnps[SnpIdx]:
            potSources.append(SrcIter)
    # If no potential sources found, raise an error (for now)
    if potSources == []:
        # This potential solution is now a dead end; skip it
        # if all potential sources have been skipped, the pool is dead and an error needs to be raised about the pool.  
        pass
    else:
        #add command line option to chane number of potential sources to use
        if len(potSources) > 3:
            centrality = FindCentral(localPotHaplotypes[potSolIter])
            ReducedSources = []
            
            modCentrality = [(x+1)/localPotHapFreqs[potSolIter][centrality[x]] if localPotHapFreqs[potSolIter][centrality[x]] > 0 else 9999 for x in xrange(len(centrality))]
            modCentIdx = sorted(range(len(modCentrality)), key = lambda x: modCentrality[x])
            centIter = 0
            while len(ReducedSources) < 3 and centIter < len(modCentIdx):
                if modCentIdx[centIter] in potSources:
                    ReducedSources.append(modCentIdx[centIter])
                centIter += 1
            potSources = ReducedSources
        for source in potSources:
            potTargets = []
            snpsInSource = [x for x in xrange(len(localPotHaplotypes[potSolIter][source])) if localPotHaplotypes[potSolIter][source][x] == 0 ]
            for CorrHapIter in xrange(len(localPotHaplotypes[potSolIter])):
                if localPotHaplotypes[potSolIter][CorrHapIter][SnpIdx] == 0:
                    # Figure out which other corrSnps are in each potential solution
                    toUse = True
                    for SnpIter2 in xrange(len(AdjPoolCorrSnps)):
                        if localPotHaplotypes[potSolIter][CorrHapIter][SnpIter2] == 0:
                            if SnpIter2 == SnpIdx:
                                pass
                            elif SnpIter2 in snpsInSource:
                                pass
                            else:
                                toUse = False
                        elif SnpIter2 in snpsInSource:
                            toUse = False
                    if toUse:
                        potTargets.append(CorrHapIter)
            if len(potTargets) > 1: 
                sys.stderr.write("More than 1 target on SNP %s: %s\n" % (SnpIdx, potTargets))
                sys.stderr.write("Used SNPs = %s\n" % usedCorrSnps)
                sys.stderr.write("Current sol: %s\n" % localPotHapFreqs[potSolIter])
                sys.stderr.write("source hap is %s\n\n" % localPotHaplotypes[potSolIter][source])
                for x in potTargets:
                    sys.stderr.write("%s\n" % localPotHaplotypes[potSolIter][x])
                raise Exception()
            if potTargets == []:
                tmpNewPotSol = localPotHapFreqs[potSolIter][:]
                tmpNewPotSol.append(AdjPoolCorrSnps[SnpIdx])
                tmpNewPotSol[source] -= AdjPoolCorrSnps[SnpIdx]
                tmpNewPotHaps = [np.copy(localPotHaplotypes[potSolIter][x]) for x in xrange(len(localPotHaplotypes[potSolIter]))]
                tmpNewPotHaps.append(np.copy(localPotHaplotypes[potSolIter][source]))
                tmpNewPotHaps[-1][SnpIdx] = 1 - tmpNewPotHaps[-1][SnpIdx]
                if tmpNewPotSol[source] == 0 and source >= numHaplotypes[potSolIter]:
                    tmpNewPotSol.pop(source)
                    tmpNewPotHaps.pop(source)
                newPotSols.append(tmpNewPotSol)
                newPotHaps.append(tmpNewPotHaps)
                newNumHaplotypes.append(numHaplotypes[potSolIter])      
            # For each potential source:
            else:
                for potTargetIter in potTargets:
                    tmpNewPotSol = localPotHapFreqs[potSolIter][:]
                    tmpNewPotSol[potTargetIter] = AdjPoolCorrSnps[SnpIdx]
                    tmpNewPotSol[source] -= AdjPoolCorrSnps[SnpIdx]
                    tmpNewPotHaps = [np.copy(localPotHaplotypes[potSolIter][x]) for x in xrange(len(localPotHaplotypes[potSolIter]))]
                    if tmpNewPotSol[source] == 0 and source >= numHaplotypes[potSolIter]:
                        tmpNewPotSol.pop(source)
                        tmpNewPotHaps.pop(source)
                    ExtantHapSetTster=[sum([numDiffs(tmpNewPotHaps[y], newPotHaps[x][y]) for y in xrange(len(newPotHaps[x]))]) if len(tmpNewPotHaps) == len(newPotHaps[x])  else 10 for x in xrange(len(newPotHaps))]
                    newPotSols.append(tmpNewPotSol)
                    newPotHaps.append(tmpNewPotHaps)
                    newNumHaplotypes.append(numHaplotypes[potSolIter])
        return(newPotSols, newPotHaps, newNumHaplotypes)

def nonParallelCorrHapProcessing(potHapsFreq, adjPoolFreqs, locPotHaps, numHaps, snpIdx, numProcesses, usedCorrSnps):
    sequence = range(len(potHapsFreq))
    
    func = partial(analyzePotHap, SnpIdx = snpIdx, localPotHapFreqs = potHapsFreq, AdjPoolCorrSnps=adjPoolFreqs, localPotHaplotypes = locPotHaps, numHaplotypes=numHaps, usedCorrSnps = usedCorrSnps)
    result = []
    for x in sequence:
        result.append(func(x))
    cleaned = [x for x in result if not x is None]
    outPotSols = []
    outPotHaps = []
    outNumHaps = []   
    for recombIter in xrange(len(cleaned)):
        outPotSols.extend(cleaned[recombIter][0])
        outPotHaps.extend(cleaned[recombIter][1])
        outNumHaps.extend(cleaned[recombIter][2])
        outHapTypes.extend(cleaned[recombIter][3])
    # not optimal but safe
    return(outPotSols, outPotHaps, outNumHaps)

def easy_parallizeCorrHapProcessing(potHapsFreq, adjPoolFreqs, locPotHaps, numHaps, snpIdx, numProcesses, usedCorrSnps):
    from multiprocessing import Pool
    pool = Pool(processes=numProcesses, maxtasksperchild=500)
    sequence = range(len(potHapsFreq))
    
    func = partial(analyzePotHap, SnpIdx = snpIdx, localPotHapFreqs = potHapsFreq, AdjPoolCorrSnps=adjPoolFreqs, localPotHaplotypes = locPotHaps, numHaplotypes=numHaps, usedCorrSnps = usedCorrSnps)
    result = pool.map(func, sequence)
    cleaned = [x for x in result if not x is None]
    pool.close()
    pool.join()
    outPotSols = []
    outPotHaps = []
    outNumHaps = []   
    for recombIter in xrange(len(cleaned)):
        outPotSols.extend(cleaned[recombIter][0])
        outPotHaps.extend(cleaned[recombIter][1])
        outNumHaps.extend(cleaned[recombIter][2])
    # not optimal but safe
    return(outPotSols, outPotHaps, outNumHaps)
    
    
def easy_parallizeLS(sequence, numProcesses, snpsFreqs, poolSize):
    from multiprocessing import Pool
    pool = Pool(processes=numProcesses, maxtasksperchild=500)
    intermediate = pool.map(easyConcat, sequence)
    cleanedIntermediate = [x for x in intermediate if not x is None]
    pool.close()
    pool.join()
    pool2 = Pool(processes=numProcesses, maxtasksperchild=500)

    func = partial(massFindFreqs, inSnpFreqs=snpsFreqs, p=poolSize)
    result = pool2.map(func, cleanedIntermediate)
    cleaned = [x for x in result if not x is None]
    # not optimal but safe
    pool2.close()
    pool2.join()
    return(cleaned)

def parallelCompareHaps(inTypeB, inTypeA, inPotHapSetsList, currType, inTypesList, inNumBaseTypes):
    if inTypesList[inTypeB] == None:
        if len(inPotHapSetsList[inTypeA]) == len(inPotHapSetsList[inTypeB]) and inNumBaseTypes[inTypeA] == inNumBaseTypes[inTypeB]:
            if comparePotHaps(inPotHapSetsList[inTypeA], inPotHapSetsList[inTypeB], inNumBaseTypes[inTypeA]):
                return(currType)
            else:
                return(-1)
        else:
            return(-1)
    else:
        return(-1)
        
def easy_parallizeCompareHaps(numProcesses, inPotHapSetsList, currType, inCompA, inHapTypes, inNumBaseTypes):
    from multiprocessing import Pool
    pool = Pool(processes=numProcesses, maxtasksperchild=500)
    sequence = range(len(inPotHapSetsList))
    func = partial(parallelCompareHaps, inTypeA = inCompA, inPotHapSetsList=inPotHapSetsList, currType=currType, inTypesList=inHapTypes, inNumBaseTypes=inNumBaseTypes)
    result = pool.map(func, sequence)
    cleaned = [x for x in result if not x is None]
    outputHapTypes = [cleaned[x] if cleaned[x] != -1 else inHapTypes[x] for x in sequence]
    # not optimal but safe
    pool.close()
    pool.join()
    return(outputHapTypes)