import numpy as np


def checkSnps(inSnps, inNumSnps, inNumPools, inPoolSize, inThreshold = 0.99):
    # Calculate correlation coefficients
    correlations = np.corrcoef(inSnps)
    corrThreshold = inThreshold
    #For each SNP, check if any other SNP is correlated to it
    snpsToQuestion = []
    for snp1 in xrange(inNumSnps):
        isCorrelated = False
        for snp2 in xrange(inNumSnps):
            if snp1 != snp2:
                if correlations[snp1, snp2] >= corrThreshold:
                    isCorrelated = True
        if isCorrelated == False:
            # figure out how many pools it is in
            containingPools = 0
            for pool in xrange(inNumPools):
                if inSnps[snp1,pool] <= 1 - (1/inPoolSize):
                    containingPools += 1
            if containingPools > 1:
                snpsToQuestion.append(snp1)
    return(snpsToQuestion)

def corrGroups(inSNP_Freqs, inThreshold):
    correlations=np.corrcoef(inSNP_Freqs)
    corrGroupings = {}
    corrThreshold = inThreshold
    usedSnps = []
    numSNPs = inSNP_Freqs.shape[0]
    for Snp1 in xrange(numSNPs):
        if Snp1 in usedSnps:
            pass
        else:
            corrGroupings[Snp1] = []
            for Snp2 in xrange(numSNPs):
                if(Snp2 in usedSnps):
                    pass
                elif correlations[Snp1,Snp2] >= corrThreshold:
                    corrGroupings[Snp1].append(Snp2)
                    usedSnps.append(Snp2)
    return([list(corrGroupings[x]) for x in corrGroupings.keys()])

def HapsToCorrHaps(SnpHaps, corGroupsList, numHaps, skip_SNPs = []):
    corrHaps = np.zeros((len(corGroupsList),numHaps))
    for corGroupIter in xrange(len(corGroupsList)):
        if corGroupIter == len(corGroupsList) - 1:
            for hap in xrange(numHaps):
                corrHaps[corGroupIter,hap] = 1
        else:
            for hap in xrange(numHaps):
                for snp in corGroupsList[corGroupIter]:
                    if SnpHaps[snp,hap] == 1 and snp not in skip_SNPs:
                        corrHaps[corGroupIter,hap] = 1
    return(corrHaps)

def SnpsToCorrSnps(SnpsFreqsMatrix, corGroupsList, numPools, snpsToSkip = []):
    corrSnpsFreqs = np.zeros((len(corGroupsList),numPools))
    corrSnpsDevs = np.zeros((len(corGroupsList),numPools))
    for corGroupIter in xrange(len(corGroupsList)):
        if corGroupIter == len(corGroupsList) - 1:
            for poolIter in xrange(numPools):
                corrSnpsFreqs[corGroupIter,poolIter] = 1
                corrSnpsDevs[corGroupIter,poolIter] = 0
        else:
            for poolIter in xrange(numPools):
                corrGroupFreqsList = []
                for snp in corGroupsList[corGroupIter]:
                    if snp not in snpsToSkip:
                        corrGroupFreqsList.append(SnpsFreqsMatrix[snp, poolIter])
                if corrGroupFreqsList != []:
                    corrSnpsFreqs[corGroupIter,poolIter] = np.mean(corrGroupFreqsList)
                    corrSnpsDevs[corGroupIter,poolIter] = np.std(corrGroupFreqsList)
    return(corrSnpsFreqs, corrSnpsDevs)

def CorrHapsToHaps(inCorrHaps, corGroupsList, numSNPs):
    numCorrHaps = inCorrHaps.shape[1]
    outSnpHaps = np.zeros((numSNPs + 1, numCorrHaps))
    for corrGroupIter in xrange(len(corGroupsList)):
        for corrHapIter in xrange(numCorrHaps):
            if inCorrHaps[corrGroupIter, corrHapIter] == 1:
                for snpIter in corGroupsList[corrGroupIter]:
                    outSnpHaps[snpIter, corrHapIter] = 1
    return(outSnpHaps)