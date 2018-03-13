import numpy as np

def  MakeHaps(inSnpSets, inPoolSize, inOldHaps, inInitialFreqs):
    # 
    possibleFreqs = [inInitialFreqs[:]]
    possibleHaps = [inOldHaps]
    initialHaps = len(inOldHaps)
    
    freqSet = 0
    testStop = len(possibleFreqs)
    loopCtr1 = 0
    
    while freqSet < testStop:
        loopCtr1 += 1
#       if loopCtr1 > 100000:
#           print(testStop)
#           raise Exception("Too many iterations at line 331")
        baseFreq = []
        for freq in xrange(len(possibleFreqs[freqSet])):
            if possibleFreqs[freqSet][freq] > 0:
                baseFreq.append(freq)
        
        newFreq = 0
        loopCtr2 = 0
        while newFreq < len(baseFreq):
            loopCtr2 += 1
            if loopCtr2 > 1000:
                raise Exception("Too many iterations at line 342 with baseFreq = %s" % len(baseFreq))
            if baseFreq[newFreq] > initialHaps:
                if newFreq == len(baseFreq) - 1:
                    # Change the origional frequency set and haplotypes set
                    possibleFreqs[freqSet].append(1)
                    possibleHaps[freqSet].append(np.copy(possibleHaps[freqSet][baseFreq[newFreq]]))
                    for iter1 in inSnpSets:
                        possibleHaps[freqSet][-1][iter1] = 1 - possibleHaps[freqSet][-1][iter1]
                else:
                    # make a copy of the origional frequency set and haplotypes set
                    possibleFreqs.append([x for x in possibleFreqs[freqSet]])
                    possibleHaps.append([np.copy(x) for x in possibleHaps[freqSet]])
                    # change the copy
                    possibleFreqs[-1].append(1)
                    possibleHaps[-1].append(np.copy(possibleHaps[freqSet][newFreq]))
                    for iter1 in inSnpSets:
                        possibleHaps[-1][-1][iter1] = 1 - possibleHaps[-1][-1][iter1]
            else:

                if newFreq == len(baseFreq) - 1:
                    # Change the origional frequency set and haplotypes set
                    possibleFreqs[freqSet].append(1)
#                                possibleFreqs[freqSet][baseFreq[newFreq]] -= key
                    possibleHaps[freqSet].append(np.copy(possibleHaps[freqSet][baseFreq[newFreq]]))
                    for iter1 in inSnpSets:
                        possibleHaps[freqSet][-1][iter1] = 1 - possibleHaps[freqSet][-1][iter1]
                    if possibleFreqs[freqSet][baseFreq[newFreq]] == 0:
                        possibleFreqs[freqSet].pop(baseFreq[newFreq])
                        possibleHaps[freqSet].pop(baseFreq[newFreq])
                else:
                    # make a copy of the origional frequency set and haplotypes set
                    possibleFreqs.append([x for x in possibleFreqs[freqSet]])
                    possibleHaps.append([np.copy(x) for x in possibleHaps[freqSet]])
                    # change the copy
                    possibleFreqs[-1].append(1)
                    possibleHaps[-1].append(np.copy(possibleHaps[freqSet][baseFreq[newFreq]]))
                    for iter1 in inSnpSets:
                        possibleHaps[-1][-1][iter1] = 1 - possibleHaps[-1][-1][iter1]
                    if possibleFreqs[freqSet][baseFreq[newFreq]] == 0:
                        possibleFreqs[freqSet].pop(baseFreq[newFreq])
                        possibleHaps[freqSet].pop(baseFreq[newFreq])
            newFreq += 1
        freqSet += 1
    return(possibleHaps)