#!/usr/bin/env python
# CallHap parallel.py
# By Brendan Kohrn
#
# This script includes multiprocessing functionality for deterministic flag in CallHap

import numpy as np
from functools import partial
from Modules.MakeHaplotypes import *
from Modules.CallHap_LeastSquares import *
from General import *
import sys
from multiprocessing import Pool
import time

def easyConcat(listHaps):
    '''One argument command for concatenating a list of arrays into a single array'''
    return(np.concatenate([x[np.newaxis].transpose() 
           for x in listHaps], axis=1))

def massFindFreqs(inHaps, inSnpFreqs, p):
    '''Find frequencies for many pools at the same time'''
    mySLSqs = []
    myFreqs = []
    for poolIter in xrange(inSnpFreqs.shape[1]):
        tmpSol = Find_Freqs(inHaps, inSnpFreqs[:, poolIter], p[poolIter])
        mySLSqs.append(tmpSol[1])
        myFreqs.append(tmpSol[0])
    myAIC = sum(mySLSqs)/len(mySLSqs)
    return(mySLSqs, myFreqs, myAIC)
        
def easy_parallizeLS(sequence, numProcesses, snpsFreqs, poolSize):
    '''Parallelization method for finding the frequencies for several potential
    haplotype sets at the same time'''
    pool = Pool(processes=numProcesses, maxtasksperchild=500)
    intermediate = pool.map(easyConcat, sequence)
    # Concatenate the haplotype sets into a single numpy array
    cleanedIntermediate = [x for x in intermediate if not x is None]
    pool.close()
    pool.join()
    pool2 = Pool(processes=numProcesses, maxtasksperchild=500)
    # Create function for single argument calling of FindFreqs
    func = partial(massFindFreqs, inSnpFreqs=snpsFreqs, p=poolSize)
    # Find haplotype frequencies for each potential haplotype set
    result = pool2.map(func, cleanedIntermediate)
    cleaned = [x for x in result if not x is None]
    pool2.close()
    pool2.join()
    return(cleaned)