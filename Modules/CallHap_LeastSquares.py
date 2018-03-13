#!/usr/bin/env python
# CallHap CallHap_LeastSquares.py
# By Brendan Kohrn
#
# The main Sum Least Squares method for CallHap_HapCallr

import numpy as np
import decimal as dec
from General import *


def Find_Freqs(A, b, p):
    '''Find the frequency of various haplotypes in a pool.  A is the haplotypes 
    matrix, b is the SNP Frequency matrix, and p is pool size'''
    # Set variables for number of haplotypes and number of SNPs
    M = A.shape[0]
    N = A.shape[1]
    # Create an empty numpy array for the current solution
    x = npDecZeros(1, N)
    # Create an empty numpy array to hold the last solution
    lastX = npDecZeros(1, N)
    # Create dummy variables to hold the last and current sum squared residuals
    currentSSR = -1
    lastSSR = -1
    
    # Run the first test to determine the best starting haplotype
    currentSSR, x = InitialTest(A, b, x, currentSSR, M, N, p)
    lastSSR = currentSSR
    lastX = np.copy(x)
    # Create finished switch and counter to check for infinite loops
    finished = False
    # Iterations:
    while not finished:
        # invoke the mail loop
        currentSSR, x = mainLoop(A, b, x, N, M, currentSSR, p)
        # If the SSR (Sum Squared Residuals; equivalent to RSS) value increases on this loop, finish
        if currentSSR >= lastSSR:
            finished = True
        else:
            lastSSR = currentSSR
            lastX = np.copy(x)

    # output frequencies are contained in lastX
    # output SSR contained in lastSSR
    return(lastX, lastSSR)

def InitialTest(A, b, xIT, curSSR, M, N, p):
    # set counter for best sum squared residuals
    bestSSR = -1
    # Check each haplotype
    for hapIndex in xrange(N):
        # Calculate the SSR if this haplotype was the only one in the pool
        testSSR = sum([resid**2 for resid in np.subtract(A[:, hapIndex], b)])
        # Check if this SSR is an improvement
        if bestSSR == -1:
            bestSSR = [hapIndex, testSSR]
        elif testSSR < bestSSR[1]:
            bestSSR[0] = hapIndex
            bestSSR[1] = testSSR
    xIT[0, bestSSR[0]] = dec.Decimal(p)
    # Return SSR value and best haplotype frequency vector
    return(bestSSR[1], np.copy(xIT))

def SSR(A, xSSR, b):
    '''Calculate the sum of squared residuals for a given solution to Ax=b'''
    if type(b) == list:
        out = sum([resid**2 for resid in np.subtract(np.sum(A * xSSR, 1), b)])
    else:
        out = sum([resid**2 for resid in np.subtract(np.sum(A * xSSR, 1), 
                                                     b.ravel())])
    return(out)

def mainLoop(A, b, xML, N, M, currSSR, p):
    # for each element of x s.t. x[x_1] > 1, subtract 1 from that element 
    # and add one to each other element (x_2) in turn; 
    bestSSR = [(-1,-1),currSSR]
    for x_1 in xrange(N):
        if xML[0, x_1] >= 1:
            for x_2 in xrange(N):
                wx = np.copy(xML)
                if x_1 != x_2:
                    wx[0, x_1] -= 1
                    wx[0, x_2] += 1
                    testSSR = SSR(A, wx / p, b)
                    if testSSR < bestSSR[1]:
                        bestSSR[0] = (x_1, x_2)
                        bestSSR[1] = testSSR
    if bestSSR[0] != (-1,-1):
        xML[0,bestSSR[0][0]] -= 1
        xML[0,bestSSR[0][1]] += 1
    return(bestSSR[1], np.copy(xML))