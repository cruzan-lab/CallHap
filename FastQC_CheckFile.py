# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:26:23 2016

@author: Brendan Kohrn
"""

import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-a', '--adapterThreshold', type=float, help="% of reads that can have an adapter and still have the file pass")

o = parser.parse_args()

mode = None
for line in sys.stdin:
    linebins = line.strip().split('\t')
    if mode == None:
        if linebins[0] == ">>Adaoter Content":
            mode = "AC"
    else:
        if linebins[0] == ">>END_MODULE":
            mode = None
    if mode == "AC":
        if float(linebins[1]) >= o.adapterThreshold:
            
    if linebins[1] == "Adapter Content":
        if linebins[0] != "PASS":
            sys.stderr.write("Warning: Sample %s failed FastQC Adapter Content test\n" % linebins[3])
    elif linebins[1] == "Per Base Sequence Quality":
        if linebins[0] != "PASS":
            sys.stderr.write("Warning: Sample %s failed FastQC Per Base Sequence Quality test\n" % linebins[3])
