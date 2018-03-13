#!/bin/python
# CallHap Config Creator
# by Brendan Kohrn


from argparse import ArgumentParser
import time

parser = ArgumentParser()
parser.add_argument(
    "--input", 
    action = 'store', 
    dest='readgroupTemplate', 
    required=True, 
    help="a csv file containing the input files and all information for readgroup creation")
parser.add_argument(
    "--adapt1", 
    action='store', 
    dest='adapt1', 
    required=True, 
    help='Forward adaptor sequence')
parser.add_argument(
    "--adapt2", 
    action='store', 
    dest='adapt2', 
    required=True, 
    help='Referse adaptor sequence')
parser.add_argument(
    "--sequencer", 
    action='store', 
    dest='RGPL', 
    required=True, 
    help="What type of sequencer did this data come from (e.g. Illumina)?")
parser.add_argument(
    "--minBaseQual", 
    action='store', 
    dest='minBaseQual', 
    default="20")
parser.add_argument(
    "--minReadQual", 
    action='store', 
    dest='minReadQual', 
    default="20")
parser.add_argument(
    "--runID", 
    action='store', 
    dest='runID', 
    required=True, help="An identifier for this run; no spaces")
o = parser.parse_args()

print("Running CallHap Config Creator on %s at %s.  " % (time.strftime("%d/%m/%Y"),time.strftime("%H:%M:%S")))
print("Command: python CallHap_ConfigCreator.py --input %s --adapt1 %s --adapt2 %s --sequencer %s --minBaseQual %s --minReadQual %s --runID %s" % (o.readgroupTemplate, o.adapt1, o.adapt2, o.RGPL, o.minBaseQual, o.minReadQual, o.runID))

rgParamsFile = open(o.readgroupTemplate, 'rb')
outFile = open("%s_config.sh" % o.runID, 'wb')
outFile.write("#!/bin/bash\n\n")
outFile.write("# Parameters:\n")
outFile.write("# Mode should be either 'pe' or 'se'\n")

firstLine = True
file1s=[]
file2s=[]
RGLBs=[]
RGSMs=[]
RGPUs=[]
Modes=[]
Refs=[]
numSamples = 0
for line in rgParamsFile:
    if firstLine:
        firstLine = False
    else:
        if '.csv' in o.readgroupTemplate:
            linebins = line.strip().split(',')
        else:
            linebins = line.strip().split()
        file1s.append(linebins[0])
        file2s.append(linebins[1])
        RGLBs.append(linebins[2])
        RGSMs.append(linebins[3])
        RGPUs.append(linebins[4])
        if linebins[5] not in ('se', 'pe'):
            raise IOerror("Mode for sample %s must be either se or pe" % RGSMs[-1])
        Modes.append(linebins[5])
        Refs.append(linebins[6])
        numSamples += 1
rgParamsFile.close()
outFile.write("# How many samples are you processing\n")
outFile.write("NumberSamples=%s\n" % numSamples)
outFile.write("# Raw .fq.gz files, as bash arrays\n")
outFile.write("inputRead1=(%s)\n" % (' '.join(file1s)))
outFile.write("inputRead2=(%s)\n" % (' '.join(file2s)))
outFile.write("# Reference genome to be used for these samples\n")
outFile.write("Refs=(%s)\n" % ' '.join(Refs))
outFile.write("# These arguments are needed for GATK use\n")
outFile.write("# Readgroup Library Names, as bash arrays\n")
outFile.write("RGLBs=(%s)\n" % (' '.join(RGLBs)))
outFile.write("# Readgroup Barcodes, as bash arrays\n")
outFile.write("RGSMs=(%s)\n" % (' '.join(RGSMs)))
outFile.write("# Readgroup Sample Names, as bash arrays\n")
outFile.write("RGPUs=(%s)\n" % (' '.join(RGPUs)))
outFile.write("# Readgroup Sequencing Platform, probably Ilumina\n")
outFile.write("inRGPL=%s\n" % o.RGPL)
outFile.write("Mode=(%s)\n" % (' '.join(Modes)))
outFile.write("minReadQuality=%s\n" % o.minReadQual)
outFile.write("minBaseQuality=%s\n" % o.minBaseQual)
outFile.write("runID=%s\n" % o.runID)
outFile.write("inAdapter1=%s\n" % o.adapt1)
outFile.write("inAdapter2=%s\n" % o.adapt2)
outFile.close()