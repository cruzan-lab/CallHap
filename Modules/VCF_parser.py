#!/usr/bin/env python
# CallHap IO.py
# By Brendan Kohrn
# 3/21/2017
#
# This is the VCF parser used by all CallHap specific programs
import numpy as np
import time

class vcfFile:
    def __init__(self, inFileName, mode, source=""):
        if mode == 'r':
            return(vcfReader(inFileName))
        elif mode == 'w':
            return(vcfWriter(inFileName,source))
            
class vcfWriter:
    '''Class to write VCF output based on an imput template.  
    Call order:
    a = vcfWriter(inName, source)
    a.writeHeader(sampNames)
    a.setFormat(formatStr)
    a.importLinesInfo(Chroms, Poss, Refs, Alts, Quals)
    for sampleName in sampNames:       
        a.importSampleValues(inValues, sampleName)
    a.writeSamples()
    a.close()
    '''
    def __init__(self, inFileName, source, commandLine, baseHead, FormatBlock):
        '''Initialize the class'''
        # Open an output file
        self.outputFile = open(inFileName, "wb")
        # Write header lines        
        self.outputFile.write("##fileformat=VCFv4.1\n")
        self.outputFile.write("".join(baseHead["headBlock"]))
        self.outputFile.write("##fileDate=%s\n" % time.strftime("%Y%m%d"))
        self.outputFile.write("##source=%s\n" % ("CallHap_VCF_parser"))
        self.outputFile.write("##commandline=\"%s\"" % commandLine)
        self.outputFile.write("\n".join(baseHead["INFO"]))
        self.outputFile.write("\n".join(FormatBlock))
        self.outputFile.write("\n".join(baseHead["contig"]))
        if len(baseHead["contig"]) != 0: self.outputFile.write("\n")
    
    def writeHeader(self, sampleNames):
        
        # Write column nanes line
        self.outputFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\t")
        self.outputFile.write("FILTER\tINFO\tFORMAT\t")
        self.outputFile.write("%s\n" % "\t".join(sampleNames))
        # Create output columns
        self.outputCols = {x:[] for x in sampleNames}
        # Save sample names
        self.sampleNames = sampleNames
        # Create counter to keep track of how many columns have been filled
        self.colsFilled = 0
        # How many columns need to be filled
        self.totalCols = len(sampleNames)
        
    def setFormat(self, formatStr):
        # Set the format string for outputs.
        self.formatStr = formatStr
        
    def importInfo(self, InfoField, InfoValues):
        '''Add text to the info field'''
        # Check that there is an info value for each row
        if len(InfoValues) != self.numRows:
            raise Exception
        # Check if there is already info data present
        if self.infosSet == True:
            # If info data is present, add to it
            for x in xrange(self.numRows):
                self.infos[x] += ";%s=%s" % (InfoField, InfoValues)
        else:
            # Create info data
            self.infos = ["%s=%s" % (InfoField, InfoValues[x]) 
                         for x in xrange(self.numRows)]
            self.infosSet = True
        
    def importLinesInfo(self, Chroms, Poss, Refs, Alts, Quals,
                        IDs = None, Filts = None, Infos = None):
        '''Add positional and quality information about the lines'''
        testLen = len(Chroms)
        # Check that all lists of values are the same length
        if len(Poss) != testLen:
            raise Exception
        elif len(Refs) != testLen:
            raise Exception
        elif len(Alts) != testLen:
            raise Exception
        elif len(Quals) != testLen:
            raise Exception
        elif IDs != None:
            if len(IDs) != testLen:
                raise Exception
        elif Filts != None:
            if len(Filts) != testLen:
                raise Exception
        elif Infos != None:
            if len(Infos) != testLen:
                raise Exception
        self.numRows = testLen
        self.chroms = Chroms
        self.pos = Poss
        if IDs == None:
            self.ID_Set = False
            self.IDs = ['.' for x in xrange(len(Chroms))]
        else:
            self.ID_Set = True
            self.IDs = IDs
        self.refs = Refs
        self.alts = [x[0] for x in Alts]
        self.quals = Quals
        if Filts == None:
            self.filts = ['.' for x in xrange(len(Chroms))]
        else:
            self.filts = Filts
        if Infos == None:
            self.infosSet = False
            self.infos = ['.' for x in xrange(len(Chroms))]
        else:            
            self.infos = Infos
    
    def importSampleValues(self, inValues, sampleName):
        '''Import cell data for one column of a VCF file'''
        # Old debugging text
        if len(inValues) != self.numRows + 1:
            print(type(inValues))
            print(len(inValues))
            print(inValues)
            print(self.numRows)
            raise Exception
        # Fill the column
        self.outputCols[sampleName] = inValues[:-1]
        self.colsFilled += 1
        
    def removeRows(self, rowsToRemove):
        '''Remove rows from the VCF output'''
        removalRows = sorted(rowsToRemove, reverse = True)
        for rowIter in removalRows:
            self.chroms.pop(rowIter)
            self.pos.pop(rowIter)
            self.IDs.pop(rowIter)
            self.refs.pop(rowIter)
            self.alts.pop(rowIter)
            self.quals.pop(rowIter)
            self.filts.pop(rowIter)
            self.infos.pop(rowIter)
            self.numRows -= 1
        self.skippedRows = rowsToRemove
    
    def writeSamples(self):
        '''Write the VCF output to file'''
        # Throw an error if not all columns have been filled
        if self.colsFilled != self.totalCols:
            raise Exception
        else:
            for lineNum in xrange(self.numRows):
                outLine = "%s\t" % self.chroms[lineNum]
                outLine += "%s\t" % self.pos[lineNum]
                outLine += "%s\t" % self.IDs[lineNum]
                outLine += "%s\t" % self.refs[lineNum]
                outLine += "%s\t" % self.alts[lineNum]
                outLine += "%s\t" % self.quals[lineNum]
                outLine += "%s\t" % self.filts[lineNum]
                outLine += "%s\t" % self.infos[lineNum]
                outLine += "%s\t" % self.formatStr
                outLine += "%s\n" % "\t".join(
                    [str(self.outputCols[self.sampleNames[x]][lineNum])
                     for x in xrange(len(self.sampleNames))]
                     )
                self.outputFile.write(outLine)
    
    def close(self):
        '''Close the output file'''
        self.outputFile.close()

class vcfReader:
    '''Method for reading VCF files'''
    def __init__(self, inFileName):
        '''Initialize the reader and read the file'''
        self.headInfo = {}
        self.headInfo["headBlock"] = []
        self.headInfo["FORMAT"] = []
        self.headInfo["contig"] = []
        self.headInfo["INFO"] = []
        self.lines = []
        # Open the file
        inFile = open(inFileName, "rb")
        for line in inFile:
            line = line.strip()
            # Check if this is a header line
            if line == "":
                pass
            elif line[0:2] == "##":
                # Parse the header line, in case that information is needed 
                # later
                
                wLine = line.strip("#").split("=", 1)
                if "INFO" in wLine[0]:
                    if "INFO" not in self.headInfo:
                        self.headInfo["INFO"] = []
                    linebins = wLine[1].strip("<>").split(",")
                    self.headInfo["INFO"].append(line)
                elif "FORMAT" in wLine[0]:
                    if "FORMAT" not in self.headInfo:
                        self.headInfo["FORMAT"] = []
                    self.headInfo["FORMAT"].append(line)
                elif "contig" in wLine[0]:
                    if "contig" not in self.headInfo:
                        self.headInfo["contig"] = []
                    self.headInfo["contig"].append(line)
                elif "INFO" not in self.headInfo:
                    self.headInfo["headBlock"].append(line)    
                elif wLine[0] not in self.headInfo:
                    self.headInfo[wLine[0]] = [line]
                elif wLine[0] in self.headInfo:
                    self.headInfo[wLine[0]].append(line)
                    
            # Check if this is the column labels line
            elif line[0] == "#":
                # Save the sample names
                self.sampNames = line.strip().split()[9:]
            elif line[0] != "":
                # Create a new VCF line with the data in this line
                self.lines.append(vcfLine(line))
        if "contig" not in self.headInfo:
            self.headInfo["contig"] = []
        # Close the input file
        inFile.close()
    
    def getData(self, target, lineTarget = None, sampTarget = None):
        '''Retrieve data about the VCF file from a specific line or information
         column'''
        if target in ("chrom", "pos", "ID", "ref", "alt", 
                      "qual", "filt","info", "form"):
            if lineTarget == 'a':
                outList = []
                for line in self.lines:
                    outList.append(line.getData(target))
                return(outList)
                
    def getNames(self):
        '''Retrieve the column names'''
        return(self.sampNames)
    
    def asNumpyArray(self, target):
        '''Output the data from this file as a numpy array'''
        # This method assumes targeting all rows and columns
        prearray = [x.toNP(target) for x in self.lines]
        return(np.array(prearray))
        
class vcfLine:
    '''Handler class for a single line of a VCF file'''
    def __init__(self, inLine):
        '''Initialize the class and read in the data'''
        linebins = inLine.split()
        self.data = {}
        self.data["chrom"] = linebins[0]
        self.data["pos"] = int(linebins[1])
        self.data["ID"] = linebins[2]
        self.data["ref"] = linebins[3]
        self.data["alt"] = linebins[4].split(",")
        self.data["qual"] = linebins[5]
        self.data["filt"] = linebins[6]
        if linebins[7] == ".":
            self.data["info"] = {}
        else:
            self.data["info"] = {a.split("=")[0]: a.split("=")[1] 
                                 for a in linebins[7].split(";")} 
        self.data["form"] = linebins[8].split(":")
        # Save the data from each column within this row as a VCF cell class
        self.data["data"] = [vcfCell(self.data["form"], a) 
                             for a in linebins[9:]]
    
    def getData(self, target, sampTarget = None):
        '''Retrieve data from this line'''
        if sampTarget == None and target in self.data.keys():
            return self.data[target]
        elif sampTarget == "a":
            return([self.data["data"][x].getData(target) 
                    for x in xrange(len(self.data["data"]))])
        elif "info" in sampTarget:
            if ":" in sampTarget:
                infoTarget = sampTarget.split(":")[1]
                if infoTarget in self.data["info"].keys():
                    return(self.data["info"][infoTarget])
                else:
                    raise
            else:
                return(self.data["info"])
        elif sampTarget < len(self.data["data"]):
            return(self.data["data"][sampTarget].getData(target))
        else:
            raise
    
    def toNP(self, target):
        '''Get the data from this row to for numpy array creation'''
        return([self.data["data"][x].toNP(target) 
                for x in xrange(len(self.data["data"]))])

    def setElmt(self, target, newValue):
        '''Set a specific value in the data from this row'''
        if target in self.data.keys():
            self.data[target] = newValue
        else:
            raise

class vcfCell:
    '''Class for holding data from a single cell of a VCF file'''
    def __init__(self, FormatList, inCellText):
        '''Create the cell'''
        cellbins = inCellText.split(":")
        if cellbins[0] == ".":
            self.data = {x: [np.nan] for x in FormatList}
        else:
            self.data = {}
            for formatIter in xrange(len(FormatList)):
                if FormatList[formatIter] == "GT":
                    self.data[FormatList[formatIter]] = [float(x) 
                        for x in cellbins[formatIter].split("/")]
                else:
                    self.data[FormatList[formatIter]] = [float(x) 
                        for x in cellbins[formatIter].split(",")]

    def getData(self, target=None):
        '''Retrieve data from the cell'''
        if target == None:
            return(self.data)
        elif target in self.data.keys():
            return(self.data[target][0])
        else:
            raise
    
    def toNP(self, target):
        '''Get data from the cell to create a numpy array'''
        return(float(self.data[target][0]))

def toNP_array(inFileName, target):
    '''Open a VCF file and create a numpy array of a specific type of data from
     that array, along with the sample names'''
    tmpVCF = vcfReader(inFileName)
    output = tmpVCF.asNumpyArray(target)
    out2 = tmpVCF.getNames()
    return(output, out2)
