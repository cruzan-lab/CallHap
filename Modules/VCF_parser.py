import numpy as np
import time

class vcfFile:
    def __init__(self, inFileName, mode, source=""):
        if mode == 'r':
            return(vcfReader(inFileName))
        elif mode == 'w':
            return(vcfWriter(inFileName,source))
            
class vcfWriter:
    '''Call order:
    a = vcfWriter(inName, source)
    a.writeHeader(sampNames)
    a.setFormat(formatStr)
    a.importLinesInfo(Chroms, Poss, Refs, Alts, Quals)
    for sampleName in sampNames:       
        a.importSampleValues(inValues, sampleName)
    a.writeSamples()
    a.close()
    '''
    def __init__(self, inFileName, source):
        self.outputFile = open(inFileName, "wb")
        self.outputFile.write("##fileformat=VCFv4.2\n")
        self.outputFile.write("##fileDate=%s\n" % time.strftime("%Y%m%d"))
        self.outputFile.write("##source=%s\n" % source)
    
    def writeHeader(self, sampleNames):
        self.outputFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(sampleNames))
        self.outputCols = {x:[] for x in sampleNames}
        self.sampleNames = sampleNames
        self.colsFilled = 0
        self.totalCols = len(sampleNames)
        
    def setFormat(self, formatStr):
        self.formatStr = formatStr
        
    def importInfo(self, InfoField, InfoValues):
        if len(InfoValues) != self.numRows:
            raise Exception
        if self.infosSet == True:
            for x in xrange(self.numRows):
                self.infos[x] += ";%s=%s" % (InfoField, InfoValues)
        else:
            self.infos = ["%s=%s" % (InfoField, InfoValues[x]) for x in xrange(self.numRows)]
            self.infosSet = True
        
    def importLinesInfo(self, Chroms, Poss, Refs, Alts, Quals,  IDs = None, Filts = None, Infos = None):
        testLen = len(Chroms)
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
        if len(inValues) != self.numRows + 1:
            print(type(inValues))
            print(len(inValues))
            print(inValues)
            print(self.numRows)
            raise Exception
        self.outputCols[sampleName] = inValues[:-1]
        self.colsFilled += 1
        
    def removeRows(self, rowsToRemove):
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
        if self.colsFilled != self.totalCols:
            raise Exception
        else:
            for lineNum in xrange(self.numRows):
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t
                self.outputFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.chroms[lineNum], self.pos[lineNum], self.IDs[lineNum], self.refs[lineNum], self.alts[lineNum], self.quals[lineNum], self.filts[lineNum], self.infos[lineNum], self.formatStr, "\t".join([str(self.outputCols[self.sampleNames[x]][lineNum]) for x in xrange(len(self.sampleNames))])))
    
    def close(self):
        self.outputFile.close()

class vcfReader:
    def __init__(self, inFileName):
        self.headInfo = {}
        self.lines = []
        inFile = open(inFileName, "rb")
        for line in inFile:
            if line[0:2] == "##":
                wLine = line.strip("#").split("=", 1)
                if "INFO" in wLine[1]:
                    if "INFO" not in self.headInfo:
                        self.headInfo["INFO"] = {}
                    linebins = wLine[1].strip("<>").split(",")
                    self.headInfo["INFO"][linebins[0].split("=")[1]] = {x.split("=")[0]: x.split("=")[1] for x in linebins}

                elif "FORMAT" in wLine[1]:
                    if "FORMAT" not in self.headInfo:
                        self.headInfo["FORMAT"] = {}
                    linebins = wLine[1].strip("<>").split(",")
                    self.headInfo["FORMAT"][linebins[0].split("=")[1]] = {x.split("=")[0]: x.split("=")[1] for x in linebins}
                else:
                    self.headInfo[wLine[0]] = wLine[1]
            elif line[0] == "#":
                self.sampNames = line.strip().split()[9:]
            else:
                self.lines.append(vcfLine(line))
        inFile.close()
    
    def getData(self, target, lineTarget = None, sampTarget = None):
        if target in ("chrom", "pos", "ID", "ref", "alt", "qual", "filt","info", "form"):
            if lineTarget == 'a':
                outList = []
                for line in self.lines:
                    outList.append(line.getData(target))
                return(outList)
    def getNames(self):
        return(self.sampNames)
    
    def asNumpyArray(self, target):
        # This method assumes targeting all rows and columns
        prearray = [x.toNP(target) for x in self.lines]
        return(np.array(prearray))
        
class vcfLine:
    def __init__(self, inLine):
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
            self.data["info"] = {a.split("=")[0]: a.split("=")[1]  for a in linebins[7].split(";")} 
        self.data["form"] = linebins[8].split(":")
        self.data["data"] = [vcfCell(self.data["form"], a) for a in linebins[9:]]
    
    def getData(self, target, sampTarget = None):
        if sampTarget == None and target in self.data.keys():
            return self.data[target]
        elif sampTarget == "a":
            return([self.data["data"][x].getData(target) for x in xrange(len(self.data["data"]))])
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
        return([self.data["data"][x].toNP(target) for x in xrange(len(self.data["data"]))])

    def setElmt(self, target, newValue):
        if target in self.data.keys():
            self.data[target] = newValue
        else:
            raise

class vcfCell:
    def __init__(self, FormatList, inCellText):
        cellbins = inCellText.split(":")
        if cellbins[0] == ".":
            self.data = {x: [np.nan] for x in FormatList}
        else:
            self.data = {}
            for formatIter in xrange(len(FormatList)):
                self.data[FormatList[formatIter]] = [float(x) for x in cellbins[formatIter].split(",")]

    def getData(self, target=None):
        if target == None:
            return(self.data)
        elif target in self.data.keys():
            return(self.data[target][0])
        else:
            raise
    
    def toNP(self, target):
        return(float(self.data[target][0]))

def toNP_array(inFileName, target):
    tmpVCF = vcfReader(inFileName)
    output = tmpVCF.asNumpyArray(target)
    out2 = tmpVCF.getNames()
    return(output, out2)