import os

def printHelp(helpType):
    if helpType == "open": # repeat for each command
        # print help for open command
        print("Open a vcf file for editing.  ")
        print("Syntax:\n open file=\"FileLocation\"")
    elif helpType == "printHeader": # repeat for each command
        # print help for printHeader command
        print("Print the header of an open VCF file.  If there is no open VCF file, raises an error.  ")
        print("Syntax:\n printHeader")
    elif helpType == "removeCols": # repeat for each command
        # print help for removeCols command
        # syntax: removeCols "0,1,2,3" 
        print("Removes columns from an open VCF file.  Takes a comma seperated list of columns to remove.  If there is no open VCF file, or the column reference is outside the number of columns, gives an error")
        print("Syntax:\n removeCols \"col1,col2,col3,...,colN\"")
    elif helpType == "write": # repeat for each command
        # write file="fileName"
        # print help for write command
        print("Write VCF output to a specified file.  ")
        print("Syntax:\n write file=\"FileLocation\"")
    elif helpType == "quit": # repeat for each command
        # print help for quit command
        print("Exit the program.  ")
        print("Syntax:\n quit")
    elif helpType == "":
        print("Welcome to Brendan's Basic Interactive VCF editor, a tool for removing VCF columns.  ")
        print("Order of operations:")
        print("Open file")
        print("Check header")
        print("Remove columns")
        print("Write file")
        print("Quit")
        print("")
        print("Commands: ")
        print("help")
        print("open")
        print("printHeader")
        print("removeCols")
        print("write")
        print("quit")
        print("Type \"help command\" for more help")
    else:
        print("Error: %s is not a valid command" % helpType)
        
def main():
    headerLines = []
    inputLines = []
    modifiedLines = []
    modified = False
    header=""
    fileOpen = False
    
    # Set up loop to wait for user input
    done = False
    print("Welcome to Brendan's Basic Interactive VCF editor, a tool for removing VCF columns.  ")
    print("Order of operations:")
    print("Open file")
    print("Check header")
    print("Remove columns")
    print("Write file")
    while not done:
        test=raw_input("> ")
        
        testbins = test.split(' ', 1)
        print(testbins)
        # Check user input
        if testbins[0] == "help": # get help
            # Check for presence of aditional option
            if len(testbins) == 1:
                # print list of commands
                printHelp("")
            else:
                printHelp(testbins[1])
        elif testbins[0] == "open": # open a new vcf
            # takes 2 arguments, as folows
            # open file="file"
            fileArgs = testbins[1].split('=', 1)
            if len(fileArgs) < 2:
                print("Error: malformed open statement")
            elif  fileArgs[0] != "file":
                print("Error: Misformatted open statement")
            else: # Correctly formated open statement
                # Check if file exists
                filePath = fileArgs[1].strip("\"")
                if os.path.isfile(filePath):
                    # Reset internal variables
                    headerLines = []
                    inputLines = []
                    modifiedLines = []
                    modified = False
                    header=""
                    # File exists; open it
                    myVcfFile = open(filePath, 'rb')
                    for line in myVcfFile:
                        if line[:2] == "##":
                            headerLines.append(line.strip())
                        elif line[0] == "#":
                            header = line.strip().split()
                        else:
                            inputLines.append(line.strip())
                    myVcfFile.close()
                    fileOpen = True
                    print("File opened")
                else:
                    print("Error: File %s not found." % filePath)
        elif testbins[0] == "printHeader":
            if not fileOpen:
                print("Error: No file open.  ")
            else:
                for colNum in xrange(len(header)):
                    print("%s: %s" % (colNum, header[colNum]))
        elif testbins[0] == "removeCols": # remove columns
            # syntax: removeCols "0,1,2,3"                
            # Check command formating
            colsToRemove = sorted([int(x) for x in testbins[1].strip('"').split(',')], reverse=True)
            # Check if file is open
            if not fileOpen:
                print("Error: No file open.  ")
            else:
                toCont = True
                for colIdx in colsToRemove:
                    if colIdx => len(header):
                        print("Error: one or more columns specified for removal not present.  Check the header.  ")
                        toCont = False
                if toCont:
                    # Check if file has already been modified
                    if modified:
                        tmpLines = modifiedLines[:]
                        modifiedLines = []
                    else:
                        tmpLines = inputLines[:]
                    
                    for line in tmpLines:
                        modLine = line.split()
                        for colIdx in colsToRemove:
                            del(modLine[colIdx])
                        modifiedLines.append("\t".join(modLine))
                    for colIdx in colsToRemove:
                        del(header[colIdx])
                    modified = True
        elif testbins[0] == "write": # write out file
            # syntax: write file="fileName"
            if not fileOpen:
                print("Error: No file open.  ")
            else:
                # check syntax
                fileArgs = testbins[1].split('=', 1)
                if len(fileArgs) < 2:
                    print("Error: malformed open statement")
                elif  fileArgs[0] != "file":
                    print("Error: Misformatted write statement")
                else:
                    continueWrite = False
                    # test if output file exists
                    filePath = fileArgs[1].strip("\"")
                    if os.path.isfile(filePath):
                        goodAnswer = False
                        while not goodAnswer:
                            replaceFile = raw_input("File exists.  Do you want to replace it? [y/n] ")
                            if replaceFile == "y" or replaceFile == "n":
                                goodAnswer = True
                            else:
                                print("Please type either \"y\" or \"n\".  ")
                        if replaceFile == "y":
                            continueWrite = True
                        else:
                            print("Error: File Exists.  ")
                    else:
                        continueWrite = True
                    if continueWrite:
                        # open output file
                        outFile = open(filePath, 'wb')
                        firstLine = True
                        for line in headerLines:
                            if firstLine:
                                outFile.write("%s" % line)
                                firstLine = False
                            else:
                                outFile.write("\n%s" % line)
                        outFile.write("\n%s" % ("\t".join(header)))
                        for line in modifiedLines:
                            outFile.write("\n%s" % line)
                        outFile.close()
        elif testbins[0] == "quit":
                done = True
        else:
            print("Error: invalid command.  ")
          
if __name__ == "__main__":
    main()