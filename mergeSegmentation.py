import glob
import sys
from shared import *

def main():
    if len(sys.argv) < 5:
        print "Usage: python mergeSegmentation.py <bed folder> <output file> <segment size> <num states>"
        exit(1)

    bedFolder = formatDir(sys.argv[1])
    oFile = open(sys.argv[2], 'w')
    segmentSize = int(sys.argv[3])
    numStates = sys.argv[4]

    offset = 0
    fileNames = glob.glob(bedFolder + '*_segments.bed')
    maxIdx = 0
    for fileName in fileNames:
        idx = int(fileName.split("_")[-3].split("cell")[-1])
        if idx > maxIdx:
            maxIdx = idx

    lastLine = ""
    for idx in range(0, maxIdx + 1):
        file = bedFolder + "cell" + str(idx) + "_" + numStates + "_segments.bed"
        print "Joining file ", file

        bedFile = open(file, 'r')
        firstLine = True
        for line in bedFile:
            splitLine = line.split("\t")
            splitLine[1] = str(int(splitLine[1]) + offset)
            splitLine[2] = str(int(splitLine[2]) + offset)
            if firstLine:
                if lastLine[0] == splitLine[0]: # if on the same chromosome
                    if lastLine[-1] == splitLine[-1]: #if in the same state
                        oFile.write(lastLine[0] + "\t" + lastLine[1] + "\t" + splitLine[2] + "\t" + splitLine[-1]) # join them
                    else:
                        oFile.write("\t".join(lastLine)) # write last line separately
                        oFile.write("\t".join(splitLine)) # write first line separately
                else:
                    oFile.write("\t").join(splitLine) # write last line separately
                    oFile.write("\t".join(splitLine)) # write first line separately
                firstLine = False
            else:
                oFile.write("\t".join(splitLine)) # if not the first line always just write the line
                lastLine = splitLine
        offset += segmentSize
    oFile.close()

main()
