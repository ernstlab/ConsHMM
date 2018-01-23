from __future__ import print_function
import glob
import sys
import gzip
import argparse
from shared import *

def main():
    parser = argparse.ArgumentParser(description='Merge ConsHMM segmentation chunks into one segmentation.', allow_abbrev=False)
    parser.add_argument('-b', '--b')
    if len(sys.argv) < 6:
        print("Usage: python mergeSegmentation.py <bed folder> <output file> <segment size> <num states> <suffix include>")
        exit(1)

    bedFolder = formatDir(sys.argv[1])
    oFile = gzip.open(sys.argv[2], 'w')
    segmentSize = int(sys.argv[3])
    numStates = sys.argv[4]
    suffix = sys.argv[5]

    offset = 0
    fileNames = glob.glob(bedFolder + '*_' + suffix + '_segments.bed')
    maxIdx = 0
    for fileName in fileNames:
        idx = int(fileName.split("_")[-4].split("cell")[-1])
        if idx > maxIdx:
            maxIdx = idx

    lastLine = ""
    for idx in range(0, maxIdx + 1):
        print("idx = ", idx)
        file = bedFolder + "cell" + str(idx) + "_" + numStates + "_" + suffix + "_segments.bed"
        print("Joining file ", file)

        bedFile = open(file, 'r')
        firstLine = True
        for line in bedFile:
            splitLine = line.split("\t")
            splitLine[1] = str(int(splitLine[1]) + offset)
            splitLine[2] = str(int(splitLine[2]) + offset)
            
            if lastLine != "":
                if lastLine[0] == splitLine[0]: # if on the same chromosome
                    if lastLine[-1] == splitLine[-1]: # if in the same state
                        lastLine[2] = splitLine[2]
                    else:
                        oFile.write("\t".join(lastLine))
                        lastLine = splitLine
                else:
                    oFile.write("\t".join(lastLine))
                    lastLine = splitLine
            else:
                lastLine = splitLine
        offset += segmentSize

    oFile.write("\t".join(lastLine))
    oFile.close()

main()
