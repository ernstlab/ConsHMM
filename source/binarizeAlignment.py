import os
import sys
import gzip
import time
import random
import numpy as np
from sklearn import linear_model, datasets, metrics
from collections import defaultdict
from shared import *

def createFeatures(splitLine, posReference):
    features = np.empty((len(splitLine) - 2) * 2) # two features for each species

    human = splitLine[posReference] # variable for keeping the reference nucleotide; does not have to be human
    idx = 0
    for i in range(1, len(splitLine)):
        if (i != posReference):
            if splitLine[i] == 'X':
                features[idx] = 0
                features[idx + 1] = 2
            else:
                features[idx] = 1
                if splitLine[i] == human:
                    features[idx + 1] = 1
                else:
                    features[idx + 1] = 0
            idx += 2
    
    return int(splitLine[0]), features

def main():
    if len(sys.argv) < 8:
        print "Usage: binarizeAlignment.py <MAF sequence file> <output directory> <chromosome> <lines per chunk> <chromosome lengths file> <reference species> <file list>"
        exit(1)

    alignmentFile = gzip.open(sys.argv[1], 'r')
    oDir = formatDir(sys.argv[2])
    chr = sys.argv[3]
    linesPerChunk = int(sys.argv[4])
    chromSizes = open(sys.argv[5], 'r')
    refSpecies = sys.argv[6]
    outputFileList = open(sys.argv[7], 'w')

    chrLength = {}
    for line in chromSizes:
        splitLine = line.strip().split()
        chrLength[splitLine[0]] = int(splitLine[1])

    totalGenomeLength = sum(chrLength.values())

    # Read file headers 
    header = alignmentFile.readline()
    headerSplit = header.strip().split(",")
    newHeader = ""
    for species in headerSplit[1:]:
        if (species != refSpecies):
            if (newHeader != ""):
                newHeader += "\t"
            newHeader = newHeader + species + "_aligned\t" + species + "_matched"
    
    posReference = headerSplit.index(refSpecies)
    numFeats = len(header.split(",")) - 2

    startTime = time.time()
    curChunk = 0

    print "Creating output file . . . ",
    ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')
    outputFileList.write(chr + "_" + str(curChunk) + "_features_binary.txt.gz\n")
    ofile.write("cell" + str(curChunk) + "\t" + chr + "\n")
    ofile.write(newHeader + "\n")
    print "Done."

    lastPos = -1
    for line in alignmentFile:
        splitLine = line.strip().split(",")
        pos, features = createFeatures(splitLine, posReference)

        if pos != (lastPos + 1):
            for i in range(lastPos + 1, pos):
                if (i != 0) and (i % linesPerChunk == 0): # check to see if we should start a new chunk
                    ofile.close()
                    curChunk += 1
                    ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')
                    ofile.write("cell" + str(curChunk) + "\t" + chr + "\n")
                    ofile.write(newHeader + "\n")

                extraFeatures = [0, 2] * numFeats
                ofile.write("\t".join([str(k) for k in extraFeatures]))
                ofile.write("\n")

        if (pos != 0) and (pos % linesPerChunk == 0):
            ofile.close()
            curChunk += 1
            ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')
            outputFileList.write(chr + "_" + str(curChunk) + "_features_binary.txt.gz\n")
            ofile.write("cell" + str(curChunk) + "\t" + chr + "\n")
            ofile.write(newHeader + "\n")

        for k in range(len(features)):
            ofile.write(str(int(features[k])) + "\t")
        ofile.write("\n")
        lastPos = pos

    if lastPos != chrLength[chr]:
        for i in range(lastPos + 1, chrLength[chr]):
            if (i % linesPerChunk == 0):
                ofile.close()
                curChunk += 1
                ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')
                ofile.write("cell" + str(curChunk) + "\t" + chr + "\n")
                ofile.write(newHeader + "\n")

            extraFeatures = [0, 2] * numFeats
            ofile.write("\t".join([str(k) for k in extraFeatures]))
            ofile.write("\n")

    outputFileList.close()
    ofile.close()

    endTime = time.time()
    print "Done.  Time: ", endTime - startTime, " seconds."

main()
