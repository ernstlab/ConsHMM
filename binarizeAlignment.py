import os
import sys
import gzip
import time
import random
import numpy as np
from sklearn import linear_model, datasets, metrics
from collections import defaultdict
from shared import *

def createFeatures(splitLine):
    features = np.empty((len(splitLine) - 2) * 2) # two features for each species

    human = splitLine[1]
    for i in range(2, len(splitLine)):
        idx = (i - 2) * 2
        if splitLine[i] == 'X':
            features[idx] = 0
            features[idx + 1] = 2
        else:
            features[idx] = 1
            if splitLine[i] == human:
                features[idx + 1] = 1
            else:
                features[idx + 1] = 0

    return int(splitLine[0]), features

def main():
    if len(sys.argv) < 6:
        print "Usage: binarizeAlignment.py <MAF features directory> <output directory> <alignment name> <chromosome> <lines per chunk>"
        exit(1)

    featuresDir = formatDir(sys.argv[1])
    oDir = formatDir(sys.argv[2])
    alignmentName = sys.argv[3]
    chr = sys.argv[4]
    linesPerChunk = int(sys.argv[5])

    chrLength = {"chr1":249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260,
            "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747,
            "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
            "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
            "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}
    totalGenomeLength = sum(chrLength.values())

    if alignmentName == "UCSC":
        numFeats = 99
    elif alignmentName == "Ensembl_EPO":
        numFeats = 36
    elif alignmentName == "Ensembl_Pecan":
        numFeats = 20
    else:
        print "Alignment name not recognized."
        exit(1)

    alignmentFile = gzip.open(featuresDir + chr + "_features.out.gz", 'r')

    # Read file headers to get them out of the way
    header = alignmentFile.readline()

    startTime = time.time()
    curChunk = 0

    print "Creating output file . . . ",
    ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')
    print "Done."

    lastPos = -1
    for line in alignmentFile:
        splitLine = line.strip().split(",")
        pos, features = createFeatures(splitLine)

        if pos != (lastPos + 1):
            for i in range(lastPos + 1, pos):
                if (i != 0) and (i % linesPerChunk == 0): # check to see if we should start a new chunk
                    ofile.close()
                    curChunk += 1
                    ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')

                ofile.write(str(i) + "\t")
                extraFeatures = [0, 2] * numFeats
                ofile.write("\t".join([str(k) for k in extraFeatures]))
                ofile.write("\n")

        if (pos != 0) and (pos % linesPerChunk == 0):
            ofile.close()
            curChunk += 1
            ofile = gzip.open(oDir + chr + "_" + str(curChunk) + "_features_binary.txt.gz", 'w')

        ofile.write(str(pos) + "\t")
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

            ofile.write(str(i) + "\t")
            extraFeatures = [0, 2] * numFeats
            ofile.write("\t".join([str(k) for k in extraFeatures]))
            ofile.write("\n")

    ofile.close()

    endTime = time.time()
    print "Done.  Time: ", endTime - startTime, " seconds."

main()
