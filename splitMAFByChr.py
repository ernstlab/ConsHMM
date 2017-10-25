from __future__ import print_function
from shared import *
import sys
import gzip

def main():
    if len(sys.argv) != 5:
        print("Usage: python splitMAFByChr.py <MAF file> <output directory> <reference species> <chromosome list>")
        exit(1)
    mafFile = gzip.open(sys.argv[1], 'r')
    oDir = formatDir(sys.argv[2])
    refSpecies = sys.argv[3]
    chromosomeList = open(sys.argv[4], 'r')

    print("Creating output files (one for each chromosome) . . .",)
    outputFiles = {}
    for line in chromosomeList:
        curChr = line.strip()
        oFile = gzip.open(oDir + curChr + ".maf.gz", "w")
        outputFiles[curChr] = oFile
    print("Done.")

    print("Splitting up MAF file per chromosome . . .")
    curBlock = []
    curChr = ""
    skippedChromosomes = set()
    for line in mafFile:
        if line[0] == "#": # comment lines will get discarded
            continue

        if line == "\n":
            if curChr not in outputFiles:
                skippedChromosomes.add(curChr)
            else:
                for curLine in curBlock:
                    outputFiles[curChr].write(curLine)
                outputFiles[curChr].write("\n")
            curBlock = []
            curChr = ""
        else:
            if line[0] == "s":
                splitLine = line.strip().split()
                curSpecies = splitLine[1].split(".")[0]
                if (curSpecies == refSpecies):
                    curChr = splitLine[1].split(".")[1]
            curBlock.append(line)
    print("Done.")

    print("Closing files . . .")
    for oFile in outputFiles:
        oFile.close()
    print("Done.")

    print("Other chromosomes found in MAF file: ")
    for chr in skippedChromosomes:
        print(chr, end = "")

main()





