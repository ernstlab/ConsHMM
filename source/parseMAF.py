import gzip
import sys
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from shared import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""parse MAF.py takes in a MAF file that contains a multiple sequence alignment for one of the chromosomes
    of the reference species in the alignment, and outputs a file containing just the sequence information. The output file will have as many rows as the
    number of bases in the input chromosome, and N + 1 columns, where N is the number of species in the multiple sequence alignment. The first column contains the
    coordinate of a base in the reference species. The following columns contain the nucleotide at that position in the alignment for each species including the
    reference. Indels and 'N' nucleotides are both encoded with 'X'.""", allow_abbrev=False)
    parser.add_argument('-m', '--MAFFile', required = True, dest='mafFile', help='.maf file containing a multi-sequence alignment.')
    parser.add_argument('-s', '--speciesFile', required=True, dest='speciesFile', help='File containing a list of the names of all the species in the alignment.')
    parser.add_argument('-o', '--outputFile', required=True, dest='outputFile', help='Output file to write the sequence information to.')
    parser.add_argument('-c', '--chromosome', required=True, dest='chr', help='Chromosome contained in the MAF file. Specify as chr*')
    parser.add_argument('-r', '--refSpecies', required=True, dest='refSpecies', help='Genome assembly name of the reference species in the alignment.')
    parser.add_argument('-cl', '--chromLengths', required=True, dest='chromLengths', help='File containing the chrosomome lengths for the reference species in the alignment.')

    args = parser.parse_args()

    mafFile = args.mafFile
    speciesFile = open(args.speciesFile, 'r')
    featFile = args.outputFile
    chr = args.chr
    refSpecies = args.refSpecies + "." + chr
    chromLengths = args.chromLengths

    f = gzip.open(mafFile, 'rt')
    o = gzip.open(featFile, 'wt')
    chromSizes = open(chromLengths, 'r')

    species = []
    for line in speciesFile:
        species.append(line.strip())
    species = sorted(species)
    
    chrLength = {}
    for line in chromSizes:
        splitLine = line.strip().split()
        chrLength[splitLine[0]] = int(splitLine[1])

    if chr not in chrLength:
        printDictKeys(chrLength)
        exit(1)

    idx = 0
    # write header to output file
    o.write("pos")
    for sp in species:
        o.write("," + sp)
    o.write("\n")

    f.readline() # read header line
    scoreLine = f.readline()
    while scoreLine[0] == '#':
        scoreLine = f.readline()

    numBases = 0
    covered = {} # dictionary of bases covered to account for blocks convering duplicate pieces
    while scoreLine != "":
        sequenceLine = f.readline()
        startHuman = -1
        seqHuman = ""
        alignedToHuman = {}

        lastAlignment = []
        reverseComplement = False
        # parse an alignment block
        while sequenceLine != "\n":
            lastAlignment.append(sequenceLine)
            splitSeq = sequenceLine.split()

            type = splitSeq[0]
            if type == 's':
                curSpecies = splitSeq[1]
                
                if curSpecies == refSpecies: # store info for human
                    trueLenHuman = int(splitSeq[3])
                    if splitSeq[4] == '-': # reverse complement
                        startHuman = chrLength[chr] - int(splitSeq[2]) - trueLenHuman
                        seqHuman = Seq(splitSeq[6].upper()).reverse_complement()
                        reverseComplement=True
                    else:
                        startHuman = int(splitSeq[2])
                        seqHuman = splitSeq[6].upper()
                    lenHumanDash = len(seqHuman)
                else:
                    curSeq = ""
                    if reverseComplement:
                        curSeq = Seq(splitSeq[6].upper()).reverse_complement()
                    else:
                        curSeq = splitSeq[6].upper()
                    alignedToHuman[curSpecies[:curSpecies.find(".")]] = curSeq # store the aligned sequence in other species
            sequenceLine = f.readline()

        # when done with an alignment block go through and output whether each species is the same as human
        pos = startHuman
        for i in range(lenHumanDash):
            if seqHuman[i] != '-':
                numBases += 1
                if not (pos in covered):
                    o.write(str(pos))
                    for sp in species:
                        if (sp == sys.argv[5]):
                            o.write("," + seqHuman[i])
                        elif (sp not in alignedToHuman) or (alignedToHuman[sp][i] == '-'):
                            o.write(",X")
                        else:
                            o.write("," + str(alignedToHuman[sp][i]))
                    o.write("\n")
                    covered[pos] = 1
                pos += 1
        scoreLine = f.readline()
    o.close()


