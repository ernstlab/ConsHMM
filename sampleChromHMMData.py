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
            features[idx + 1] = 0
        else:
            features[idx] = 1
            if splitLine[i] == human:
                features[idx + 1] = 1
            else:
                features[idx + 1] = 0

    return int(splitLine[0]), features

def getChrSamples(chr, chrLength, numSamples, sampleLength):
    pickedStarts = defaultdict(bool)
    numPicked = 0
    while (numPicked < numSamples):
        start = random.randint(0, chrLength - sampleLength)
        end = start + sampleLength

        # check if any other sample picked already contains this sample
        alreadyPicked = False
        for otherStart in pickedStarts:
            otherEnd = otherStart + sampleLength
            if ((otherStart <= start) and (start <= otherStart + sampleLength)) or ((otherStart <= end) and (end <= otherEnd)):
                alreadyPicked = True
                break

        if alreadyPicked:
            continue

        # if it doens't overlap anything add it to the samples set
        pickedStarts[start] = True
        numPicked += 1
    return pickedStarts.keys()

def main():
    if len(sys.argv) < 6:
        print "Usage: sampleChromHMMData.py <binarized features directory> <output directory> <total Bases> <sample length> <chromosome>"
        exit(1)

    featuresDir = formatDir(sys.argv[1])
    oDir = formatDir(sys.argv[2])
    totalBases = int(sys.argv[3])
    sampleLength = int(sys.argv[4])
    chr = sys.argv[5]

    chrLength = {"chr1":249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260,
            "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747,
            "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
            "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
            "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}
    totalGenomeLength = sum(chrLength.values())

    allSpecies = sorted(['amaVit1', 'latCha1', 'criGri1', 'punNye1', 'odoRosDiv1', 'lepWed1', 'anoCar2', 'chiLan1', 'loxAfr3', 'rn5', 'tetNig2', 'oreNil2', 'astMex1', 'eptFus1', 'vicPac2', 'geoFor1', 'oryLat2', 'sorAra2', 'gadMor1', 'eriEur2', 'oryCun2', 'micOch1', 'araMac1', 'octDeg1', 'dasNov3', 'taeGut2', 'mayZeb1', 'macEug2', 'musFur1', 'ficAlb2', 'myoLuc2', 'gorGor3', 'pteVam1', 'galGal4', 'pelSin1', 'xipMac1', 'anaPla1', 'ornAna1', 'sarHar1', 'takFla1', 'equCab2', 'panHod1', 'rheMac3', 'camFer1', 'cerSim1', 'jacJac1', 'macFas5', 'calJac3', 'melUnd1', 'neoBri1', 'conCri1', 'monDom5', 'pteAle1', 'oviAri3', 'saiBol1', 'ailMel1', 'susScr3', 'canFam3', 'echTel2', 'papHam1', 'capHir1', 'ponAbe2', 'pseHum1', 'hapBur1', 'mm10', 'xenTro7', 'falPer1', 'fr3', 'colLiv1', 'gasAcu1', 'nomLeu3', 'zonAlb1', 'tupChi1', 'chrAsi1', 'cavPor3', 'panTro4', 'cheMyd1', 'triMan1', 'lepOcu1', 'apaSpi1', 'mesAur1', 'eleEdw1', 'bosTau7', 'ochPri3', 'melGal1', 'petMar2', 'oryAfe1', 'hetGla2', 'turTru2', 'speTri2', 'allMis1', 'chlSab1', 'otoGar3', 'myoDav1', 'falChe1', 'danRer7', 'orcOrc1', 'felCat5', 'chrPic1'])  # list of all the possible species in the alignment

    numFeats = 99 * 2

    createDir(oDir)

    numTotalSamples = totalBases / sampleLength # total number of samples
    numSamples = (numTotalSamples * chrLength[chr]) / totalGenomeLength # how many samples will correspond to this chromosome
    sampleStarts = sorted(getChrSamples(chr, chrLength[chr], numSamples, sampleLength))

    print "Extracting ", numSamples, " samples from ", chr, " . . ."
    for idx in range(len(sampleStarts) - 1):
        if (sampleStarts[idx] + sampleLength) < (sampleStarts[idx] + 1):
            print "Error: samples overlap."
            exit(1)
    startTime = time.time()

    featuresFile = gzip.open(featuresDir + chr + "_features.out.gz", 'r')

    # Read file headers to get them out of the way
    header = featuresFile.readline()

    featuresPos = -1

    for j in range(len(sampleStarts)):
        print "Creating file for sample ", j, " . . .",
        start = sampleStarts[j]

        ofile = gzip.open(oDir + chr + "_" + str(j) + "_marks_binary.txt.gz", 'w')
        ofile.write("cell1\t" + chr + "\n")
        for species in allSpecies:
            ofile.write(species + "\t")
        ofile.write("\n")

        for genomePos in range(start, start + sampleLength):
            while featuresPos < genomePos:
                splitLine = featuresFile.readline().strip().split(",")
                if splitLine[0] == '':
                    featuresPos = 300000000 # mark that there is nothing more in UCSC
                    break
                featuresPos, featuresFile = createFeatures(splitLine)

            if (featuresPos > genomePos):
                ConsHMMfeatures = np.array([0] * numFeats)

            for k in range(len(ConsHMMfeatures)):
                ofile.write(str(int(ConsHMMfeatures[k])) + "\t")
            ofile.write("\n")

        print "Done. Output file: ", oDir + chr + "_" + str(j) + "_marks_binary.txt.gz"
        ofile.close()

    endTime = time.time()
    print "Done.  Time: ", endTime - startTime, " seconds."

main()
