'''
createMAFFeatures.py takes in an alignment file in MAF format and extracts just the sequence information
creates a feature file in which
each line corresponds to a nucleotide in the human genome, and contains columns explaining whether
the other species have the same nucleotide at that position or not.
----The stuff below is not accurate----
The output format is:
pos,species1,species2,...,species99
pos = position in the human genome
speciesX = 0 if there is an aligned nucleotide in that particular species
         = 1 if there is no aligned nucleotide in that particular species
The order of the species is the one specified in the first line of the file.
'''

__author__ = 'adriana'
'''
createFeatures.py takes in an alignment file in MAF format and creates a feature file in which
each line corresponds to a nucleotide in the human genome, and contains columns explaining whether
the other species have the same nucleotide at that position or not.
The output format is:
pos,species1,species2,...,species99
pos = position in the human genome
speciesX = 0 if there is an aligned nucleotide in that particular species
         = 1 if there is no aligned nucleotide in that particular species
The order of the species is the one specified in the first line of the file.
'''

import gzip
import sys
from collections import defaultdict
from Bio.Seq import Seq

def main():
    if len(sys.argv) < 6:
        print "Usage: python createMAFFeatures.py <alignment file> <alignment source> <output file> <chromosome> <reference species> <start position (optional)"
        exit(1)
    mafFile = sys.argv[1]
    source = sys.argv[2]
    featFile = sys.argv[3]
    chr = sys.argv[4]
    refSpecies = sys.argv[5]
    if len(sys.argv) == 7:
        startPos = int(sys.argv[7])
    f = gzip.open(mafFile, 'rb')
    o = gzip.open(featFile, 'w')

    if source == "UCSC_hg19":
        species = sorted(['amaVit1', 'latCha1', 'criGri1', 'punNye1', 'odoRosDiv1', 'lepWed1', 'anoCar2', 'chiLan1', 'loxAfr3', 'rn5', 'tetNig2', 'oreNil2', 'astMex1', 'eptFus1', 'vicPac2', 'geoFor1', 'oryLat2', 'sorAra2', 'gadMor1', 'eriEur2', 'oryCun2', 'micOch1', 'araMac1', 'octDeg1', 'dasNov3', 'taeGut2', 'mayZeb1', 'macEug2', 'musFur1', 'ficAlb2', 'myoLuc2', 'gorGor3', 'pteVam1', 'galGal4', 'pelSin1', 'xipMac1', 'anaPla1', 'ornAna1', 'sarHar1', 'takFla1', 'equCab2', 'panHod1', 'rheMac3', 'camFer1', 'cerSim1', 'jacJac1', 'macFas5', 'calJac3', 'melUnd1', 'neoBri1', 'conCri1', 'monDom5', 'pteAle1', 'oviAri3', 'saiBol1', 'ailMel1', 'susScr3', 'canFam3', 'echTel2', 'papHam1', 'capHir1', 'ponAbe2', 'pseHum1', 'hapBur1', 'mm10', 'xenTro7', 'falPer1', 'fr3', 'colLiv1', 'gasAcu1', 'nomLeu3', 'zonAlb1', 'tupChi1', 'chrAsi1', 'cavPor3', 'panTro4', 'cheMyd1', 'triMan1', 'lepOcu1', 'apaSpi1', 'mesAur1', 'eleEdw1', 'bosTau7', 'ochPri3', 'melGal1', 'petMar2', 'oryAfe1', 'hetGla2', 'turTru2', 'speTri2', 'allMis1', 'chlSab1', 'otoGar3', 'myoDav1', 'falChe1', 'danRer7', 'orcOrc1', 'felCat5', 'chrPic1'])  # list of all the possible species in the alignment
        humanSpecies = "hg19." + chr
    elif source == "UCSC_hg38":
        species = sorted(["panTro4","gorGor3","ponAbe2","nomLeu3","rheMac3","macFas5","papAnu2","chlSab2","calJac3","saiBol1","otoGar3","tupChi1","speTri2","jacJac1","micOch1","criGri1","mesAur1","mm10","rn6","hetGla2","cavPor3","chiLan1","octDeg1","oryCun2","ochPri3","susScr3","vicPac2","camFer1","turTru2","orcOrc1","panHod1","bosTau8","oviAri3","capHir1","equCab2","cerSim1","felCat8","canFam3","musFur1","ailMel1","odoRosDiv1","lepWed1","pteAle1","pteVam1","myoDav1","myoLuc2","eptFus1","eriEur2","sorAra2","conCri1","loxAfr3","eleEdw1","triMan1","chrAsi1","echTel2","oryAfe1","dasNov3","monDom5","sarHar1","macEug2","ornAna1","falChe1","falPer1","ficAlb2","zonAlb1","geoFor1","taeGut2","pseHum1","melUnd1","amaVit1","araMac1","colLiv1","anaPla1","galGal4","melGal1","allMis1","cheMyd1","chrPic2","pelSin1","apaSpi1","anoCar2","xenTro7","latCha1","tetNig2","fr3","takFla1","oreNil2","neoBri1","hapBur1","mayZeb1","punNye1","oryLat2","xipMac1","gasAcu1","gadMor1","danRer10","astMex1","lepOcu1","petMar2"])
        humanSpecies = "hg38." + chr
        chrLength = {"chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,"chr6":170805979,"chr7":159345973,"chrX":156040895,"chr8":145138636,"chr9":138394717,"chr11":135086622,"chr10":133797422,"chr12":133275309,"chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285,"chr20":64444167,"chr19":58617616,"chrY":57227415,"chr22":50818468,"chr21":46709983}
    elif source == "Ensembl_EPO":
        species = sorted(['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'macaca_mulatta', 'callithrix_jacchus', 'mus_musculus', 'rattus_norvegicus', 'oryctolagus_cuniculus', 'canis_familiaris', 'felis_catus', 'equus_caballus', 'bos_taurus', 'ovis_aries', 'sus_scrofa', 'nomascus_leucogenys', 'tarsius_syrichta', 'microcebus_murinus', 'otolemur_garnettii', 'tupaia_belangeri', 'dipodomys_ordii', 'ictidomys_tridecemlineatus', 'cavia_porcellus', 'ochotona_princeps', 'ailuropoda_melanoleuca', 'mustela_putorius_furo', 'myotis_lucifugus', 'pteropus_vampyrus', 'tursiops_truncatus', 'vicugna_pacos', 'erinaceus_europaeus', 'sorex_araneus', 'loxodonta_africana', 'procavia_capensis', 'echinops_telfairi', 'dasypus_novemcinctus', 'choloepus_hoffmanni'])  # list of all the possible species in the alignment
        humanSpecies = "homo_sapiens." + chr[3:]
    elif source == "Ensembl_Pecan":
        species = sorted(['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'macaca_mulatta', 'callithrix_jacchus', 'mus_musculus', 'rattus_norvegicus', 'oryctolagus_cuniculus', 'canis_familiaris', 'felis_catus', 'equus_caballus', 'bos_taurus', 'ovis_aries', 'sus_scrofa', 'monodelphis_domestica', 'ornithorhynchus_anatinus', 'gallus_gallus', 'meleagris_gallopavo', 'taeniopygia_guttata', 'anolis_carolinensis'])
        humanSpecies = "homo_sapiens." + chr[3:]
    elif source == "UCSC_mm10":
        species = sorted(['mm10', 'cavPor3', 'dipOrd1', 'hetGla2', 'ochPri2', 'oryCun2', 'rn5', 'speTri2', 'tupBel1', 'calJac3', 'gorGor3', 'hg19', 'micMur1', 'nomLeu2', 'otoGar3', 'panTro4', 'papHam1', 'ponAbe2', 'rheMac3', 'saiBol1', 'tarSyr1', 'ailMel1', 'bosTau7', 'canFam3', 'choHof1', 'dasNov3', 'echTel1', 'equCab2', 'eriEur1', 'felCat5', 'loxAfr3', 'myoLuc2', 'oviAri1', 'proCap1', 'pteVam1', 'sorAra1', 'susScr3', 'triMan1', 'turTru2', 'vicPac1', 'anoCar2', 'chrPic1', 'danRer7', 'fr3', 'gadMor1', 'galGal4', 'gasAcu1', 'latCha1', 'macEug2', 'melGal1', 'melUnd1', 'monDom5', 'oreNil2', 'ornAna1', 'oryLat2', 'petMar1', 'sarHar1', 'taeGut1', 'tetNig2', 'xenTro3'])
        humanSpecies = refSpecies + "." + chr
    else:
        print source
        print "Alignment source not recognized. Please input UCSC or Ensembl."
        exit(1)

    chrLength = {"chr1":249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520, "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}
    if refSpecies == "mm10":
        chrLength = {"chr1":195471971,"chr2":182113224,"chrX":171031299,"chr3":160039680,"chr4":156508116,"chr5":151834684,"chr6":149736546,"chr7":145441459,"chr10":130694993,"chr8":129401213,"chr14":124902244,"chr9":124595110,"chr11":122082543,"chr13":120421639,"chr12":120129022,"chr15":104043685,"chr16":98207768,"chr17":94987271,"chrY":91744698,"chr18":90702639,"chr19":61431566}
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

                if curSpecies == humanSpecies: # store info for human
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
                    o.write(str(pos) + "," + seqHuman[i])
                    for sp in species:
                        if (sp not in alignedToHuman) or (alignedToHuman[sp][i] == '-'):
                            o.write(",X")
                        else:
                            o.write("," + str(alignedToHuman[sp][i]))
                    o.write("\n")
                    covered[pos] = 1
                pos += 1

        scoreLine = f.readline()
    print numBases, " of the genome covered."
    o.close()
main()

panTro4	gorGor3	ponAbe2	nomLeu3	rheMac3	macFas5	papHam1	chlSab1	calJac3	saiBol1	otoGar3	tupChi1	speTri2	jacJac1	micOch1	criGri1	mesAur1	mm10	rn5	hetGla2	cavPor3	chiLan1	octDeg1	oryCun2	ochPri3	susScr3	vicPac2	camFer1	turTru2	orcOrc1	panHod1	bosTau7	oviAri3	capHir1	equCab2	cerSim1	felCat5	canFam3	musFur1	ailMel1	odoRosDiv1	lepWed1	pteAle1	pteVam1	myoDav1	myoLuc2	eptFus1	eriEur2	sorAra2	conCri1	loxAfr3	eleEdw1	triMan1	chrAsi1	echTel2	oryAfe1	dasNov3	monDom5	sarHar1	macEug2	ornAna1	falChe1	falPer1	ficAlb2	zonAlb1	geoFor1	taeGut2	pseHum1	melUnd1	amaVit1	araMac1	colLiv1	anaPla1	galGal4	melGal1	allMis1	cheMyd1	chrPic1	pelSin1	apaSpi1	anoCar2	xenTro7	latCha1	tetNig2	fr3	takFla1	oreNil2	neoBri1	hapBur1	mayZeb1	punNye1	oryLat2	xipMac1	gasAcu1	gadMor1	danRer7	astMex1	lepOcu1	petMar2

