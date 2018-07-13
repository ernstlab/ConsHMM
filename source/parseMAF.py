import gzip
import argparse
from Bio.Seq import Seq
import importlib
from shared import *
import sys


def parse_maf(main_args):
    maf_file = main_args.maf_file
    species_file = open(main_args.species_file, 'r')
    output_file = main_args.output_file
    chromosome = main_args.chromosome
    reference_species = main_args.reference_species + "." + chromosome
    chromosome_sizes_file_name = main_args.chromosome_sizes_file_name

    f = gzip.open(maf_file, 'rt')
    o = gzip.open(output_file, 'wt')
    chromosome_sizes_file = open(chromosome_sizes_file_name, 'r')

    species = []
    for line in species_file:
        species.append(line.strip())
    species = sorted(species)

    chromosome_sizes_dict = {}
    for line in chromosome_sizes_file:
        split_line = line.strip().split()
        chromosome_sizes_dict[split_line[0]] = int(split_line[1])

    if chromosome not in chromosome_sizes_dict:
        print_dict_keys(chromosome_sizes_dict)
        exit(1)

    # write header to output file
    o.write("pos")
    for sp in species:
        o.write("," + sp)
    o.write("\n")

    score_line = f.readline()  # read first line

    num_bases = 0
    covered = {}  # dictionary of bases covered to account for blocks convering duplicate pieces
    while score_line != "":
        if score_line[0] == "#":  # skip comment lines
            score_line = f.readline()
            continue

        sequence_line = f.readline()
        start_human = -1
        seq_human = ""
        aligned_to_human = {}

        last_alignment = []
        reverse_complement = False

        # parse an alignment block
        len_human_dash = -1
        while sequence_line != "\n":
            last_alignment.append(sequence_line)
            split_seq = sequence_line.split()

            line_type = split_seq[0]
            if line_type == 's':
                cur_species = split_seq[1]

                if cur_species == reference_species:  # store info for human
                    true_len_human = int(split_seq[3])
                    if split_seq[4] == '-':  # reverse complement
                        start_human = chromosome_sizes_dict[chromosome] - int(split_seq[2]) - true_len_human
                        seq_human = Seq(split_seq[6].upper()).reverse_complement()
                        reverse_complement = True
                    else:
                        start_human = int(split_seq[2])
                        seq_human = split_seq[6].upper()
                    len_human_dash = len(seq_human)
                else:
                    if reverse_complement:
                        cur_seq = Seq(split_seq[6].upper()).reverse_complement()
                    else:
                        cur_seq = split_seq[6].upper()
                    # store the aligned sequence in other species
                    aligned_to_human[cur_species[:cur_species.find(".")]] = cur_seq
            sequence_line = f.readline()

        # when done with an alignment block go through and output whether each species is the same as human
        if len_human_dash == -1:
            print("Reference species not found in alignment block:")
            print(last_alignment)
            exit(1)

        pos = start_human
        for i in range(len_human_dash):
            if seq_human[i] != '-':
                num_bases += 1
                if not (pos in covered):
                    o.write(str(pos))
                    for sp in species:
                        if sp == reference_species[:reference_species.find(".")]:
                            o.write("," + seq_human[i])
                        elif (sp not in aligned_to_human) or (aligned_to_human[sp][i] == '-'):
                            o.write(",X")
                        else:
                            o.write("," + str(aligned_to_human[sp][i]))
                    o.write("\n")
                    covered[pos] = 1
                pos += 1
        score_line = f.readline()
    o.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""parse MAF.py takes in a MAF file that contains a multiple sequence
        alignment for one of the chromosomes of the reference species in the alignment, and outputs a file containing
        just the sequence information. The output file will have as many rows as the number of bases in the input
        chromosome, and N + 1 columns, where N is the number of species in the multiple sequence alignment. The first
        column contains the coordinate of a base in the reference species. The following columns contain the nucleotide
        at that position in the alignment for each species including the reference. Indels and 'N' nucleotides are
        both encoded with 'X'.""")

    parser.add_argument('-m', '--MAFFile', required=True, dest='maf_file',
                        help='.maf file containing a multi-sequence alignment.')
    parser.add_argument('-s', '--speciesFile', required=True, dest='species_file',
                        help='File containing a list of the names of all the species in the alignment.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file to write the sequence information to.')
    parser.add_argument('-c', '--chromosome', required=True, dest='chromosome',
                        help='Chromosome contained in the MAF file. Specify as chr*')
    parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')
    parser.add_argument('-cl', '--chromLengths', required=True, dest='chromosome_sizes_file_name',
                        help='File containing the chrosomome lengths for the reference species in the alignment.')

    args = parser.parse_args()

    parse_maf(args)
