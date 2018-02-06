import gzip
import time
import argparse
from shared import *
import numpy as np


def create_features(split_line, pos_reference):
    line_features = np.empty((len(split_line) - 2) * 2)  # two features for each species

    human = split_line[pos_reference]  # variable for keeping the reference nucleotide; does not have to be human
    idx = 0
    for i in range(1, len(split_line)):
        if i != pos_reference:
            if split_line[i] == 'X':
                line_features[idx] = 0
                line_features[idx + 1] = 2
            else:
                line_features[idx] = 1
                if split_line[i] == human:
                    line_features[idx + 1] = 1
                else:
                    line_features[idx + 1] = 0
            idx += 2
    
    return int(split_line[0]), line_features


def binarize_alignment(main_args):
    alignment_file = gzip.open(main_args.seq_file, 'rt')
    output_directory = format_dir(main_args.output_directory)
    chromosome = main_args.chromosome
    num_bases = main_args.num_bases
    chromosome_sizes_file = open(main_args.chromosome_sizes_file_name, 'r')
    reference_species = main_args.reference_species
    output_file_list = open(output_directory + chromosome + "_binary_list.txt", 'w')

    chromosome_sizes_dict = {}
    for line in chromosome_sizes_file:
        split_line = line.strip().split()
        chromosome_sizes_dict[split_line[0]] = int(split_line[1])

    # Read file headers
    header = alignment_file.readline()
    header_split = header.strip().split(",")
    new_header = ""
    for species in header_split[1:]:
        if species != reference_species:
            if new_header != "":
                new_header += "\t"
            new_header = new_header + species + "_aligned\t" + species + "_matched"

    pos_reference = header_split.index(reference_species)
    num_feats = len(header.split(",")) - 2

    start_time = time.time()
    cur_chunk = 0

    print("Creating file", output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz", ". . . ",)
    output_file = gzip.open(output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz",
                            'wt')
    output_file_list.write(chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz\n")
    output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
    output_file.write(new_header + "\n")
    print("Done.")

    last_pos = -1
    for line in alignment_file:
        split_line = line.strip().split(",")
        pos, features = create_features(split_line, pos_reference)

        if pos != (last_pos + 1):
            for i in range(last_pos + 1, pos):
                if (i != 0) and (i % num_bases == 0):  # check to see if we should start a new chunk
                    output_file.close()
                    cur_chunk += 1
                    print("Creating file", output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz", ". . . ",)
                    output_file = gzip.open(output_directory + chromosome + "_" + str(cur_chunk) +
                                            "_features_binary.txt.gz", 'wt')
                    print("Done.")
                    output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                    output_file.write(new_header + "\n")

                extra_features = [0, 2] * num_feats
                output_file.write("\t".join([str(k) for k in extra_features]))
                output_file.write("\n")

        if (pos != 0) and (pos % num_bases == 0):
            output_file.close()
            cur_chunk += 1
            print("Creating file", output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz", ". . . ",)
            output_file = gzip.open(output_directory + chromosome + "_" + str(cur_chunk) +
                                    "_features_binary.txt.gz", 'wt')
            print("Done.")
            output_file_list.write(chromosome + "_" + str(cur_chunk) + "_binary.txt.gz\n")
            output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
            output_file.write(new_header + "\n")

        for k in range(len(features)):
            output_file.write(str(int(features[k])) + "\t")
        output_file.write("\n")
        last_pos = pos

    if last_pos != chromosome_sizes_dict[chromosome]:
        for i in range(last_pos + 1, chromosome_sizes_dict[chromosome]):
            if i % num_bases == 0:
                output_file.close()
                cur_chunk += 1
                print("Creating file", output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz", ". . . ",)
                output_file = gzip.open(output_directory + chromosome + "_" + str(cur_chunk) +
                                        "_features_binary.txt.gz", 'wt')
                print("Done.")
                output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                output_file.write(new_header + "\n")

            extra_features = [0, 2] * num_feats
            output_file.write("\t".join([str(k) for k in extra_features]))
            output_file.write("\n")

    output_file_list.close()
    output_file.close()

    end_time = time.time()
    print("Done.  Time: ", end_time - start_time, " seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""binarizeAlignment.py takes in a file generated by parseMAF.py
        and splits it up into segments of a certain length, then binarizes said segments into the format needed by
        ChromHMM.""")

    parser.add_argument('-s', '--seqFile', required=True, dest='seq_file',
                        help='Sequence file generated by parseMAF.py.')
    parser.add_argument('-o', '--outputDirectory', required=True, dest='output_directory',
                        help='Output directory in which to put the binarized segments.')
    parser.add_argument('-c', '--chr', required=True, dest='chromosome',
                        help='Chromosome of the reference species contained in the input file.')
    parser.add_argument('-n', '--numBases', required=True, dest='num_bases', type=int,
                        help='Number of bases in each segment.')
    parser.add_argument('-cl', '--chromLengths', required=True, dest='chromosome_sizes_file_name',
                        help='File containing the chrosomome lengths for the reference species in the alignment.')
    parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')

    args = parser.parse_args()

    binarize_alignment(args)
