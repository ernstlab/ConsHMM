import gzip
import time
import argparse
from shared import *
import os
import sys
from collections import deque

def write_empty_line(output_file, num_feats, dummy = False):
    extra_features = [0, 2] * num_feats
    output_file.write("\t".join([str(k) for k in extra_features]))
    if dummy:
        output_file.write("\t0")
    output_file.write("\n")


def write_dummy_line(output_file, num_feats):
    extra_features = [0] * (2 * num_feats)
    output_file.write("\t".join([str(k) for k in extra_features]))
    output_file.write("\t1\n")


def output_current_block(output_file, output_meta_file, cur_block_lines, cur_block_feats, pos_reference, num_feats,
                         middle_pos):
    variant_line = cur_block_lines[10]
    for variant in ['A', 'C', 'T', 'G']:
        if variant != variant_line[pos_reference]:
            var_pos, var_features = create_features(variant_line, pos_reference, variant)

            for i in range(0, 10):  # write block in front of variant
                for k in range(len(cur_block_feats[i])):
                    output_file.write(str(int(cur_block_feats[i][k])) + "\t")
                output_file.write("0\n")  # add dummy observation at the end
            for k in range(len(var_features)):  # write variant line
                output_file.write(str(int(var_features[k])) + "\t")
            output_file.write("0\n") # add dumy observation at the end
            for i in range(11, 21):  # write block after variant
                for k in range(len(cur_block_feats[i])):
                    output_file.write(str(int(cur_block_feats[i][k])) + "\t")
                output_file.write("0\n")  # add dummy observation at the end
            write_dummy_line(output_file, num_feats)

            output_meta_file.write(str(middle_pos) + "\t" + variant_line[pos_reference] + "\t" + variant + "\n")


def create_file(file_name, output_file_list = None, file_type = "normal"):
    print("Creating file", file_name, ". . . ",)
    if file_type == "normal":
        output_file = open(file_name, 'w')
    elif file_type == "gzip":
        output_file = gzip.open(file_name, 'wt')
    else:
        sys.exit("Output file type must be either normal or gzip. You provided " + file_type + ".")

    if output_file_list is not None:
        slash_pos = file_name.rfind("/")
        output_file_list.write(file_name[(slash_pos + 1):] + "\n")
    print("Done.")

    return output_file


def binarize_reference(main_args):
    alignment_file = gzip.open(main_args.seq_file, 'rt')
    output_directory = format_dir(main_args.output_directory) + "Binarized_Input/"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    chromosome = main_args.chromosome
    num_bases = main_args.num_bases
    chromosome_sizes_file = open(main_args.chromosome_sizes_file_name, 'r')
    reference_species = main_args.reference_species
    binary_list_directory = format_dir(main_args.output_directory) + "Binarized_Input_File_Lists/"
    if not os.path.exists(binary_list_directory):
        os.makedirs(binary_list_directory)
    output_file_list = open(binary_list_directory + chromosome + "_binarized_files_list.txt", 'w')

    chromosome_sizes_dict = {}
    for line in chromosome_sizes_file:
        split_line = line.strip().split()
        chromosome_sizes_dict[split_line[0]] = int(split_line[1])

    # Read file header
    header = alignment_file.readline()
    header_split = header.strip().split(",")
    new_header = ""
    for species in header_split[2:]:
        if species != reference_species:
            if new_header != "":
                new_header += "\t"
            new_header = new_header + species + "_aligned\t" + species + "_matched"

    if reference_species not in header_split:
        sys.exit("Reference species provided " + reference_species + " is not found in header of file provided to -s "
                                                                    "flag.")

    pos_reference = header_split.index(reference_species)
    num_feats = len(header.split(",")) - 3

    start_time = time.time()
    cur_chunk = 0

    print("Creating binarized files of length ", num_bases, " bases.")
    output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz",
                              output_file_list, file_type="gzip")
    output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
    output_file.write(new_header + "\n")

    last_pos = -1
    for line in alignment_file:
        split_line = line.strip().split(",")
        pos = int(split_line[1])

        pos, features = create_features(split_line, pos_reference)

        if pos != (last_pos + 1):
            for i in range(last_pos + 1, pos):
                if (i != 0) and (i % num_bases == 0):  # check to see if we should start a new chunk
                    output_file.close()
                    cur_chunk += 1

                    output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                              "_features_binary.txt.gz", output_file_list, file_type = "gzip")
                    output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                    output_file.write(new_header + "\n")

                extra_features = [0, 2] * num_feats
                output_file.write("\t".join([str(k) for k in extra_features]))
                output_file.write("\n")

        if (pos != 0) and (pos % num_bases == 0):
            output_file.close()
            cur_chunk += 1

            output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                      "_features_binary.txt.gz", output_file_list, file_type="gzip")
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

                output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                          "_features_binary.txt.gz", output_file_list, file_type="gzip")
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
    parent_parser = argparse.ArgumentParser(description='binarizeAlignment.py takes in a file generated by parseMAF.py'
        'and splits it up into segments of a certain length, then binarizes said segments into the format needed by'
        'ChromHMM.', add_help=False)

    parent_parser.add_argument('-s', '--seqFile', required=True, dest='seq_file',
                        help='Sequence file generated by parseMAF.py.')
    parent_parser.add_argument('-o', '--outputDirectory', required=True, dest='output_directory',
                        help='Output directory in which to put the binarized segments.')
    parent_parser.add_argument('-c', '--chr', required=True, dest='chromosome',
                        help='Chromosome of the reference species contained in the input file.')
    parent_parser.add_argument('-n', '--numBases', required=True, dest='num_bases', type=int,
                        help='Number of bases in each segment if doing no or full variation or number of variants in '
                             'each file if doing local variation.')
    parent_parser.add_argument('-cl', '--chromLengths', required=True, dest='chromosome_sizes_file_name',
                        help='File containing the chrosomome lengths for the reference species in the alignment.')
    parent_parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')

    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help="Available modes", dest='command', required=True)

    parser_reference = subparsers.add_parser('reference', parents=[parent_parser], help='Default version binarizing '
                                             'an alignment with respect to the reference genome of the reference '
                                             'species in it, ignoring any potential variation.')

    main_args = main_parser.parse_args()
    binarize_reference(main_args)
    # parser_variation = subparsers.add_parser('variation', parents=[parent_parser], help='BETA VERSION; NOT STABLE; '
    #                                          'Mode for binarizing alignments including variants either for certain '
    #                                          'positions or for generating all possible variants in the genome.')
    #
    # parser_variation.add_argument('-v', '--variationFile', dest='variation_file', help='File containing variants in the'
    #                     ' reference genome (Optional). Format required is BED: chr, start, end (start + 1), rsID, '
    #                     'reference allele, alternate allele (if multiple separate by comma, but only first will be '
    #                     'kept')
    # parser_variation.add_argument('-lv', '--localVariationBinarize', type=int, dest='lv_size', help='If included only '
    #                     'locations within this many bases around the variants will be binarized')
    # parser_variation.add_argument('-d', '--dummyState', action='store_true', dest='dummy_state', help='Include if dummy'
    #                     'state should be used to separate variants.')
    # parser_variation.add_argument('-fv', '--fullVariation', action='store_true', dest='full_variation', help='Include '
    #                     'togenerate binarized file for every single possible variant on this chromosome. Do not use '
    #                     'this flag lightly -- it will generate a lot of files. By default this uses a 10bp window '
    #                     'around each variant and the dummy state to separate variants. Will cause -lv and and -v '
    #                     'flags to be ignored')
    # parser_variation.add_argument('-fi', '--fullIndex', type=int, dest='full_index', help='The index of the current '
    #                     'chunk of maf sequence file for which to create all possible variants. Since the genome is so '
    #                     'big the output of parseMAF.py was split into manageable chunks, and this script requires the '
    #                     'index of the current file for some housekeeping.')
    #if (args.lv_size) and (not args.variation_file):
    #    parser_main.error('A variation file is required if the -lv flag is specified.')
    #if (args.full_variation) and (args.full_index is None):
    #    parser_main.error('The -fi flag is required when running full variation (-fv flag)')
