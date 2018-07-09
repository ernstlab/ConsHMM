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


def create_file(file_name, output_file_list = None, type = "normal"):
    print("Creating file", file_name, ". . . ",)
    if type == "normal":
        output_file = open(file_name, 'w')
    elif type == "gzip":
        output_file = gzip.open(file_name, 'wt')
    else:
        sys.exit("Output file type must be either normal or gzip. You provided ", type, ".")

    if output_file_list is not None:
        slash_pos = file_name.rfind("/")
        output_file_list.write(file_name[(slash_pos + 1):] + "\n")
    print("Done.")

    return output_file


def binarize_alignment(main_args):
    # need to use chunk index to update proper behavior of this script. Should skip first 10 lines if not first chunk
    #  etc. last chunk should be different etc.
    alignment_file = gzip.open(main_args.seq_file, 'rt')
    output_directory = format_dir(main_args.output_directory)
    chromosome = main_args.chromosome
    num_bases = main_args.num_bases
    chromosome_sizes_file = open(main_args.chromosome_sizes_file_name, 'r')
    reference_species = main_args.reference_species
    output_file_list = open(output_directory + chromosome + "_binary_list.txt", 'w')
    dummy_state = main_args.dummy_state
    full_variation = main_args.full_variation
    full_index = main_args.full_index

    chromosome_sizes_dict = {}
    for line in chromosome_sizes_file:
        split_line = line.strip().split()
        chromosome_sizes_dict[split_line[0]] = int(split_line[1])

    # If provided read in SNPs
    SNP_positions = {}
    if main_args.variation_file is not None:
        print("Reading in variation file ", main_args.variation_file, " . . .", end="", flush=True)
        lv_size = main_args.lv_size
        SNP_positions, SNP_lines = parse_variation_file(main_args.variation_file, chromosome)
        print("Done.")

    # Read file header
    header = alignment_file.readline()
    header_split = header.strip().split(",")
    new_header = ""
    for species in header_split[1:]:
        if species != reference_species:
            if new_header != "":
                new_header += "\t"
            new_header = new_header + species + "_aligned\t" + species + "_matched"
    if dummy_state or full_variation:
        new_header += "\tdummy"

    pos_reference = header_split.index(reference_species)
    num_feats = len(header.split(",")) - 2

    start_time = time.time()
    cur_chunk = 0
    num_SNPs_done = 0

    if (main_args.lv_size is not None) and (main_args.full_variation is None):
        lv_size = main_args.lv_size
        print("Creating binarized files with a window of ", lv_size, " around the SNPs provided in the variation "
                                                                     "file . . .", flush=True)

        sorted_SNPs = sorted(SNP_positions.keys())
        pos = -1 

        for SNP_pos in sorted_SNPs:
            # open a file for current batch of SNPs if needed
            if (num_SNPs_done % num_bases) == 0:
                output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                          "_features_binary.txt.gz", output_file_list, type="gzip")
                output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                output_file.write(new_header + "\n")

            binary_features_dict = {}
            while (pos <= (SNP_pos[1] + lv_size)) and (pos < chromosome_sizes_dict[chromosome]):
                # keep reading from alignment file until get to window of SNP
                split_line = alignment_file.readline().strip().split(",")
                if not split_line[0]:
                    break
                pos = int(split_line[0])
                if (pos >= SNP_pos[1] - lv_size) and (pos <= SNP_pos[1] + lv_size): # if we are in the window
                    if (chromosome, pos) in SNP_positions:
                        pos, features = create_features(split_line, pos_reference, SNP_positions[(chromosome, pos)])
                    else:
                        pos, features = create_features(split_line, pos_reference)
                    binary_features_dict[pos] = features

            for i in range(SNP_pos[1] - lv_size, SNP_pos[1] + lv_size + 1):
                if i in binary_features_dict: # if we have alignment at this position
                    for k in range(len(binary_features_dict[i])):
                        output_file.write(str(int(binary_features_dict[i][k])) + "\t")
                    if dummy_state:
                        output_file.write("0\n") # add dummy observation at the end
                else:
                    write_empty_line(output_file, num_feats, dummy=dummy_state)

            if dummy_state:
                write_dummy_line(output_file, num_feats)
            num_SNPs_done += 1

            if (num_SNPs_done % num_bases) == 0:
                output_file.close()
                cur_chunk += 1

        if (num_SNPs_done % num_bases) != 0:
            output_file.close()
    elif main_args.full_variation:
        print('Creating binarized files for every possible variant on this chromosome, with a window of 10bp on both '
              'sides of each variant, 3 possible variants at each position, 3000 bases per file')
        #dir_index = 0
        cur_chunk = 0
        #cur_output_dir = output_directory + chromosome + "_" + str(full_index) + "_binarized_files/"
        # try:
        #     os.stat(cur_output_dir)
        # except:
        #     print('Creating', cur_output_dir, '. . .')
        #     os.mkdir(cur_output_dir)
        #     print('Done')

        output_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" + str(cur_chunk) +
                                  "_features_binary.txt.gz", output_file_list, type = "gzip")
        output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
        output_file.write(new_header + "\n")

        output_meta_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" + str(cur_chunk) +
                                       "_meta_binary.txt.gz", type = "gzip")

        last_pos = -1
        num_bases_done = 0
        line_index = 0
        cur_block_lines = deque()
        cur_block_feats = deque()
        for line in alignment_file:
            split_line = line.strip().split(",")
            pos = int(split_line[0])
            line_index += 1

            # fill in the gap if there is one unless we are in the first 10 lines on a non 0 chunk
            if (pos != last_pos + 1) and (line_index > 10 or full_index == 0):
                for gap_pos in range(last_pos + 1, pos):
                    extra_features = [0, 2] * num_feats
                    cur_block_lines.append([str(gap_pos)] + ['X'] * (num_feats + 1))
                    cur_block_feats.append(extra_features)

                    if len(cur_block_feats) > 21:
                        cur_block_feats.popleft()
                        cur_block_lines.popleft()

                    if len(cur_block_feats) == 21:
                        if num_bases_done == num_bases:
                            output_file.close()
                            output_meta_file.close()
                            cur_chunk += 1

                            output_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" +
                                                      str(cur_chunk) + "_features_binary.txt.gz", output_file_list,
                                                      type = "gzip")
                            output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                            output_file.write(new_header + "\n")

                            output_meta_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" +
                                                           str(cur_chunk) + "_meta_binary.txt.gz", type = "gzip")

                            num_bases_done = 0
                        output_current_block(output_file, output_meta_file, cur_block_lines, cur_block_feats,
                                             pos_reference, num_feats, gap_pos - 10)
                        num_bases_done += 1

            pos, features = create_features(split_line, pos_reference)
            last_pos = pos
            cur_block_feats.append(features)
            cur_block_lines.append(split_line)

            if len(cur_block_feats) > 21:
                cur_block_feats.popleft()
                cur_block_lines.popleft()

            if len(cur_block_feats) == 21:
                if num_bases_done == num_bases:
                    output_file.close()
                    output_meta_file.close()
                    cur_chunk += 1

                    output_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" +
                                              str(cur_chunk)+ "_features_binary.txt.gz", output_file_list, type="gzip")
                    output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                    output_file.write(new_header + "\n")

                    output_meta_file = create_file(output_directory + chromosome + "_" + str(full_index) + "_" +
                                                   str(cur_chunk) + "_meta_binary.txt.gz",
                                                   type = "gzip")

                    num_bases_done = 0
                output_current_block(output_file, output_meta_file, cur_block_lines, cur_block_feats, pos_reference,
                                     num_feats, pos - 10)
                num_bases_done += 1
        output_file.close()
        output_meta_file.close()
    else:
        print("Creating binarized files of length ", num_bases, " bases.")
        output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) + "_features_binary.txt.gz",
                                  output_file_list, type="gzip")
        output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
        output_file.write(new_header + "\n")

        last_pos = -1
        for line in alignment_file:
            split_line = line.strip().split(",")
            pos = int(split_line[0])

            if (chromosome, pos) in SNP_positions:
                pos, features = create_features(split_line, pos_reference, SNP_positions[(chromosome, pos)])
            else:
                pos, features = create_features(split_line, pos_reference)

            if pos != (last_pos + 1):
                for i in range(last_pos + 1, pos):
                    if (i != 0) and (i % num_bases == 0):  # check to see if we should start a new chunk
                        output_file.close()
                        cur_chunk += 1

                        output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                                  "_features_binary.txt.gz", output_file_list, type = "gzip")
                        output_file.write("cell" + str(cur_chunk) + "\t" + chromosome + "\n")
                        output_file.write(new_header + "\n")

                    extra_features = [0, 2] * num_feats
                    output_file.write("\t".join([str(k) for k in extra_features]))
                    output_file.write("\n")

            if (pos != 0) and (pos % num_bases == 0):
                output_file.close()
                cur_chunk += 1

                output_file = create_file(output_directory + chromosome + "_" + str(cur_chunk) +
                                          "_features_binary.txt.gz", output_file_list, type="gzip")
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
                                              "_features_binary.txt.gz", output_file_list, type="gzip")
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
    parser.add_argument('-n', '--numBases', dest='num_bases', type=int,
                        help='Number of bases in each segment if doing no or full variation or number of variants in '
                             'each file if doing local variation.')
    parser.add_argument('-cl', '--chromLengths', required=True, dest='chromosome_sizes_file_name',
                        help='File containing the chrosomome lengths for the reference species in the alignment.')
    parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')
    parser.add_argument('-v' '--variationFile', dest='variation_file', help='File containing variants in the '
                        'reference genome (Optional). Format required is BED: chr, start, end (start + 1), rsID, '
                        'reference allele, alternate allele (if multiple separate by comma, but only first will be '
                        'kept')
    parser.add_argument('-lv', '--localVariationBinarize', type=int, dest='lv_size', help='If included only locations '
                        'within this many bases around the variants will be binarized')
    parser.add_argument('-d', '--dummyState', action='store_true', dest='dummy_state', help='Include if dummy state '
                        'should be used to separate variants.')
    parser.add_argument('-fv', '--fullVariation', action='store_true', dest='full_variation', help='Include to '
                        'generate binarized file for every single possible variant on this chromosome. Do not use '
                        'this flag lightly -- it will generate a lot of files. By default this uses a 10bp window '
                        'around each variant and the dummy state to separate variants. Will cause -lv and and -v '
                        'flags to be ignored')
    parser.add_argument('-fi', '--fullIndex', type=int, dest='full_index', help='The index of the current chunk of '
                        'maf sequence file for which to create all possible variants. Since the genome is so big '
                        'the output of parseMAF.py was split into manageable chunks, and this script requires the '
                        'index of the current file for some housekeeping.')

    args = parser.parse_args()
    if not args.num_bases:
        parser.error('The -n flag must be specified.')
    if (args.lv_size) and (not args.variation_file):
        parser.error('A variation file is required if the -lv flag is specified.')
    if (args.full_variation) and (args.full_index is None):
        parser.error('The -fi flag is required when running full variation (-fv flag)')

    binarize_alignment(args)
