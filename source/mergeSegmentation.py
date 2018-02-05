import glob
import gzip
import argparse
from shared import *


def merge_segmentation(main_args):
    bed_directory = format_dir(main_args.bed_directory)
    output_file = gzip.open(main_args.output_file, 'wt')
    num_bases = main_args.num_bases
    num_states = main_args.num_states
    chromosome = main_args.chromosome

    offset = 0
    file_names = glob.glob(bed_directory + '*_' + chromosome + '_segments.bed')
    max_idx = 0
    for file_name in file_names:
        idx = int(file_name.split("_")[-4].split("cell")[-1])
        if idx > max_idx:
            max_idx = idx

    last_line = ""
    for idx in range(0, max_idx + 1):
        cur_file = bed_directory + "cell" + str(idx) + "_" + str(num_states) + "_" + chromosome + "_segments.bed"
        print("Joining file ", cur_file)

        try:
            bed_file = open(cur_file, 'r')
        except IOError:
            print(cur_file, "not found in the directory you provided for the '-b, --bedDirectory flag'.")
            exit(1)

        for line in bed_file:
            split_line = line.split("\t")
            split_line[1] = str(int(split_line[1]) + offset)
            split_line[2] = str(int(split_line[2]) + offset)

            if last_line != "":
                if last_line[0] == split_line[0]:  # if on the same chromosome
                    if last_line[-1] == split_line[-1]:  # if in the same state
                        last_line[2] = split_line[2]
                    else:
                        output_file.write("\t".join(last_line))
                        last_line = split_line
                else:
                    output_file.write("\t".join(last_line))
                    last_line = split_line
            else:
                last_line = split_line
        offset += num_bases

    output_file.write("\t".join(last_line))
    output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merge ConsHMM segmentation chunks into one segmentation for one
        chromosome. Chunk files are expected to be named cell*_num_states_chr*_segments.bed. This is the naming
        MakeSegmentation will output if run with -i chr* parameter.""")

    parser.add_argument('-b', '--bedDirectory', required=True, dest='bed_directory',
                        help='Directory containing the bed files that need to be merged.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file to contain the merged segmentation.')
    parser.add_argument('-n', '--numBases', required=True, dest='num_bases', type=int,
                        help='Number of bases in each segment.')
    parser.add_argument('-s', '--numStates', required=True, dest='num_states', type=int,
                        help='Number of states in the segmentation.')
    parser.add_argument('-c', '--chr', required=True, dest='chromosome',
                        help='Chromosome of the reference species contained in the input file.')

    args = parser.parse_args()

    merge_segmentation(args)
