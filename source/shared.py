import os
import numpy as np
import gzip

def parse_variation_file(variation_file_name, chromosome = None):
    # Expected file format is BED: chr, start, end (just start + 1), rsID, refAllele, altAllele (if more than 1 only
    # 1 will be kept.
    if variation_file_name[-3:] == ".gz":
        variation_file = gzip.open(variation_file_name, "rt")
    else:
        variation_file = open(variation_file_name, "r")

    SNP_positions = {}
    SNP_lines = {}

    for line in variation_file:
        split_line = line.strip().split("\t")
        cur_chr = split_line[0]
        pos = int(split_line[1])
        alt_allele = split_line[5].split(",")[0]
        if (not chromosome) or (cur_chr == chromosome):
            if (cur_chr, pos) in SNP_positions:
                print(cur_chr, pos)
            SNP_positions[(cur_chr, pos)] = alt_allele
            SNP_lines[(cur_chr, pos)] = line.strip()

    return SNP_positions, SNP_lines

def create_features(split_line, pos_reference, alt_allele = None):
    line_features = np.empty((len(split_line) - 3) * 2)  # two features for each species
    pos = int(split_line[1])
    # Compare other species to the nucleotide in the reference genome, or to the variant of that nucleotide
    if alt_allele:
        reference = alt_allele
    else:
        reference = split_line[pos_reference]

    idx = 0
    for i in range(2, len(split_line)):
        if i != pos_reference:
            if split_line[i] == '-':
                line_features[idx] = 0
                line_features[idx + 1] = 2
            else:
                line_features[idx] = 1
                if split_line[i] == reference:
                    line_features[idx + 1] = 1
                else:
                    line_features[idx + 1] = 0
            idx += 2

    return pos, line_features

def format_dir(directory):
    if not os.path.isdir(directory):
        print(directory, "directory not found. Attempting to create . . .")
        try:
            os.makedirs(directory, exist_ok=False)
        except OSError:
            print("Could not create ", directory, ".")
            exit(1)
        print("Done.")
    if directory[-1] != '/':
        return directory + '/'
    return directory


def print_dict_keys(dictionary):
    for key in dictionary.keys():
        print(key + " ")

if __name__ == '__main__':
    print ("Loaded shared.py module.")
