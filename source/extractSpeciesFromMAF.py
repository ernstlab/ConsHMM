## Authors: Brooke Felsheim, Adriana Arneson

import re
import gzip
import argparse

def extract_species(main_args):
    maf_file = main_args.maf_file
    num_species = main_args.num_species
    output = main_args.output

    f = gzip.open(maf_file, 'rt')

    species_set = set()

    line = f.readline()
    while line != "":
        if len(species_set) == num_species:
            break
        if line[0] != "s":
            line = f.readline()
            continue
        split_seq = line.split()
        species_name = re.sub('\..*', '', split_seq[1])
        if species_name != 'ancestral_sequence' and species_name != 'ancestral_sequences':
            species_set.add(re.sub('\..*', '', split_seq[1]))
        line = f.readline()

    with open(output, 'w') as o:
        for n in species_set:
            print(n, file=o)

    if len(species_set) != num_species:
        print("Warning: number of species found in MAF file is less than the anticipated number of species in the "
              "alignment")
        print("Number of expected species:", num_species, "Number of species found:", len(species_set))

    if len(species_set) == num_species:
        print(num_species, "species were expected and" , len(species_set) , "species were extracted from the MAF "
              "file.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extracts species list from a MAF file, given the expected number of'
                                    'species in the alignment and the name of the output file. Prints a warning message'
                                    'if less than the number of expected species are found. Note that if the expected '
                                    'number of species in the alignment is smaller than the number of species in the '
                                    'MAF file, it will extract only up to the number of expected species for '
                                    'efficiency.')
    parser.add_argument('-m', '--MAFFile', required=True, dest='maf_file', help='.maf.gz file containing multiple '
                        'species alignment.')
    parser.add_argument('-n', '--numSpecies', required=True, dest='num_species', type=int, help='Number of species in '
                        'the multiple species alignment.')
    parser.add_argument('-o', '--outputFile', required=True, dest= 'output', help='name of output file to write species'
                        'list to.')

    args = parser.parse_args()

    extract_species(args)
