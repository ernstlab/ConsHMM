from shared import *
import gzip
import argparse


def split_maf_by_chr(main_args):
    maf_file = gzip.open(main_args.maf_file, 'rt')
    output_directory = format_dir(main_args.output_directory)
    reference_species = main_args.reference_species
    chromosome_list_file = open(main_args.chromosome_list_file_name, 'r')

    print("Creating output files (one for each chromosome) . . .",)
    output_files = {}
    for line in chromosome_list_file:
        cur_chr = line.split()
        output_file = gzip.open(output_directory + cur_chr[0] + ".maf.gz", "wt")
        output_files[cur_chr] = output_file
    print("Done.")

    print("Splitting up MAF file per chromosome . . .")
    cur_block = []
    cur_chr = ""
    skipped_chromosomes = set()
    for line in maf_file:
        if line[0] == "#":  # comment lines will get discarded
            continue

        if line == "\n":
            if cur_chr not in output_files:
                skipped_chromosomes.add(cur_chr)
            else:
                for curLine in cur_block:
                    output_files[cur_chr].write(curLine)
                output_files[cur_chr].write("\n")
            cur_block = []
            cur_chr = ""
        else:
            if line[0] == "s":
                split_line = line.strip().split()
                cur_species = split_line[1].split(".")[0]
                if cur_species == reference_species:
                    cur_chr = split_line[1].split(".")[1]
            cur_block.append(line)
    print("Done.")

    print("Closing files . . .")
    for output_file in output_files:
        output_file.close()
    print("Done.")

    print("Other chromosomes found in MAF file: ")
    for chromosome in skipped_chromosomes:
        print(chromosome, end="")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""splitMafByChr.py splits multiple sequence alignmenst that are
        provided in one single MAF files into multiple MAF files -- one for each chromosome in the reference genome.""")

    parser.add_argument('-m', '--MAFFile', required=True, dest='maf_file',
                        help='.maf file containing a multi-sequence alignment.')
    parser.add_argument('-o', '--outputDirectory', required=True, dest='output_directory',
                        help='Output directory in which to put the chromosome files.')
    parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')
    parser.add_argument('-cl', '--chromList', required=True, dest='chromosome_list_file_name',
                        help='File containing a list of all the chromosomes in the reference species genome for which '
                             'you want to generate a separate MAF file.')

    args = parser.parse_args()
