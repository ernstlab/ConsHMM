import gzip
import argparse
from Bio.Seq import Seq
import importlib
#from shared import *
import sys
import os.path
from collections import defaultdict

def open_file(file_name, mode):
    if file_name[-3:] == ".gz":
        file_handle = gzip.open(file_name, mode)
    else:
        file_handle = open(file_name, mode)
    return file_handle


def write_header(output_file, print_pos, species):
    if not print_pos: # if only interested in position in reference species
        output_file.write('chr_reference,pos_reference')
        for sp in species:
            output_file.write("," + sp)
        output_file.write("\n")
    else: # if interested in positions in all other species in the alignment
        num_species = 0
        for species_name in species:
            num_species += 1
            output_file.write('pos_' + species_name + ',strand_' + species_name + ',nuc_' + species_name)
            if num_species < len(species):
                output_file.write(',')
        output_file.write('\n')


def parse_maf(main_args):
    debug_mode = False

    species_file = open_file(main_args.species_file, 'r')
    reference_species = main_args.reference_species
    print_pos = main_args.print_pos
    allow_insertions = main_args.allow_insertions
    output_dir = main_args.output_dir
    if main_args.debug_file is not None:
        debug_file = open_file(main_args.debug_file, 'wt')
        debug_mode = True

    # Read species list
    species = []
    for line in species_file:
        species.append(line.strip())
    species = sorted(species)

    maf_files = []
    output_files = {}

    if main_args.command == 'singleChromosome': # Only reading one file and writing one file
        num_chr = 1
        maf_files = [open_file(main_args.maf_file, 'rt')]
    elif main_args.command == 'fullGenome': # Have to iterate over the whole genome
        chrom_needed = main_args.cur_chrom
        maf_dir = main_args.maf_dir
        output_dir = main_args.output_dir

        # read list of chromosome names
        chrom_list = open_file(main_args.chrom_list, 'rt')
        for chrom in chrom_list:
            maf_files.append(open_file(maf_dir + '/' + chrom.strip() + '.maf.gz', 'rt'))
    else:
        sys.exit('This should never happen. Argparse must have had an anneurysm.')

    count_files_done = 0
    count_blocks_done = 0
    for f in maf_files:
        # Parse MAF file
        if debug_mode:
            debug_file.write('Starting MAF file ' + str(count_files_done + 1) + '.\n')
            debug_file.flush()

        score_line = f.readline()
        while score_line[0] != "a":  # skip lines until we get to the start of an alignment block
            score_line = f.readline()
        
        covered = {}  # dictionary of bases covered to account for blocks convering duplicate pieces
        count_blocks_skipped = 0
        while score_line != "": # should be on an alignment line when we first get here
            sequence_line = f.readline()
            # skip all other lines at the beginning of block until the first sequence line
            
            while sequence_line[0] != 's':
                sequence_line = f.readline()
            cur_alignment_block = {}
            cur_alignment_raw = []

            reverse_complement_ref = False

            # parse an alignment block
            chrom_species = {} # will store chromosome for each species in the alignment
            start_species_seq = {} # will store start position in every species on positive strand
            end_species_seq = {} # will store end position in every species on positive strand
            strand_species = {} # will store on which strand the aligned sequence is
            cur_pos_per_species = {} # will store where we are in each species when outputting position

            # keep reading until we get to the next block or to the end of the file
            while sequence_line != "" and sequence_line[0] != 'a': # should be on a sequence line when we first get here
                cur_alignment_raw.append(sequence_line)
                split_seq = sequence_line.split()

                if len(split_seq) == 0:
                    line_type = 'e'
                else:
                    line_type = split_seq[0]

                if line_type != 's':
                    sequence_line = f.readline()
                    continue

                cur_species = split_seq[1].split('.')[0]
                cur_chromosome = '.'.join(split_seq[1].split('.')[1:])
                if cur_chromosome[0] != 'c':
                    cur_chromosome = 'chr' + cur_chromosome

                true_len_species = int(split_seq[3]) # number of nucleotides excluding '-' character for gaps
                chrom_size = int(split_seq[5])
                chrom_species[cur_species] = cur_chromosome # store chromosome
                start_species_seq[cur_species] = int(split_seq[2]) # store start position
                end_species_seq[cur_species] = int(split_seq[2]) + true_len_species - 1
                strand_species[cur_species] = '+' # store strand
                cur_alignment_block[cur_species] = split_seq[6].upper() # store actual sequence

                if split_seq[4] == '-':  # reverse complement so we have to flip coordinates to + strand
                    start_species_seq[cur_species] = chrom_size - int(split_seq[2]) - true_len_species
                    end_species_seq[cur_species] = chrom_size - int(split_seq[2]) - 1
                    strand_species[cur_species] = '-'
                    if cur_species == reference_species:
                        reverse_complement_ref = True

                sequence_line = f.readline()
            
            score_line = sequence_line
            count_blocks_done += 1

            # when done with an alignment block go through and output whether each species is the same as human
            # skiping blocks that don't contain the species or chromosome of interest
            if main_args.command == 'singleChromosome':
                chrom_needed = chrom_species[reference_species]
            if (reference_species not in cur_alignment_block) or (chrom_species[reference_species] != chrom_needed):
                count_blocks_skipped += 1
                continue

            pos = start_species_seq[reference_species]

            rev = []
            # Reverse complement sequences and adjust start and end coordinates accordingly
            for sp in species:
                if sp in cur_alignment_block:
                    if sp == reference_species:
                        cur_pos_per_species[sp] = start_species_seq[sp]
                        if reverse_complement_ref:
                            cur_alignment_block[sp] = Seq(cur_alignment_block[sp].upper()).reverse_complement() 
                    else:
                        if reverse_complement_ref: # if negative strand in reference
                            rev.append(sp) 
                            # if positive strand in other species go from start to end on positive strand for coordinates
                             # but reverse complement the sequence
                            cur_alignment_block[sp] = Seq(cur_alignment_block[sp].upper()).reverse_complement()
                            if strand_species[sp] == '+':
                                cur_pos_per_species[sp] = end_species_seq[sp]
                                strand_species[sp] = '-'
                            else:
                                # otherwise go from start to end because they are both on negative strand
                                cur_pos_per_species[sp] = start_species_seq[sp]
                        else:
                            if strand_species[sp] == '+':
                                cur_pos_per_species[sp] = start_species_seq[sp]
                            else:
                                cur_pos_per_species[sp] = end_species_seq[sp]

            for i in range(len(cur_alignment_block[reference_species])):
                cur_nucleotide = cur_alignment_block[reference_species][i]
                if (cur_nucleotide == '-') or (cur_nucleotide == '.'):
                    if not allow_insertions:
                        continue

                if chrom_species[reference_species] not in output_files:
                    output_files[chrom_species[reference_species]] = open_file(output_dir + '/' + chrom_species[
                        reference_species] + '_maf_sequence.csv.gz', 'wt')
                    write_header(output_files[chrom_species[reference_species]], print_pos, species)

                o = output_files[chrom_species[reference_species]]
                if pos not in covered:
                    if not print_pos: # if only interested in position in reference species
                        # write position in reference species
                        o.write(chrom_species[reference_species] + ',')
                        if cur_alignment_block[reference_species][i] in ['-', '.']:
                            o.write('d')
                        else:
                            o.write(str(start_species_seq[reference_species]))
                            start_species_seq[reference_species] += 1
                            covered[pos] = 1
                            pos += 1
                        if (start_species_seq[reference_species] - 1) == 230245683:
                            print(len(species))
                            for sp in species:
                                print(sp)
                        for sp in species:
                            if (sp not in cur_alignment_block) or (cur_alignment_block[sp][i] in ['-', '.']):
                                o.write(',-')
                            else:
                                o.write(',' + str(cur_alignment_block[sp][i]))
                        o.write('\n')
                    else: # if we want to print position and strand for every other species
                        num_species = 0

                        for sp in species:
                            num_species += 1
                            if sp not in cur_alignment_block:
                                o.write('X,X,X') # this means species does not have an alignment at all in this block
                            else:
                                if cur_alignment_block[sp][i] in ['-', '.']: # species has a gap
                                    # all gaps converted to '-' even though can also be '.'
                                    o.write(chrom_species[sp] + ':d,X,-')
                                else: # species has alignment
                                    #if sp == reference_species: print("normal")
                                    o.write(chrom_species[sp] + ':' + str(cur_pos_per_species[sp]) + \
                                            ',' + str(strand_species[sp] + ',' + str(cur_alignment_block[sp][i])))
                                    if reverse_complement_ref:
                                        if strand_species[sp] == '+':
                                            cur_pos_per_species[sp] -= 1
                                        else:
                                            cur_pos_per_species[sp] += 1
                                    else:
                                        if strand_species[sp] == '+':
                                            cur_pos_per_species[sp] += 1
                                        else:
                                            cur_pos_per_species[sp] -= 1
                                    if sp == reference_species:
                                        covered[pos] = 1
                                        pos += 1

                            if num_species < len(species):
                                o.write(',')
                        o.write('\n')

        count_files_done += 1
        if debug_mode:
            debug_file.write('Done with ' + str(count_files_done) + ' MAF files.\n')
            debug_file.flush()

    # Close output files
    for i in output_files:
        output_files[i].close()

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(description="""parse MAF.py takes in a MAF file that contains a multiple 
    sequence alignment for one of the chromosomes of the reference species in the alignment, and outputs a file 
    containing just the sequence information. The output file will have as many rows as the number of bases in the input
    chromosome, and N + 1 columns, where N is the number of species in the multiple sequence alignment. The first
    column contains the coordinate of a base in the reference species. The following columns contain the nucleotide
    at that position in the alignment for each species including the reference. Indels and 'N' nucleotides are
    both encoded with 'X'.""", add_help=False)

    parent_parser.add_argument('-s', '--speciesFile', required=True, dest='species_file',
                        help='File containing a list of the names of all the species in the alignment.')
    parent_parser.add_argument('-r', '--refSpecies', required=True, dest='reference_species',
                        help='Genome assembly name of the reference species in the alignment.')
    parent_parser.add_argument('-p', '--printPos', action='store_true', dest='print_pos',
                        help='If included, positions of each base in every species will be printed, instead of '
                             'just the positions in the reference species')
    parent_parser.add_argument('-i', '--allowInsertion', dest='allow_insertions', action='store_true',
                        help='If included, insertions in other species will not be dropped in the output file, '
                             'so there will be positions where the reference species will have a gap. Do not '
                             'use this option if running for ConsHMM.')
    parent_parser.add_argument('-o', '--outputDirectory', required=True, dest='output_dir',
                        help='Output directory to write the parsed files separated by chromosome to.')
    parent_parser.add_argument('-d', '--debugMode', dest='debug_file', help='If a file is specified, the code will '
                              'print debug output to the specified file ')


    main_parser = argparse.ArgumentParser()
    subparsers = main_parser.add_subparsers(help='Mode in which to run the script.', dest='command', required=True)

    # For when running on files that are already split by chromosome of the reference species
    single_chromosome = subparsers.add_parser('singleChromosome', help='Process a single file that contains the '
                                              'alignment for a full chromosome of the reference species',
                                              parents = [parent_parser])
    single_chromosome.add_argument('-m', '--MAFFile', required=True, dest='maf_file',
                        help='.maf file containing a multi-sequence alignment.')

    full_genome = subparsers.add_parser('fullGenome', help='Process multiple files that cover the whole genome of a '
                                        'reference species, but the files are split up by a different species in the '
                                        'alignment. Note that the resulting files in this case will not be sorted by '
                                        'the position in the reference species, but they will be separated by '
                                        'chromosomes', parents = [parent_parser])
    full_genome.add_argument('-m', '--MAFDirectory', required=True, dest='maf_dir',
                        help='Directory containing .maf files named chr*.maf.gz, where the * is the number '
                             'of a chromosome in the reference species, containing a multi-sequence alignment.')
    full_genome.add_argument('-c', '--chromList', required=True, dest='chrom_list',
                        help='File containing the names of the chromosomes in the reference of the multi-sequence '
                             'alignment, one per line.')
    full_genome.add_argument('-cc', '--curChr', required=True, dest='cur_chrom',
                             help='Which chromosome of the species of interest to look for. This way the script will '
                                  'skip blocks that contain a different chromosome.')


    main_args = main_parser.parse_args()
    parse_maf(main_args)
