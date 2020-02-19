'''
This script is necessary because of the overlap at the end of chromosome chunks and existing gaps in the MAF files,
so we can't just concatenate all the pieces.
'''

import argparse
import gzip

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge variant segmentation files for one chromosome.',
                                     allow_abbrev=False)
    parser.add_argument('-c', '--chromosome', required=True, dest='chromosome', help='Current chromosome.')
    parser.add_argument('-d', '--directory', required=True, dest='directory', help='Directory where the variant '
                                                                                   'segmentations are stored.')
    parser.add_argument('-n', '--numFiles', required=True, dest='numFiles', type=int, help='Number of files for the '
                                                                                           'current chromosome.')
    parser.add_argument('-o', '--output', required=True, dest='output_file', help='Output file for the merged '
                                                                                  'segmentation. Must be .gz')

    args = parser.parse_args()
    chromosome = args.chromosome
    directory = args.directory
    numFiles = args.numFiles
    output = gzip.open(args.output_file, 'w')

    seen_bases = {}

    for i in range(numFiles):
        cur_file = gzip.open(directory + '/' + chromosome + '_' + format(i, '02') + '_variant_segmentation.txt.gz')
        for line in cur_file:
            split_line = line.strip().split()
            cur_pos = split_line[0]

            if cur_pos not in seen_bases:
                output.write(line)
                seen_bases[cur_pos] = 1
        cur_file.close()

    output.close()
