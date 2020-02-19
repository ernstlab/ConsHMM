import argparse
import gzip
from collections import defaultdict

def parse_segmentation(segmentation_file):
    state_freqs = defaultdict(int)

    for line in segmentation_file:
        split_line = line.strip().split('\t')
        start = int(split_line[1])
        end = int(split_line[2])
        state = split_line[3][1:]

        state_freqs[state] += (end - start)

    return state_freqs

def update_model(model_file, state_freqs, total_bases, output_file):
    for state in state_freqs:
        state_freqs[state] /= float(total_bases)

    for line in model_file:
        split_line = line.strip().split("\t")
        if split_line[0] == "probinit":
            state = int(split_line[1])
            split_line[-1] = str(state_freqs[state])
            output_file.write("\t".join(split_line) + "\n")
        else:
            output_file.write(line)

    model_file.close()
    output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes in a ConsHMM model and a genome wide ConsHMM segmentation '
                                     'and outputs an updated parameter set where the inital state parameters are '
                                     'replaced by the genome wide frequency of each state in the segmentation.',
                                     allow_abbrev=False)
    parser.add_argument('-m', '--modelFile', required=True, dest='model_file_name', help='ConsHMM model file')
    parser.add_argument('-s', '--segmentationFile', required=True, dest='segmentation_file', help='File containing'
                        ' the genome wide ConsHMM segmentation.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file_name', help='Output file where the '
                        'model with the updated parameters will be written')

    args = parser.parse_args()
    model_file = gzip.open(args.model_file_name, 'rt')
    segmentation_file = gzip.open(args.segmentation_file_name, 'rt')
    output_file = gzip.open(args.output_file_name, 'wt')

    state_freqs = parse_segmentation(segmentation_file)
    total_bases = sum(state_freqs.values())

    update_model(model_file, state_freqs, total_bases, output_file)

