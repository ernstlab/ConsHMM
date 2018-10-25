import argparse


def extract_species_names(main_args):
    nh_file = open(main_args.nh_file, 'r')
    output_file = open(main_args.output_file, 'w')
    for line in nh_file:
        parsed_line = line.replace("(", "")
        parsed_line = parsed_line.replace(")", "")
        parsed_line = parsed_line.replace('\n', "")
        parsed_line = parsed_line.split(",")

        for i in parsed_line:
            parsed_element = i.split(":")
            for j in parsed_element:
                if (len(j) > 0) and (not j[0].isdigit()):
                    output_file.write(j.rstrip().lstrip() + "\n")

    output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""extractSpeciesNames.py attempts to extract the species names from a
        .nh file describing a phylogenetic tree. Output should output should be checked manually, as .nh files are
        not well standardized""")

    parser.add_argument('-n', '--nhFile', required=True, dest='nh_file',
                        help='.nh file describing a phylogenetic tree.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file to write the list of species to.')

    args = parser.parse_args()

    extract_species_names(args)
