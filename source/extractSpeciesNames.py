import sys
import argparse

def main():

    if len(sys.argv) < 3:
        print("Usage: python extractSpeciesNames.py <.nh species file> <output file>")
        exit(1)

    speciesFile = open(sys.argv[1], 'r')
    oFile = open(sys.argv[2], 'w')
    for line in speciesFile:
        parsedLine = line.replace("(", "")
        parsedLine = parsedLine.replace(")", "")
        parsedLine = parsedLine.replace('\n', "")
        parsedLine = parsedLine.split(",")

        for i in parsedLine:
            parsed_element = i.split(":")
            for j in parsed_element:
                if (len(j) > 0) and (not j[0].isdigit()):
                    oFile.write(j.rstrip().lstrip() + "\n")
        
    oFile.close()

main()
