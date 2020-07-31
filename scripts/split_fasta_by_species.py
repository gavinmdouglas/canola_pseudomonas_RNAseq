#!/usr/bin/python3

import argparse
import sys
import re
import gzip
import os


def parse_info_source(line_input):
    '''Parses string matching [Source.+] at end of fasta headerlines.
    If there is a match will return the string and return the headerline
    after removing this string from the headerline. If there is no match then
    will return empty string and unaltered headerline.'''

    # Initialize info_source as empty.
    info_source = ""

    # First parse out "Source" string based on regular expression.
    info_source_match = re.match(r'.+( \[Source:.+\])', line_input)

    # If there was a match above then parse with below commands.
    if info_source_match:

        info_source = info_source_match.group(1)

        # Remove this string from line (replace with empty string).
        line_input = line_input.replace(info_source, "")

        # Check that last element (info_source) starts and ends with
        # closed brackets (ignoring first character, which should be a space).

        if info_source[1] != "[" or info_source[-1] != "]":
            sys.exit("Stopping: source info needs to start and end with "
                     "closed brackets, instead got:" + info_source)

        elif "Source:" not in info_source:
            sys.exit("Stopping: source info does not contain the string "
                     "\"Source:\", instead got:" + info_source)

        # Remove leading space, closed brackets, and "Source:" from
        # info_source string.
        info_source = info_source[9:]
        info_source = info_source[0:-1]

    return(info_source, line_input)


def main():

    parser = argparse.ArgumentParser(description="Reads in gzipped FASTA file and \
prints out each different species' sequences of a specified genus to \
different files. NOTE: This script assumes the genus and species will be the second\
and third fields in the headerline after splitting by the specified delimiter",

epilog='''Usage example:

python3 split_fasta_by_species.py -f \
current_Bacteria_aligned.fa -o out_folder -g Pseudomonas

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to GZIPPED FASTA file", required=True)
    parser.add_argument("-o", "--output", metavar="OUT FOLDER", type=str,
                        help="Path to output folder", required=True)
    parser.add_argument("-g", "--genus", metavar="STRING", type=str,
                        help="Genus to split into fastas for each species",
                        required=True)
    parser.add_argument("-d", "--delimiter", metavar="STRING", type=str,
                        help="String to separate header fields. Splits by \
                        whitespace by default", default="")

    args = parser.parse_args()

    # Check if outfolder exists and create if not.
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Initialize empty dictionary to contain filehandles.
    fh = {}
    last_species = None

    # Inititalize flag for whether seq should be printed.
    print_seq = False

    # Keep track of number of 16S per species
    counts = {}

    # Read in FASTA line-by-line.
    with gzip.open(args.fasta, "rt") as fasta:

        for line in fasta:

            # Remove newline characters from end.
            line = line.rstrip("\n\r")

            # Check if header.
            if line[0] == ">":

                # Split line on delimiter if set.
                if args.delimiter:
                    line_split = line.split(args.delimiter)
                else:
                    line_split = line.split()

                # Do not print seq and skip if not genus.
                if args.genus != line_split[1]:
                    print_seq = False
                    continue

                species = line_split[2]

                # If species ends in ";" then remove.
                if species[-1] == ";":
                    species = species[0:len(species)-1]

                last_species = species

                # Check if species is "sp." or "sp"
                if species == "sp." or species == "sp":
                    print_seq = False
                    continue
                else:
                    print_seq = True

                if species not in fh:
                    outfile = args.output + "/" + species + ".fa"
                    fh[species] = open(outfile, "w")
                    counts[species] = 0

                counts[species] += 1
                print(">" + args.genus + "_" + species + "_" + 
                      str(counts[species]), file=fh[species])
                

            # Otherwise print out if sequence line and header was genus.
            elif print_seq:
                print(line, file=fh[last_species])

    # Loop through and close all filehandles.
    for file in fh.values():
        file.close()


if __name__ == '__main__':
    main()
