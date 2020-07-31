#!/usr/bin/python3

import argparse
import sys
import re


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


def parse_description(line_input):
    '''Parses out elements of fasta headerline that correspond to
    "description:". If there is a match will return the string and
    return the headerline after removing this string from the headerline.
    If there is no match then will return empty string and unaltered
    headerline.'''

    # Initialize descsription as empty.
    description = ""

    # First parse out "description" string based on regular expression.
    description_match = re.match(r'.+( description:.+$)', line_input)

    # If there was a match above then parse with below commands.
    if description_match:

        description = description_match.group(1)

        # Remove this string from line (replace with empty string).
        line_input = line_input.replace(description, "")

        # Remove " description:" from start of string.
        description = description[13:]

    return(description, line_input)


def parse_gene_symbol(line_input):
    '''Parses out elements of fasta headerline that correspond to
    "gene_symbol:". If there is a match will return the string and
    return the headerline after removing this string from the headerline.
    If there is no match then will return empty string and unaltered
    headerline.'''

    # Initialize gene_symbol as empty.
    gene_symbol = ""

    # First parse out "gene_symbol" string based on regular expression.
    gene_symbol_match = re.match(r'.+( gene_symbol:.+$)', line_input)

    # If there was a match above then parse with below commands.
    if gene_symbol_match:

        gene_symbol = gene_symbol_match.group(1)

        # Remove this string from line (replace with empty string).
        line_input = line_input.replace(gene_symbol, "")

        # Remove " gene_symbol:" from start of string.
        gene_symbol = gene_symbol[13:]

    return(gene_symbol, line_input)


def main():

    parser = argparse.ArgumentParser(description="Reads in ENSEMBL FASTA file \
and parses out all information into a table. This script assumes that \
\n\n\
Here is an example fasta file header: \
\n\n\
>CDY60769 cds supercontig:AST_PRJEB5043_v1:LK033913:35160:37571:-1 \
gene:GSBRNA2T00030031001 gene_biotype:protein_coding \
transcript_biotype:protein_coding gene_symbol:BnaAnng17200D \
description:BnaAnng17200D protein \
[Source:UniProtKB/TrEMBL;Acc:A0A078J9M5]",

                                     epilog='''Usage example:

python3 parse_ENSEMBL_fasta_header.py -f \
Brassica_napus.AST_PRJEB5043_v1.cds.all.fa -o parsed_headers.tab

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)
    parser.add_argument("-o", "--output", metavar="OUTFILE", type=str,
                        help="Path to output file", required=True)

    parser.add_argument("--one_transcript", action="store_true",
                        help="Flag to indicate that only one line should \
be output per gene (assuming that numbers after \".\" in gene names are \
alternative transcripts.")

    args = parser.parse_args()

    # Make list with ordered columns of itnerest.
    info_columns = ["gene_name", "annotation_type", "position", "gene",
                    "gene_biotype", "gene_symbol", "description", "source"]

    # Make set of categories that should be ignored.
    columns2skip = set(["transcript_biotype"])

    # Open output file.
    outfile = open(args.output, "w")

    # Print out headerline of output table.
    print("\t".join(info_columns), file=outfile)

    # Create set that will contain all past genes if one_transcript flag set.
    if args.one_transcript:
        past_genes = set([])

    # Read in FASTA line-by-line.
    with open(args.fasta, "r") as fasta:

        for line in fasta:

            # Remove newline characters from end.
            line = line.rstrip("\n\r")

            # If not header line then skip.
            if line[0] != ">":
                continue

            # Otherwise need to parse header.

            # This is what an example headerline looks like (no newlines):
            # >CDY69013 cds supercontig:AST_PRJEB5043_v1:LK038450:1714:4174:-1
            # gene:GSBRNA2T00082019001 gene_biotype:protein_coding
            # transcript_biotype:protein_coding gene_symbol:BnaCnng61390D
            # description:BnaCnng61390D protein
            # [Source:UniProtKB/TrEMBL;Acc:A0A078JNA5]

            # Inititalize dictionary with columns of interest as keys.
            # Missing values will be output as empty strings.
            info = {}
            for column in info_columns:
                info[column] = ""

            # Parse out Source, description and gene_symbol based on regular
            # expressions. If the regex match the line then they will be
            # removed from the original line.
            info_source, line = parse_info_source(line)
            description, line = parse_description(line)
            gene_symbol, line = parse_gene_symbol(line)

            # Add this info to dictionary.
            info["source"] = info_source
            info["description"] = description
            info["gene_symbol"] = gene_symbol

            # Then split entire line on spaces.
            line_split = line.split(" ")

            # Remove first element (gene name), but ignore the > character.
            gene_name = line_split.pop(0)[1:]

            # Remove transcript number if --one_transcript set and check if
            # already seen.
            if args.one_transcript:
                gene_name = gene_name.split(".")[0]

                # Skip line if already seen.
                if gene_name in past_genes:
                    continue

                # Otherwise add to set.
                past_genes.add(gene_name)

            # Add gene name to info dictionary.
            info["gene_name"] = gene_name

            # Remove annotation type as well (next element)
            info["annotation_type"] = line_split.pop(0)

            # Also remove position (next element).
            info["position"] = line_split.pop(0)

            # Loop over remaining info elements.
            for info_element in line_split:

                # Split on ":" and get column name from first element.
                info_split = info_element.split(":")

                # Check if one of the categories that should be skipped.
                if info_split[0] in columns2skip:
                    continue

                # Otherwise if known column name then add to dictionary.
                elif info_split[0] in info:
                    info[info_split[0]] = ":".join(info_split[1:])

                # If not in this known set then throw warning and skip.
                else:
                    print("Warning: element " + info_split[0] +
                             " is not " +
                             "one of the expected columns: " + 
                             " ".join(info_columns) + "." +
                             " This is the full line:\n\n" + line, 
                             file=sys.stderr)

            # Print out info from headerline to output table.
            outline = []
            for column in info_columns:
                outline = outline + [info[column]]
            
            print("\t".join(outline), file=outfile)

if __name__ == '__main__':
    main()
