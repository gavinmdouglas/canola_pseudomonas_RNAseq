#!/usr/bin/python3

import argparse


def main():

    parser = argparse.ArgumentParser(
                 description="Script to read through "
                             "mmquant files sum abundances by AT genes.",

                 formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="FILE", type=str,
                        help="Path to input mmquant file", required=True)

    parser.add_argument("-m", "--mapfile", metavar="FILE", type=str,
                        help="Path to mapfile", required=True)

    parser.add_argument("-b", "--background", metavar="FILE", type=str,
                        help="Path to background list of genes", required=True)

    parser.add_argument("-o", "--output", metavar="FILE", type=str,
                        help="Path to output file", required=True)

    args = parser.parse_args()

    # Initialize dictionary to map B. napus ids to A. thaliana
    Bn2At = {}

    lc = 0
    with open(args.mapfile, "r") as map_in:
        for line in map_in:
            if lc == 0:
                lc += 1
                continue

            line = line.rstrip("\n\r")
            line_split = line.split("\t")

            # Skip if At gene is "NA"
            if line_split[5] == "NA":
                continue

            Bn2At[line_split[1]] = line_split[5]

    # Initialize dictionary that contains counts per At homolog.
    At_counts = {}

    with open(args.background, "r") as genes_in:
        for line in genes_in:
            At_counts[line.rstrip("\n\r")] = 0

    lc = 0

    with open(args.input, "r") as mmquant_in:
        for line in mmquant_in:
            if lc == 0:
                lc += 1
                continue

            line = line.rstrip("\n\r")
            line_split = line.split("\t")

            # Split gene ids if there are multiple
            Bn_ids = [line_split[0]]

            if "--" in Bn_ids[0]:
                Bn_ids = Bn_ids[0].split("--")

            # Get unique At ids.
            At_hits = set()
            for Bn_id in Bn_ids:
                if Bn_id in Bn2At:
                    At_hits.add(Bn2At[Bn_id])

            for At_hit in At_hits:
                At_counts[At_hit] += int(line_split[1])

    mmquant_out_fh = open(args.output, "w")

    for At, counts in At_counts.items():
        print(At + "    " + str(counts), file=mmquant_out_fh)

    mmquant_out_fh.close()


if __name__ == '__main__':
    main()
