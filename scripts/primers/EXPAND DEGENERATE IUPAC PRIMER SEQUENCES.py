#!/usr/bin/env python3

############################################################
# EXPAND DEGENERATE IUPAC PRIMER SEQUENCES
############################################################

# This script reads a FASTA file containing primer sequences with
# degenerate IUPAC nucleotide codes and generates all possible
# non-degenerate sequence combinations.
#
# For each expanded sequence, the FASTA header is modified by adding
# a tag corresponding to the bases selected at degenerate positions.
#
# Usage:
# python3 expand_iupac_primers.py <INPUT_PRIMER_FASTA>
#
# Output:
# expanded_primers.fasta


import sys
from itertools import product


############################################################
# IUPAC nucleotide degeneracy dictionary
############################################################

iupac_degenerate_bases = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],

    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],

    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],

    "N": ["A", "C", "G", "T"]
}


############################################################
# Expand one IUPAC sequence
############################################################

def expand_iupac_sequence(sequence):
    """
    Generate all possible non-degenerate sequences from a sequence
    containing IUPAC degenerate nucleotide codes.

    Returns:
    - tag: bases selected at the degenerate positions only
    - expanded_sequence: complete expanded nucleotide sequence
    """

    sequence = sequence.upper()

    base_options = [
        iupac_degenerate_bases.get(base, [base])
        for base in sequence
    ]

    for combination in product(*base_options):

        expanded_sequence = "".join(combination)

        tag_characters = [
            combination[index]
            for index, base in enumerate(sequence)
            if (
                base in iupac_degenerate_bases
                and len(iupac_degenerate_bases[base]) > 1
            )
        ]

        tag = "".join(tag_characters)

        yield tag, expanded_sequence


############################################################
# Parse FASTA file
############################################################

def parse_fasta(fasta_path):
    """
    Read a FASTA file and return a list of records.

    Each record is returned as:
    - header
    - sequence
    """

    records = []

    with open(fasta_path, "r") as infile:

        header = None
        sequence_lines = []

        for line in infile:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records.append((header, "".join(sequence_lines)))

                header = line[1:]
                sequence_lines = []

            else:
                sequence_lines.append(line)

        if header is not None:
            records.append((header, "".join(sequence_lines)))

    return records


############################################################
# Main function
############################################################

def main():

    if len(sys.argv) != 2:
        print("Usage: python3 expand_iupac_primers.py <INPUT_PRIMER_FASTA>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = "expanded_primers.fasta"

    records = parse_fasta(input_fasta)

    with open(output_fasta, "w") as outfile:

        for header, sequence in records:

            for tag, expanded_sequence in expand_iupac_sequence(sequence):

                if tag:
                    expanded_header = f"{header}_{tag}"
                else:
                    expanded_header = header

                outfile.write(f">{expanded_header}\n")
                outfile.write(f"{expanded_sequence}\n")

    print(f"Completed. File saved as: {output_fasta}")


if __name__ == "__main__":
    main()