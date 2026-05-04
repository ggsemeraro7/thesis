#!/usr/bin/env python3

############################################################
# EXTRACT AMPLICON SEQUENCES FROM FASTPCR OUTPUT
############################################################

# This script parses a FastPCR in silico primer search output file
# and extracts all detected amplicon sequences into a FASTA file.
#
# The FASTA header includes:
# - target sequence identifier
# - forward primer name and binding position
# - reverse primer name and binding position
# - amplicon length


import re


############################################################
# Input and output files
############################################################

input_file = "<FASTPCR_OUTPUT_TXT>"
output_fasta = "<AMPLICON_OUTPUT_FASTA>"


############################################################
# Initialize counter
############################################################

count_amplicon_hits = 0


############################################################
# Read FastPCR output
############################################################

with open(input_file, "r") as infile:
    lines = infile.readlines()


############################################################
# Initialize variables
############################################################

current_target = None

forward_primer = None
forward_position = None

reverse_primer = None
reverse_position = None

amplicon_length = None
amplicon_sequence = []

inside_amplicon = False


############################################################
# Parse FastPCR output and write amplicons to FASTA
############################################################

with open(output_fasta, "w") as outfile:

    for line in lines:
        line = line.rstrip()


        ####################################################
        # Detect a new target sequence section
        ####################################################

        match = re.match(r"In silico Primer\(s\) search for:\s+(.+)", line)

        if match:

            # Save the previous amplicon, if present
            if inside_amplicon and amplicon_sequence:
                header = (
                    f">{current_target}"
                    f"|{forward_primer}:{forward_position}"
                    f"|{reverse_primer}:{reverse_position}"
                    f"|{amplicon_length}bp"
                )

                outfile.write(header + "\n")
                outfile.write("".join(amplicon_sequence) + "\n")

                count_amplicon_hits += 1

            # Reset variables for the new target sequence
            current_target = match.group(1).strip()

            forward_primer = None
            forward_position = None

            reverse_primer = None
            reverse_position = None

            amplicon_length = None
            amplicon_sequence = []

            inside_amplicon = False

            continue


        ####################################################
        # Detect forward primer binding site
        ####################################################

        match = re.match(r">(\S+)\s+(\d+)->(\d+)", line)

        if match and forward_primer is None:
            forward_primer = match.group(1)
            forward_position = f"{match.group(2)}-{match.group(3)}"

            continue


        ####################################################
        # Detect reverse primer binding site
        ####################################################

        match = re.match(r">(\S+)\s+(\d+)<-(\d+)", line)

        if match:
            reverse_primer = match.group(1)
            reverse_position = f"{match.group(3)}-{match.group(2)}"

            continue


        ####################################################
        # Detect amplicon header
        ####################################################

        match = re.match(r">(\d+)-(\d+)\s+Amplicon size:\s+(\d+)bp", line)

        if match:
            amplicon_length = match.group(3)
            amplicon_sequence = []
            inside_amplicon = True

            continue


        ####################################################
        # Read amplicon sequence lines
        ####################################################

        if inside_amplicon:

            # End of amplicon block
            if line == "" or line.startswith("___"):

                if amplicon_sequence:
                    header = (
                        f">{current_target}"
                        f"|{forward_primer}:{forward_position}"
                        f"|{reverse_primer}:{reverse_position}"
                        f"|{amplicon_length}bp"
                    )

                    outfile.write(header + "\n")
                    outfile.write("".join(amplicon_sequence) + "\n")

                    count_amplicon_hits += 1

                amplicon_sequence = []
                inside_amplicon = False

            else:
                amplicon_sequence.append(line)


    ########################################################
    # Save the final amplicon if the file ends inside a block
    ########################################################

    if inside_amplicon and amplicon_sequence:
        header = (
            f">{current_target}"
            f"|{forward_primer}:{forward_position}"
            f"|{reverse_primer}:{reverse_position}"
            f"|{amplicon_length}bp"
        )

        outfile.write(header + "\n")
        outfile.write("".join(amplicon_sequence) + "\n")

        count_amplicon_hits += 1


############################################################
# Print summary
############################################################

print(f"count_amplicon_hits: {count_amplicon_hits}")