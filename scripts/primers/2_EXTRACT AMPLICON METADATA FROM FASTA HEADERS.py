#!/usr/bin/env python3

############################################################
# EXTRACT AMPLICON METADATA FROM FASTA HEADERS
############################################################

# This script parses a FASTA file containing amplicons extracted from
# FastPCR output and generates a tab-separated summary table.
#
# Expected FASTA header format:
#
# >sequence_id|forward_primer:forward_start-forward_end|reverse_primer:reverse_start-reverse_end|amplicon_lengthbp
#
# Output table columns:
# - seqID
# - forward primer name
# - reverse primer name
# - amplicon start
# - amplicon end
# - amplicon length
# - amplicon sequence
# - number of hits per sequence ID
# - input file identifier


import csv
import os


############################################################
# Input and output files
############################################################

fasta_file = "<AMPLICON_FASTA>"
output_file = "<AMPLICON_SUMMARY_TSV>"


############################################################
# Generate input file identifier
############################################################

# Use the input FASTA filename without extension as file identifier.

input_file_id = os.path.splitext(os.path.basename(fasta_file))[0]


############################################################
# Check whether output file already exists
############################################################

file_exists = os.path.exists(output_file)


############################################################
# Parse FASTA file
############################################################

records = {}
header = None
sequence_lines = []


def process_record(header, sequence_lines, records):
    """
    Parse one FASTA record and update the records dictionary.
    """

    if not header or not sequence_lines:
        return

    amplicon_sequence = "".join(sequence_lines)

    header_parts = header.split("|")

    forward_field = header_parts[-3].strip()
    reverse_field = header_parts[-2].strip()
    length_field = header_parts[-1].strip()

    sequence_id = "|".join(header_parts[:-3]).strip()

    forward_name, forward_position = forward_field.split(":")
    forward_start, forward_end = map(int, forward_position.split("-"))

    reverse_name, reverse_position = reverse_field.split(":")
    reverse_pos_1, reverse_pos_2 = map(int, reverse_position.split("-"))

    reverse_start = min(reverse_pos_1, reverse_pos_2)
    reverse_end = max(reverse_pos_1, reverse_pos_2)

    amplicon_start = forward_start
    amplicon_end = reverse_end

    amplicon_length_bp = length_field.replace("bp", "").strip()

    if sequence_id not in records:
        records[sequence_id] = [
            sequence_id,
            forward_name,
            reverse_name,
            amplicon_start,
            amplicon_end,
            amplicon_length_bp,
            amplicon_sequence,
            1
        ]
    else:
        records[sequence_id][-1] += 1


with open(fasta_file, "r") as infile:

    for line in infile:
        line = line.rstrip()

        if line.startswith(">"):

            process_record(header, sequence_lines, records)

            header = line[1:]
            sequence_lines = []

        else:
            sequence_lines.append(line)


    ########################################################
    # Process final FASTA record
    ########################################################

    process_record(header, sequence_lines, records)


############################################################
# Write or append summary table
############################################################

rows_added = 0

with open(output_file, "a", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")

    if not file_exists:
        writer.writerow([
            "seqID",
            "fwd",
            "rev",
            "start",
            "end",
            "length_bp",
            "amplicon_seq",
            "notes",
            "ID_file"
        ])

    for record in records.values():
        writer.writerow(record + [input_file_id])
        rows_added += 1


############################################################
# Print summary
############################################################

print(f"Done. Analysis added. Rows added: {rows_added}")