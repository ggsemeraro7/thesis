#!/usr/bin/env python3

############################################################
# BUILD PRIMER SPECIFICITY TABLE
############################################################

# This script combines:
# 1. an amplicon summary table
# 2. an ORF annotation table
#
# It checks whether each amplicon overlaps annotated ORFs
# and adds the corresponding functional annotation fields
# to the amplicon table.
#
# Expected input:
# - Amplicon table generated from extracted FastPCR amplicons
# - ORF table containing contig IDs, ORF IDs, and annotation columns
#
# Expected output:
# - Amplicon specificity table with ORF and functional annotations


import pandas as pd


############################################################
# Input and output files
############################################################

amplicon_table = "<AMPLICON_SUMMARY_TSV>"
orf_table = "<ORF_ANNOTATION_TABLE>"
output_table = "<SPECIFICITY_OUTPUT_TSV>"


############################################################
# Helper functions
############################################################

def ranges_overlap(start1, end1, start2, end2):
    """
    Return True if two coordinate ranges overlap at least partially.
    """

    return (start1 <= end2) and (start2 <= end1)


def parse_orf_id(orf_id):
    """
    Parse an ORF identifier containing genomic coordinates.

    Expected ORF ID format:
    <contig_id>_<start>-<end>

    Example:
    sample_contig_001_25000-25878

    Returns:
    - contig identifier
    - ORF start coordinate
    - ORF end coordinate
    """

    parts = str(orf_id).rsplit("_", 1)

    if len(parts) == 2 and "-" in parts[1]:
        coord_start, coord_end = parts[1].split("-")

        try:
            return parts[0], int(coord_start), int(coord_end)
        except ValueError:
            return parts[0], None, None

    return str(orf_id), None, None


def aggregate_matches(field, orf_matches):
    """
    Safely aggregate values from matched ORFs using semicolon separation.
    Empty values and missing values are ignored.
    """

    values = []

    for match in orf_matches:
        value = match.get(field)

        if value is None:
            continue

        value_string = str(value).strip()

        if value_string == "" or value_string.lower() == "nan":
            continue

        values.append(value_string)

    return ";".join(values)


############################################################
# Read input tables
############################################################

amplicons = pd.read_csv(
    amplicon_table,
    sep="\t",
    dtype=str
)

orfs = pd.read_csv(
    orf_table,
    sep="\t",
    dtype=str,
    comment="#"
)


############################################################
# Clean column names
############################################################

amplicons.columns = amplicons.columns.str.strip()
orfs.columns = orfs.columns.str.strip()


############################################################
# Convert amplicon coordinates to numeric values
############################################################

amplicons["start"] = pd.to_numeric(amplicons["start"], errors="coerce")
amplicons["end"] = pd.to_numeric(amplicons["end"], errors="coerce")


############################################################
# Add output annotation columns
############################################################

columns_to_add = [
    "ORF_ID_found",
    "Gene_name",
    "Tax",
    "KEGG_ID",
    "KEGGFUN",
    "KEGGPATH",
    "COG_ID",
    "COGFUN",
    "COGPATH",
    "PFAM",
    "Hits"
]

for column in columns_to_add:
    amplicons[column] = ""


############################################################
# Match amplicons to overlapping ORFs
############################################################

for index, row in amplicons.iterrows():

    input_file_id = str(row.get("ID_file", "")).lower()

    # Process only rows generated from specificity analysis files.
    # Modify or remove this filter if all rows should be processed.
    if "magspec" not in input_file_id:
        continue

    sequence_id = row["seqID"]
    amplicon_start = row["start"]
    amplicon_end = row["end"]

    if pd.isna(amplicon_start) or pd.isna(amplicon_end):
        continue

    # Select ORFs from the same contig or sequence.
    matching_orfs = orfs[orfs["Contig ID"] == sequence_id]

    if matching_orfs.empty:
        continue

    orf_matches = []

    for _, orf_row in matching_orfs.iterrows():

        orf_full_id = orf_row["ORF ID"]

        _, orf_start, orf_end = parse_orf_id(orf_full_id)

        if orf_start is None or orf_end is None:
            continue

        if ranges_overlap(amplicon_start, amplicon_end, orf_start, orf_end):
            orf_matches.append({
                "ORF_ID": orf_full_id,
                "Gene_name": orf_row.get("Gene name", ""),
                "Tax": orf_row.get("Tax", ""),
                "KEGG_ID": orf_row.get("KEGG ID", ""),
                "KEGGFUN": orf_row.get("KEGGFUN", ""),
                "KEGGPATH": orf_row.get("KEGGPATH", ""),
                "COG_ID": orf_row.get("COG ID", ""),
                "COGFUN": orf_row.get("COGFUN", ""),
                "COGPATH": orf_row.get("COGPATH", ""),
                "PFAM": orf_row.get("PFAM", ""),
                "Hits": orf_row.get("Hits", "")
            })

    if orf_matches:
        amplicons.at[index, "ORF_ID_found"] = aggregate_matches("ORF_ID", orf_matches)
        amplicons.at[index, "Gene_name"] = aggregate_matches("Gene_name", orf_matches)
        amplicons.at[index, "Tax"] = aggregate_matches("Tax", orf_matches)
        amplicons.at[index, "KEGG_ID"] = aggregate_matches("KEGG_ID", orf_matches)
        amplicons.at[index, "KEGGFUN"] = aggregate_matches("KEGGFUN", orf_matches)
        amplicons.at[index, "KEGGPATH"] = aggregate_matches("KEGGPATH", orf_matches)
        amplicons.at[index, "COG_ID"] = aggregate_matches("COG_ID", orf_matches)
        amplicons.at[index, "COGFUN"] = aggregate_matches("COGFUN", orf_matches)
        amplicons.at[index, "COGPATH"] = aggregate_matches("COGPATH", orf_matches)
        amplicons.at[index, "PFAM"] = aggregate_matches("PFAM", orf_matches)
        amplicons.at[index, "Hits"] = aggregate_matches("Hits", orf_matches)


############################################################
# Save output table
############################################################

amplicons.to_csv(
    output_table,
    sep="\t",
    index=False
)

print(f"Done. Specificity table saved to: {output_table}")