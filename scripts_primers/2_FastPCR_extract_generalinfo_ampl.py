#Attiva ambiente venv: 
#source "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/PRIMERS/tesiR/MAGspecificity/venv/bin/activate"

#Usa tramite
#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/2_FastPCR_extract_generalinfo_ampl.py"

#Tailored script
import csv
import os
from datetime import date

fasta_file = "linb12_mis2_MAGspec_amplicons.fasta"
output_file = "linb12_amplicons.tsv"

# ID_file = nome file senza estensione .txt
ID_file = os.path.basename(fasta_file).replace(".txt", "")

file_exists = os.path.exists(output_file)

with open(fasta_file) as f:

    records = {}
    header = None
    seq_lines = []

    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if header and seq_lines:
                fasta_seq = "".join(seq_lines)

                parts = header.split("|")

                fwd_field = parts[-3].strip()
                rev_field = parts[-2].strip()
                length_part = parts[-1].strip()
                seqID = "|".join(parts[:-3]).strip()

                fwd_name, fwd_pos = fwd_field.split(":")
                fwd_start, fwd_end = map(int, fwd_pos.split("-"))

                rev_name, rev_pos = rev_field.split(":")
                r1, r2 = map(int, rev_pos.split("-"))
                rev_start = min(r1, r2)
                rev_end = max(r1, r2)

                start = fwd_start
                end = rev_end
                length_bp = length_part.replace("bp", "").strip()

                if seqID not in records:
                    records[seqID] = [
                        seqID, fwd_name, rev_name,
                        start, end,
                        length_bp, fasta_seq,
                        1  # notes = numero hit iniziali
                    ]
                else:
                    records[seqID][-1] += 1

            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)

    # ultima entry
    if header and seq_lines:
        fasta_seq = "".join(seq_lines)
        parts = header.split("|")

        fwd_field = parts[-3].strip()
        rev_field = parts[-2].strip()
        length_part = parts[-1].strip()
        seqID = "|".join(parts[:-3]).strip()

        fwd_name, fwd_pos = fwd_field.split(":")
        fwd_start, fwd_end = map(int, fwd_pos.split("-"))

        rev_name, rev_pos = rev_field.split(":")
        r1, r2 = map(int, rev_pos.split("-"))
        rev_start = min(r1, r2)
        rev_end = max(r1, r2)

        start = fwd_start
        end = rev_end
        length_bp = length_part.replace("bp", "").strip()

        if seqID not in records:
            records[seqID] = [
                seqID, fwd_name, rev_name,
                start, end,
                length_bp, fasta_seq,
                1
            ]
        else:
            records[seqID][-1] += 1


# 🔢 contatore righe aggiunte
righe_aggiunte = 0

with open(output_file, "a", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    if not file_exists:
        writer.writerow([
            "seqID", "fwd", "rev",
            "start", "end",
            "length(bp)", "amplicon_seq",
            "notes", "ID_file"
        ])

    for rec in records.values():
        writer.writerow(rec + [ID_file])
        righe_aggiunte += 1


print(f"Done! Analysis added. Righe aggiunte: {righe_aggiunte}")

