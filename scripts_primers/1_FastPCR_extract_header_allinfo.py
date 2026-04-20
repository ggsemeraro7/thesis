#Attiva ambiente venv: 
#source "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/PRIMERS/tesiR/MAGspecificity/venv/bin/activate"

#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/1_FastPCR_extract_header_allinfo.py"

import re

input_file = "linb12_mis2_MAGspec.txt"
output_fasta = "linb12_mis2_MAGspec_amplicons.fasta"

# contatore ampliconi trovati
count_hits = 0

with open(input_file) as f:
    lines = f.readlines()

current_seq_full = None
primer_f = None
pos_f = None
primer_r = None
pos_r = None
amplicon_len = None
amplicon_seq = []
in_amplicon = False

with open(output_fasta, "w") as out:

    for line in lines:
        line = line.rstrip()

        # --- Nuova sezione target ---
        m = re.match(r"In silico Primer\(s\) search for:\s+(.+)", line)
        if m:
            # salva amplicone precedente (se esiste)
            if in_amplicon and amplicon_seq:
                header = (
                    f">{current_seq_full}|{primer_f}:{pos_f}|{primer_r}:{pos_r}|{amplicon_len}bp"
                )
                out.write(header + "\n")
                out.write("".join(amplicon_seq) + "\n")
                count_hits += 1  # ++ contatore

            # reset
            current_seq_full = m.group(1).strip()
            primer_f = None
            pos_f = None
            primer_r = None
            pos_r = None
            amplicon_len = None
            amplicon_seq = []
            in_amplicon = False
            continue

        # --- Binding site forward ---
        m = re.match(r">(\S+)\s+(\d+)->(\d+)", line)
        if m and primer_f is None:
            primer_f = m.group(1)
            pos_f = f"{m.group(2)}-{m.group(3)}"
            continue

        # --- Binding site reverse ---
        m = re.match(r">(\S+)\s+(\d+)<-(\d+)", line)
        if m:
            primer_r = m.group(1)
            pos_r = f"{m.group(3)}-{m.group(2)}"
            continue

        # --- Amplicon header ---
        m = re.match(r">(\d+)-(\d+)\s+Amplicon size:\s+(\d+)bp", line)
        if m:
            amplicon_len = m.group(3)
            amplicon_seq = []
            in_amplicon = True
            continue

        # --- Amplicon sequence lines ---
        if in_amplicon:
            if line == "" or line.startswith("___"):
                if amplicon_seq:
                    header = (
                        f">{current_seq_full}|{primer_f}:{pos_f}|{primer_r}:{pos_r}|{amplicon_len}bp"
                    )
                    out.write(header + "\n")
                    out.write("".join(amplicon_seq) + "\n")
                    count_hits += 1  # ++ contatore

                amplicon_seq = []
                in_amplicon = False
            else:
                amplicon_seq.append(line)

    # --- ultimo amplicone se fine file ---
    if in_amplicon and amplicon_seq:
        header = (
            f">{current_seq_full}|{primer_f}:{pos_f}|{primer_r}:{pos_r}|{amplicon_len}bp"
        )
        out.write(header + "\n")
        out.write("".join(amplicon_seq) + "\n")
        count_hits += 1  # ++

# alla fine stampa il numero di hits
print(f"count_amplicon_hits: {count_hits}")
