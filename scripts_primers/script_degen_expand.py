##degeneration explicited 
#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/script_degen_expand.py" linb_our_raw.txt


import sys
from itertools import product

# dizionario di degenerazione IUPAC
degenerate = {
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

def expand_iupac(seq):
    """
    Genera tutte le possibili espansioni per una sequenza con
    basi IUPAC degeneranti.

    Ritorna tuple (tag, expanded_seq) dove
    tag = solo le basi esplicitate nei punti degenerati
    """
    # costruisce le liste di possibili basi per ogni posizione
    pools = [degenerate.get(base.upper(), [base]) for base in seq]
    # genera tutte le combinazioni
    for combo in product(*pools):
        expanded_seq = "".join(combo)
        # il tag è concatenazione solo delle basi esplicitate
        # nei punti degenerati originali
        tag_chars = [
            combo[i]
            for i, base in enumerate(seq)
            if base.upper() in degenerate and len(degenerate[base.upper()]) > 1
        ]
        tag = "".join(tag_chars)
        yield tag, expanded_seq

def parse_fasta(fasta_path):
    """
    Legge un file FASTA e restituisce una lista di
    (header, sequence) tuple.
    """
    records = []
    with open(fasta_path) as f:
        header = None
        seq_lines = []

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records

def main():
    if len(sys.argv) != 2:
        print("Uso: python3 expand_iupac.py <file_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = "linb_expanded_primers_our.fasta"

    records = parse_fasta(fasta_file)

    with open(output_file, "w") as out:
        for header, seq in records:
            for tag, expanded_seq in expand_iupac(seq):
                # se non ci sono basi degenerati, tag può essere vuoto
                # allora metti "_" + tag solo se tag non è vuoto
                if tag:
                    new_header = f"{header}_{tag}"
                else:
                    new_header = header
                out.write(f">{new_header}\n")
                out.write(f"{expanded_seq}\n")

    print(f"Completed! File saved as: {output_file}")

if __name__ == "__main__":
    main()
