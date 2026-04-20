##FASTA TO CSV
#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/fasta_to_csv.py" linb_expanded_primers_our.fasta linb_expanded_primers_our.csv


import csv
import sys

def fasta_to_csv(fasta_file, csv_file):
    with open(fasta_file, "r") as fa, open(csv_file, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["primer_name", "sequence"])

        header = None
        seq = ""

        for line in fa:
            line = line.strip()
            if line.startswith(">"):
                # scrivi la riga precedente se esistente
                if header is not None:
                    writer.writerow([header, seq])
                # imposta nuovo header senza '>'
                header = line[1:]
                seq = ""
            else:
                seq += line

        if header is not None:
            writer.writerow([header, seq])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 fasta_to_csv.py <input_fasta> <output_csv>")
    else:
        fasta_to_csv(sys.argv[1], sys.argv[2])
        print(f"CSV done: {sys.argv[2]}")
