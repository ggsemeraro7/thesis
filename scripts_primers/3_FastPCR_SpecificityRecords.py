#Attiva ambiente venv: 
#source "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/PRIMERS/tesiR/MAGspecificity/venv/bin/activate"

#PIPELINE PER CREARE TABELLA SPECIFICITà (da tabella ampliconi e da tabella ortable) 
#Usa tramite:
#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/3_FastPCR_SpecificityRecords.py"


#Tailored script

import pandas as pd

# --- FUNZIONE DI SOVRAPPOSIZIONE DI RANGE ---
def ranges_overlap(start1, end1, start2, end2):
    """
    Ritorna True se i due range [start1,end1] e [start2,end2] si sovrappongono anche parzialmente.
    """
    return (start1 <= end2) and (start2 <= end1)

# --- FUNZIONE PER PARSARE ORF ID ---
def parse_orf_id(orf_id):
    """
    Esegue il parsing di ORF ID tipo
    COL_SA_mg_mt_15701_25000-25878

    Ritorna:
    - contig_part = 'COL_SA_mg_mt_15701'
    - start2, end2 (coordinate)
    """
    parts = orf_id.rsplit("_", 1)
    if len(parts) == 2 and "-" in parts[1]:
        coord_start, coord_end = parts[1].split("-")
        try:
            return parts[0], int(coord_start), int(coord_end)
        except:
            return parts[0], None, None
    return parts[0], None, None

# --- FUNZIONE DI AGGREGAZIONE SICURA ---
def agg_matches(field, ORF_matches):
    vals = []
    for d in ORF_matches:
        v = d.get(field)
        if v is None:
            continue
        v_str = str(v).strip()
        if v_str == "" or v_str.lower() == "nan":
            continue
        vals.append(v_str)
    return ";".join(vals)

# --- LETTURA DELLE TABELLE ---  
tab1 = pd.read_csv("linb12_amplicons.tsv", sep="\t", dtype=str)

tab2 = pd.read_csv(
    "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/PRIMERS/tesiR/MAGspecificity/13.COL_SA_mg_mt.orftable",
    sep="\t",
    dtype=str,
    comment="#"
)

# Pulizia header tab2 (rimuove spazi indesiderati)
tab2.columns = tab2.columns.str.strip()

# convertire start/end numerici (Tabella1)
tab1["start"] = pd.to_numeric(tab1["start"], errors="coerce")
tab1["end"] = pd.to_numeric(tab1["end"], errors="coerce")

# colonne da aggiungere
cols_to_add = [
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

for c in cols_to_add:
    tab1[c] = ""

# --- LOOP SU TABELLA1 ---
for idx, row in tab1.iterrows():

    # match case-insensitive su ID_file per 'MAGspecificity'
    idfile_val = str(row.get("ID_file", "")).lower()
    if "magspec" not in idfile_val:
        continue

    seqid = row["seqID"]
    s1 = row["start"]
    e1 = row["end"]

    # subfilter Tabella2 per stesso contig
    sub2 = tab2[tab2["Contig ID"] == seqid]
    if sub2.empty:
        continue

    ORF_matches = []

    for i2, r2 in sub2.iterrows():
        orf_full = r2["ORF ID"]

        contig2, s2, e2 = parse_orf_id(orf_full)
        if s2 is None or e2 is None:
            continue

        # controllo di sovrapposizione range
        if ranges_overlap(s1, e1, s2, e2):
            ORF_matches.append({
                "ORF_ID": orf_full,
                "Gene_name": r2.get("Gene name", ""),
                "Tax": r2.get("Tax", ""),
                "KEGG_ID": r2.get("KEGG ID", ""),
                "KEGGFUN": r2.get("KEGGFUN", ""),
                "KEGGPATH": r2.get("KEGGPATH", ""),
                "COG_ID": r2.get("COG ID", ""),
                "COGFUN": r2.get("COGFUN", ""),
                "COGPATH": r2.get("COGPATH", ""),
                "PFAM": r2.get("PFAM", ""),
                "Hits": r2.get("Hits", "")
            })

    # se ci sono match, aggrega con ';'
    if ORF_matches:
        tab1.at[idx, "ORF_ID_found"] = agg_matches("ORF_ID", ORF_matches)
        tab1.at[idx, "Gene_name"]    = agg_matches("Gene_name", ORF_matches)
        tab1.at[idx, "Tax"]          = agg_matches("Tax", ORF_matches)
        tab1.at[idx, "KEGG_ID"]      = agg_matches("KEGG_ID", ORF_matches)
        tab1.at[idx, "KEGGFUN"]      = agg_matches("KEGGFUN", ORF_matches)
        tab1.at[idx, "KEGGPATH"]     = agg_matches("KEGGPATH", ORF_matches)
        tab1.at[idx, "COG_ID"]       = agg_matches("COG_ID", ORF_matches)
        tab1.at[idx, "COGFUN"]       = agg_matches("COGFUN", ORF_matches)
        tab1.at[idx, "COGPATH"]      = agg_matches("COGPATH", ORF_matches)
        tab1.at[idx, "PFAM"]         = agg_matches("PFAM", ORF_matches)
        tab1.at[idx, "Hits"]         = agg_matches("Hits", ORF_matches)

# --- SALVATAGGIO ---
tab1.to_csv("linb12_amplicons_specs.tsv", sep="\t", index=False)

print("Done. File saved")

