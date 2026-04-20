##script_final tables
#nomi_ID_primers 
#python3 "/mnt/c/Users/gseme/OneDrive - University of Pisa/UniPi/Tesi/scripts_primers/script_p_fusionID.py"

import pandas as pd

# legge il file intero senza header
df_full = pd.read_excel("1_primer_design_table_info.xlsx", header=None)

# salva le prime due righe così come sono (metadati)
header_rows = df_full.iloc[:2].copy()

# i dati veri partono da riga 3
data = df_full.iloc[2:].copy()

# la riga 3 diventa header dei dati
data.columns = data.iloc[0]
data = data[1:].reset_index(drop=True)

# pulisce gli spazi nei nomi delle colonne
data.columns = [col.strip() for col in data.columns]

# funzione per creare p_final_ID
def make_final_id(row):
    origin = str(row["db_origin"])
    origin_initial = origin[0] if origin else ""
    base = f"{row['p_ID_new']}{row['p_direction']}"
    return f"{base}_{row['p_dvlp_program']}_{origin_initial}_{row['db_seq']}"

# aggiunge la nuova colonna
data["p_final_ID"] = data.apply(make_final_id, axis=1)

# ricombina prime due righe + dati con colonna aggiunta
result = pd.concat([header_rows, data], ignore_index=True)

# salva *senza header e senza index* così come è nel file originale
result.to_excel("primer_design_table_info.xlsx", index=False, header=False)

print("Job completed!")
