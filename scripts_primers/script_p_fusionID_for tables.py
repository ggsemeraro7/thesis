import pandas as pd

df = pd.read_excel("1_primer_design_table_info.xlsx", header=2)

def make_final_id(row):
    origin_initial = row["db_origin"][0] if pd.notna(row["db_origin"]) else ""
    base = f"{row['p_ID_new']}{row['p_direction']}"
    return f"{base}_{row['p_dvlp_program']}_{origin_initial}_{row['db_seq']}"

df["p_final_ID_deg"] = df.apply(make_final_id, axis=1)

df.to_excel("primer_design_table_info.xlsx", index=False)

print("Job completed!")