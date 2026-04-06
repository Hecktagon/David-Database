import pandas as pd
from pathlib import Path
from functools import reduce

def get_column_name_from_file(filename, prefixes_to_remove = ("UP_KW_", "UP_")):
    last_part = filename.split("2")[-1].split(".")[0]

    for prefix in prefixes_to_remove:
        if last_part.startswith(prefix):
            last_part = last_part[len(prefix):]
            break  # only remove the first matching prefix

    return last_part.lower()


def combine_raw_david_knowledgebase(folder_path, output_file = "databases/combined_knowledgebase.tsv"):
    dfs = []

    for file in Path(folder_path).glob("*.txt"):
        col_name = get_column_name_from_file(file.name)
        df = pd.read_csv(file, sep="\t", header=None, names=["uniprot_accession", col_name])
        dfs.append(df)

    # inner join all files on accession
    combined_df = reduce(lambda left, right: pd.merge(left, right, on="uniprot_accession", how="inner"), dfs)
    if output_file:
        combined_df.to_csv(output_file, sep="\t")
    return combined_df

# combined_knowledgebase = raw_david_to_pandas("knowledgebase")
# print(combined_knowledgebase)

# df = pd.read_csv("databases/combined_knowledgebase.tsv", sep ="\t")