import pandas as pd
from pathlib import Path


def get_enrichments(gene_list, file_or_folder="knowledgebase"):
    """
    for file in knowledgebase:
        which attributes from file are most common in my gene list?
    """
    path = Path(file_or_folder)
    enrichments = []

    # if the path is a file, load it
    if path.is_file() and path.suffix == ".txt":
        col_name = get_column_name_from_file(path.name)
        df = pd.read_csv(path, sep="\t", header=None, names=["uniprot_accession", col_name])
        enrichments.append(enrichment_for_df(gene_list, df))
    # if the path is a directory, load each file within it.
    elif path.is_dir():
        for file in path.glob("*.txt"):
            col_name = get_column_name_from_file(file.name)
            df = pd.read_csv(file, sep="\t", header=None, names=["uniprot_accession", col_name])
            enrichments.append(enrichment_for_df(gene_list, df))

    return enrichments


def get_column_name_from_file(filename, prefixes_to_remove=("UP_KW_", "UP_")):
    last_part = filename.split("2")[-1].split(".")[0]

    for prefix in prefixes_to_remove:
        if last_part.startswith(prefix):
            last_part = last_part[len(prefix):]
            break  # only remove the first matching prefix

    return last_part.lower()


def enrichment_for_df(gene_list, df):
    annotation_col = df.columns[1]

    # df with only rows that contain genes from gene_list
    gene_annotations_df = df[df["uniprot_accession"].isin(gene_list)]

    enrichment_df = (
        gene_annotations_df.groupby(annotation_col)
        .agg(
            uniprot_accessions=("uniprot_accession", list),
            count=("uniprot_accession", "count")
        )
        .reset_index()
    )

    enrichment_df = enrichment_df.sort_values(by="count", ascending=False)
    return enrichment_df


if __name__ == "__main__":
    gene_accessions = [
        "Q80TN5",
        "B3KPC1",
        "E9PLF1",
        "F5H056",
        "D3Z7Q9",
        "P21673",
        "Q08AN1",
        "J3QSG9",
        "Q923F5",
        "C0LUL2",
        "J3KQR4",
        "fake_id"
    ]
    enrichment_dfs = get_enrichments(gene_accessions, "knowledgebase")
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_colwidth", None)
    pd.set_option("display.width", None)
    for enrichment in enrichment_dfs:
        print()
        print(enrichment.head(10))
