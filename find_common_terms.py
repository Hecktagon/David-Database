import pandas as pd


def top_shared_by_count(df1, df2, index_col, cols_to_keep, n):
    top_dict = {
        index_col: [],
        "Shared-Count": []
    }

    for col in cols_to_keep:
        top_dict[col] = []

    # index both dfs by index_col
    d1 = df1.set_index(index_col)
    d2 = df2.set_index(index_col)

    # find shared indices
    shared_idx = d1.index.intersection(d2.index)

    # compute shared min counts
    shared_counts = {
        idx: min(d1.loc[idx, "Count"], d2.loc[idx, "Count"])
        for idx in shared_idx
    }

    # get top n indices by shared count
    top_items = sorted(shared_counts.items(), key=lambda x: x[1], reverse=True)[:n]

    for idx, count in top_items:
        top_dict[index_col].append(idx)
        top_dict["Shared-Count"].append(count)

        for col in cols_to_keep:
            val1 = d1.loc[idx, col]
            val2 = d2.loc[idx, col]
            top_dict[col].append([val1, val2])

    return pd.DataFrame(top_dict)


def keep_common_genes(df):
    df = df.copy()

    def intersect(pair):
        g1 = pair[0].split(",")
        g2 = pair[1].split(",")
        return [list(set(g1) & set(g2))]

    df["Genes"] = df["Genes"].apply(intersect)
    return df


def split_gene_lists(df):
    df = df.copy()

    def split_pair(pair):
        g1 = pair[0].split(",")
        g2 = pair[1].split(",")
        return [g1, g2]

    df["Genes"] = df["Genes"].apply(split_pair)
    return df


def find_common_terms(df1, df2, outfile, only_common_genes=True, n_top_annotations=10, index_col="Category-Term", cols_to_keep=None):
    if cols_to_keep is None:
        cols_to_keep = set()
    cols_to_keep.add("Genes")

    df_result = top_shared_by_count(df1, df2, index_col=index_col, cols_to_keep=cols_to_keep, n = n_top_annotations)

    if only_common_genes:
        df_result = keep_common_genes(df_result)
    else:
        df_result = split_gene_lists(df_result)

    if outfile:
        sep = "\t" if outfile.endswith(".tsv") else ","
        df_result.to_csv(outfile, sep=sep, index=False)
        print(f"Common term results written to {outfile}")

    return df_result