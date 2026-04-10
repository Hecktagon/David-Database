import pandas as pd
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind


def file_to_numeric_df(file, numeric_col, genes, gene_col="analyte_id"):
    # Takes a deuterater gene file and turns it into a df with only the specified genes, and only the column of interest.
    sep = "\t" if file.endswith(".tsv") else ","
    df = pd.read_csv(file, sep=sep)
    filtered_df = df[df[gene_col].isin(genes)]
    numeric_df = filtered_df[[gene_col, numeric_col]]
    numeric_df = numeric_df.rename(columns={gene_col: 'Genes'}).set_index("Genes")
    return numeric_df


def paired_t_test_dfs_by_col(df1, df2, col):
    # runs a T-test between two dfs for a column of interest.

    combined_df = df1.join(df2, on="Genes", how="inner", lsuffix="_1", rsuffix="_2")

    combined_df[col + '_1'] = pd.to_numeric(combined_df[col + '_1'], errors='coerce')
    combined_df[col + '_2'] = pd.to_numeric(combined_df[col + '_2'], errors='coerce')
    combined_df = combined_df.dropna(subset=[col + "_1", col + "_2"])

    t_stat, p_val = ttest_rel(combined_df[col + "_1"], combined_df[col + "_2"])

    return t_stat, p_val, len(combined_df)


def compare_by_genes_paired(file1, file2, genes, numeric_col, gene_col):
    df1 = file_to_numeric_df(file1, numeric_col, genes, gene_col=gene_col)
    df2 = file_to_numeric_df(file2, numeric_col, genes, gene_col=gene_col)

    return paired_t_test_dfs_by_col(df1, df2, col=numeric_col)


def compare_by_genes_independent(file1, file2, genes1, genes2, numeric_col, gene_col):
    """
    Runs an independent t-test (Welch's t-test) between two files for the numeric_col.
    genes1 and genes2 are lists of genes for file1 and file2 respectively.
    Works even if the gene lists are different and lengths don't match.
    """
    # Load numeric data for each file
    df1 = file_to_numeric_df(file1, numeric_col, genes1, gene_col=gene_col)
    df2 = file_to_numeric_df(file2, numeric_col, genes2, gene_col=gene_col)

    # Convert to numeric and drop NA
    col1 = pd.to_numeric(df1[numeric_col], errors='coerce')
    col2 = pd.to_numeric(df2[numeric_col], errors='coerce')
    col1 = col1.dropna()
    col2 = col2.dropna()

    # Run independent t-test (Welch's t-test)
    t_stat, p_val = ttest_ind(col1, col2, equal_var=False)

    # returns t-test results and the length of the smaller gene list as the count.
    return t_stat, p_val, min([len(col1), len(col2)])


def compare_by_genes(common_annotations_df, numeric_col, gene_col, file1, file2, only_common_genes = True):
    t_results = {
        "Category-Term": [],
        "Count-After-DropNA": [],
        "T_stat": [],
        "P_value": []
    }
    for i, row in common_annotations_df.iterrows():
        if only_common_genes:
            genes = row["Genes"][0]
            t_stat, p_val, count = compare_by_genes_paired(file1, file2, genes, numeric_col, gene_col)
        else:
            genes1, genes2 = row["Genes"]
            t_stat, p_val, count = compare_by_genes_independent(file1, file2, genes1, genes2, numeric_col, gene_col)

        t_results["Category-Term"].append(row["Category-Term"])
        t_results["Count-After-DropNA"].append(count)
        t_results["T_stat"].append(t_stat)
        t_results["P_value"].append(p_val)

    return pd.DataFrame(t_results)

