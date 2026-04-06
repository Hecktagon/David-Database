import pandas as pd
import numpy as np

def diff_csv(input1, input2, output, delimiter=","):
    sep1 = "\t" if input1.endswith(".tsv") else ","
    df1 = pd.read_csv(input1, sep=sep1, header=0)
    sep2 = "\t" if input2.endswith(".tsv") else ","
    df2 = pd.read_csv(input2, sep=sep2, header=0)
    print(df1.columns)
    print(df2.head())

    # create combined key
    df1["Category_term"] = df1["Category"] + ":" + df1["Term"]
    df2["Category_term"] = df2["Category"] + ":" + df2["Term"]

    # set as index
    df1 = df1.set_index("Category_term")
    df2 = df2.set_index("Category_term")

    # align rows by index
    df1, df2 = df1.align(df2, join="inner")

    # numeric columns only
    numeric_cols = df1.select_dtypes(include=np.number).columns

    diff_df = df1.copy()
    diff_df[numeric_cols] = (df1[numeric_cols] - df2[numeric_cols]).abs()

    diff_df.to_csv(output, sep=delimiter)

    print(f"Diffs written to {output}")


if __name__ == "__main__":
    file1 = "functional_annotation_results_old.tsv"
    file2 = "functional_annotation_results_old2.tsv"
    outfile = "func_annotation_diffs.tsv"

    diff_csv(file1, file2, outfile, delimiter="\t")