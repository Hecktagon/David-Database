#!/usr/bin/env python3

import os
import glob
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def load_gene_list(filepath):
    with open(filepath) as f:
        genes = set(line.strip() for line in f if line.strip())
    return genes


def load_mapping_file(filepath):
    """
    Returns:
        term_to_genes: dict {term: set(uniprot_ids)}
        all_genes: set of all uniprot_ids in file
    """
    df = pd.read_csv(filepath, sep="\t", header=None, names=["gene", "term"])
    df = df.dropna()

    term_to_genes = {}
    for term, group in df.groupby("term"):
        term_to_genes[term] = set(group["gene"])

    all_genes = set(df["gene"])

    return term_to_genes, all_genes


def compute_enrichment(term_to_genes, input_genes, background_genes, category_name, min_term_size):
    """
    Performs hypergeometric enrichment for one annotation file.
    """
    results = []

    N = len(background_genes)
    n = len(input_genes)

    for term, genes_in_term in term_to_genes.items():

        genes_in_term_bg = genes_in_term & background_genes
        K = len(genes_in_term_bg)

        if K < min_term_size:
            continue

        genes_in_list = genes_in_term_bg & input_genes
        k = len(genes_in_list)

        if k == 0:
            continue

        # Hypergeometric test (P[X >= k])
        pval = hypergeom.sf(k - 1, N, K, n)

        # Fold enrichment
        fold_enrichment = (k / n) / (K / N)

        results.append({
            "Category-Term": f"{category_name}:{term}",
            "Count": k,
            "Percent": 100 * k / n,
            "PValue": pval,
            "Genes": ",".join(sorted(genes_in_list)),
            "List Total": n,
            "Pop Hits": K,
            "Pop Total": N,
            "Fold Enrichment": fold_enrichment
        })

    return results


def do_functional_annotation(input_genes, knowledgebase_directory, outfile, min_background_term_size = 5):

    all_results = []

    knowledgebase_files = glob.glob(os.path.join(knowledgebase_directory, "*.txt"))

    for filepath in knowledgebase_files:

        category_name = os.path.basename(filepath).replace(".txt", "")

        print(f"\tProcessing {category_name}")

        term_to_genes, background_genes = load_mapping_file(filepath)

        # Restrict input genes to those present in background
        input_genes_filtered = input_genes & background_genes

        results = compute_enrichment(
            term_to_genes,
            input_genes_filtered,
            background_genes,
            category_name,
            min_term_size=min_background_term_size
        )

        all_results.extend(results)

    if not all_results:
        print("\tWARNING: No enrichment found.")
        return

    results_df = pd.DataFrame(all_results)

    # Multiple testing correction
    pvals = results_df["PValue"].values

    bonf = multipletests(pvals, method="bonferroni")[1]
    bh = multipletests(pvals, method="fdr_bh")[1]

    results_df["Bonferroni"] = bonf
    results_df["Benjamini"] = bh

    # Sort by p-value
    results_df = results_df.sort_values("PValue")

    sep = "\t" if outfile.endswith(".tsv") else ","
    results_df.to_csv(outfile, sep=sep, index=False)

    print(f"\n\tResults written to {outfile}")

    return results_df


