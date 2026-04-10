from functional_annotation.extract_gene_list import extract_id_column
from functional_annotation.functional_annotation import do_functional_annotation
from functional_annotation.find_common_terms import find_common_terms
from functional_annotation.compare_by_genes import compare_by_genes
from ui.inputs_ui import run_gui

def main(infile1, infile2, id_col, numeric_col, only_common_genes = True, n_top_annotations = 10, knowledgebase_directory = "knowledgebase", annotation_outfile1 = "outputs/functional_annotation_results1.tsv", annotation_outfile2 = "outputs/functional_annotation_results_2.tsv", shared_outfile = "outputs/shared_annotation_genes.tsv", final_outfile = "outputs/t_test_results.tsv"):
    """
    in_file1: A filepath to a file containing uniprot IDs and other information.
    in_file2: Another filepath to compare.
    id_col: Name of uniprot ID column in the input files.
    knowledgebase_directory: Path to DAVID knowled gebase directory.
    numeric_col: The numerical column of interest from the deuterater inputs to be used in T-testing.
    only_common_genes: keeps only shared annotation genes and runs a pairwise T-test if true, else does an independent (Welch's) T-test.
    n_top_annotations: The number of most-common annotations between the two input files to keep.
    """

    # extract the uniprot accession IDs as sets from the input files:
    gene_set1 = extract_id_column(infile1, id_col)
    gene_set2 = extract_id_column(infile2, id_col)

    # calculate the functional annotations for each gene set and write them to a file.
    annotations_df1 = do_functional_annotation(gene_set1, knowledgebase_directory, annotation_outfile1)
    annotations_df2 = do_functional_annotation(gene_set2, knowledgebase_directory, annotation_outfile2)

    # find the categories that had the most genes in common between the two annotations, and extract the genes in common.
    common_annotations_df = find_common_terms(annotations_df1, annotations_df2, shared_outfile, only_common_genes=only_common_genes, n_top_annotations=n_top_annotations)
    print(common_annotations_df.head())

    # run a t-test (pairwise if only_common_genes, else Welch's) on the two input files for the list of genes in each of the top n shared annotations.
    t_results_df = compare_by_genes(common_annotations_df, numeric_col, id_col, infile1, infile2, only_common_genes=only_common_genes)
    t_results_df.to_csv(final_outfile, sep="\t")


if __name__ == "__main__":
    inputs = run_gui()

    main(
        # Make sure these inputs match your input files.
        infile1=inputs["infile1"],
        infile2=inputs["infile2"],
        id_col=inputs["id_col"],
        numeric_col=inputs["numeric_col"],

        # Optionally, change these settings as necessary.
        only_common_genes=inputs["only_common_genes"],
        n_top_annotations=inputs["n_top_annotations"],
        knowledgebase_directory=inputs["knowledgebase_directory"],
        annotation_outfile1=inputs["annotation_outfile1"],
        annotation_outfile2=inputs["annotation_outfile2"],
        shared_outfile=inputs["shared_outfile"],
        final_outfile=inputs["final_outfile"]
    )