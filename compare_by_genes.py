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
    # df1, df2 = df1.align(df2, join="inner")

    combined_df = df1.join(df2, on="Genes", how="inner", lsuffix="_1", rsuffix="_2")
    print(combined_df.head())

    combined_df[col + '_1'] = pd.to_numeric(combined_df[col + '_1'], errors='coerce')
    combined_df[col + '_2'] = pd.to_numeric(combined_df[col + '_2'], errors='coerce')
    combined_df = combined_df.dropna(subset=[col + "_1", col + "_2"])

    t_stat, p_val = ttest_rel(combined_df[col + "_1"], combined_df[col + "_2"])

    print("Pairwise t-test results:")
    print("t-statistic:", t_stat)
    print("p-value:", p_val)
    return t_stat, p_val


def compare_by_genes_paired(file1, file2, genes, numeric_col):
    df1 = file_to_numeric_df(file1, numeric_col, genes)
    df2 = file_to_numeric_df(file2, numeric_col, genes)

    return paired_t_test_dfs_by_col(df1, df2, col=numeric_col)


def compare_by_genes_independent(file1, file2, genes1, genes2, numeric_col, gene_col="analyte_id"):
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

    print("Independent t-test (Welch's) results:")
    print("t-statistic:", t_stat)
    print("p-value:", p_val)
    return t_stat, p_val


def compare_by_genes(column_of_interest, file1, file2, genes1, genes2=None):
    if not genes2:
        return compare_by_genes_paired(file1, file2, genes1, column_of_interest)
    else:
        return compare_by_genes_independent(infile1, infile2, genes1, genes2, col_of_interest)


if __name__ == "__main__":
    infile1 = "sample_rates_from_deuterater1.csv"
    infile2 = "sample_rates_from_deuterater2.csv"
    geneset1 = {'P17182', 'P41216', 'P99029', 'A2AQP0', 'P20029', 'P45376', 'P97457', 'P61982', 'P14824', 'P47857',
             'Q99LC3', 'P28474', 'Q8CGP0', 'Q8VHX6', 'P62259', 'Q8CI51', 'Q03265', 'Q8BMS1', 'Q9ET01', 'Q3V1D3',
             'P09542', 'P56480', 'Q9JI91', 'Q9R062', 'Q9QZ47', 'P68369', 'P06151', 'O70209', 'Q9DCS9', 'P05063',
             'Q9Z2I9', 'Q9CZU6', 'P05201', 'Q9JKS4', 'Q9D6J6', 'P11499', 'Q8VEM8', 'Q5SX39', 'P99024', 'P05977',
             'P08249', 'Q8BFR5', 'O70250', 'P16858', 'P21550', 'Q9D051', 'Q07417', 'Q9WUZ5', 'Q04447', 'Q9D0K2',
             'P02088', 'Q8QZT1', 'P16045', 'C0HKE3', 'P32848', 'Q9JKB3', 'P17751', 'P17879', 'P04104', 'P26443',
             'P58771', 'Q8CI94', 'O88492', 'P05202', 'Q02053', 'Q8CC35', 'P50544', 'P17156', 'O35129', 'Q9QXX4',
             'Q9QWL7', 'Q9JJW5', 'P51885', 'P15626', 'P45952', 'P50171', 'Q9CRB9', 'P62897', 'P07356', 'Q6PIE5',
             'P35486', 'Q7TQ48', 'Q99LX0', 'Q922F4', 'Q9D2G2', 'Q5EBG6', 'Q9WUB3', 'A2AAJ9', 'P13412', 'P11798',
             'P07724', 'P15532', 'P14152', 'O88456', 'O55143', 'P51174', 'O09165', 'Q06185', 'Q921I1', 'P09103',
             'P18760', 'P68368', 'Q921G7', 'Q8K2B3', 'Q5SX40', 'P62869', 'Q9D6Y9', 'P58252', 'P68510', 'O54724',
             'O08539', 'P68372', 'P70404', 'Q60854', 'Q61171', 'P11404', 'P24527', 'P58774', 'Q60932', 'E9Q557',
             'Q8R429', 'P42232', 'P01942', 'P18572', 'E9PZQ0', 'A2AUC9', 'Q8BTM8', 'P23927', 'Q91VR2', 'P70670',
             'Q9DC70', 'O88346', 'P08228', 'Q64522', 'Q9JK37', 'P16015', 'Q91Z83', 'Q9D0F9', 'Q9DCW4', 'P19123',
             'Q8CHT0', 'P62631', 'P63038', 'P04247', 'Q9R0Y5', 'P14602', 'Q9Z2I0', 'P57776', 'Q60930', 'P56375',
             'P35700', 'O08553', 'Q99MN9', 'P10649', 'Q8VDN2', 'Q9QXS1', 'Q6P8J7', 'Q9CWF2', 'O70622', 'Q9Z2U0',
             'P07901', 'Q9D6R2', 'P62962', 'P60843', 'P15864', 'Q922U2', 'P13707', 'P70349', 'P67778', 'P70296',
             'Q8BWT1', 'Q91ZJ5', 'A2AMM0', 'P0CG50', 'P06745', 'Q99LC5', 'Q9CZ13', 'Q62234', 'Q8CAQ8', 'Q63836',
             'P09411', 'A2ASS6', 'P21107', 'Q9D6F9', 'Q9EQ20', 'P97807', 'Q60597', 'Q9ESD7', 'Q3TXS7', 'P17742',
             'P09041', 'P51667', 'P48962', 'Q64105', 'Q01853', 'Q9DB77', 'O08749', 'P19783', 'P38647', 'P04117',
             'O08709', 'P20152', 'P63017', 'P05064', 'Q9CXZ1'}
    geneset2 = {'P17182', 'P41216', 'P99029', 'A2AQP0', 'P20029', 'P45376', 'P97457', 'P61982', 'P14824', 'P47857',
             'Q99LC3', 'P28474', 'Q8CGP0', 'Q8VHX6', 'P62259', 'Q8CI51', 'Q03265', 'Q8BMS1', 'Q9ET01', 'Q3V1D3',
             'P09542', 'P56480', 'Q9JI91', 'Q9R062', 'Q9QZ47', 'P68369', 'P06151', 'O70209', 'Q9DCS9', 'P05063',
             'Q9Z2I9', 'Q9CZU6', 'P05201', 'Q9JKS4', 'Q9D6J6', 'P11499', 'Q8VEM8', 'Q5SX39', 'P99024', 'P05977',
            'P62869', 'Q9D6Y9', 'P58252', 'P68510', 'O54724',
             'O08539', 'P68372', 'P70404', 'Q60854', 'Q61171', 'P11404', 'P24527', 'P58774', 'Q60932', 'E9Q557',
             'Q8R429', 'P42232', 'P01942', 'P18572', 'E9PZQ0', 'A2AUC9', 'Q8BTM8', 'P23927', 'Q91VR2', 'P70670',
             'Q9DC70', 'O88346', 'P08228', 'Q64522', 'Q9JK37', 'P16015', 'Q91Z83', 'Q9D0F9', 'Q9DCW4', 'P19123',
             'Q8CHT0', 'P62631', 'P63038', 'P04247', 'Q9R0Y5', 'P14602', 'Q9Z2I0', 'P57776', 'Q60930', 'P56375',
             'P35700', 'O08553', 'Q99MN9', 'P10649', 'Q8VDN2', 'Q9QXS1', 'Q6P8J7', 'Q9CWF2', 'O70622', 'Q9Z2U0',
             'P07901', 'Q9D6R2', 'P62962', 'P60843', 'P15864', 'Q922U2', 'P13707', 'P70349', 'P67778', 'P70296',
             'Q8BWT1', 'Q91ZJ5', 'A2AMM0', 'P0CG50', 'P06745', 'Q99LC5', 'Q9CZ13', 'Q62234', 'Q8CAQ8', 'Q63836',
             'P09411', 'A2ASS6', 'P21107', 'Q9D6F9', 'Q9EQ20', 'P97807', 'Q60597', 'Q9ESD7', 'Q3TXS7', 'P17742',
             'P09041', 'P51667', 'P48962', 'Q64105', 'Q01853', 'Q9DB77', 'O08749', 'P19783', 'P38647', 'P04117',
             'O08709', 'P20152', 'P63017', 'P05064', 'Q9CXZ1'}

    col_of_interest = "Abundance rate"
    compare_by_genes_paired(infile1, infile2, geneset1, col_of_interest)
    compare_by_genes_independent(infile1, infile2, geneset1, geneset2, col_of_interest)

