# **Local David Database**



## **Quick Use Guide:**
### 1. Download DAVID Knowledgebase Files 

- Go to the link below and download any knowledgebase files you are interested in using.  
[Download DAVID Knowledgebase Files Here](https://davidbioinformatics.nih.gov/knowledgebase/knowledgebaseRequest.jsp)  
- Put those files in the folder named "knowledgebase".

### 2. Select Gene Files to Compare

- Add two .csv or .tsv files to the `inputs` folder. These files must at least contain an ID column and some numeric column that you want to compare. The names of these columns must match between the two files.

### 3. Run main.py

- Run main.py.
- Fill in the top 4 values (seen below) in the UI menu to match your input files. Optionally change the other settings as needed.
  1. *Input File 1: One input file containing at least an ID column and a numeric column to compare.*
  2. *Input File 2: A second input file for comparison. ID column and numeric column headers must match.* 
  3. *ID Column: The column header for the ID column in both input files.* 
  4. *Numeric Column: The column header for the numeric column of interest in both input files.*
- Outputs will be written to the `outputs` folder.



## **Pipeline:**  
This tool compares two proteomics datasets by functional annotation category, identifying which biological terms show statistically significant differences between the two samples.  

### Stage 1 — Input Collection  

A Tkinter GUI collects the two input files, the ID column name, and the numeric column of interest, along with optional settings. Parameters are passed directly into main().  


### Stage 2 — Gene ID Extraction (extract_id_column)  

The specified ID column (e.g. UniProt accession IDs) is read from each input file into a Python set. Both .csv and .tsv formats are supported.  


### Stage 3 — Functional Annotation (do_functional_annotation)  

Each gene set is independently annotated against a local DAVID-style knowledgebase. For every .txt mapping file in the knowledgebase directory (each representing one annotation category such as GO terms or KEGG pathways), a hypergeometric enrichment test is run to determine whether genes from the input are overrepresented in each term compared to the background. Results are corrected for multiple testing using both Bonferroni and Benjamini-Hochberg (FDR) methods, then written to a .tsv output file. The result is a dataframe of enriched annotation terms with p-values, fold enrichment, and the list of genes driving each term.  


### Stage 4 — Finding Common Terms (find_common_terms)

The two annotation dataframes are compared to find terms that appear in both. Terms are ranked by the minimum gene count shared between the two datasets, and the top N are kept. Depending on the only_common_genes setting:  

True — only genes present in both datasets for a given term are kept (intersection)  
False — the full gene lists for each dataset are kept separately for independent testing  

Results are written to a shared annotations output file.  


### Stage 5 — Statistical Comparison (compare_by_genes)  

For each of the top shared annotation terms, a t-test is run on the numeric column of interest (e.g. abundance rate) between the two input files, using only the genes belonging to that term:

only_common_genes=True → paired t-test on the shared gene set
only_common_genes=False → Welch's independent t-test on the separate gene lists

The final output is a .tsv file with the t-statistic, p-value, and gene count for each annotation term.



## Credit:  
B.T. Sherman, M. Hao, J. Qiu, X. Jiao, M.W. Baseler, H.C. Lane, T. Imamichi and W. Chang. DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update). Nucleic Acids Research. 23 March 2022. doi:10.1093/nar/gkac194.[PubMed]  

Huang DW, Sherman BT, Lempicki RA. Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources. Nature Protoc. 2009;4(1):44-57.  [PubMed]

