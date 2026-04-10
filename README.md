# **Local David Database**

## **Quick Use Guide:**
### 1. Download DAVID Knowledgebase Files 

- Go to the link below and download any knowledgebase files you are interested in using.  
[Download DAVID Knowledgebase Files Here](https://davidbioinformatics.nih.gov/knowledgebase/knowledgebaseRequest.jsp)  
- Put those files in the folder named "knowledgebase".

### 2. Select Gene Files to Compare

- Add two .csv or .tsv files to the inputs folder. These files must at least contain an ID column and some numeric column that you want to compare. The names of these columns must match between the two files.

### 3. Run main.py

- Update the inputs to the main function call at the bottom of main.py.  
- Run main.py, outputs will go to the outputs folder.



### Credit:  
B.T. Sherman, M. Hao, J. Qiu, X. Jiao, M.W. Baseler, H.C. Lane, T. Imamichi and W. Chang. DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update). Nucleic Acids Research. 23 March 2022. doi:10.1093/nar/gkac194.[PubMed]  

Huang DW, Sherman BT, Lempicki RA. Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources. Nature Protoc. 2009;4(1):44-57.  [PubMed]

