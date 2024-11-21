# G9-HTG


## Supplementary Resource to [PMID 39566464](https://pubmed.ncbi.nlm.nih.gov/39566464/):  

Denkert C et al. Cell Rep Med. 2024 Nov 19;5(11):101825.

Molecular adaptation to neoadjuvant immunotherapy in triple-negative breast cancer.

[DOI: 10.1016/j.xcrm.2024.101825](https://doi.org/10.1016/j.xcrm.2024.101825)

************************************************************

## This resource contains the following data regarding the analyses described in the paper:


1. [*2024-10-18-G9-HTG-Suppl.pdf*](https://github.com/tkarn/G9-HTG/blob/main/2024-10-18-G9-HTG-Suppl.pdf):  An R Markdown file containing the analyses and generation of figures in the paper based on OR and HR data for 2549 genes from the HTG EdgeSeq Oncology Biomarker Panel measured on the HTG-Molecular EdgeSeq platform. This file also includes validation analyses and figures in three additional TNBC cohorts with RNA-Seq data for the respective genes.
2. [*2024-10-18-G9-HTG-Suppl.R*](https://github.com/tkarn/G9-HTG/blob/main/2024-10-18-G9-HTG-Suppl.R):  The R-script that generates this R-Markdown file and the output below.


## Directory [*input*](https://github.com/tkarn/G9-HTG/tree/main/input/) contains:
1. [*G9_combined_results_all_genes_24JUL2023.xlsx*](https://github.com/tkarn/G9-HTG/tree/main/input/G9_combined_results_all_genes_24JUL2023.xlsx):  An MS Excel file containing response (pCR) and survival (DDFS) data separately for each GeparNuevo treatment arm for all 2549 HTG-genes, as well as longditudinal changes in expression during treatment along with annotation of gene signalling pathways. This file contains the input data for all analyses on GeparNuevo samples in the R script above.
2. [*valid_datasets_HTG.RData*](https://github.com/tkarn/G9-HTG/tree/main/input/valid_datasets_HTG.RData):  An RData file containing both clinical data and RNA-seq data for the 2549 HTG-genes from three validation cohorts. This file provides the input data for the validation analyses in the R script above. 


## Directory [*out*](https://github.com/tkarn/G9-HTG/tree/main/out/) contains: 
Figures and data generated by the script above and presented in the paper.

************************************************************

 
