# Pathway over-representation analysis

## Rationale
The goal is to perform pathway over-representation analysis with DOSE on various lists of differential genes.

## Steps

The script is ``enrichment_analysis.R`` and has some command line parameters :
- __dir_tables__ : directory to the differential genes table
- __outpath__ : output path
- __threshold__ : LogFC threshold (genes below that threshold ar ignores for the analysis)
- __pval__  : differential analysis adjusted pvalue threshold
- __isg__ : logical, keep ISG genes (TRUE) for the analysis or exclude them (FALSE)

This script will create directories inside the output path. Inside will be for each gene set the results in tsv format, a R object with the result and the logFC of the genes, and various plots summarizing the results.

## Usage example

``Rscript enrichment_analysis.R ../output/DE_tables ../output 1 5e-2 TRUE``
