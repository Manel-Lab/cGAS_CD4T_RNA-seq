# Transcription factor enrichment analysis

## Rationale
The goal is to perform transcription factor enrichment analysis using the TFEA.ChiP package on various gene sets.

## Steps

The script is ``TF_enrichment.R`` and has some command line parameters :
- __list_DE__ : path to the differential gene tables
- __outpath__ : output path
- __isg__ : logical, keep ISG genes (TRUE) for the analysis or exclude them (FALSE)
- __list_TFs__ : list of TF of interest to highlight in plot or in tables

This script will create directories inside the output path :
- in ``ranking``, result of the GSEA ranking of the Fisher over-representation analysis, in tsv format aswell as in a GSEA ranking plot.
- in ``over-representation``, result of the over-representation analysis in table format or in a plot showing the Odd Ratio and the pvalue of the analysis for each ChIP experiment in the database.
- in ``Genes_associated_with_TF``, tables showing the DE genes associated with ChiP with our TFs of interest
- in ``IFN_Encode``, tables showing IFN with peaks in all the experiments in the Encode database.
- in ``TF_summarize``, same tables as before but summarized to TF instead of showing all experiments.

## Usage example

``Rscript TF_enrichment ../output/DE_tables ../output TRUE "c('RELA', 'IRF3')"``
