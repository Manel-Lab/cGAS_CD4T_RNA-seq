# Differential analysis

## Rationale
The goal is to do a differential analysis using the DESeq model created. The comparisons of interest are :
- genes modified by cGAMP treatment in baseline
- genes modified by TSA treatment in baseline
- genes modified by TSA treatment when already treated by cGAMP
- genes that have a different effect in TSA treatment

## Steps

The script is ``DE_analysis.R`` has some command line parameters :
- __Robj__ : path to DESeq model R object
- __previous_tsv__ : Differential gene table of a previous experiment
- __outpath__ : Output path
- __logFC_t__ : LogFC threshold

This script will create the differential tables in a _"DE_tables"_ directory in the output path aswell as heatmaps of gene expression for the DE genes, heatmaps with the logFC of the genes in various comparisons, volcano plots and venn diagram intersecting the genes in common.

## Usage example

``Rscript DE_analysis.R ../output/DESeq_Model.R ../data/DEGS.tsv ../output 1``
