# Downstream analysis of cGAMP/TSA RNA-seq

## Introduction

This repository contains R scripts analyzing a RNA-seq of CD4T+ cells treated with cGAMP, TSA, both, or neither.

This downstream analysis is comprised of several steps than can be executed separately. However, for reproducibility purposes, a R script called ``downstream_pipeline.R`` is available to re-run the full pipeline from  start to finish with parameters used for the published results. It takes as arguments :

1. Path to the raw counts csv table
2. Path to TPM normalized counts csv table
3. Path to gene annotation csv table
4. Path to metadata tsv table
5. Path to txt file with list of ISG (one per line)
6. Output path for the whole analysis
7. Differential gene table of a previous experiment
8. Json file for the figures part

## Steps

1. [Creating DESeq model](docs/model.md)
2. [Differential analysis](docs/de_analysis.md)
3. [Pathway over-representation analysis](docs/pathway.md)
4. [Transcription factor enrichment analysis](docs/tfactor.md)
5. [Figures](docs/plots.md)
