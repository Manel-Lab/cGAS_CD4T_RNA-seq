# Create DESeq2 model

## Rationale
The goal is to add to a create a DESeq2 model for the project. It also creates several figures to assess the quality of the model aswell as the data's.

## Steps

The script is ``create_DESeq_model.R`` has some command line parameters :
- __csvRaw__ : Path to the raw counts csv table
- __csvTpm__ : Path to TPM normalized counts csv table
- __csvAnnot__ :  Path to gene annotation csv table
- __tsvMetadata__ : Path to metadata tsv table
- __dir_ISG__ : Path to directory with ISG list, with a text file per list inside
- __outpath__ : Output path for the whole analysis

It will create in the output path qqplot and distribution plot with the count data, aswell as PCA, and also a Robj containing the DESeq model.

## Usage example

``Rscript create_DESeq_model.R ../data/tablecounts_raw.csv ../data/tablecounts_tpm.csv ../data/tableannot.csv ../data/metadata.tsv ../data/dir_ISGs ../output``
