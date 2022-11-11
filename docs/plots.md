# Plots

## Rationale
This script will create various additional plots.

## Steps

The script is ``drafts_results.R`` and has some command line parameters :
- __tpm_csv__ : path to csv TPM counts file
- __isg_file__ : directory with ISG txt lists
- __annot_csv__ : path to gene annotation csv
- __metadata_file__ : path to samples metadata tsv
- __filter__ : TPM filter
- __jfile__ : json of arguments for individual figures
- __outpath__ : output path

Additional arguments found in the json are :

- __HM_associated__ :
  - _make :_ wether to make a heatmap showing the Z-scores of the TPM values of genes associated with a TF of interest (TRUE or FALSE)
  - _title :_ plot title
  - _file_int :_ path to the file of interest containing the genes associated with a certain TF
  - _filename :_ name of the file
- __TF__ :
  - _make :_ wether to make a dotplot summarazing the results of the GSEA ranking for our differential genes sets (TRUE or FALSE)
  - _filepaths :_ path to file of TFEA.ChIP GSEA ranking results
  - _names :_ name list of the sets, in the order they are in the directory
  - _orders :_ integer list orders in which to display them in the plot
  - _title :_ plot title
  - _filename :_ name of the file
- __Scatter_plot__ :
  - _make :_ wether to make a scatter plot showing the logFC of genes in two different conditions (TRUE or FALSE)
  - _list_DE_file :_ list with the two differential genes table
  - _title :_ title of the plot
  - _filename :_ name of the file
- __Recap_HM__ :
  - _make :_ wether to make a a Z-scores TPM values of heatmap of the union of DE genes across sets (TRUE or FALSE)
  - _list_DE_file :_ list of differential gene files of interest
  - _number_of_plots :_ number of plots to create : it will do so by iterating over the list_DE_file with a step equal to this number. So your list need to be a multiple of this number. If you have four tables in your list and ask for 2 plots, it will make one plot with the first two files, and one with the last two.
  - _titles :_ titles of each plot, must be of length ``number_of_plots``
  - _filenames :_ name of the files, must be of length ``number_of_plots``
  - _list_genes :_ list of genes to display in the HM, empty if you don't want to show gene names
  - _n :_ maximum DE genes to display
- __Pathway_plots__ :
  - _make :_ wether to make barplots depicting the adjusted pvalue for a pathway (TRUE or FALSE)
  - _list_pathways :_ list to pathways results in Robj format
  - _list_DE_file :_ corresponding DE tables
  - _titles :_ plot titles, as many as elements in the former lists
  - _filenames :_ file names
- __Venn_DE_sets__ :
  - _make :_ wether to make venn diagram of three sets of genes intersecting (TRUE or FALSE)
  - _list_DE_files :_ list of DE tables to intersect, length must be a multiple of 3
  - _filenames :_ names of the files
  - _categories :_ list of list of category names for each plot
- __Boxplot_IFN__ :
  - _make :_ wether to make boxplots depicting the log2(TPM + 1) values of inteferon genes (TRUE or FALSE)

## Usage example

``Rscript drafts_results.R ../data/tablecounts_tpm.csv ../data/ISGs/ ../data/tableannot.csv ../data/metadata.tsv 1 figures_filter1.json ../output/Draft_figures/``
