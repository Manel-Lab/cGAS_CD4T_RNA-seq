## ---------------------------
##
## Purpose of script: Utiliy functions for DE analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-02-19
##
## ---------------------------

## Packages

library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)

## ------ FUNCTIONS ------- ##

# Import counts frrom a csv table
# csv: file path to the csv file
import_counts <- function(csv){
  counts <- read_csv(csv) %>%
    dplyr::rename(Gene = X1)
  return(counts)
}

# Write a count table with annotations to tsv format
# tibble: count tibble
# annot: annotation table
# filepath: output file path
tsv_counts <- function(tibble, annot, filepath) {
  my_tsv <- tibble %>%
    left_join(annot) %>%
    dplyr::select(-Haplotype)
  
  write_tsv(my_tsv, filepath)
}

# Get the gene set for a DESeq result object
# resDESeq = results DESeq object
# annot = csv with symbol/ENSEMBL correspondance
# n = number of the set
# pval = pvalue cutoff
# logFC = log fold-change cutoff
get_set <- function(resDESeq, annot, n, pval = NA, logFC = 0.5) {
  set <- as_tibble(resDESeq) %>% 
    add_column(Gene = rownames(resDESeq), .before = 1) %>%
    # Put the set number as a column
    add_column(Set = n) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    # Filter on adjusted pvalue and logFC
    filter(abs(log2FoldChange) >= logFC) %>%
    left_join(annot) %>%
    arrange(desc(abs(log2FoldChange)))
  
  if(!is.na(pval)){
    set <- set %>% filter(padj < pval)
  }
  
  return(set)
}

# Return a list of sets 
# Plus write in a tsv UP and DOWN regulated genes for a set
# list_sets = empty list
# set = a gene set for "get_set" function
# outpath = output path 
# i = number of the set
# logFC_t = logFC cutoff
list_sets <- function(list_sets, set, outpath = NA, i, name, logFC_t = 0.5) {
  up <- set %>%
    filter(log2FoldChange > logFC_t)
  down <- set %>%
    filter(log2FoldChange < (-1*logFC_t))
  # Add that particular set of genes to the
  list_sets[[i]] <- set
  if(!is.na(outpath)){
    write_tsv(up, file.path(outpath, paste0("UP_", name, ".tsv")))
    write_tsv(down, file.path(outpath, paste0("DOWN_", name, ".tsv")))
    write_tsv(set, file.path(outpath, paste0("FULL_", name, ".tsv")))
  }
  return(list_sets)
}

# Create a dataframe recapitulating the genes differentially
# expressed across set and with what fold change
# list_of_sets = list of genes sets
# outputpath = output path
# namedlist = correspondance namedlist between sets and conditions
data_RecapHM <- function(list_of_sets, outputpath, namedlist){
  my_table <- dplyr::bind_rows(list_of_sets) 
  for(i in unique(my_table$Set)){
    # Column for set number i
    column1 <- namedlist[[i]]
    # Column for logFC value in set number i
    column2 <- paste0("LogFC_", column1)
    my_table <- my_table %>% 
      # Boolean column : 1 if differential in this set
      mutate(!!column1 := ifelse(Set == i, 1, 0)) %>%
      # Value is LogFC if differential, else NA
      mutate(!!column2 := ifelse(Set == i, log2FoldChange, 0))
  }
  
  # Collapse table
  my_table <- my_table %>% 
    dplyr::select(-baseMean, -log2FoldChange, -lfcSE, -stat, -pvalue, -padj, -Set, -Gene) %>% 
    group_by(Symbol) %>% 
    # There is i value of logFC and set for each gene
    # all in i different columns 
    # therefore we can sum while collapsing to keep one line per gene
    # Without modyfing the logFC for a set
    summarize_if(is.numeric, sum) %>%
    mutate_if(is.numeric,  list(~na_if(., 0)))
  
  write_tsv(my_table, file.path(outputpath, "resume_diffgenes.tsv"))
  return(my_table)
}

# Filter the recapitulating dataframe from data_RecapHM function
# with a logFC based threshold
# my_table = result of data_RecapHM function 
# threshold = logFC threshold
filter_table <- function(my_table, threshold = 2) {
  # Apply a threshold of LogFC
  my_table <- dplyr::filter_if(my_table, is.numeric, any_vars(abs(.) > threshold))
  # Rename some variable for Heatmap
  my_table <- dplyr::select(my_table, starts_with(c("Symbol", "Log")))
  
  names(my_table) <- names(my_table) %>%
    gsub("logFC_", "", .)
  
  return(my_table)
}

# Create a Heatmap of logFC across sets
# mat = matrix with genes in columns and sets in rows
# colors = color vector
# annotation_col = dataframe annotating the genes
# ann_color = color for annotations
# myBreaks = vector of breaks
recap_heatmap <- function(mat, colors, cluster_rows = FALSE, cluster_cols = TRUE, display_numbers = TRUE,
                          annotation_col = NULL, ann_color = NULL, myBreaks = NULL, title = NA) {
  
  if(is.null(myBreaks)){
    myBreaks <- c(seq(min(mat), 0, length.out=ceiling(11/2) + 1), 
                  seq(max(mat)/11, max(mat), 
                      length.out=floor(11/2)
                  )
    )
  }
  
  logFCheatmap <- pheatmap::pheatmap(mat, color=colors, display_numbers = display_numbers,
                                     cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                                     breaks=myBreaks, main = title,
                                     annotation_col = annotation_col, ann_color = ann_color,
                                     cellwidth = 22, fontsize = 12, cellheight = 22, fontsize_number = 8)
  
  return(logFCheatmap)
}

# Heatmap of differential genes
# mat = matrix with genes in columns and samples in rows
# row_annot = annotation of samples
# ann_color = color for annotations
# cols_annotation = annotation for genes
# cluster_rows = cluster samples ?
# cluster_cols = cluster_genes ?
# myBreaks = vector of breaks
# ... = additional arguments passed to pheatmap::pheatmap function
heatmap_DEG <- function(mat, colors, row_annot = NULL, ann_color = NULL, cols_annotation = NULL,
                        myBreaks = NULL, title=NA,
                        ...) {
  
  HM <- pheatmap::pheatmap(mat, color=colors, breaks = myBreaks,
                           annotation_col = cols_annotation, annotation_row = row_annot, annotation_colors = ann_color,
                           fontsize = 7, scale="none", clustering_method="complete", main=title,
                           ...)
  return(HM)
}

# Create all volcano plots
# list_res = list of results from DE analysis
# outpath = output path of the volcano plots
# annot = annotation table
# include_isg = wether to filter out ISGs (FALSE) or not (TRUE)
create_volcanos <- function(list_res, list_names, outpath, annot, include_isg = TRUE) {
  
  for(i in seq_along(list_res)){
    
    if(!include_isg){
      # Filter out ISG
      my_res <- list_res[[i]][rownames(list_res[[i]]) %in% annot[annot$ISG == FALSE,]$Gene,]
    }
    else{
      my_res <- list_res[[i]]
    }
    
    vlcplot <- EnhancedVolcano(my_res,
                               lab =  annot[match(rownames(my_res), annot$Gene),]$Symbol,
                               x = 'log2FoldChange',
                               y = 'padj',
                               xlim = c(-8, 8),
                               title = list_names[[i]],
                               pCutoff = 0.05,
                               FCcutoff = 0.5)
    pdf(file.path(outpath, paste0(list_names[[i]], "_volcanoplot.pdf")), height=15, width=18)
    print(vlcplot)
    dev.off()
  }
  
}

## --------- RUN ---------- ##
