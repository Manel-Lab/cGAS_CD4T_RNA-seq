## ---------------------------
##
## Purpose of script: Pathway enrichment analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-03-04
##
## ---------------------------

## Arguments

if (sys.nframe() == 0){
  .libPaths(c("~/R_packages",.libPaths()))
  ## Arguments
  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  dir_tables <- args[1] # Directory of the DE tables
  outpath <- args[2] # Output path
  threshold <- as.numeric(args[3]) # LogFC threshold
  pval <- as.numeric(args[4]) # pval threshold
  isg <- as.logical(args[5]) # keep ISG (TRUE) or exclude them (FALSE)
}

## ---------------------------

## Packages

library(clusterProfiler)
library("org.Hs.eg.db")
library(enrichplot)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(doParallel)

## ------ FUNCTIONS ------- ##

# Enrichment analysis with GO annotations
# ids = list of genes
# ont = type of ontologie (CC, BP, MF)
# keyType = annotation of genes in list
do_GO <- function(ids, ont, keyType = "ENTREZID"){
  enrich <- enrichGO(ids, OrgDb = org.Hs.eg.db, ont=ont, keyType=keyType)
  enrich <- clusterProfiler::simplify(enrich, cutoff=0.7, by="p.adjust", select_fun=min)
  return(enrich)
}

# Enrichment analysis with KEGG annotations
# ids = list of genes
do_KEGG<- function(ids){
  return(enrichKEGG(ids))
}

# Personal heatplot
# enrich = result from enrichment analysis
# my_tibble = DE genes
# positive = boolean (TRUE if there are over expressed pathways)
# negative = boolean (TRUE if there are under expressed pathways)
# analysis = top_pathways or top_genes
# n = number of pathways/gene
# ... = additional arguments passed to pheatmap::pheatmap
myHeatplot <- function(enrich, my_tibble, positive, negative,
                       analysis = "top_pathways", n = 50,
                       ...) {
  
  # Different color schemes if plotting under or over expressed pathways
  if(positive && negative){
    colors <- rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu"))
  } else if(positive){
    colors <- RColorBrewer::brewer.pal(n = 7, name ="Reds")
  } else{
    colors <- rev(RColorBrewer::brewer.pal(n = 7, name ="Blues"))
  }
  
  
  # Filter the enrichment result
  logtib <- my_tibble[c("Symbol", "log2FoldChange")]
  # Extract pathways and symbols from enrichment results
  pathtib <- as_tibble(na.omit(enrich@result))
  
  if(analysis == "top_pathways"){
    # Take the top 20 pathways (pvalue ordered)
    top_path <- pathtib[1:n,]$Description
  }
  if(analysis == "top_genes"){
    # Top 50 genes (logFC ordered)
    top_genes <- logtib[1:n,]$Symbol
  }
  
  pathtib <- pathtib[c("Description", "geneID")] %>%
    # Split geneID columns at each /
    separate_rows(geneID,sep="/") %>% 
    # Join
    inner_join(logtib, by=c("geneID" = "Symbol"))
  
  if(analysis == "top_genes"){
    pathtib <- pathtib %>% filter(geneID %in% top_genes)
    main = "Top DE genes and their enriched pathways"
    if(dim(pathtib[1]) < 1){
      print("No top DE genes are over-enriched in pathways")
      return(NULL)
    }
  } else if(analysis == "top_pathways"){
    main = "Top enriched pathways and their genes"
    pathtib <- pathtib %>% filter(Description %in% top_path)
  } else{
    print("analysis argument is not valid.")
    return(NULL)
  }

  # Filter
  pathtib <- pathtib %>% 
    # Keep genes appearing in at least two pathways
    #filter(geneID %in% pathtib$geneID[duplicated(pathtib$geneID)]) %>%
    # Keep pathways with at least two genes
    #filter(Description %in% pathtib$Description[duplicated(pathtib$Description)]) %>%
    pivot_wider(names_from="Description", values_from="log2FoldChange")
  
  # Plus clustering
  # Get enriched pathways and create a boolean clustering
  clust <-  as_tibble(enrich@result)[c("Description", "geneID")] %>%
    # Split geneID columns at each /
    separate_rows(geneID,sep="/") %>%
    add_column(Value = 1)
  # Create a tibble of boolean
  clust <- pivot_wider(clust, names_from = "Description", values_from="Value") %>%
    replace(is.na(.), 0) %>%
    replace(NULL, 0)
  # Filter gene IDs to those in the results of pathways enrichment analysis
  clust <- clust %>% filter(geneID %in% pathtib$geneID) %>%
    arrange(match(geneID, pathtib$geneID))
  
  # Create a matrix of boolean
  clustmat <- as.matrix(clust[-1])
  rownames(clustmat) <- clust$geneID
  
  if(dim(pathtib)[1] > 2 & dim(pathtib)[2] > 2){
    #Creating the plot itself
    mat <- t(as.matrix(pathtib[-1]))
    colnames(mat) <- pathtib$geneID
    HM <- pheatmap::pheatmap(mat, color = colors, main=main, cluster_rows=FALSE, 
                             cluster_cols=hclust(dist(clustmat, "binary")), na_col="#FFFFFF",
                             ...)
    return(list(HM, pathtib))
    }
  return(NULL)
}


# Enrichment analysis function
# df = dataframe with differential genes in rownames
# pdf_file = name of the pdf_file with output plots 
# FUN = which enrichment function to use (do_KEGG or do_GO)
# + additional arguments passed to FUN
sc_enrichment <- function(my_tibble, outpath, filename, FUN, ...){
  # Get Entrez Id for Ensembl and join with symbols
  entrez_ids <- bitr(my_tibble$Gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
  entrez_ids <- entrez_ids %>% left_join(my_tibble[c("Gene", "Symbol")], by=c("ENSEMBL" = "Gene"))
  # Do the enrichment analysis
  enrich <- FUN(entrez_ids$ENTREZID, ...)
  
  # If the enrichment analysis is conclusive
  if(!is.null(enrich)){
    for(i in seq_along(enrich@result$geneID)){
      # Get the gene of each pathway
      init <- ""
      ids <- unlist(strsplit(enrich@result$geneID[i], split = "/"))
      for (j in seq_along(ids)){
        j <- entrez_ids[(match(ids[j], entrez_ids$ENTREZID)), "Symbol"]
        init <- paste(init, j, sep = "/")
      }
      init <- sub('/', '', init)
      enrich@result$geneID[i] <- init
    }
    
    # Get the Fold Change to plot
    geneList <- my_tibble[my_tibble$Symbol %in% entrez_ids[!is.na(entrez_ids$ENTREZID),]$Symbol,]
    foldChange <- geneList$log2FoldChange
    names(foldChange) <- geneList$Symbol
    
    if(!is.na(outpath)){
      save(enrich, foldChange, list=c("enrich", "foldChange"), 
           file=file.path(outpath, paste0(filename, ".Robj")))
    }
    return(list(enrich, foldChange))
  }
  return(NULL)
}


# Create the plots for the pathways
# enrich = result from enrichment pathways
# foldChange = sorted gene list with foldChange
# my_tibble = DE genes
# outpath = output path
# filename = name of the files
sc_plotting <- function(enrich, foldChange, my_tibble, outpath, filename) {
  # Positive is true if there are overexpressed genes
  positive <- any(foldChange > 0)
  # Negative is true if there are underexpessed genes
  negative <- any(foldChange < 0)
  # Filter enrichment value
  enrich <- enrich[enrich$Count > 5, asis=T]
  enrich <- enrich[enrich$pvalue < 0.05, asis=T]
  rowN <- dim(data.frame(enrich))[1]
  
  # Plot if there is at least two pathways resulting from the filtering
  if(rowN > 1){
    
    dotP <- dotplot(enrich, showCategory = min(rowN, 20), font.size = 10) + 
      ggtitle(paste0("Dotplot for ", filename))
    #heatP_gene <- myHeatplot(enrich, my_tibble, positive, negative, "top_genes")
    #heatP_path <- myHeatplot(enrich, my_tibble, positive, negative, "top_pathways")
    upsetP <- upsetplot(enrich, n = min(rowN, 20))
    cnetP <- cnetplot(enrich, showCategory = min(rowN, 20), font.size = 8, categorySize = "pvalue", foldChange = foldChange) + 
      ggtitle(paste0("Cnetplot for ", filename))
    emaP <- emapplot(enrich, showCategory = min(rowN, 20), font.size = 8) +
      ggtitle(paste0("Emaplot for ", filename))
    
    #Dotplot
    pdf(file.path(outpath, paste0("Dotplot_", paste0(filename, ".pdf"))), height=8, width=18)
    print(dotP)
    dev.off()
    
    #Heatplot
    #pdf(file.path(outpath, paste0("Heatmap_topGenes_", paste0(filename, ".pdf"))), height=7, width=25)
    #print(heatP_gene[[1]])
    #dev.off()
    #pdf(file.path(outpath, paste0("Heatmap_topPath_", paste0(filename, ".pdf"))), height=7, width=25)
    #print(heatP_path[[2]])
    #dev.off()
    
    #Upsetplot
    pdf(file.path(outpath, paste0("Upsetplot_", paste0(filename, ".pdf"))), height=14, width=25)
    print(upsetP)
    dev.off()
    
    # Cnetplot & Emaplot
    pdf(file.path(outpath, paste0("Cnetplot_", paste0(filename, ".pdf"))), height=16, width=16)
    print(cnetP)
    print(emaP)
    dev.off()
  }
}

# Combine pdfs with pdftk
# path = path where the single pdfs are
# output_pdf = full path of the output pdf
combine_pdfs = function(path, output_pdf){
  system(sprintf("pdftk %s/*pdf cat output %s", path, output_pdf))
}

# Limit genes for pathway analysis
# to top hits 
# based on a LogFC threshold
# my_tibble = DE genes table
# threshold = logFC threshold of detection
top_hits <- function(my_tibble, threshold, pval){
  my_tibble <- my_tibble %>% 
    dplyr::filter(abs(log2FoldChange) > threshold) %>%
    dplyr::filter(padj < pval)
  return(my_tibble)
}

# Limit the number of pathways considered
# based on a top X p.adjust threshold
# Or a gene Count threshold
# df = dataframe of pathways
# pcutoff = pvalue cutoff
# count = minimum number of genes for this pathway
top_pathways <- function(df, pcutoff = NULL, count = NULL){
  if(!is.null(p.adjust)){
    df <- df %>% arrange(dplyr::desc(p.adjust)) %>% top_n(pcutoff)
  }
  if(!is.null(count)){
    df <- df %>% dplyr::filter(Count > count)
  }
  return(df)
}

# Function to cluster genes according to the
# pathways they belong to
# enrich = enrichment results
cluster_pathways <- function(my_tibble, enrich){
  
  # Get enriched pathways and create a boolean clustering
  pathtib <-  as_tibble(enrich@result)[c("Description", "geneID")] %>%
    # Split geneID columns at each /
    separate_rows(geneID,sep="/") %>%
    add_column(Value = 1)
  # Create a tibble of boolean
  pathtib <- pivot_wider(pathtib, names_from = "Description", values_from="Value") %>%
    replace(is.na(.), 0)
  
  # Need at least 3 genes to cluster
  if(nrow(pathib) <= 2){
    return(NULL)
  }
  
  # Create a matrix of boolean
  clust <- as.matrix(pathtib[-1])
  rownames(clust) <- pathtib$geneID
  return(clust)
  
}

# Create enrichment analysis plots and table
# dir_tables = directory with DE tables
# outpath = output path
# threshold = logFC threshold of input DE genes
# pval = adjusted pvalue threshold of input DE genes
# concat = for a gene set, concatenate all plots in one PDF (pdftk need to be installed)
# isg = logical, wether to filter ISG genes out (FALSE) or not (TRUE)
main_enrichment <- function(dir_tables, outpath, threshold, pval, isg) {
  # Set seed
  set.seed(0)
  # Create the output path in case it doesn't exist
  dir.create(outpath)
  dir.create(file.path(outpath, "UP"))
  dir.create(file.path(outpath, "DOWN"))
  dir.create(file.path(outpath, "FULL"))
  
  # List of comparisons
  list_names = list("+cGAMP-TSA_+cGAMP+TSA", "-cGAMP-TSA_+cGAMP-TSA",
                    "-cGAMP-TSA_-cGAMP+TSA", "Diff_effect_TSA")
  
  # UP & DOWN & FULL
  directory <- list.files(dir_tables)
  print(length(directory))
  for(i in 1:length(directory)) {
    print(directory[i])
    if(!grepl("set5", directory[i], fixed = TRUE)){
      # Construct the name of the file with different elements
      strfile <- strsplit(strsplit(directory[i], ".tsv")[[1]][1], "_")[[1]]
      # The name of the set is one element
      set <- paste(strfile[-1], collapse="_")
      # The sense (UP, DOWN) is another
      sense <- strfile[1]
      my_tibble <- read_tsv(file.path(dir_tables, directory[i]))
      my_tibble <- top_hits(my_tibble, threshold, pval)
      print("DE file read")
      if(!isg){
        my_tibble <- my_tibble %>% dplyr::filter(ISG == FALSE)
      }
      
      # Over representation analysis
      list_dbs <- c("CC", "MF", "BP", "Kegg")
      for(dbs in c("CC", "MF", "BP", "Kegg")){
        if(dbs == "Kegg"){
          list_sc <- sc_enrichment(my_tibble, file.path(outpath, sense), 
                                   paste0(set, "_", dbs, "_", sense),
                                   do_KEGG)
          print("Kegg enrichment done")
          enrich <- list_sc[[1]]
          foldChange <- list_sc[[2]]
          if(!is.null(enrich)){
            sc_plotting(enrich, foldChange, my_tibble, file.path(outpath, sense), 
                        paste0(set, "_", dbs, "_", sense))
            print("Kegg plotting done")
          }
        }
        else{
          list_sc <- sc_enrichment(my_tibble, file.path(outpath, sense), 
                                   paste0(set, "_", dbs, "_", sense),
                                   do_GO, dbs)
          print("GO enrichment done")
          enrich <- list_sc[[1]]
          foldChange <- list_sc[[2]]
          if(!is.null(enrich)){
            sc_plotting(enrich, foldChange, my_tibble, file.path(outpath, sense), 
                        paste0(set, "_", dbs, "_", sense))
            print("GO plotting done")
          }
        }
        if(!is.null(enrich)){
          write_tsv(as.tibble(enrich), file.path(outpath, sense, 
                                                 paste0("Results_", set, "_", dbs, "_", sense, ".tsv")))
          print("Enrichment results written to file")
        }
      }
    }
  }
}

## --------- RUN ---------- ##

if(sys.nframe() == 0){
  main_enrichment(dir_tables, outpath, threshold, pval, isg)
}
