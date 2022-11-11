## ---------------------------
##
## Purpose of script: TF enrichment analysis
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
  
  list_DE <- args[1]
  outpath <- args[2]
  ISG <- as.logical(args[3]) # TRUE
  list_TFs <- eval(parse(text=args[4]))
}

## ---------------------------

## Packages

library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(TFEA.ChIP)
library(clusterProfiler)
library(org.Hs.eg.db)

## ------ FUNCTIONS ------- ##

# Main enrichment function
# DE_file = DE table tsv file
# outpath = output path
# ISG = logical, wether to filter ISG genes out (FALSE) or not (TRUE)
# list_TFs = list of TF of interest to highlight in over-representation plot
# trad = wether genes are already in Entrez IDs (FALSE) or not (TRUE)
main_TF_enrichment <- function(DE_file, outpath, ISG, list_TFs, trad){
  # Set seed
  set.seed(0)
  # Get name of file
  list_path <- strsplit(DE_file, "/")
  filename <- list_path[[1]][length(list_path[[1]])]
  title <- strsplit(filename, ".tsv")[[1]][1]
  
  # Import data
  DE_genes <- read_tsv(DE_file)
  # Filter out ISG genes
  if(!ISG){
    DE_genes <- DE_genes[DE_genes$ISG == FALSE,]
  }
  # If need to translate to Entrez ID
  if(!trad){
    colnames(DE_genes)[1] <- "Genes"
    colnames(DE_genes)[7] <- "pval.adj"
    # Get Entrez IDs
    DE_genes <- preprocessInputData(as.data.frame(DE_genes))
    # Define genes tested
    Genes_tested <- DE_genes$Genes
  } else {
    Genes_tested <- na.omit(DE_genes$Entrez)
  }
  
  # Create contigency matrix and test  
  mat <- contingency_matrix(Genes_tested)
  pval_mat <- getCMstats(mat)
  
  ####################
  ### GSEA RANKING ###
  ####################
  
  # Apply GSEA to data
  dir.create(file.path(outpath, "ranking"))
  pval_mat <- na.omit(pval_mat)
  if(dim(pval_mat)[1] > 2 & length(unique(pval_mat$TF)) > 1){
    TF_ranking <- rankTFs(pval_mat, rankMethod = "gsea", makePlot = TRUE,
                          plotTitle = paste0("Ranking of TFs enriched in DE genes of ", title))
    print(as_tibble(TF_ranking[["TF_ranking"]]))
    rank <- as_tibble(TF_ranking[["TF_ranking"]]) %>% dplyr::select(-`arg.ES`)
    p <- TF_ranking[["TFranking_plot"]]  
    # Plot the ranking
    pdf(file.path(outpath, "ranking", paste0(title, "_ranking_TFs.pdf")))
    print(p)
    dev.off()
    # Write the result of the GSEA analysis
    write_tsv(rank, file.path(outpath, "ranking", paste0(title, "_ranking_TFs.tsv")))
  
    # Plot the ranking results in an alternate way
    # Create a color column
    rank <- rank %>%
      mutate(Color = ifelse((abs(ES)>0.9 & pVal < 0.05) | TF %in% c("RELA", "IRF3"), "red", "black"))
  
      # Plot
    p <- ggplot(rank, aes(pVal, ES, color=Color)) +
      geom_point() +
      scale_color_identity() +
      geom_text_repel(aes(label=ifelse(Color=="red",as.character(TF),'')), hjust = 1) +
      labs(title=paste0("GSEA results of ", title, "\nGenes represented : ES > 0.9")) +
      geom_vline(xintercept=0.05, linetype="dashed", 
                 color = "red", size=0.5)
    pdf(file.path(outpath, "ranking", paste0(title, "_summary_GSEA_TF.pdf")), width=10)
    print(p)
    dev.off()
  }
  
  
  # Make matrix as tibble
  pval_tib <- as_tibble(pval_mat) %>%
    # Join with metadata
    left_join(MetaData, by=c("Accession", "Cell", "Treatment", "TF"))
  
  #########################
  ## OVER REPRESENTATION ##
  ### ANALYSIS RESULTS ####
  #########################
  
  # Create an OR limit
  OR_lim <- floor(max(pval_tib$OR)/1.5)
  if(OR_lim == Inf){OR_lim <- 4}
  # Create a color column
  pval_tib <- pval_tib %>%
    mutate(Color = ifelse(OR>OR_lim | TF %in% c("NFKB1", "NFKB2", "RELA", "RELB", "REL", "IRF3"), "red", "black"))
  
  p <- ggplot(pval_tib, aes_string("log10.adj.pVal", "log2.OR", color="Color")) +
    geom_point() +
    scale_color_identity() +
    geom_text_repel(aes(label=ifelse(Color=="red",as.character(TF),'')), hjust = 1) +
    labs(title=paste0("TFs analysis of ", title, "\nGenes represented : OR > ", OR_lim)) +
    geom_vline(xintercept=log10(0.05)*-1, linetype="dashed", 
               color = "red", size=0.5) +
    xlim(c(-0.5, max(pval_tib$`log10.adj.pVal`)+2))
  
  # Plot a "volcano" plot of the results
  dir.create(file.path(outpath, "over-representation"))
  pdf(file.path(outpath, "over-representation", paste0(title, "_plots_TFs.pdf")), width=10)
  print(p)
  dev.off()
  
  # Write the results to table
  pval_tib <- pval_tib %>%
    dplyr::select(-Color)
  pval_tib <- pval_tib[,c("Accession", "Cell", "Treatment","Cell.Type", "Antibody",	"TF",	"p.value",	
                          "OR",	"log2.OR", "adj.p.value",	"log10.adj.pVal",	"distance")]
  write_tsv(pval_tib, file.path(outpath, "over-representation", paste0(title, "_TF_table.tsv")))

  #########################
  ##### CHIP METADATA #####
  #########################
  
  ### Get the genes associated with TF of interest ###
  
  # Get Symbol of genes from our data
  chip <- as_tibble(Mat01) %>%
    add_column(ENTREZID = rownames(Mat01), .before=1)
  symbols <- bitr(chip$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb = org.Hs.eg.db, drop=FALSE)
  chip <- left_join(chip, as_tibble(symbols))
  chip <- chip[c(1, dim(chip)[2], 2:(dim(chip)[2]-1))]
  
  
  ## 1) For our TFs of interest
  ## Get the genes with binding sites
  ## According to the ChIP metadata
  dir.create(file.path(outpath, "Genes_associated_with_TF"))
  for(my_TF in list_TFs){
    # Get TF of interest encode annotations
    interest <- pval_tib %>% filter(TF %in% c(my_TF))
    # Get genes associated
    chip_int <- chip[c("ENTREZID", "SYMBOL", interest$Accession)]
    chip_int <- chip_int[rowSums(chip_int[-1:-2]) != 0,]
    # Get genes in common with our DE genes
    chip_int_cross <- filter(chip_int, SYMBOL %in% DE_genes$Symbol)
    write_tsv(chip_int_cross, file.path(outpath, "Genes_associated_with_TF", paste(title, my_TF, "associated_genes.tsv", sep="_")))
  }
  
  ### 2) Table of peaks detected for IFN genes
  dir.create(file.path(outpath, "IFN_Encode"))
  chip_genes <- chip %>% 
    filter(startsWith(SYMBOL, "IFN"))
  chip_genes <- chip_genes[, colSums(chip_genes != 0) > 0]
  write_tsv(chip_genes, file.path(outpath, "IFN_Encode", "IFN_genes_TFEA_all.tsv"))

  ### 3) Table of peaks in TF enriched in our experiment
  # Keep statistically significant TF
  pval_tib_filt <- pval_tib %>% filter(`adj.p.value` < 5.e-02)
  # Keep ChIP data if the TF is in our analysis result
  chip_genes <- chip_genes[,c(TRUE, TRUE, c(colnames(chip_genes) %in% pval_tib_filt$Accession)[-1:-2])]
  write_tsv(chip_genes, file.path(outpath, "IFN_Encode", paste(title, "IFN_genes_ENCODE-Chip.tsv", sep="_")))
  
  ### 4) Summarized long version of the previous table
  # Get the TF name without the experiment identifier
  TF_names <- lapply(colnames(chip_genes)[-1:-2], function(x) gsub(".*[.]([^.]+)[.].*", "\\1", x))
  colnames(chip_genes) <- c("ENTREZID", "SYMBOL", unlist(TF_names))
  # Create a long version of the previous table
  # But summarized to TFs instead of putting the peaks found
  # for every single experiment in the ChIP database
  # Verify first if there are results (more than 2 columns)
  dir.create(file.path(outpath, "TF_summarize"))
  if(dim(chip_genes)[2] > 2){
    chip_long <- pivot_longer(chip_genes, cols=c(-ENTREZID, -SYMBOL), names_to = "TF") %>%
      group_by(ENTREZID, SYMBOL, TF) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      mutate(value=ifelse(.$value==0, 0, 1)) %>%
      dplyr::rename(peak_found = value)
    write_tsv(chip_long, file.path(outpath, "TF_summarize", paste(title, "IFN_genes_TF.tsv", sep="_")))
  }
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  # Create outpaths if they don't exist
  dir.create(outpath)
  # Load metadata
  data("MetaData", package = "TFEA.ChIP")
  # Load Chip Data
  data("Mat01")
  for(DE_file in list.files(list_DE)){
    print(DE_file)
    if(!grepl("FULL", DE_file, fixed = TRUE)){
      if(grepl("set5", DE_file, fixed = TRUE)){trad <- TRUE}else{trad <- FALSE}
      main_TF_enrichment(file.path(list_DE, DE_file), outpath, ISG, list_TFs, trad)
    }
  }
}
