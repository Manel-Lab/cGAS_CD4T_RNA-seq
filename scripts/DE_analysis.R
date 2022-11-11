## ---------------------------
##
## Purpose of script: Differential expression analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-03-03
##
## ---------------------------

if (sys.nframe() == 0){
  ## Arguments
  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 4) {
    print("Not enough arguments provided")
    print("Usage : Rscript DE_analysis Robj previous_tsv outpath logFC")
    q()
  } else if (length(args) > 4) {
    print("Too many arguments provided")
    print("Usage : Rscript DE_analysis Robj previous_tsv outpath logFC")  
    q()
  }
  
  Robj <- args[1]
  previous_tsv <- args[2]
  outpath <- args[3]
  logFC_t <- as.numeric(args[4])
  
  # Set seed
  set.seed(0)
}

## ---------------------------

## Packages

source("utils/DESeq_functions.R")
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

## ------ FUNCTIONS ------- ##

# Create an expression matrix from a DE result table
# set = result from DESeq2
# annot = gene annotation tibble
# ddsEstim = DESeq result object
# include_isg = logical, wether to filter out ISG
# n = top number of DE genes to display
create_expression_matrix <- function(set, annot, ddsEstim,
                                     include_isg = TRUE, n = 50) {
  # Sort by descreasing absolute FoldChange (highest FC top of the tibble)
  set <- set %>% arrange(dplyr::desc(abs(log2FoldChange)))
   if(!include_isg){
     # Filter out ISG genes
     set <- set %>% dplyr::filter(ISG == FALSE)
   }
   # Get n genes with highest logFC
   topDE <- set[1:n, "Gene"]$Gene
   # Get TPM values of the top n
   tab <- as.tibble(ddsEstim@assays@data[["tpm"]]) %>%
     add_column(Gene=rownames(ddsEstim@assays@data[["tpm"]]), .before=1) %>% 
     left_join(annot)
   tab <- tab[tab$Gene %in% topDE,]
   ncols <- dim(tab)[2]
   mat <- as.matrix(tab[c(-1, -(ncols-6):-ncols)])
   rownames(mat) <- tab$Symbol
   mat <- t(mat)
   # Scale each feature
   mat <- scale(mat)
   # Give a particular order
   mat <- mat[c("T01", "T05", "T09", "T13",
                "T02", "T06", "T10", "T14",
                "T03", "T07", "T11", "T15",
                "T04", "T08", "T12", "T16"),
              ]
   return(mat)
}

# Create all HM, ISG or no ISG
# list_of_sets = list of results from DE analysis
# ddsEstim = DESeq object
# outpath = output path for the HM
# annot = annotation table
# list_names = list of set names
# include_ISG = wether to filter out ISG (FALSE) or not (TRUE)
# n = top number of DE to display
# ... = additional arguments passed to the
create_heatmaps <- function(list_of_sets, ddsEstim, outpath, annot, 
                            list_names, include_isg = TRUE, n = 50,
                            ...) {
  for(i in seq_along(list_of_sets)){
    mat <- create_expression_matrix(list_of_sets[[i]], annot, ddsEstim,
                                    include_isg, n)
    heatM <- pheatmap::pheatmap(mat, ...)
    pdf(file.path(outpath, paste0(list_names[[i]], "_heatmap.pdf")), height=7, width=18)
    print(heatM)
    dev.off()
  }
}

# Create all HM, ISG or no ISG, vertical format
# list_of_sets = list of results from DE analysis
# ddsEstim = DESeq object
# outpath = output path for the HM
# annot = annotation table
# list_names = list of set names
# include_ISG = wether to filter out ISG (FALSE) or not (TRUE)
# n = top number of DE genes to display
# ... = additional argument passed to the heatmap_DEG() function
create_heatmaps2 <- function(list_of_sets, ddsEstim, outpath, annot, list_names,
                             include_isg = TRUE, n = 50,
                             ...) {
  for(i in seq_along(list_of_sets)){
    mat <- create_expression_matrix(list_of_sets[[i]], annot, ddsEstim,
                                    include_isg, n)
    rownames(mat) <- c(rep("Control", 4), rep("cGAMP", 4), rep("TSA", 4), rep("cGAMP+TSA", 4))
    heatM <- pheatmap::pheatmap(t(mat), ...)
    pdf(file.path(outpath, paste0(list_names[[i]], "_heatmap_vertical.pdf")), height=15, width=7)
    print(heatM)
    dev.off()
  }
}

# Create all direction heatmaps
# table_recap = recapitulating table of the results
# ddsEstim = DESeq object
# outpath = output path for the HM
# annot = annotation table
# include_ISG = wether to filter out ISG (FALSE) or not (TRUE)
# ... = additional arguments passed to heatmap_DEG() function
create_diff_direction_heatmaps <- function(table_recap, ddsEstim, outpath, annot, include_isg = TRUE, 
                                           ...) {
  if(!include_isg){
    # Filter out ISG
    table_recap <- table_recap %>% 
      filter(Symbol %in% annot[annot$ISG == FALSE,]$Symbol)
  }
  

  keep <- table_recap %>% dplyr::filter(`+cGAMP-TSA_+cGAMP+TSA` == 1, 
                                        `-cGAMP-TSA_-cGAMP+TSA` == 1, 
                                        `Diff_effect_TSA` == 1) %>%
    # Get genes under-expressed in Set 4
    dplyr::filter(LogFC_Diff_effect_TSA < 0) %>%
    # And expressed/repressed in only Set 3 or Set 1 (different directions)
    dplyr::mutate(Keep = ifelse(`LogFC_-cGAMP-TSA_-cGAMP+TSA`/`LogFC_+cGAMP-TSA_+cGAMP+TSA` < 0, "TRUE", "FALSE")) %>%
    dplyr::filter(Keep == TRUE) %>%
    left_join(annot)
  tab <- as.tibble(ddsEstim@assays@data[["tpm"]]) %>%
    add_column(Gene=rownames(ddsEstim@assays@data[["tpm"]]), .before=1) %>% 
    left_join(annot) %>%
    mutate(Symbol = ifelse(Haplotype, paste0(Symbol, " (", Gene, ")"), Symbol))
  tab <- tab[tab$Gene %in% keep$Gene,]
  ncols <- dim(tab)[2]
  mat <- as.matrix(tab[c(-1, -(ncols-6):-ncols)])
  rownames(mat) <- tab$Symbol
  mat <- t(mat)
  mat <- scale(mat)
  heatM <- heatmap_DEG(mat, ...)
  pdf(file.path(outpath, "Opposite_directions_genes.pdf"), height=7, width=18)
  print(heatM)
  dev.off()
}

# Create all recap heatmaps
# list_of_sets = list of results from DE analysis
# outpath = output path for the HM
# annot = annotation table
# include_ISG = wether to filter out ISG (FALSE) or not (TRUE)
# ... = additional arguments for the recap_heatmap() function
create_recaps <- function(list_of_sets, outpath, annot, namedlist, include_isg = TRUE,
                          ...){
  
  table_recap <- data_RecapHM(list_of_sets, outpath, namedlist) %>% 
    dplyr::filter(Diff_effect_TSA == 1)
  
  if(!include_isg){
    table_recap <- table_recap %>% 
      filter(Symbol %in% annot[annot$ISG == FALSE,]$Symbol)
  }
  
  
  # Keep the top 50 up-regulated
  recap <- table_recap  %>%
    arrange(dplyr::desc(LogFC_Diff_effect_TSA)) %>%
    dplyr::slice(1:50) 
  
  recap <- filter_table(recap, threshold = 0)
  mat <- as.matrix(recap[-1])
  rownames(mat) <- recap$Symbol
  mat <- t(mat)
  recapHM <- recap_heatmap(mat, ...)
  pdf(file.path(outpath, "LogFC_UP_Heatmap.pdf"), height=5, width=22)
  print(recapHM)
  dev.off()
  
  
  # Keep the top 50 down-regulated
  recap <- table_recap %>%
    arrange(LogFC_Diff_effect_TSA) %>%
    dplyr::slice(1:50) 
  
  recap <- filter_table(recap, threshold = 0)
  mat <- as.matrix(recap[-1])
  rownames(mat) <- recap$Symbol
  mat <- t(mat)
  recapHM <- recap_heatmap(mat, ...)
  pdf(file.path(outpath, "LogFC_DOWN_Heatmap.pdf"), height=5, width=22)
  print(recapHM)
  dev.off()
}

# Create Vienn Diagram comparing two lists of DE genes
# list1 = list of DE genes from first experiment
# list2 = list of DE genes from second experiment
# print = wether to save Venn Diagram or not
# filename = name of the file
# title = title of the plot
create_venn_diagram <- function(list1, list2, print, filename, title){
  venn.diagram(
    list("Experiment 1" = list1, "Experiment 2" = list2),
    filename = filename,
    imagetype="png" ,
    print.mode = print,
    main = title,
    main.cex = 0.3,
    height = 500, 
    width = 500, 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col=c("#440154ff", '#21908dff'),
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = c("#440154ff", '#21908dff')
  )
}

# Create a new Set of genes
# table1 = DE table
# table2 = DE table
# filename = name of the file
create_set5 <- function(table1, table2, filename) {
  set5 <- table1 %>% inner_join(table2, by=c("Gene", "Symbol", "Description", "Haplotype", "ISG"), suffix=c("_set1", "_set2"))
  trad <- as_tibble(bitr(set5$Symbol, fromType="SYMBOL", toType="ENTREZID", drop=FALSE, OrgDb="org.Hs.eg.db"))
  names(trad) = c("Symbol", "Entrez")
  set5 <- set5 %>% left_join(trad)
  write_tsv(set5, file.path(outpath,"DE_tables", filename))
  return(set5)
}

# Create the DE tables and plots (heatmaps, volcano plots etc.)
# Robj = R object of the DESeq model
# previous_tsv = DE table of a previous analysis
# outpath = output path
# logFC_t = logFC threshold (under this threshold, gene is not a DEG)
main_DE_analysis <- function(Robj, previous_tsv, outpath, logFC_t){
  # Create output path in case it doesn't exist
  dir.create(outpath)
  dir_diffgenes <- file.path(outpath, "DE_tables")
  dir.create(dir_diffgenes)
  # Load model
  load(Robj)
  # Create all the results
  list_of_sets <- list()
  # List of names of analysis
  list_names = list("+cGAMP-TSA_+cGAMP+TSA", "-cGAMP-TSA_+cGAMP-TSA",
                    "-cGAMP-TSA_-cGAMP+TSA", "Diff_effect_TSA")
  # 1) The TSA effect in the cGAMP condition
  resDESeq <- results(ddsEstim, contrast=list(c("TSA_treated_vs_untreated", "cGAMPtreated.TSAtreated")))
  # 2) What is the effect on genes of adding cGAMP ?
  resDESeq2 <- results(ddsEstim, contrast=c("cGAMP", "treated", "buffer"))
  # 3) What is the effect on the genes of adding TSA ?
  resDESeq3 <- results(ddsEstim, contrast=c("TSA", "treated", "untreated"))
  # 4) The difference of TSA effect between cGAMP/nocGAMP (interaction term)
  resDESeq4 <- results(ddsEstim, name="cGAMPtreated.TSAtreated")
  list_res <- list(resDESeq, resDESeq2, resDESeq3, resDESeq4)
  # Creating list of sets (DE genes with padj < 0.05)
  for(i in seq_along(list_res)){
    set <- get_set(list_res[[i]], annot, i, 5e-2, logFC_t)
    list_of_sets <- list_sets(list_of_sets, set, dir_diffgenes, i, 
                              list_names[[i]], logFC = logFC_t)
  }
  
  # Annotation for samples
  row_annot <- as.data.frame(colData(ddsEstim)[c("TSA","cGAMP")])
  # Colors for annotation
  ann_color <- list(
    # Seurat colors
    TSA = setNames(c("#CC79A7", "#56B4E9"), levels(row_annot$TSA)),
    cGAMP = setNames(c("#009E73", "#D55E00"), levels(row_annot$cGAMP))
  )
  # Colors for gene expression heatmaps
  colors <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))
  
  # Create outpaths if they do not exist
  outpath_isg <- file.path(outpath, "with_ISG")
  outpath_no_isg <- file.path(outpath, "without_ISG")
  dir.create(outpath_isg)
  dir.create(outpath_no_isg)
  
  # Heatmap of top 100 for each set
  myBreaks <- seq(-2.5,2.5,by=0.5)
  create_heatmaps(list_of_sets, ddsEstim, outpath_isg, annot, list_names, TRUE, 100,
                  color = colors, annotation_row = row_annot, annotation_colors = ann_color, breaks = myBreaks,
                  main = "Gene expression heatmaps of Top DE Genes. ISG genes included")
  create_heatmaps(list_of_sets, ddsEstim, outpath_no_isg, annot, list_names, FALSE, 100,
                  color = colors, annotation_row = row_annot, annotation_colors = ann_color, breaks = myBreaks,
                  main = "Gene expression heatmaps of Top DE Genes. ISG genes not included")
  
  # Vertical heatmaps
  create_heatmaps2(list_of_sets, ddsEstim, outpath_isg, annot, list_names, TRUE, 100,
                   color = colors, breaks = myBreaks,
                   main  =  "Gene expression heatmaps of Top DE Genes, Z-scores of TPM values.",
                   cluster_rows = TRUE, cluster_cols = FALSE, 
                   angle_col = 315, cellheight = 8, cellwidth = 13)
  create_heatmaps2(list_of_sets, ddsEstim, outpath_no_isg, annot, list_names, FALSE, 100,
                   color = colors, breaks = myBreaks,
                   main  =  "Gene expression heatmaps of Top DE Genes, Z-scores of TPM values.",
                   cluster_rows = TRUE, cluster_cols = FALSE, 
                   angle_col = 315, cellheight = 8, cellwidth = 13)
  
  
  # Volcano Plots
  create_volcanos(list_res, list_names, outpath_isg, annot, TRUE)
  create_volcanos(list_res, list_names, outpath_no_isg, annot, FALSE)
  
  # Recap heatmap : keep genes with a different TSA effect in cGAMP
  myBreaks <- c(-8, -6, -4, -2, -0.5, 0.5, 2, 4, 6, 8)
  # Tweak colors to have white color between -0.5 and 0.5 logFC
  colors_logFC <- rev(RColorBrewer::brewer.pal(11, "RdBu")[-1])[-1]
  create_recaps(list_of_sets, outpath_isg, annot, list_names, TRUE, 
                colors = colors_logFC, myBreaks = myBreaks)
  create_recaps(list_of_sets, outpath_no_isg, annot, list_names, FALSE, 
                colors = colors_logFC, myBreaks = myBreaks)
  
  # Compare with previous analysis
  previous <- read_tsv(previous_tsv,
                       col_names = c("ENSEMBL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
                       skip = 1) %>%
    separate(ENSEMBL, c("ENSEMBL", "Symbol"), "[\\|]")
  # Separate UP and DOWN
  previous <- previous %>% mutate(ENSEMBL = gsub("\\..*","",ENSEMBL))
  previous_DOWN <- previous %>% filter(logFC < (-1*logFC_t)) %>%
    filter(adj.P.Val < 0.05)
  previous_UP <- previous %>% filter(logFC > logFC_t) %>%
    filter(adj.P.Val < 0.05)
  # Same for the set we want to compare with
  set2_UP <- list_of_sets[[2]] %>% filter(log2FoldChange > logFC_t)
  set2_DOWN <- list_of_sets[[2]] %>% filter(log2FoldChange < (-1*logFC_t))
  
  # Get Common Genes
  common_DOWN <- set2_DOWN %>% filter(set2_DOWN$Gene %in% previous_DOWN$ENSEMBL)
  write_tsv(common_DOWN, file.path(outpath, "common_DOWN_0.05.tsv"))
  common_UP <- set2_UP %>% filter(set2_UP$Gene %in% previous_UP$ENSEMBL)
  write_tsv(common_UP, file.path(outpath, "common_UP_0.05.tsv"))
  # Venn Diagrams
  create_venn_diagram(previous_DOWN$Symbol, set2_DOWN$Symbol, "percent",
                      file.path(outpath, "venn_DOWN_0.05_RAW.png"), 
                      "DOWN-regulated genes, adj pval < 0.05")
  create_venn_diagram(previous_UP$Symbol, set2_UP$Symbol, "percent",
                      file.path(outpath, "venn_UP_0.05_RAW.png"), 
                      "UP-regulated genes, adj pval < 0.05")
  
  # Create Set 5
  up_set1 <- read_tsv(file.path(outpath, "DE_tables", "UP_+cGAMP-TSA_+cGAMP+TSA.tsv"))
  up_set2 <- read_tsv(file.path(outpath, "DE_tables", "UP_-cGAMP-TSA_+cGAMP-TSA.tsv"))
  down_set1 <- read_tsv(file.path(outpath, "DE_tables", "DOWN_+cGAMP-TSA_+cGAMP+TSA.tsv"))
  down_set2 <- read_tsv(file.path(outpath, "DE_tables", "DOWN_-cGAMP-TSA_+cGAMP-TSA.tsv"))
  
  up_set5 <- create_set5(down_set1, up_set2, "UP_set5.tsv")
  down_set5 <- create_set5(down_set2, up_set1, "DOWN_set5.tsv")
  
  # Add set 5 to list of results
  list_res[[5]] <- up_set5
  list_res[[6]] <- down_set5
  list_names[[5]] <- "UP_Set5"
  list_names[[6]] <- "DOWN_Set5"
  
  # Save DESeq model, annotations, and list of results
  save(ddsEstim, annot, list_res, list_names,
       list = c("ddsEstim", "annot", "list_res", "list_names"), 
       file=file.path(outpath, "DESeqResults.RData")
  )
  
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_DE_analysis(Robj, previous_tsv, outpath, logFC_t)
}

