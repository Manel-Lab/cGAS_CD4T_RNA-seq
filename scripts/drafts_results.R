## ---------------------------
##
## Purpose of script: Creates the figures for the "drafts" of the results
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-08-10
##
## ---------------------------

## Arguments

if (sys.nframe() == 0){
  ## Arguments
  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  DESeqres <- args[1]
  jfile <- args[2]
  outpath <- args[3]
  TPM_file <- args[4]
  
  # Read from json files
  figures <- rjson::fromJSON(file = jfile)
}

## ---------------------------

## Packages

library(readr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
source("utils/DESeq_functions.R")
source("enrichment_analysis.R")
source("create_DESeq_model.R")

## ------ FUNCTIONS ------- ##

#TODO:Comment
PCA <- function(ddsEstim, returnData = TRUE) {
  # Vsd transform the data
  vsd <- varianceStabilizingTransformation(ddsEstim, blind = TRUE,
                                           fitType = "parametric")
  
  # PCA
  colData(vsd)$cGAMP <- recode(colData(vsd)$cGAMP, treated="cGAMP")
  colData(vsd)$TSA <- recode(colData(vsd)$TSA, treated="TSA")
  PCAplot <- DESeq2::plotPCA(vsd, intgroup=c("cGAMP", "TSA"), ntop=1000, returnData=returnData)
  return(PCAplot)
}

# Add a column with the logFC in a set for a gene
# scatter_data = tibble
# set = tibble of a DE table from DESeq2
# logName = name of the new column
add_logFC <- function(scatter_data, Set, logName) {
  scatter_data <- left_join(scatter_data, Set) %>% 
    replace_na(list(log2FoldChange = 0)) %>%
    dplyr::rename(!!logName := log2FoldChange)
  return(scatter_data)
}


# Apply hierarchical clustering on genes of a matrix
# mat = expression matrix
# annotation_cluster = dataframe, cluster annotation
cluster_genes_HM <- function(mat, annotation_gene) {
  # Cluster cells INSIDE function clusters (Uninfected, 0, 1, 2)
  mat_ISG <- t(mat[,colnames(mat) %in% annotation_gene[annotation_gene$ISG==TRUE,]$Symbol])
  dISG <- dist(mat_ISG)
  clustISG <- hclust(dISG)
  mat_nonISG <- t(mat[,colnames(mat) %in% annotation_gene[annotation_gene$ISG==FALSE,]$Symbol])
  dnonISG <- dist(mat_nonISG)
  clustnonISG <- hclust(dnonISG)
  # Take order of these clusters and create one vector for ordering the matrix
  odISG <- clustISG$order
  odnonISG <- clustnonISG$order
  od = c(rownames(mat_nonISG[odnonISG,]), 
         rownames(mat_ISG[odISG,]))
  return(od)
}

# Create an expression heatmap
# expression = expression matrix
# genelist = vector of genes to limit the matrix to
# annot = gene annotation tibble
# gene_interest = vector of genes to display the names of 
# order = logical, wether to order by ISG or not
# clustering = wether to cluster genes, separately for ISG and non ISG
# .. = additional arguments for other functions
HM_draft <- function(expression, genelist, annot, gene_interest, order = TRUE, clustering = TRUE,
                    ...) {
  # Add gene symbol annotations
  tab <- expression %>% 
    filter(Gene %in% genelist) %>%
    left_join(annot[c("Gene", "Symbol")])
  # Create the matrix for the heatmap
  ncols <- dim(tab)[2]
  mat <- as.matrix(tab[c(-1, -ncols)])
  rownames(mat) <- tab$Symbol
  mat <- t(mat)
  # Scale each feature
  mat <- scale(mat)
  
  # Give a particular order of samples
  mat <- mat[c("T01", "T05", "T09", "T13",
               "T02", "T06", "T10", "T14"),
             ]
  # Give a particular order of genes
  if(clustering){
    order_genes <- cluster_genes_HM(mat, annot)
    order_genes <-  order_genes[order_genes %in% colnames(mat)]
    mat <- mat[,order_genes]
  }else{
    order_genes <- arrange(annot, ISG)$Symbol
    order_genes <-  order_genes[order_genes %in% colnames(mat)]
    mat <- mat[,order_genes]
  }
  # Change rownames to the condition
  rownames(mat) <- c(rep("Control", 4), rep("cGAMP", 4))
  # Different row labels
  custom_label <- colnames(mat)
  custom_label[!(custom_label %in% gene_interest)] = ""
  # Heatmap
  heatM <- pheatmap::pheatmap(t(mat), labels_row = custom_label,
                              ...)
  return(heatM)
}

# Create a scatter plot of LogFC in two different sets
# figures = json of figure annotation
# annot = gene annotation tibble
# outpath = output path 
# ylim = limits of y axis
create_scatter_plot <- function(list_res, figures, annot, outpath, ylim = c(-10, 10)) {
  # Import DE results for scatter plot
  idx <- figures$Scatter_plot$indexes
  title <- figures$Scatter_plot$title
  filename <- figures$Scatter_plot$filename
  # Import DE results for scatter plot
  Res1 <- add_column(as_tibble(list_res[[idx[1]]]), Gene = rownames(list_res[[idx[1]]])) %>%
    filter(padj < 5e-2, abs(log2FoldChange) > 1)
  Res2 <- add_column(as_tibble(list_res[[idx[2]]]), Gene = rownames(list_res[[idx[2]]])) %>%
    filter(padj < 5e-2, abs(log2FoldChange) > 1)
  
  # Keep relevant columns for the scatter_data
  Set1 <- dplyr::select(Res1, Gene, log2FoldChange)
  Set2 <- dplyr::select(Res2, Gene, log2FoldChange)
  
  # Create the tibble for the scatter plot
  # By adding the logFC and the set name
  scatter_data <- add_logFC(annot, Set1,  "log2FoldChange_1")
  scatter_data <- add_logFC(scatter_data, Set2, "log2FoldChange_2")
  
  # Add annotation of ISG (for color)
  scatter_data <- mutate(scatter_data, 
                         Annotation = case_when(ISG == FALSE ~ "Not an ISG",
                                                BAYLOR == TRUE & Shoggins == FALSE & MX1 == FALSE ~ "Baylor",
                                                BAYLOR == FALSE & Shoggins == TRUE & MX1 == FALSE ~ "Shoggins",
                                                BAYLOR == TRUE & Shoggins == TRUE & MX1 == FALSE ~ "Baylor & Shoggins")
  )
  # Change the name of ISG column
  scatter_data <- mutate(scatter_data, ISG = ifelse(ISG == TRUE, "ISG", "Not an ISG"))
  scatter_data <- mutate(scatter_data, ISG = ifelse(startsWith(Symbol, "IFN"), "IFN", ISG))
  print(scatter_data)

  
  # Create sub tibble for the names on the scatter plot
  subName <- scatter_data %>% filter(Symbol %in% c("IFNL1", "IFNL2", "IFNL3", "IFNB1", "IFNW1", "IRF7", "IFNGR2"))
  
  # Create sub tibble for the names on the scatter plot
  l1 <- append(arrange(subName, log2FoldChange_1)[1:5,]$Symbol, arrange(subName, desc(log2FoldChange_1))[1:5,]$Symbol)
  l2 <- append(arrange(subName, log2FoldChange_2)[1:5,]$Symbol, arrange(subName, desc(log2FoldChange_2))[1:5,]$Symbol)
  list_DEs <- append(l1, l2)
  subName <- subName %>% filter(Symbol %in% list_DEs)
  # Personal colors
  colors = c("Not an ISG" = "#999999",
             "ISG" = "#0000FF",
             "IFN" = "#FF0000")
  
  # LogFC scatter plot
  logFC_scatter <- ggplot(scatter_data,  aes(log2FoldChange_2, log2FoldChange_1, color=ISG)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = colors) +
    geom_text_repel(data = subName, mapping=aes(label=Symbol), nudge_x=-1, nudge_y=-3, direction = "x") +
    theme_bw() +
    labs(title = title) +
    xlab("Log2FoldChange, +cGAMP-TSA vs -cGAMP-TSA") + ylab("log2FoldChange, +cGAMP+TSA vs +cGAMP-TSA") + ylim(c(-8,8)) +
    theme(text = element_text(size=10))
  
  # Write to PDF
  pdf(file.path(outpath, filename), width=8, height = 4, compress = TRUE)
  print(logFC_scatter)
  dev.off()
  
  # Write table
  write_tsv(scatter_data, file.path(outpath, gsub(".pdf", "_raw_data.tsv", filename)))
}


# Main function to create expression heatmaps
# tpmC_f = tibble of tpm counts
# metadata = tibble of samples metadata
# annot = tibble of gene annotations
# figures = json of figures args
create_draft_heatmaps <- function(tpmC_f, metadata, annot, figures,
                                  only_ISG = FALSE, height = 15, width = 7, 
                                  ...) {
  if(only_ISG){
    annot <- filter(annot, (BAYLOR == TRUE | Shoggins == TRUE))
  }
  
  # Annotation for samples
  sample_annot <- as.data.frame(dplyr::select(metadata, TSA, cGAMP))
  rownames(sample_annot) <- metadata$Sample
  # Colors for annotation
  sample_color <- list(
    # Seurat colors
    TSA = setNames(c("#CC79A7", "#56B4E9"), levels(sample_annot$TSA)),
    cGAMP = setNames(c("#009E73", "#D55E00"), levels(sample_annot$cGAMP))
  )
  
  if(!only_ISG){
    # Annotation for genes
    genes_annot <- as.data.frame(annot[c("ISG")])
    rownames(genes_annot) <- annot$Symbol
    genes_annot$ISG <- as.factor(genes_annot$ISG)
    # Colors for annotation
    gene_color <- list(
      ISG = setNames(c("#CC79A7", "#56B4E9"), levels(genes_annot$ISG))
    )
  } else{
    genes_annot <- NULL
    genes_color <- NULL
  }

  
  # Colors for gene expression heatmaps
  colors <- rev(RColorBrewer::brewer.pal(8, "RdYlBu"))
  
  # Arguments from the json file
  list_DE_file <- figures$list_DE_file
  nplots <- figures$number_of_plots
  n <- figures$n
  filenames <- figures$filenames
  titles <- figures$titles
  logic <- figures$logic
  if(only_ISG){
    list_genes <- annot$Symbol
  } else{
    list_genes <- figures$list_genes
  }
  
  
  # Vector that will contain the list of DE
  list_DE <- c()
  # Gaps to separate the conditions
  gaps=c(4,8)
  # Breaks for the colors, to be symettrical around 0
  myBreaks <- seq(-2,2,by=0.5)
  
  # f = counter for the list of titles and filenames
  f <- 1
  # Iterate over the list of files but increment by the number of plots to produce
  for(i in seq(1, length(list_DE_file), nplots)){
    # For each file of DE we want to plot together
    files <- list_DE_file[i:(i+nplots-1)]
    for(j in seq_along(files)){
      # Extract the list of differential genes
      set <- read_tsv(files[j]) %>% arrange(desc(abs(log2FoldChange)))
      logic_i <- logic[[f]]
      if(logic_i == "intersect" && j > 1){
        list_sets <- lapply(files, function(x) arrange(read_tsv(x), desc(abs(log2FoldChange))))
        list_columns <- lapply(list_sets, function(x) dplyr::select(x, Gene)$Gene)
        print(list_columns)
        list_DE <- Reduce(intersect, list_columns)
      } else{
        list_DE <- append(list_DE, na.omit(set[1:n,])$Gene)
      }
    }
    list_DE <- unique(list_DE)
    # Create the heatmap and write it to file
    my_HM <- HM_draft(tpmC_f, list_DE, annot, list_genes, breaks=myBreaks, color=colors, 
                      main = titles[f], annotation_row = genes_annot, annotation_colors = gene_color, gaps_col=gaps,
                      ...)
    pdf(file.path(outpath, filenames[f]), height=height, width=width)
    print(my_HM)
    dev.off()
    # Increment by one the counter for filenames
    f <- f+1
    # Empty the list of DE to start anew
    list_DE <- c()
  }
}

# Main function to create pathways plots (for now, just barplot)
# outpath = output path
# annot = gene annotation tibble
# figures = json of args
create_pathway_plots <- function(outpath, annot, figures, n = 20) {
  list_pathways <- figures$Pathway_plots$list_pathways
  list_DE_file <- figures$Pathway_plots$list_DE_file
  titles <- figures$Pathway_plots$titles
  filenames <- figures$Pathway_plots$filenames
  
  # Annotation for genes
  genes_annot <- as.data.frame(annot[c("ISG")])
  rownames(genes_annot) <- annot$Symbol
  genes_annot$ISG <- as.factor(genes_annot$ISG)
  # Colors for annotation
  gene_color <- list(
    # Seurat colors
    ISG = setNames(c("#CC79A7", "#56B4E9"), levels(genes_annot$ISG))
  )
  
  # Iterate over the list of files but increment by the number of plots to produce
  for(i in seq_along(list_pathways)){
    # Extract the list of differential genes
    set <- read_tsv(list_DE_file[i])
    print(set$symbol)
    pathways <- get(load(list_pathways[i]))
    # Positive is true if there are overexpressed genes
    positive <- any(set$log2FoldChange > 0)
    # Negative is true if there are underexpessed genes
    negative <- any(set$log2FoldChange < 0)
    # Create the heatmap and write it to file
    if(dim(as.data.frame(pathways))[1] > 0){
      print(as.data.frame(pathways)$geneID)
      #heatP_path <- myHeatplot(pathways, set, positive, negative, "top_pathways", 
      #                         annotation_col = genes_annot, annotation_colors = gene_color)
      
      # Create barplots of the pathways and their p.adjust values
      bar_data <- as_tibble(pathways[,c("Description", "p.adjust")]) 
      #bar_data <- bar_data %>% filter(Description %in% colnames(heatP_path[[2]]))
      bar_data$Description <- factor(bar_data$Description, levels=rev(unique(bar_data$Description)))
      # Replace long pathways in bar_data
      bar_data <- arrange(bar_data, p.adjust)[1:n,]
      bar_data <- na.omit(bar_data)
      bar_data$Description <- replace(bar_data$Description, 
                                      bar_data$Description=="DNA damage response, signal transduction by p53 class mediator resulting in transcription of p21 class mediator",
                                      c("Transduction by p53 class mediator resulting in transcription of p21 class mediator"))
      
      barplot_pathway <- ggplot(data=bar_data, mapping=aes(x=Description, y=-log10(`p.adjust`))) +
        geom_bar(stat="identity", width=0.5, fill="#A8E0DF", color="black")  + 
        coord_flip() +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(size = 13),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        labs(title = titles[i])
      
      pdf(file.path(outpath, filenames[i]), height=10, width=10)
      #print(heatP_path[[1]])
      print(barplot_pathway)
      dev.off()
      
      write_tsv(bar_data, file.path(outpath, gsub(".pdf", "_Data.tsv", filenames[i])))
    }

  }
}

# Main function to create venn diagrams of three sets
# outpath = output path
# figures = json of figures args
create_draft_venn <- function(outpath, figures) {
  list_DE_file <- figures$Venn_DE_sets$list_DE_file
  filenames <- figures$Venn_DE_sets$filenames
  categories <- figures$Venn_DE_sets$categories
  
  # f = counter for the list of titles and filenames
  f <- 1
  # Iterate over the list of files but increment by the number of plots to produce
  for(i in seq(1, length(list_DE_file), 3)){
    set1 <- read_tsv(list_DE_file[i])$Gene
    set2 <- read_tsv(list_DE_file[i+1])$Gene
    set3 <- read_tsv(list_DE_file[i+2])$Gene
    vennplot <- draw.triple.venn(length(set1), length(set2), length(set3),
                                 length(intersect(set1, set2)), length(intersect(set2, set3)),
                                 length(intersect(set1, set3)), length(intersect(intersect(set2, set3), set1)),
                                 category = categories[[f]],
                                 fill = c("pink", "orange", "blue"),
                                 lty="blank",
                                 cat.cex = rep(1, 3),
                                 cat.col = c("pink", "orange", "blue"),
                                 print.mode = c("raw", "percent"))
    pdf(file.path(outpath, filenames[f]), height = 10, width = 10)
    grid.draw(vennplot)
    dev.off()
    # Increment by one the counter for filenames
    f <- f+1
  }
}




# Create a tibble for the boxplots of IFN
# tpm = tibble of TPM values
# annot = gene annotations tibble
# metadata = samples metadata
tibble_ifn <- function(tpm, annot, metadata,
                       list_genes = c("IFNB1", "IFNW1", "IFNG", "IFNG−AS1", "IFNL1", 
                                      "IFNL2", "IFNL3", "IFNAR1",
                                      "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1")) {
  tpm <- left_join(tpm, annot[c("Gene", "Symbol")])
  tpm <- tpm %>% 
    # Keep relevant genes
    filter(Symbol %in% list_genes) %>%
    # Change to long format for boxplots
    pivot_longer(c(-Gene, -Symbol), names_to = "Sample", values_to = "TPM") %>% 
    left_join(metadata)
  # Change some names for the boxplots
  tpm <- tpm %>% mutate(cGAMP = ifelse(cGAMP == "buffer", "buffer", "cGAMP"))
  # Factorize and set levels to have a particular order for the boxplot
  tpm$cGAMP <- factor(tpm$cGAMP, levels = c("buffer", "cGAMP"))
  tpm$TSA <- factor(tpm$TSA, levels = c("untreated", "treated"))
  tpm$Symbol <- factor(tpm$Symbol, levels = list_genes)
  tpm <- mutate(tpm, Condition = case_when(Condition == "BUF" ~ "Control",
                                           Condition == "cGAMP" ~ "cGAMP",
                                           Condition == "BUF-TSA" ~ "TSA",
                                           Condition == "cGAMP-TSA" ~ "cGAMP + TSA")
  )
  tpm$Condition <- factor(tpm$Condition, levels = c("Control", "TSA", "cGAMP", "cGAMP + TSA"))
  return(tpm)
}


# Main function for expression boxplots
# tpm = from tibble_ifn() function
# outpath = output path
create_boxplot_ifn <- function(tpm, outpath, 
                               title = "log2(TPM + 1) value of interferon genes",
                               filename = "Boxplot_IFN.pdf"){
  p <- ggplot(tpm, aes(x = Condition, y = log2(TPM + 1), fill = Condition)) +
    # Boxplot are brought closer together
    geom_boxplot(aes(x = Condition, y = log2(TPM + 1), fill = Condition), outlier.shape = 20, width = 0.9, position = position_dodge2(width = 0.9)) +
    # Make a grid for each gene, to have all genes one after another
    facet_grid(~ Symbol, space="free_x", scales="free_x", switch="x") +
    scale_colour_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    xlab("") +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.border = element_rect(colour="grey70"),
          panel.spacing=unit(0,"cm"),
          text = element_text(size=14),
          legend.position = 'bottom',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    ggtitle(title)

  pdf(file.path(outpath, filename), height = 7, width = 13)
  print(p)
  dev.off()
}

# Create tibble for plotting
# df = empty tibble
# filepaths = name of the comparisons
# orderspath = order of sorting for plotting
TF_tibble <- function(df, filepaths, orderspath){
  for(i in seq_along(filepaths)){
    # Read the GSEA ranking result of TF
    tab <- read_tsv(filepaths[i])
    tab <- arrange(tab, ES)
    # Give a set name for identification
    tab$set <- names(filepaths[i])
    # Give the order of plotting
    tab$order <- names(orderspath[i])
    # Bind with previous df (or empty df if first iteration)
    df <- bind_rows(df, tab)
  }
  
  # Keep pvalue < 0.05 results
  df <- df[df$pVal <= 0.05,]
  df <- na.omit(df)
  colnames(df) <- c("TF", "ES", "pvalue", "nChip", "set", "order")
  df$order <- as.numeric(df$order)
  
  # Sort by Condition and ES for better vizualisation
  df <- df %>%
    group_by(order) %>%
    arrange(order, desc(ES)) %>%
    ungroup() %>%
    mutate(`Transcription factor` = factor(TF, levels=unique(TF))) %>%
    mutate(Condition = factor(set, levels=unique(set)))
  return(df)
}

# Dotplot TF graph function
# df = tibble with data
# title = name of the graph
recap_TF <- function(df, title){
  g <- ggplot(df, aes(factor(Condition, levels=rev(levels(Condition))), `Transcription factor`)) + 
    geom_point(aes(size=-pvalue, color=ES)) +
    scale_color_gradient2(low="blue", mid="white", high="red") +
    scale_size_continuous(range = c(0.5, 3.5)) +
    labs(title = "Enriched transcription factors") +
    ylab("Transcription factors") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
}

# Create heatmap of associated TF
# tpmC_f = tpm tibble
# metadata = samples metadata
# annot = annotation of genes
# figures = json of args
# outpath = output path
HM_associated_TF <- function(tpmC_f, metadata, annot, figures, outpath) {
  
  # Arguments from json
  title <- figures$HM_associated$title
  filename <- figures$HM_associated$filename
  file_int <- figures$HM_associated$file_int
  
  genes_int <- read_tsv(file_int)
  genes <- annot %>% filter(Symbol %in% genes_int$SYMBOL)
  tpmC_f <- filter(tpmC_f, Gene %in% genes$Gene)
  newtpm <- tpmC_f %>% left_join(genes[c("Gene", "Symbol")])
  mat <- as.matrix(newtpm[c(-1,-dim(newtpm)[2])])
  rownames(mat) <- newtpm$Symbol
  colors <- rev(RColorBrewer::brewer.pal(8, "RdYlBu"))
  # Annotation for samples
  row_annot <- metadata[c("TSA","cGAMP")]
  # Colors for annotation
  ann_color <- list(
    # Seurat colors
    TSA = setNames(c("#CC79A7", "#56B4E9"), levels(row_annot$TSA)),
    cGAMP = setNames(c("#009E73", "#D55E00"), levels(row_annot$cGAMP))
  )
  myBreaks <- seq(-2,2,by=0.5)
  mat <- t(mat)
  mat <- scale(mat)
  HM <- pheatmap::pheatmap(mat, color = colors, annotation_rows = row_annot, annottion_colors = ann_color, 
                           breaks=myBreaks, main=title)
  pdf(file.path(outpath, filename), width=25, height=5)
  print(HM)
  dev.off()
}

# Main function for the script
# tpm_csv = file to TPM counts
# isg_file = ISG directory
# annot_tsv = gene annotations
# metadata_file = samples metadata
# filter = TPM filter
# outpath = output path
# figures = json of args
main_draft_figures <- function(DESeqres, outpath, figures, TPM_file) {
  #### Load data ####
  
  # Create outpath
  dir.create(outpath)
  
  load(DESeqres)
  
  tpmC <- read_csv(TPM_file)
  tpmC <- tpmC[1:17]
  colnames(tpmC)[1] <- "Gene"
  tpmC_f <- as_tibble(ddsEstim@assays@data[["tpm"]]) %>% 
    add_column(Gene = rownames(ddsEstim@assays@data[["tpm"]]), .before=1)
  metadata <- as_tibble(colData(ddsEstim)) %>%
    add_column(Sample = rownames(colData(ddsEstim)), .before=1) %>%
    dplyr::select(-Day, -sizeFactor)

  # Write expression tsv
  write_tsv(left_join(tpmC_f, annot), file.path(outpath, "TPM_expression_data.tsv"))
  
  # Write md5sum hash of files used
  md5ash <- unlist(lapply(c(DESeqres, jfile), function(x) tools::md5sum(x)))
  write.table(md5ash, file=file.path(outpath, "drafts_files.csv"), 
              sep=",", col.names=c("md5")
  )
  
  # Write args used to file
  my_args <- rjson::toJSON(figures)
  write(my_args, file=file.path(outpath, "drafts_args_used.json"))
  
  ###### PCA plot ######
  
  pdf(file.path(outpath, "PCA_plot.pdf"))
  print(PCA(ddsEstim, FALSE))
  dev.off()
  
  PCAdata <- PCA(ddsEstim, TRUE)
  write.table(PCAdata, file.path(outpath, "PCA_raw_data.tsv"), sep = "\t")
  
  
  ###### LogFC scatter plot #######
  make <- as.logical(figures$Scatter_plot$make)
  if(make){
    create_scatter_plot(list_res, figures, annot, outpath, c(-10, 10))
  }
  
  #### HM ####
  make <- as.logical(figures$Recap_HM$make)
  if(make){
    this_tpmC_f <- tpmC_f %>% dplyr::select(-T03, -T07, -T11, -T15, -T04,-T08,-T12,-T16)
    this_metadata <- metadata %>% dplyr::filter(!(Sample %in% c("T03", "T07", "T11", "T15", "T04","T08","T12","T16")))
    create_draft_heatmaps(this_tpmC_f, this_metadata, annot, figures$Recap_HM, clustering = FALSE,
                          cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, 
                          show_colnames = TRUE, angle_col = 315, 
                          cellheight = 0.3, cellwidth = 15,
                          width=5, height=4, fontsize_col=9)
  }
  
  #### HM ONLY ISG #### We'll have to do it differently
  make <- as.logical(figures$Recap_HM_ISG_only$make)
  if(make){
    this_tpmC_f <- tpmC_f %>% dplyr::select(-T03, -T07, -T11, -T15, -T04,-T08,-T12,-T16)
    this_metadata <- metadata %>% dplyr::filter(!(Sample %in% c("T03", "T07", "T11", "T15", "T04","T08","T12","T16")))
    create_draft_heatmaps(tpmC_f, metadata, annot, figures$Recap_HM_ISG_only, TRUE, 22, 9, clustering = FALSE,
                          cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, 
                          show_colnames = TRUE , angle_col = 315, 
                          cellheight = 9, cellwidth = 12)
  }
  
  ### Pathways ###
  make <- as.logical(figures$Pathway_plots$make)
  if(make){
    create_pathway_plots(outpath, annot, figures)
  }
  
  ## Venn Diagram
  make <- as.logical(figures$Venn_DE_sets$make)
  if(make){
    create_draft_venn(outpath, figures)
  }
  
  ## Boxplot IFN
  make <- as.logical(figures$Boxplot_IFN$make)
  if(make){
    tpm <- tibble_ifn(tpmC_f, annot, metadata)
    to_write <- pivot_wider(dplyr::select(tpm, -Gene), values_from=c("TPM"), names_from = c("Symbol"))
    write_tsv(to_write, file.path(outpath, "Data_Boxplot_IFN.tsv"))
    create_boxplot_ifn(tpm, outpath)
  }
  
  ## Boxplot IFN all IFN
  if(make){
    tpm <- tibble_ifn(tpmC, annot, metadata, 
                      list_genes = c("IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21"))
    to_write <- pivot_wider(dplyr::select(tpm, -Gene), values_from=c("TPM"), names_from = c("Symbol"))
    write_tsv(to_write, file.path(outpath, "Data_Boxplot_IFNA.tsv"))
    create_boxplot_ifn(tpm, outpath, filename = "Boxplot_IFNA.pdf")
  }
  
  if(make){
    tpm <- tibble_ifn(tpmC, annot, metadata, 
                      list_genes = c("IFNB1", "IFNE", "IFNK", "IFNW1", "IFNG", "IL10R2", "IFNL1", "IFNL2", "IFNL4", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1"))
    to_write <- pivot_wider(dplyr::select(tpm, -Gene), values_from=c("TPM"), names_from = c("Symbol"))
    write_tsv(to_write, file.path(outpath, "Data_Boxplot_IFN_non_IFNA.tsv"))
    create_boxplot_ifn(tpm, outpath, filename = "Boxplot_IFN_except_A.pdf")
  }
  
  ## Boxplot genes interest
  if(make){
    tpm<- tibble_ifn(tpmC_f, annot, metadata,
                     list_genes = c("IRF7", "RELA", "IRF3", "STING1", "CGAS", "TBK1"))
    to_write <- pivot_wider(dplyr::select(tpm, -Gene), values_from=c("TPM"), names_from = c("Symbol"))
    write_tsv(to_write, file.path(outpath, "Data_Boxplot_Genes_Interest.tsv"))
    create_boxplot_ifn(tpm, outpath, "log2(TPM + 1) value", "Boxplot_Genes_interest.pdf")
  }
  
  ## Recap of TF enriched
  make <- as.logical(figures$TF$make)
  if(make){
    title <- figures$TF$title
    filepaths <- figures$TF$filepaths
    filename <- figures$TF$filename
    names <- figures$TF$names 
    orders <- figures$TF$orders
    df <- tibble()
    orderspath <- filepaths
    # Give names for better handling later
    names(filepaths) <- names
    names(orderspath) <- orders
    df <- TF_tibble(df, filepaths, orderspath)
    recap <- recap_TF(df, title)
    pdf(file.path(outpath, filename), width=4, height=10)
    print(recap)
    dev.off()
  }
  
  # Heatmap of associated genes #
  make <- as.logical(figures$HM_associated$make)
  if(make){
    HM_associated_TF(tpmC_f, metadata, annot, figures, outpath) 
  }
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_draft_figures(DESeqres, outpath, figures, TPM_file)
}

