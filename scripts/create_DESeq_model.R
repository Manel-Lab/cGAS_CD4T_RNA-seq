## ---------------------------
##
## Purpose of script: Create a DESeq model
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-02-19
##
## ---------------------------

if (sys.nframe() == 0){
  ## Arguments
  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 7) {
    print("Not enough arguments provided")
    print("Usage : Rscript csvRaw csvTpm csvAnnot tsvMetadata fileISG outpath")
    q()
  } else if (length(args) > 7) {
    print("Too many arguments provided")
    print("Usage : Rscript csvRaw csvTpm csvAnnot tsvMetadata fileISG outpath")  
    q()
  }
  
  csvRaw <- args[1] # csv of the raw counts table
  csvTpm <- args[2] # csv of the TPM count table
  csvAnnot <- args[3] # gene annotation csv
  tsvMetadata <- args[4] # tsv of the samples metadata
  dir_ISGs <- args[5] # txt file with the ISG genes, one per line
  outpath <- args[6] # output path
  filter <- as.numeric(args[7]) # filter of TPM
}

## ---------------------------

## Packages

source("utils/DESeq_functions.R")
library(biomaRt)

## ------ FUNCTIONS ------- ##

# Import gene annotation table
# csv = the csv annotation table
import_annotations <- function(csv, dir_ISGs, exclude = "None") {
  annot <- read_csv(csv) %>%
    dplyr::select(-X1) %>%
    dplyr::rename(Gene = gene_id, Symbol = gene_name)
  
  # Get gene IDs and HUGO gene description from ENSEMBL API
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id", values = annot$Gene, mart = ensembl
  )
  mapping <- as_tibble(mapping)
  colnames(mapping) <- c("Gene", "Name", "Description")
  
  # Add IDs and Gene description
  annot <- annot %>%
    left_join(mapping) %>%
    mutate(Symbol = ifelse(is.na(Name), Symbol, Name)) %>%
    dplyr::select(-Name)
  # Add a column indicating if gene has haplotypes
  # (if we find the same gene more than once in the table)
  annot <- annot %>%
    dplyr::mutate(Haplotype = ifelse(Symbol %in% annot$Symbol[duplicated(annot$Symbol)], "TRUE", "FALSE")) %>%
    mutate(Symbol = ifelse(Haplotype, paste0(Symbol, " (", Gene, ")"), Symbol))
  
  # Prefill ISG column with "FALSE"
  annot$ISG <- "FALSE"
  # Add ISG annotation
  list_names <- c()
  for(file in list.files(dir_ISGs)){
    # Get name of the ISG list
    name <- regmatches(file, regexpr("^.*(?=_)", file, perl = TRUE))
    list_names <- append(list_names, name)
    # Read ISG file
    isg <- read_tsv(file.path(dir_ISGs, file))
    # Add column indicating wether it is an ISG according to this list
    annot <- mutate(annot, 
                    !!name := ifelse(Symbol %in% isg$official, "TRUE", "FALSE"))
    # And if it is, mark them as "TRUE" in ISG column
    # On the condition it's not from a particular ISG list
    # that we want to "exclude"
    if(name != exclude){
      annot <- mutate(annot, ISG = ifelse(Symbol %in% isg$official, "TRUE", ISG))
    }
  }
  
  return(annot)
}

# Creates the DESeq model and save it as a R object
# csvRaw = csv of the raw counts
# csvTpm = csv of the TPM counts
# csvAnnot = csv of gene annotations
# tsvMetadata = tsv of sample metadata
# fileISG = file with ISG genes
# outpath = output path where objects are saved
# filter = at least one sample must have this much TPM value for a gene to not be filtered out
main_create_model <- function(csvRaw, csvTpm, csvAnnot, 
                              tsvMetadata, dir_ISGs, outpath, filter){
   
  # Create the outpath in case it doesn't exist
  dir.create(outpath)
  # Write md5sum hash of files used
  md5ash <- unlist(lapply(c(csvRaw, csvTpm, csvAnnot, tsvMetadata), function(x) tools::md5sum(x)))
  write.table(md5ash, file=file.path(outpath, "create_DESeq_model_arguments.csv"), 
              sep=",", col.names=c("md5")
  )

  # Import counts data
  rawC <- import_counts(csvRaw)
  tpmC <- import_counts(csvTpm)
  # Import annotations
  annot <- import_annotations(csvAnnot, dir_ISGs, exclude = "MX1")
  # Import metadata
  metadata <- read_tsv(tsvMetadata)
  
  # Rename samples
  colnames(rawC) <- sub("D286-D281", "", colnames(rawC))
  colnames(tpmC) <- sub("D286-D281", "", colnames(tpmC))
  
  # Filter genes with low counts
  tpmC_f <- tpmC %>%
    filter_if(is.numeric, any_vars(. > filter))
  rawC_f <- rawC %>% 
    filter(Gene %in% tpmC_f$Gene)
  
  # Write a count tables with annotations to tsv format
  tsv_counts(rawC, annot, file.path(outpath, "Counts_raw_ALL.tsv"))
  tsv_counts(tpmC, annot, file.path(outpath, "Counts_tpm_ALL.tsv"))
  tsv_counts(rawC_f, annot, file.path(outpath, "Counts_raw_filer.tsv"))
  tsv_counts(tpmC_f, annot, file.path(outpath, "Counts_tpm_filter.tsv"))
  
  # Create a long version of the raw counts table
  cnts_long <- pivot_longer(rawC_f, cols=-Gene, names_to="Sample", values_to="Counts") %>%
    left_join(metadata)
  cnts_long$Donor <- as.factor(cnts_long$Donor)
  cnts_long$Condition <- as.factor(cnts_long$Condition)
  x <- sample(log(cnts_long$Counts+1), 5000)
  
  # Plot the distribution 
  distPlot <- ggplot(cnts_long, aes(Counts, color = Donor)) + 
    geom_freqpoly(binwidth = 1) + 
    scale_x_continuous(trans='log', labels = function(x) round(x)) +
    facet_grid(~Condition) +
    labs(title = "Log transformed gene distribution")
  
  # QQplot to see if the counts are log normal
  qqPlot <- ggqqplot(x) + labs(title="Log transformed gene distribution vs normal distribution")
  
  # Save to pdf file
  pdf(file.path(outpath, "Distribution.pdf"))
  print(distPlot)
  print(qqPlot)
  dev.off()
  
  # Prepare metadata
  meta_df <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df$Sample
  meta_df <- meta_df[-1]
  matrix_cnt <- as.matrix(rawC_f[-1])
  rownames(matrix_cnt) <- rawC_f$Gene
  # Reorder colnames
  matrix_cnt <- matrix_cnt[,rownames(meta_df)]
  # Convert variables to factor
  meta_df$Condition <- as.factor(meta_df$Condition)
  meta_df$Donor <- as.factor(meta_df$Donor)
  meta_df$TSA <- as.factor(meta_df$TSA)
  meta_df$cGAMP <- as.factor(meta_df$cGAMP)
  # Set reference level
  meta_df$TSA <- relevel(meta_df$TSA, "untreated")
  meta_df$cGAMP <- relevel(meta_df$cGAMP, "buffer")
  
  # Create the DESeq model
  dds <- DESeqDataSetFromMatrix(countData = matrix_cnt,
                                colData = meta_df,
                                design = ~ Donor + cGAMP + TSA + cGAMP:TSA)
  ddsObjet <- estimateSizeFactors(dds)
  ddsEstim <- DESeq(ddsObjet, parallel=FALSE)
  
  # Add TPM assay before saving DESeqDataSet object
  tpm <- as.matrix(tpmC_f[-1])
  rownames(tpm) <- tpmC_f$Gene
  ddsEstim@assays@data[["tpm"]] <- tpm
  
  # Save model and custom annotations
  save(ddsEstim, annot, list = c("ddsEstim", "annot"), file=file.path(outpath, "DESeqModel.RData"))
  
  # Vsd transform the data
  vsd <- varianceStabilizingTransformation(ddsEstim, blind = TRUE,
                                           fitType = "parametric")
  
  # PCA
  colData(vsd)$cGAMP <- recode(colData(vsd)$cGAMP, treated="cGAMP")
  colData(vsd)$TSA <- recode(colData(vsd)$TSA, treated="TSA")
  PCAplot <- DESeq2::plotPCA(vsd, intgroup=c("cGAMP", "TSA"), ntop=1000, returnData=FALSE)
  
  # Save to pdf file
  pdf(file.path(outpath, "Count_model.pdf"))
  plotDispEsts(ddsEstim, main="Estimated dispersion of the data")
  plotMA(ddsEstim, alpha = 0.01, ylim=c(-2,2), main="log2FC scatter plot")
  print(PCAplot)
  dev.off()
}

## --------- RUN ---------- ##

if (sys.nframe() == 0){
  main_create_model(csvRaw, csvTpm, csvAnnot, tsvMetadata, dir_ISGs, outpath, filter)
}
