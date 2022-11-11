## ---------------------------
##
## Purpose of script: Pipeline to recreate our analysis
##
## Author: Kévin De Azevedo (kdeazevedo@github)
## Email: kevin.de-azevedo@curie.fr
##
## Date Created: 2020-07-10
##
## ---------------------------

## Arguments
.libPaths(c("~/R_packages", .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  print("Not enough arguments provided")
  print("Usage : Rscript csvRaw csvTpm csvAnnot tsvMetadata fileISG outpath previous_tsv jfile")
  q()
} else if (length(args) > 8) {
  print("Too many arguments provided")
  print("Usage : Rscript csvRaw csvTpm csvAnnot tsvMetadata fileISG outpath previous_tsv jfile")
  q()
}


csvRaw <- args[1]       # csv of the raw counts table
csvTpm <- args[2]       # csv of the TPM count table
csvAnnot <- args[3]     # gene annotation csv
tsvMetadata <- args[4]  # tsv of the samples metadata
fileISG <- args[5]      # txt file with the ISG genes, one per line
outpath <- args[6]      # output path
previous_tsv <- args[7] # Previous DE genes from other experiment
jfile <- args[8]        # Json argument file

## Global variables

# Read from json files
figures <- rjson::fromJSON(file = jfile)

## ---------------------------

## Packages

source("create_DESeq_model.R")
source("DE_analysis.R")
source("enrichment_analysis.R")
source("TF_enrichment.R")
source("drafts_results.R")

## ------ FUNCTIONS ------- ##

## --------- RUN ---------- ##

print("1) Creating model")
main_create_model(csvRaw, csvTpm, csvAnnot,
                  tsvMetadata, fileISG, outpath, 1)
print("2) Differential analysis")
main_DE_analysis(file.path(outpath, "DESeqModel.RData"), previous_tsv, outpath, 0.5)
print("3) Pathway over representation analysis")
main_enrichment(file.path(outpath, "DE_tables"), file.path(outpath, "with_ISG", "Pathway_analysis"), 
                0.5, 0.05, TRUE)
main_enrichment(file.path(outpath, "DE_tables"), file.path(outpath, "without_ISG", "Pathway_analysis"),  
                0.5, 0.05, FALSE)
print("4) Transcription factor enrichment analysis")
# Create outpaths if they don't exist
dir.create(file.path(outpath, "with_ISG", "TF_analysis"))
dir.create(file.path(outpath, "without_ISG", "TF_analysis"))
# Load metadata
data("MetaData", package = "TFEA.ChIP")
# Load Chip Data
data("Mat01")
for(DE_file in list.files(file.path(outpath, "DE_tables"))){
  print(DE_file)
  if(!grepl("FULL", DE_file, fixed = TRUE)){
    if(grepl("set5", DE_file, fixed = TRUE)){trad <- TRUE}else{trad <- FALSE}
    main_TF_enrichment(file.path(file.path(outpath, "DE_tables"), DE_file), 
                       file.path(outpath, "with_ISG", "TF_analysis"), TRUE, c("RELA", "IRF3"), trad)
    main_TF_enrichment(file.path(file.path(outpath, "DE_tables"), DE_file), 
                       file.path(outpath, "without_ISG", "TF_analysis"), FALSE, c("RELA", "IRF3"), trad)
  }
}
print("5) Figures")
#With ISG only
main_draft_figures(csvTpm, fileISG, csvAnnot, tsvMetadata,
                   1, file.path(outpath, "Draft_figures"), figures)
