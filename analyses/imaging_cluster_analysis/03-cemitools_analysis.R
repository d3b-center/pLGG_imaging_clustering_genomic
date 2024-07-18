# Author: Komal S. Rathi, Run Jin
# Perform CEMiTools analysis for pathways 
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(CEMiTool)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "imaging_cluster_analysis")

# source function
source(file.path(analysis_dir, "utils", "run_cemitools.R"))

# parse command line options
option_list <- list(
  make_option(c("--mat_file"), type = "character",
              help = "gene matrix data, preferably count matrix (.rds) "),
  make_option(c("--anno_file"), type = "character",
              help = "path to cluster annotation file"),
  make_option(c("--gtf_file"), type = "character",
              help = "path to gtf annotation file"),
  make_option(c("--msigdb"), type = "character", 
              help = "reactome, kegg or canonical"),
  make_option(c("--cor_method"), type = "character", default = "pearson",
              help = "correlation method"),
  make_option(c("--network_type"), type = "character", default = "unsigned",
              help = "network type"),
  make_option(c("--tom_type"), type = "character", default = "signed",
              help = "TOM type"),
  make_option(c("--plots_dir"), type = "character",
              help = "path to plots directory"),
  make_option(c("--results_dir"), type = "character",
              help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
msigdb <- opt$msigdb
anno_file <- opt$anno_file
plots_dir <- opt$plots_dir
results_dir <- opt$results_dir

# output directories
dir.create(plots_dir, showWarnings = F, recursive = T)
dir.create(results_dir, showWarnings = F, recursive = T)

# read gene count file
mat_matrix <- readRDS(opt$mat_file)

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = opt$gtf_file)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter gene matrix file to contain only protein coding gene
mat_matrix_coding <- mat_matrix[rownames(mat_matrix) %in% gencode_gtf$gene_name,]

# read in the cluster annotation file 
cluster_anno <- readr::read_tsv(anno_file) %>% 
  # also make rownames kids first biospecimen id
  dplyr::mutate(tmp = Kids_First_Biospecimen_ID) %>% 
  tibble::column_to_rownames("tmp")

################# run CEMiTools on each cluster 
# filter gene matrix to only contain samples from this disease
mat_of_interest <- mat_matrix_coding %>%
  dplyr::select(any_of(cluster_anno$Kids_First_Biospecimen_ID))
cluster_anno <- cluster_anno %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(mat_of_interest))

# create input gmt file on the fly
if(msigdb == "kegg"){
  gmt_file <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", "CP:KEGG")
} else if(msigdb == "reactome"){
  gmt_file <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", "CP:REACTOME")
} else if(msigdb == "canonical"){
  gmt_file <- msigdbr::msigdbr(category = "C2")
  gmt_file <- gmt_file %>%
    filter(gs_subcat %in% c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS"))
}
# modify column per the requirement of CEMiTool
gmt_file <- gmt_file %>%
  dplyr::rename("term" = "gs_name", 
                "gene" = "gene_symbol") %>%
  dplyr::select(term, gene)

# run the tool on particular cluster
set.seed(100)
run_cemitools_functions(expr_df = mat_of_interest, 
                        annot_df = cluster_anno, 
                        cor_method = opt$cor_method,
                        network_type = opt$network_type,
                        tom_type = opt$tom_type,
                        n = 100, 
                        output_dir = results_dir, 
                        plots_dir = plots_dir,
                        gmt_file = gmt_file)
