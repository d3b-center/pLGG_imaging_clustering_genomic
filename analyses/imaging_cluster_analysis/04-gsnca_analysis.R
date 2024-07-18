# Author: Komal S. Rathi, Run Jin
# Perform GSNCA analysis across different clusters
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(CEMiTool)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "imaging_cluster_analysis")

# source function
source(file.path(analysis_dir, "utils", "gsnca_calc.R"))

# parse command line options
option_list <- list(
  make_option(c("--mat_file"), type = "character",
              help = "gene matrix data, preferably tpm matrix (.rds) "),
  make_option(c("--anno_file"), type = "character",
              help = "path to cluster annotation file"),
  make_option(c("--gtf_file"), type = "character",
              help = "path to gtf annotation file"),
  make_option(c("--msigdb"), type = "character", 
              help = "reactome, kegg or canonical"),
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
mat <- readRDS(opt$mat_file)

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = opt$gtf_file)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter gene matrix to contain only protein coding gene
mat_coding <- mat[rownames(mat) %in% gencode_gtf$gene_name,]

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
# modify column per the requirement of GSNCA
gmt_file <- gmt_file %>%
  dplyr::rename("term" = "gs_name", 
                "gene" = "gene_symbol") %>%
  dplyr::select(term, gene)

# now perform analysis for each disease of interest 
# read in the cluster annotation file 
cluster_anno <- readr::read_tsv(anno_file) %>% 
  # also make rownames kids first biospecimen id
  dplyr::mutate(tmp = Kids_First_Biospecimen_ID) %>% 
  tibble::column_to_rownames("tmp")

# filter to samples in cancer group
mat_coding_cg_filtered <-  mat_coding %>%
  dplyr::select(any_of(cluster_anno$Kids_First_Biospecimen_ID)) %>%
  filter_low_expr_df()

# iterate through all cluster pair combination to get results 
cluster_n <- length(unique(cluster_anno$cluster_assigned))

# get all combinations
for(p in 1:(cluster_n - 1)){
  for(q in (p+1):cluster_n){
    group1_cluster <- p
    group2_cluster <- q
    
    # generate annotation file with only cluster of interest
    cluster_anno_each <- cluster_anno %>%
      dplyr::filter(cluster_assigned %in% c(p, q)) %>%
      dplyr::mutate(group = case_when(
        cluster_assigned == p ~ "1",
        cluster_assigned == q ~ "2"
      )) %>% 
      dplyr::select(Kids_First_Biospecimen_ID, group)
    
    # filter to these samples 
    mat_coding_cg_filtered_each <- mat_coding_cg_filtered %>%
      dplyr::select(any_of(cluster_anno_each$Kids_First_Biospecimen_ID))
    cluster_anno_each <- cluster_anno_each %>%
      filter(Kids_First_Biospecimen_ID %in% colnames(mat_coding_cg_filtered_each))
    
    # run the analysis 
    set.seed(100)
    gsnca_analysis_plot(tpm_matrix = as.matrix(mat_coding_cg_filtered_each), 
                        cluster_anno = cluster_anno_each, 
                        pathway_df = gmt_file, 
                        comparison = paste0("cluster", p, "_vs_cluster", q),
                        output_file_dir = results_dir, 
                        output_plot_dir = plots_dir, 
                        top_bar = 20, 
                        top_net = 5)
  }
}
