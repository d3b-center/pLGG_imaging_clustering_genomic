# Author: Komal S. Rathi
# differential genes and pathway enrichment analyses
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(DGCA)
  library(DESeq2)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "imaging_cluster_analysis")

# source function 
source(file.path(analysis_dir, "utils", "perform_enrichment_gsea.R"))

# parse command line options 
option_list <- list(
  make_option(c("--mat_file"), type = "character",
              help = "gene matrix data, preferably count matrix (.rds) "),
  make_option(c("--anno_file"), type = "character",
              help = "path to cluster annotation file"),
  make_option(c("--gtf_file"), type = "character",
              help = "path to gtf annotation file"),
  make_option(c("--plots_dir"), type = "character",
              help = "path to plots directory"),
  make_option(c("--results_dir"), type = "character",
              help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
plots_dir <- opt$plots_dir
results_dir <- opt$results_dir

# output directories
dir.create(plots_dir, showWarnings = F, recursive = T)
dir.create(results_dir, showWarnings = F, recursive = T)

# read count data
count_mat <- readRDS(opt$mat_file)

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = opt$gtf_file)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
count_mat <- count_mat %>%
  filter(rownames(count_mat) %in% gencode_gtf$gene_name)

# read annotation 
anno_file <- read_tsv(opt$anno_file)
anno_file$cluster_assigned <- as.factor(anno_file$cluster_assigned)
anno_file <- anno_file %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(count_mat))

# subset expression matrix to annotation file
count_mat <- count_mat %>%
  dplyr::select(any_of(anno_file$Kids_First_Biospecimen_ID))

# use DCGA to filter out low count, low variance features
count_mat <- DGCA::filterGenes(inputMat = count_mat, 
                               filterTypes = c("central", "dispersion"),
                               filterDispersionType = "cv", 
                               filterDispersionPercentile = 0.2,
                               sequential = TRUE)

# run DESeq2 to perform differential expression 
anno_file$cluster_assigned <- as.factor(anno_file$cluster_assigned)
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = anno_file,
                              design = ~ cluster_assigned)
dds <- DESeq(dds)
saveRDS(object = dds, file = file.path(results_dir, "dds_output.rds"))

# save tsv file of differentially expressed genes across all comparisons
deseq_output <- results(dds) %>% 
  as.data.frame() %>%
  rownames_to_column("gene") 
write_tsv(x = deseq_output, file = file.path(results_dir, "deseq2_output.tsv"))

# differential expression
comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
for(i in 1:length(comparisons)){
  
  # prefix for output files and plots
  prefix <- paste0("cluster", comparisons[[i]][1], "_vs_", "cluster", comparisons[[i]][2])
  
  # pull differentially expressed genes
  res <- results(dds, contrast = c("cluster_assigned", comparisons[[i]][1], comparisons[[i]][2]))
  res <- res %>%
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    dplyr::rename("log2FC" = "log2FoldChange") %>%
    rownames_to_column("gene")
  write_tsv(res, file.path(results_dir, paste0(prefix, "_deseq_output.tsv")))
  
  # filter to top 50 genes arranged by adjusted p-value
  top50_genes <- res %>%
    as.data.frame() %>%
    arrange(padj) %>%
    head(n = 50)
  
  # generate heatmap for top 50 differentially expressed genes for each comparison
  anno_file_sub <- anno_file %>%
    filter(cluster_assigned %in% comparisons[[i]]) %>%
    arrange(cluster_assigned, molecular_subtype)
  
  # ranked pathway enrichment for each comparison using reactome pathways
  reactome_pathways <- msigdbr::msigdbr(category = "C2", subcategory = "CP:REACTOME")
  reactome_pathways <- reactome_pathways %>% 
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename("term" = "gs_name",
                  "gene" = "gene_symbol")
  perform_enrichment_gsea(diffexpr_res = res,
                           pathways = reactome_pathways,
                           prefix = prefix,
                           plots_dir = file.path(plots_dir, "reactome"),
                           results_dir = file.path(results_dir, "reactome"))
}
