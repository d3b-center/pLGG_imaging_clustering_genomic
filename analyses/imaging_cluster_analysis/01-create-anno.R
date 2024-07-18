# Author: Komal S. Rathi
# create lgg annotation file from imaging excel file 
suppressPackageStartupMessages({
  library(tidyverse)
  library(gplots)
  library(corrplot)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "imaging_cluster_analysis")

# source function for chisq test + balloon plots
source(file.path(analysis_dir, "utils", "corr_plots.R"))

# output directories
output_dir <- file.path(analysis_dir, "results", "cluster_anno")
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- file.path(analysis_dir, "plots", "cluster_anno")
dir.create(plots_dir, showWarnings = F, recursive = T)

# read excel file 
lgg_clusters <- readxl::read_xlsx(file.path(data_dir, "ImagingClusterAssignment_Aug2023.xlsx"))
lgg_clusters <- lgg_clusters %>%
  dplyr::select(SubjectID, `Cluster Assignment`, Cohort)

# if minimum value is 0, then add +1
if(min(lgg_clusters$`Cluster Assignment`) == 0){
  lgg_clusters$`Cluster Assignment` <- lgg_clusters$`Cluster Assignment` + 1
}

# 1) molecular subtypes from 20230826_release histologies vs imaging clusters 
# filter histologies from OT to RNA-Seq and participants in the excel file
histology_file = file.path(data_dir, "20230826_release.annotated_histologies.tsv")
histology_file <- read_tsv(histology_file)
histology_file <- histology_file %>%
  filter(experimental_strategy == "RNA-Seq",
         cohort_participant_id %in% lgg_clusters$SubjectID) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cohort_participant_id, sample_id, experimental_strategy, pathology_diagnosis, cancer_group, short_histology, broad_histology, molecular_subtype)

# only use samples with actual RNA-seq data 
count_mat = file.path(data_dir, "20230826_release-gene-counts-rsem-expected_count.collapsed.rds")
count_mat <- readRDS(count_mat)
histology_file <- histology_file %>%
  filter(Kids_First_Biospecimen_ID %in% colnames(count_mat))

# combine both 
lgg_cluster_file <- histology_file %>%
  inner_join(lgg_clusters, by = c("cohort_participant_id" = "SubjectID")) %>%
  dplyr::rename("cluster_assigned" = "Cluster Assignment",
                "cohort" = "Cohort") 

# distribution before collapsing
plyr::count(lgg_cluster_file, "molecular_subtype")

# reduce molecular subtypes by grouping them
lgg_cluster_file$molecular_subtype <- as.character(lgg_cluster_file$molecular_subtype)
lgg_cluster_file$molecular_subtype_summarized <- lgg_cluster_file$molecular_subtype
# Group all samples containing 'BRAF V600E' in the subtype string into an 'LGG, BRAF V600E' group
lgg_cluster_file$molecular_subtype_summarized[grep("BRAF V600E", lgg_cluster_file$molecular_subtype)] <- "LGG, BRAF V600E"
# Group all samples containing 'BRAF-KIAA1549' in the subtype string into an 'LGG, BRAF-KIAA1549' group
lgg_cluster_file$molecular_subtype_summarized[grep("BRAF-KIAA1549", lgg_cluster_file$molecular_subtype)] <- "LGG, BRAF-KIAA1549"
# Group all remaining samples containing 'NF1' in the subtype string into an 'LGG_NF1' group.
lgg_cluster_file$molecular_subtype_summarized[grep("NF1", lgg_cluster_file$molecular_subtype)] <- "LGG_NF1"
# Group all remaining samples containing 'CDKN2A/B' in the subtype string into an 'LGG_CDKN2A/B' group.
lgg_cluster_file$molecular_subtype_summarized[grep("CDKN2A/B", lgg_cluster_file$molecular_subtype)] <- "LGG_CDKN2A/B"
# Group all remaining samples containing 'IDH' in the subtype string into an 'LGG_IDH' group.
lgg_cluster_file$molecular_subtype_summarized[grep("IDH", lgg_cluster_file$molecular_subtype)] <- "LGG_IDH"
# Group all remaining samples containing 'RTK' or 'FGFR' in the subtype string into an 'LGG_RTK' subgroup
lgg_cluster_file$molecular_subtype_summarized[grep("RTK|FGFR", lgg_cluster_file$molecular_subtype)] <- "LGG_RTK"
# Modify MYB other MAPK
lgg_cluster_file$molecular_subtype_summarized[grep("LGG, MYB/MYBL1, other MAPK", lgg_cluster_file$molecular_subtype)] <- "LGG, MYB/MYBL1"

# cluster distribution after collapsing
plyr::count(lgg_cluster_file, "molecular_subtype_summarized")
lgg_cluster_file <- lgg_cluster_file %>%
  dplyr::select(-c(molecular_subtype)) %>%
  dplyr::rename("molecular_subtype" = "molecular_subtype_summarized")

# save full file for downstream analyses
write_tsv(lgg_cluster_file, file = file.path(output_dir, 'LGG_cluster_info.tsv'))

# chisq and corrplots
corr_plots(lgg_cluster_file = lgg_cluster_file, 
           prefix = "imaging_clusters_vs_molsubtype")
