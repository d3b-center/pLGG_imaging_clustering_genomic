#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Define directory and input files
data_dir="../../data"
count_file="${data_dir}/20230826_release-gene-counts-rsem-expected_count.collapsed.rds"
tpm_file="${data_dir}/20230826_release-gene-expression-rsem-tpm.collapsed.rds"
gtf_file="${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz"
msigdb="reactome"

# Run script to generate cluster annotation file
Rscript --vanilla 01-create-anno.R

# Run differential expression analysis
Rscript --vanilla 02-diffexpr-per-cluster-deseq.R \
--mat_file $count_file \
--anno_file "results/cluster_anno/LGG_cluster_info.tsv" \
--gtf_file $gtf_file \
--plots_dir "plots/deseq" \
--results_dir "results/deseq"

# Run CEMiTools network analysis using REACTOME
# cor_method = "spearman", network_type = "unsigned", and tom_type = "signed"
Rscript --vanilla 03-cemitools_analysis.R \
--mat_file $count_file \
--anno_file "results/cluster_anno/LGG_cluster_info.tsv" \
--gtf_file $gtf_file \
--msigdb $msigdb \
--cor_method "spearman" \
--network_type "unsigned" \
--tom_type "signed" \
--plots_dir "plots/network_analysis_spearman_unsigned_signed" \
--results_dir "results/network_analysis_spearman_unsigned_signed"

# Run GSNCA co-expression analysis using REACTOME
Rscript --vanilla 04-gsnca_analysis.R \
--mat_file $tpm_file \
--anno_file "results/cluster_anno/LGG_cluster_info.tsv" \
--gtf_file $gtf_file \
--msigdb $msigdb \
--plots_dir "plots/gsnca_analysis" \
--results_dir "results/gsnca_analysis"

# Run ssGSEA using the GSVA implementation to generate Reactome pathway scores (curated c2 pathways from msigdb) 
# on the rnaseq data derived from LGG participant list 
Rscript 05-ssgsea.R

# ML-based prediction models on ssgsea matrix generated above
Rscript 06-ML-prediction_models.R
