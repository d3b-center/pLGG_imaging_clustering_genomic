## Author: Komal S. Rathi
 
## Purpose

This module performs the following analyses:

1) Differential gene expression and identification of significantly deregulated pathways. 
2) CEMITool Network and GSNCA Co-expression analysis for LGG Imaging clusters using REACTOME database.
3) Compute and plot Tumor inflammation signature scores.

### Data version

[Release 20230826](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230826_release) of master genomics.

Imaging clusters: `ImagingClusterAssignment_Aug2023.xlsx`

###  Run Analysis

The full analysis is run on `172 unique biospecimens` from `201 cohort participants` from the imaging clusters file. 

```
# run full analysis
bash run_analysis.sh
```
***
### 1. Annotation file

`01-create-anno.R`: The purpose of this script is to create an annotation file by mapping important fields like `Kids_First_Biospecimen_ID` , `molecular_subtype` etc from the histologies file to the imaging clusters file. This is required for running the downstream scripts. This script also performs a chisq test between the RNA-derived molecular subtypes and the imaging clusters to see the concordance between the two groups. 

#### Input

```
../data
└── ImagingClusterAssignment_Aug2023.xlsx
```

#### Output

Results:

```
results/cluster_anno
├── LGG_cluster_info.tsv # full cluster info
└── imaging_clusters_vs_molsubtype_chisq_test.txt # chisq test results for imaging clusters vs RNA-derived molecular subtypes
```

Plots:

```
plots/cluster_anno
├── imaging_clusters_vs_molsubtype_balloonplot.pdf # frequency plot of imaging clusters vs RNA-derived molecular subtypes
└── imaging_clusters_vs_molsubtype_corrplot.pdf # correlation plot of imaging clusters vs RNA-derived molecular subtypes
```

#### Summary

Breakdown of molecular subtypes before collapsing, majority of samples belong to the `LGG, KIAA1549-BRAF` subtype:

```
                  molecular_subtype freq
1                   GNG, BRAF V600E    4
2                GNG, KIAA1549-BRAF    3
3                    GNG, MYB/MYBL1    1
4                   GNG, other MAPK    2
5                          GNG, RTK    1
6             GNG, To be classified    1
7                     GNG, wildtype    5
8            GNT, FGFR, NF1-somatic    1
9                   LGG, BRAF V600E   10
10        LGG, BRAF V600E, CDKN2A/B    5
11            LGG, BRAF V600E, FGFR    1
12   LGG, BRAF V600E, KIAA1549-BRAF    2
13    LGG, BRAF V600E, NF1-germline    1
14        LGG, CDKN2A/B, other MAPK    1
15                        LGG, FGFR    3
16         LGG, FGFR, KIAA1549-BRAF    2
17          LGG, FGFR, NF1-germline    2
18                         LGG, IDH    1
19               LGG, KIAA1549-BRAF   74
20 LGG, KIAA1549-BRAF, NF1-germline    1
21                   LGG, MYB/MYBL1    1
22                LGG, NF1-germline    1
23   LGG, NF1-germline, NF1-somatic    1
24    LGG, NF1-germline, other MAPK    1
25                 LGG, NF1-somatic    3
26                  LGG, other MAPK   11
27                         LGG, RTK    6
28            LGG, To be classified    1
29                    LGG, wildtype   25
30              SEGA, CDKN2A/B, RTK    1
```

After collapsing to larger groups:

```
   molecular_subtype_summarized freq
1            GNG, KIAA1549-BRAF    3
2                GNG, MYB/MYBL1    1
3               GNG, other MAPK    2
4         GNG, To be classified    1
5                 GNG, wildtype    5
6                  LGG_CDKN2A/B    6
7                       LGG_IDH    1
8                       LGG_NF1    8
9                       LGG_RTK   17
10              LGG, BRAF V600E   16
11           LGG, KIAA1549-BRAF   74
12               LGG, MYB/MYBL1    1
13              LGG, other MAPK   11
14        LGG, To be classified    1
15                LGG, wildtype   25
```
***
### 2. Differential expression and pathway enrichment 

#### 2.1 DESeq2

`02-diffexpr-per-cluster-deseq.R`: This script performs differential expression between each pair of imaging cluster using the `DESeq2` R package. Then, using the differentially expressed genes at an `FDR of < 0.05`, it performs a pre-ranked enrichment test on `REACTOME` , `KEGG MEDICUS`  as well as `Canonical (BIOCARTA, KEGG, PID, REACTOME, WIKIPATHWAYS)` pathways using the `clusterProfiler::GSEA` function. 

Using the enrichment output, a barplot of top 10 upregulated and downregulated pathways, a network plot and a dotplot is plotted. This script can be run using the following command: 

```
# run differential expression and identify significant pathways
Rscript 02-diffexpr-per-cluster-deseq.R
```

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
├── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output

The results are written under  `results/deseq` and plots are generated under `plots/deseq` folder.  

For each pairwise comparison, the DESeq2 output with FDR adjusted p-value < 0.05 is saved under `_deseq_output.tsv` and the GSEA output with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.

```
results/deseq
├── c2_canonical
│   ├── cluster1_vs_cluster2_gsea.tsv
│   ├── cluster1_vs_cluster3_gsea.tsv
│   └── cluster2_vs_cluster3_gsea.tsv
├── cluster1_vs_cluster2_deseq_output.tsv
├── cluster1_vs_cluster3_deseq_output.tsv
├── cluster2_vs_cluster3_deseq_output.tsv
├── dds_output.rds # full DESeq2 output
├── deseq2_output.tsv # differential genes across all clusters
├── kegg_medicus
│   ├── cluster1_vs_cluster2_gsea.tsv
│   ├── cluster1_vs_cluster3_gsea.tsv
│   └── cluster2_vs_cluster3_gsea.tsv
└── reactome
    ├── cluster1_vs_cluster2_gsea.tsv
    ├── cluster1_vs_cluster3_gsea.tsv
    └── cluster2_vs_cluster3_gsea.tsv
```

A heatmap of top 50 differentially expressed genes across all comparisons (`top50_genes_heatmap.pdf`) and for all pairwise comparisons (`*_heatmap.pdf`) is created. 

Additionally, a barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by GSEA at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`, dotplot is generated under `*_gsea_dotplot.pdf` and network under `*_gsea_cnet.pdf`.

```
plots/deseq
├── c2_canonical
│   ├── cluster1_vs_cluster2_gsea_barplot.pdf
│   ├── cluster1_vs_cluster2_gsea_cnet.pdf
│   ├── cluster1_vs_cluster2_gsea_dotplot.pdf
│   ├── cluster1_vs_cluster3_gsea_barplot.pdf
│   ├── cluster1_vs_cluster3_gsea_cnet.pdf
│   ├── cluster1_vs_cluster3_gsea_dotplot.pdf
│   ├── cluster2_vs_cluster3_gsea_barplot.pdf
│   ├── cluster2_vs_cluster3_gsea_cnet.pdf
│   └── cluster2_vs_cluster3_gsea_dotplot.pdf
├── cluster1_vs_cluster2_heatmap.pdf
├── cluster1_vs_cluster3_heatmap.pdf
├── cluster2_vs_cluster3_heatmap.pdf
├── kegg_medicus
│   ├── cluster1_vs_cluster2_gsea_barplot.pdf
│   ├── cluster1_vs_cluster2_gsea_cnet.pdf
│   ├── cluster1_vs_cluster2_gsea_dotplot.pdf
│   ├── cluster1_vs_cluster3_gsea_barplot.pdf
│   ├── cluster1_vs_cluster3_gsea_cnet.pdf
│   └── cluster1_vs_cluster3_gsea_dotplot.pdf
├── reactome
│   ├── cluster1_vs_cluster2_gsea_barplot.pdf
│   ├── cluster1_vs_cluster2_gsea_cnet.pdf
│   ├── cluster1_vs_cluster2_gsea_dotplot.pdf
│   ├── cluster1_vs_cluster3_gsea_barplot.pdf
│   ├── cluster1_vs_cluster3_gsea_cnet.pdf
│   ├── cluster1_vs_cluster3_gsea_dotplot.pdf
│   ├── cluster2_vs_cluster3_gsea_barplot.pdf
│   ├── cluster2_vs_cluster3_gsea_cnet.pdf
│   └── cluster2_vs_cluster3_gsea_dotplot.pdf
└── top50_genes_heatmap.pdf
```

#### 2.2 NOISeq


`02-diffexpr-per-cluster-noiseq.R`: This script performs differential expression between each pair of imaging cluster using the `NOISeq` R package. Then, using the differentially expressed genes at an `FDR of < 0.05`, it performs a pre-ranked enrichment test on `REACTOME` , `KEGG MEDICUS`  as well as `Canonical (BIOCARTA, KEGG, PID, REACTOME, WIKIPATHWAYS)` pathways using the `clusterProfiler::GSEA` function. 

Using the enrichment output, a barplot of top 10 upregulated and downregulated pathways, a network plot and a dotplot is plotted. This script can be run using the following command: 

```
# run differential expression and identify significant pathways
Rscript 02-diffexpr-per-cluster-noiseq.R
```

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
├── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output


The results are written under  `results/noiseq` and plots are generated under `plots/noiseq` folder.  

For each pairwise comparison, the NOISeq output with FDR adjusted p-value < 0.05 is saved under `_noiseq_output.tsv` and the GSEA output with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.

```
results/noiseq
├── c2_canonical
│   ├── cluster1_vs_cluster2_gsea.tsv
│   ├── cluster1_vs_cluster3_gsea.tsv
│   └── cluster2_vs_cluster3_gsea.tsv
├── cluster1_vs_cluster2_noiseq_output.tsv
├── cluster1_vs_cluster3_noiseq_output.tsv
├── cluster2_vs_cluster3_noiseq_output.tsv
├── kegg_medicus
│   ├── cluster1_vs_cluster2_gsea.tsv
│   ├── cluster1_vs_cluster3_gsea.tsv
│   └── cluster2_vs_cluster3_gsea.tsv
└── reactome
    ├── cluster1_vs_cluster2_gsea.tsv
    ├── cluster1_vs_cluster3_gsea.tsv
    └── cluster2_vs_cluster3_gsea.tsv
```

A heatmap of top 50 differentially expressed genes for all pairwise comparisons (`*_heatmap.pdf`) is created. 

Additionally, a barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by GSEA at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`, dotplot is generated under `*_gsea_dotplot.pdf` and network under `*_gsea_cnet.pdf`.

```
plots/noiseq
├── c2_canonical
│   ├── cluster1_vs_cluster2_gsea_barplot.pdf
│   ├── cluster1_vs_cluster2_gsea_cnet.pdf
│   ├── cluster1_vs_cluster2_gsea_dotplot.pdf
│   ├── cluster1_vs_cluster3_gsea_barplot.pdf
│   ├── cluster1_vs_cluster3_gsea_cnet.pdf
│   ├── cluster1_vs_cluster3_gsea_dotplot.pdf
│   ├── cluster2_vs_cluster3_gsea_barplot.pdf
│   ├── cluster2_vs_cluster3_gsea_cnet.pdf
│   └── cluster2_vs_cluster3_gsea_dotplot.pdf
├── cluster1_vs_cluster2_heatmap.pdf
├── cluster1_vs_cluster3_heatmap.pdf
├── cluster2_vs_cluster3_heatmap.pdf
├── kegg_medicus
│   ├── cluster1_vs_cluster2_gsea_barplot.pdf
│   ├── cluster1_vs_cluster2_gsea_cnet.pdf
│   ├── cluster1_vs_cluster2_gsea_dotplot.pdf
│   ├── cluster1_vs_cluster3_gsea_barplot.pdf
│   ├── cluster1_vs_cluster3_gsea_cnet.pdf
│   └── cluster1_vs_cluster3_gsea_dotplot.pdf
├── noiseq_pca.pdf # PCA plot of before and after batch correction
└── reactome
    ├── cluster1_vs_cluster2_gsea_barplot.pdf
    ├── cluster1_vs_cluster2_gsea_cnet.pdf
    ├── cluster1_vs_cluster2_gsea_dotplot.pdf
    ├── cluster1_vs_cluster3_gsea_barplot.pdf
    ├── cluster1_vs_cluster3_gsea_cnet.pdf
    ├── cluster1_vs_cluster3_gsea_dotplot.pdf
    ├── cluster2_vs_cluster3_gsea_barplot.pdf
    ├── cluster2_vs_cluster3_gsea_cnet.pdf
    └── cluster2_vs_cluster3_gsea_dotplot.pdf
```
***
### 3. Heatmap of Dark kinases

`03-dark-kinases-heatmap.R`: This script takes top 50 differentially expressed dark kinases from the differential expression output generated in script 2 and plots a heatmap of log2 TPM (z-scored).  

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
└── 20230826_release-gene-expression-rsem-tpm.collapsed.rds

# deseq2 full output from script 2
results/deseq
└── deseq2_output.tsv
```

#### Output

There were 38 differentially expressed kinases identified for the LGG cohort:

```
plots/deseq
└── top50_darkkinases_heatmap.pdf
```
***
## The following analyses are from the `../rnaseq_analysis`  module


### 4. CEMITool Network analysis

`../rnaseq_analysis/01-cemitools_analysis.R`: Perform CEMITool network analysis using the `REACTOME` pathways. This script creates output files and plots under `results/network_analysis` and plots under `plots/network_analysis`.

We used two approaches with different input parameters (changing the correlation method while keeping network_type constant did not result in any difference):

1. `cor_method = "pearson"`, `network_type = "signed"`, `tom_type = "signed"` 
2. `cor_method = "spearman"`, `network_type = "unsigned"`, `tom_type = "signed"` 

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
├── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output

Output files:

```
results/network_analysis_{pearson_signed_signed, spearman_unsigned_signed}
├── diagnostics.html
├── enrichment_es.tsv
├── enrichment_nes.tsv
├── enrichment_padj.tsv
├── hubs.rds
├── interactions.tsv
├── module.tsv
├── modules_genes.gmt
├── ora.tsv
├── parameters.tsv
├── report.html
├── selected_genes.txt
├── summary_eigengene.tsv
├── summary_mean.tsv
└── summary_median.tsv
```

Plots:

```
plots/network_analysis_{pearson_signed_signed, spearman_unsigned_signed}
├── beta_r2.pdf
├── gsea.pdf
├── hist.pdf
├── interaction.pdf
├── mean_k.pdf
├── mean_var.pdf
├── ora.pdf
├── profile.pdf
├── qq.pdf
└── sample_tree.pdf
```
***
### 5. GSNCA co-expression analysis

`../rnaseq_analysis/02-gsnca_analysis.R`: Perform between-cluster GSNCA co-expression network analysis, using the `REACTOME` pathways. This script creates output files and plots under `results/gsnca_analysis` and plots under `plots/gsnca_analysis`.
 
#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
├── 20230826_release-gene-expression-rsem-tpm.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output

Output files:

```
results/gsnca_analysis
├── cluster1_vs_cluster2_GSNCA_anaysis.tsv
├── cluster1_vs_cluster3_GSNCA_anaysis.tsv
└── cluster2_vs_cluster3_GSNCA_anaysis.tsv
```

Plots:

```
plots/gsnca_analysis
├── cluster1_vs_cluster2_GSNCA_plots.pdf
├── cluster1_vs_cluster3_GSNCA_plots.pdf
└── cluster2_vs_cluster3_GSNCA_plots.pdf
```
***
### 6. Tumor inflammation signature

`../rnaseq_analysis/03-tis-profiling.R`: This script uses Tumor inflammation signature (TIS) genes to calculate a TIS score per imaging cluster. This script creates output plots under `plots/tis_analysis`. 

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../data
└── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
```

#### Output

Output files:

```
results/tis_analysis
├── tis_anova_results.tsv
└── tis_tukey_results.tsv
```

For each imaging cluster, a Tumor inflammation score is calculated and a boxplot of scores by imaging cluster is generated. Further, a z-scored heatmap of individual TIS genes X sample, ordering and annotating by imaging cluster and secondarily ordering and annotating by molecular subtype is created.

```
plots/tis_analysis
├── TIS_boxplot.pdf
└── TIS_heatmap.pdf
```

#### Summary

No TIS enrichment was observed in any specific cluster vs the other clusters.
***

### 7. Clustering Analyses

`05-clustering_analysis.R`: This script pulls samples from Imaging clusters and applies multiple clustering algorithms like 1) `IntNMF`, 2) `Consensus Clustering` , 3) `MClust` and 4) `HDBSCAN` on the corresponding count matrix via [ClusTarIDseq](https://github.com/d3b-center/ClusTarIDseq) package. For each clustering method the following transformation types and feature selection were evaluated:

1. Transformation types: 
    - None (no transformation)
    - Log2 
    - Rank 
2. Feature selection: 
    - Top 10% variable features (RNA-seq dataset), 100% input features (PCA, Filtered dataset)
    - Features selected using dip.test `pvalue < 0.05`

`06-compare_clusters.R`: This script pulls derived clusters from `05-clustering_analysis.R` and compares them to clusters derived from imaging analysis as well as RNA-derived molecular subtypes. 

For each of the methods evaluated (i.e. hdbscan, best fit CCP, IntNMF), it generates summary stats `{method}_{transformation_type}_{feature_selection}_summary_stats.tsv` that contains comparisons between the two classes like `Adjusted Rand Index`, `Chisq test statistic` and `Chisq test p-value`. Using the stats, it also generates corresponding balloon plots `{method}_{transformation_type}_{feature_selection}_balloonplot.pdf` and correlation plots `{method}_{transformation_type}_{feature_selection}_corrplot.pdf`.  

`06-compare_clusters_ccp.R`: This script is for manual evaluation of CCP clusters that did not come up as the best fit using the `final_composite_score.R script`. It allows users to input specific clustering output from CCP and compares them to clusters derived from imaging analysis as well as RNA-derived molecular subtypes. 

#### Input

```
# 20230826_release data subsetted to samples from the Imaging clusters file
../data
├── ImagingClusterAssignment_Aug2023.xlsx
├── 20230826_release.annotated_histologies_subset.tsv
└── 20230826_release-gene-counts-rsem-expected_count.collapsed_subset.rds
```

#### Output

Example of log2 dip-test output directory:

```
results/image_clustering/log2/dip.test
├── ccp_output
│   ├── hc_euclidean_3_log2_dip.test_balloonplot.pdf
│   ├── ccp_optimal_clusters.tsv
│   ├── hc_euclidean_3_log2_dip.test_corrplot.pdf
│   ├── {hc, pam, km}_{binary, canberra, euclidean, maximum, pearson, spearman}_dip_test.pdf
│   ├── {hc, pam, km}_{binary, canberra, euclidean, maximum, pearson, spearman}_dip_test.rds
│   ├── {hc, pam, km}_{binary, canberra, euclidean, maximum, pearson, spearman}_dip_test.tsv
│   ├── lspline_output.tsv
│   └── hc_euclidean_3_log2_dip.test_summary_stats.tsv
├── final_score
│   └── final_clustering_output.tsv
├── hdbscan_output
│   ├── hdbscan_log2_dip.test_balloonplot.pdf
│   ├── hdbscan_log2_dip.test_corrplot.pdf
│   ├── hdbscan_optimal_clusters.tsv
│   ├── hdbscan_plot.pdf
│   └── hdbscan_log2_dip.test_summary_stats.tsv
└── intnmf_output
    ├── IntNMF_log2_dip.test_balloonplot.pdf
    ├── IntNMF_log2_dip.test_corrplot.pdf
    ├── intnmf_best_fit.rds
    ├── intnmf_clusterstats.tsv
    ├── intnmf_consensus_plot.pdf
    ├── intnmf_fit_all.rds
    ├── intnmf_optimal_clusters.tsv
    ├── intnmf_silhouette_plot.pdf
    └── IntNMF_log2_dip.test_summary_stats.tsv
```

#### Summary

The summary stats from each evaluated method + transformation + feature selection combination is combined into one file for easy navigation of relevant results:

```
results/image_clustering
└── combined_summary_stats.tsv
```

`07-ari_df_global.R`: This script pulls all relevant clustering outputs (CCP, MClust, IntNMF and HDBSCAN) and generates a global dataframe containing the Adjusted Rank Index between the derived clusters vs Imaging clusters and Molecular subtypes.

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# Imaging clusters and 20230826_release data subsetted to samples from the Imaging clusters file
../data
├── ImagingClusterAssignment_Aug2023.xlsx
└── 20230826_release.annotated_histologies_subset.tsv
```

#### Output

```
results
└── ccp_ari_index.tsv
```

***

#### Run script

The clustering analyses can be run using a single bash script:

```
# for RNA-seq clustering
bash run_analysis_clustering.sh

# for filtered dataset clustering
bash run_analysis_clustering_filtered.sh

# for PCA dataset
bash run_analysis_clustering_pca.sh
```
