## Authors: Komal S. Rathi, Adam Kraya
 
### Data version

[Release 20230826](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230826_release) of master genomics.

Imaging clusters: `ImagingClusterAssignment_Aug2023.xlsx`

Risk groups: `RiskScores_Grouping.xlsx`

###  Run Analysis

The full analysis is run on `161 unique RNA-seq biospecimens` mapped to `150 unique cohort participants` out of a total of `201 unique cohort participants` from the imaging clusters file. 

```
# run full analysis
bash run_analysis.sh
```
***
### 1. Annotation file

`01-create-anno.R`: The purpose of this script is to create an annotation file by mapping relevant fields like `Kids_First_Biospecimen_ID`, `sample_id`, `experimental_strategy`, `pathology_diagnosis`, `cancer_group`, `short_histology`, `broad_histology`, `molecular_subtype`  from the histologies file to the imaging clusters file using the `cohort_participant_id` field. This is required for running the downstream scripts.

Using this, we obtained `161 unique RNA-seq biospecimens` mapped to `150 unique cohort ids`.  

This script also performs a `chisq test` between the `RNA-derived molecular subtypes` and the `imaging clusters` to check the concordance between the two groups. 

#### Input

```
../../data
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
8                   LGG, BRAF V600E    9
9         LGG, BRAF V600E, CDKN2A/B    4
10            LGG, BRAF V600E, FGFR    1
11   LGG, BRAF V600E, KIAA1549-BRAF    2
12    LGG, BRAF V600E, NF1-germline    1
13        LGG, CDKN2A/B, other MAPK    1
14                        LGG, FGFR    3
15         LGG, FGFR, KIAA1549-BRAF    2
16          LGG, FGFR, NF1-germline    2
17                         LGG, IDH    1
18               LGG, KIAA1549-BRAF   73
19 LGG, KIAA1549-BRAF, NF1-germline    1
20                   LGG, MYB/MYBL1    1
21                LGG, NF1-germline    1
22   LGG, NF1-germline, NF1-somatic    1
23    LGG, NF1-germline, other MAPK    1
24                 LGG, NF1-somatic    3
25                  LGG, other MAPK   11
26                         LGG, RTK    5
27            LGG, To be classified    1
28                    LGG, wildtype   19
29              SEGA, CDKN2A/B, RTK    1
```

After collapsing to larger groups:

```
   molecular_subtype_summarized freq
1            GNG, KIAA1549-BRAF    3
2                GNG, MYB/MYBL1    1
3               GNG, other MAPK    2
4         GNG, To be classified    1
5                 GNG, wildtype    5
6                  LGG_CDKN2A/B    5
7                       LGG_IDH    1
8                       LGG_NF1    8
9                       LGG_RTK   15
10              LGG, BRAF V600E   15
11           LGG, KIAA1549-BRAF   73
12               LGG, MYB/MYBL1    1
13              LGG, other MAPK   11
14        LGG, To be classified    1
15                LGG, wildtype   19
```
***
### 2. Differential expression and pathway enrichment 

#### 2.1 DESeq2

`02-diffexpr-per-cluster-deseq.R`: This script performs differential expression between each pair of imaging cluster using the `DESeq2` R package. Then, using the differentially expressed genes at an `FDR of < 0.05`, it performs a pre-ranked enrichment test on `REACTOME` pathways using the `clusterProfiler::GSEA` function. 

Using the enrichment output, a barplot of top 10 upregulated and downregulated pathways is generated. This script can be run using the following command: 

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
../../data
├── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output

The results are written under  `results/deseq` and plots are generated under `plots/deseq` folder.  

For each pairwise comparison, the DESeq2 output with FDR adjusted p-value < 0.05 is saved under `_deseq_output.tsv` and the GSEA output with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.

```
results/deseq
├── cluster1_vs_cluster2_deseq_output.tsv
├── cluster1_vs_cluster3_deseq_output.tsv
├── cluster2_vs_cluster3_deseq_output.tsv
├── dds_output.rds # full DESeq2 output
├── deseq2_output.tsv # differential genes across all clusters
└── reactome
    ├── cluster1_vs_cluster2_gsea.tsv
    ├── cluster1_vs_cluster3_gsea.tsv
    └── cluster2_vs_cluster3_gsea.tsv
```

A barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by GSEA at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`.

```
plots/deseq
└── reactome
    ├── cluster1_vs_cluster2_gsea_barplot.pdf
    ├── cluster1_vs_cluster3_gsea_barplot.pdf
    └── cluster2_vs_cluster3_gsea_barplot.pdf
```

***

### 3. CEMITool Network analysis

`03-cemitools_analysis.R`: Perform CEMITool network analysis using the `REACTOME` pathways. This script creates output files and plots under `results/network_analysis_spearman_unsigned_signed` and plots under `plots/network_analysis_spearman_unsigned_signed`.

We used two approaches with different input parameters (changing the correlation method while keeping `network_type` constant did not result in any difference):

1. `cor_method = "pearson"`, `network_type = "signed"`, `tom_type = "signed"` 
2. `cor_method = "spearman"`, `network_type = "unsigned"`, `tom_type = "signed"` 

Based on our manual inspection, we chose the second method to be optimal for our analyses.

#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../../data
├── 20230826_release-gene-counts-rsem-expected_count.collapsed.rds
└── gencode.v39.primary_assembly.annotation.gtf.gz
```

#### Output

Output files:

```
results/network_analysis_spearman_unsigned_signed
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
plots/network_analysis_spearman_unsigned_signed
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
### 4. GSNCA co-expression analysis

`04-gsnca_analysis.R`: Perform between-cluster GSNCA co-expression network analysis, using the `REACTOME` pathways. This script creates output files and plots under `results/gsnca_analysis` and plots under `plots/gsnca_analysis`.
 
#### Input

```
# input cluster annotation file generated in the first step
results/cluster_anno
└── LGG_cluster_info.tsv

# 20230826_release data
../../data
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
### 5. SSGSEA

`05-ssgsea.R`: `161 unique RNA-seq biospecimens` corresponding to `150 of the 201 cohort participant identifiers` in the risk group file were pulled from the  `20230826_release`  histology file. 

GSVA using  `ssgsea`  method was performed on  `REACTOME`  gene sets with minimum set size of 10 and maximum set size of 500. The resulting matrix of 1292 gene sets across 161 biospecimens was saved under  `ssgsea_matrix.rds`. 

#### Input

```
# risk scores file
../../data
└── RiskScores_Grouping.xlsx

# 20230826_release data
../../data
├── 20230826_release.annotated_histologies.tsv
└── 20230826_release-gene-expression-rsem-tpm.collapsed.rds
```

#### Output

Output files:

```
# SSGSEA matrix
results/ssgsea_output
└── ssgsea_matrix.rds
```
***

### 6. Prediction Models

`06-ML-prediction_models.R`:  **TBD**


#### Input

```
# Imaging clustering assignments
../../data
└── ImagingClusterAssignment_Aug2023.xlsx

# SSGSEA matrix
results/ssgsea_output
└── ssgsea_matrix.rds

# 20230826_release data
../../data
└── 20230826_release.annotated_histologies.tsv
```

#### Output

Output files:

```
results/imaging_cluster_prediction
└── gbm
    ├── influence_score_one_v_rest.txt
    ├── influence_score_three_v_rest.txt
    └── influence_score_two_v_rest.txt
```

#### Plots

```
plots/imaging_cluster_prediction
└── gbm
    ├── gbm_roc_curves_is_cluster1.pdf
    ├── gbm_roc_curves_is_cluster2.pdf
    ├── gbm_roc_curves_is_cluster3.pdf
    ├── gbm_top_20_features_one_v_rest.pdf
    ├── gbm_top_20_features_three_v_rest.pdf
    └── gbm_top_20_features_two_v_rest.pdf
```
