# BRAF fusions and survival by cluster

Module authors: Ryan Corbett (@rjcorb), Jo Lynne Rokita (@jharenza)

The purpose of this module is to assess survival in imaging clusters +/- BRAF fusion type

## Usage
### script to run analysis
<br>**Run shell script to make final tables to be used for plotting below**
```
bash run-analysis.sh
```
Input files:
```
../../data/fusion-dgd.tsv.gz
../../data/histologies.tsv
20231003_merged_imaging_clustering_risk_stratification.xlsx
imaging_clustering_cohort.csv
```

## Folder content
* `run-analysis.sh` shell script to run analysis
* `01-fusion-breakpoints-by-imaging-cluster.Rmd` Assess prevalence of BRAF fusion by cluster 
* `02-imaging-cluster-survival.Rmd` Assess survival based on imaging cluster and/or BRAF fusion type

## Directory structure
```
.
├── 01-fusion-breakpoints-by-imaging-cluster.Rmd
├── 01-fusion-breakpoints-by-imaging-cluster.html
├── 02-imaging-cluster-survival.Rmd
├── 02-imaging-cluster-survival.html
├── README.md
├── input
│   ├── 20231003_merged_imaging_clustering_risk_stratification.xlsx
│   └── imaging_clustering_cohort.csv
├── plots
│   ├── breakpoint-group-by-cluster-ht.pdf
│   ├── breakpoint-type-by-cluster-ht.pdf
│   ├── imaging_forest_efs_braf_add_resection_cluster.pdf
│   ├── imaging_forest_efs_braf_add_resection_cluster_group.pdf
│   ├── imaging_forest_efs_braf_add_resection_cluster_type.pdf
│   ├── imaging_forest_efs_int_resection_fusion_cluster.pdf
│   ├── imaging_forest_efs_int_resection_subtype_cluster.pdf
│   ├── imaging_forest_efs_nonbraf_add_resection_cluster.pdf
│   ├── km_efs_imaging_braf_breakpoint_group.pdf
│   ├── km_efs_imaging_braf_breakpoint_type.pdf
│   ├── km_efs_imaging_braf_cluster.pdf
│   ├── km_efs_imaging_braf_molecular_subtype.pdf
│   ├── km_efs_imaging_cluster.pdf
│   └── km_efs_imaging_nonbraf_cluster.pdf
├── results
│   ├── braf-fusion-breakpoints-by-patient.tsv
│   ├── coxph_imaging_braf_efs_add_resection_cluster.RDS
│   ├── coxph_imaging_efs_int_resection_fusion_cluster.RDS
│   ├── coxph_imaging_efs_int_resection_subtype_cluster.RDS
│   ├── coxph_imaging_nonbraf_efs_add_resection_cluster.RDS
│   ├── coxph_imaging_nonbraf_efs_add_resection_cluster_breakpoint_group.RDS
│   ├── coxph_imaging_nonbraf_efs_add_resection_cluster_breakpoint_type.RDS
│   ├── imaging_clusters_plus_subtypes_breakpoints.tsv
│   ├── imaging_forest_efs_braf_add_resection_cluster.pdf
│   ├── imaging_forest_efs_braf_add_resection_cluster_type.pdf
│   ├── imaging_forest_efs_int_resection_fusion_cluster.pdf
│   ├── imaging_forest_efs_nonbraf_add_resection_cluster.pdf
│   ├── lgg-braf-fusion-common-breakpoint-freq-image-clusters.tsv
│   └── lgg-braf-fusion-common-breakpoint-type-image-clusters.tsv
├── run-analysis.sh
└── util
    ├── heatmap_function.R
    └── survival_models_interaction.R
```