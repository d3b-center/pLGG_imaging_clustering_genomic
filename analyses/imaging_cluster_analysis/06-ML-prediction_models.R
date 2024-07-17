# Authors: Adam Kraya
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(plyr)
  library(readxl)
  library(readr)
  library(purrr)
  library(tidyr)
  library(caret)
  library(pROC)
  library(hdi)
  library(ggplot2)
  library(effects)
  library(hydroGOF)
  library(MLmetrics)
  library(doParallel)
  library(mboost)
  library(import)
  library(gbm)
  library(RColorBrewer)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "imaging_cluster_analysis")

# output directory
output_dir <-
  file.path(analysis_dir, "results", "imaging_cluster_prediction")
gbm_output_dir <- file.path(output_dir, "gbm")
dir.create(gbm_output_dir,
           showWarnings = F,
           recursive = T)

# plots dirs
plots_dir <-
  file.path(analysis_dir, "plots", "imaging_cluster_prediction")
gbm_plots_dir <- file.path(plots_dir, "gbm")
dir.create(gbm_plots_dir, showWarnings = F, recursive = T)

# read histology file
histology_file <-
  file.path(data_dir, "20230826_release.annotated_histologies.tsv") %>%
  fread()

# read imaging risk file and pull corresponding bs identifiers
lgg_clusters <-
  file.path(data_dir, "ImagingClusterAssignment_Aug2023.xlsx") %>%
  read_excel()

histology_file <- histology_file %>%
  dplyr::filter(
    cohort_participant_id %in% lgg_clusters$SubjectID,
    experimental_strategy == "RNA-Seq"
  ) %>%
  dplyr::select(
    cohort_participant_id,
    Kids_First_Biospecimen_ID,
    molecular_subtype,
    CNS_region,
    reported_gender,
    race,
    age_at_diagnosis_days
  )

cluster_labels <- merge(histology_file,
                        lgg_clusters,
                        by.x = 'cohort_participant_id',
                        by.y = 'SubjectID') %>%
  dplyr::select(
    cluster_assignment = "Cluster Assignment",
    Kids_First_Biospecimen_ID,
    molecular_subtype,
    CNS_region,
    reported_gender,
    race,
    age_at_diagnosis_days
  )

# read ssgsea_matrix file
ssgsea_scores <-
  readRDS(file.path(analysis_dir, "results", "ssgsea_output", "ssgsea_matrix.rds")) %>% t() %>% as.data.frame()

# remove zero variance or near zero variance predictors
mads <- apply(ssgsea_scores, 2, mad)
ssgsea_scores <- ssgsea_scores[, rev(order(mads))[1:100]]

# Join with cluster labels to create a matrix
ssgsea_scores <- ssgsea_scores %>%
  rownames_to_column('Kids_First_Biospecimen_ID')

joint_table <- ssgsea_scores %>%
  dplyr::inner_join(cluster_labels, by = 'Kids_First_Biospecimen_ID') %>%
  dplyr::mutate(cluster_assignment = cluster_assignment + 1) %>%
  dplyr::mutate(cluster_assignment = as.factor(cluster_assignment)) %>%
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) %>%
  dplyr::mutate(molecular_subtype = gsub('-', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub('\\, ', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub(' ', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub('\\/', '_', molecular_subtype)) %>%
  dplyr::mutate(CNS_region = gsub('-', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub('\\, ', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub(' ', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub('\\/', '_', CNS_region)) %>%
  dplyr::mutate(race = gsub('-', '_', race)) %>%
  dplyr::mutate(race = gsub('\\, ', '_', race)) %>%
  dplyr::mutate(race = gsub(' ', '_', race)) %>%
  dplyr::mutate(race = gsub('\\/', '_', race)) %>%
  dplyr::mutate_if(is.character, as.factor)



# Check for any missing data or unacceptable format

joint_table <-
  joint_table %>% dplyr::select(!where(~ any(grepl("e-", .))))

# Code for one vs rest analyses

joint_table <- joint_table %>%
  dplyr::mutate(
    'is_cluster1' = ifelse(cluster_assignment == 1, 'cluster1', 'rest'),
    'is_cluster2' = ifelse(cluster_assignment == 2, 'cluster2', 'rest'),
    'is_cluster3' = ifelse(cluster_assignment == 3, 'cluster3', 'rest'),
  ) %>%
  dplyr::mutate_if(is.character, as.factor)

# Split the dataset into train and test for model

set.seed(3456)

train_index <-
  createDataPartition(
    joint_table$cluster_assignment,
    p = 0.7,
    list = FALSE,
    times = 1
  )

train_data <- joint_table[train_index,]
test_data <- joint_table[-train_index,]

covariates <- c(
  'molecular_subtype',
  'CNS_region',
  'reported_gender',
  'race',
  'age_at_diagnosis_days'
)

pathway_names <-
  colnames(joint_table)[!names(joint_table) %in% c(
    "Kids_First_Biospecimen_ID",
    "cluster_assignment",
    "is_cluster1",
    "is_cluster2",
    "is_cluster3",
    covariates
  )]
train_formula <-
  paste0("cluster_assignment ~ ", paste(c(pathway_names, covariates), collapse = " + "))
formulas <- list(
  'one_v_rest' =
    paste0("is_cluster1 ~ ", paste(c(
      pathway_names, covariates
    ), collapse = " + ")),
  'two_v_rest' =
    paste0("is_cluster2 ~ ", paste(c(
      pathway_names, covariates
    ), collapse = " + ")),
  'three_v_rest' =
    paste0("is_cluster3 ~ ", paste(c(
      pathway_names, covariates
    ), collapse = " + "))
)

responses <- c('is_cluster1', 'is_cluster2', 'is_cluster3')

# GBMboost
train_data <- train_data %>%
  dplyr::mutate(
    is_cluster1 = case_when(is_cluster1 == 'rest' ~ 0,
                            .default = 1),
    is_cluster2 = case_when(is_cluster2 == 'rest' ~ 0,
                            .default = 1),
    is_cluster3 = case_when(is_cluster3 == 'rest' ~ 0,
                            .default = 1)
  )

test_data <- test_data %>%
  dplyr::mutate(
    is_cluster1 = case_when(is_cluster1 == 'rest' ~ 0,
                            .default = 1),
    is_cluster2 = case_when(is_cluster2 == 'rest' ~ 0,
                            .default = 1),
    is_cluster3 = case_when(is_cluster3 == 'rest' ~ 0,
                            .default = 1)
  )

for (i in 1:length(formulas)) {
  foi <- formulas[[i]]
  gbm_fit <-
    gbm(
      formula = as.formula(foi),
      data = as.data.frame(train_data),
      distribution = 'bernoulli',
      n.trees = 500,
      cv.folds = 10
    )
  
  # Train and test ROC_AUC
  train_predict <-
    predict(
      gbm_fit,
      newdata = dplyr::select(
        train_data,
        -cluster_assignment,
        -Kids_First_Biospecimen_ID,
        -is_cluster1,
        -is_cluster2,
        -is_cluster3
      ),
      type = 'response'
    )
  
  train_predict <- as.data.frame(train_predict)
  train_predict <- train_predict %>%
    dplyr::mutate('observed_labels' = train_data[[responses[i]]],
                  'group' = 'train') %>%
    dplyr::rename('probabilities' = 'train_predict')
  
  test_predict <-
    predict(
      gbm_fit,
      newdata = dplyr::select(
        test_data,
        -cluster_assignment,
        -Kids_First_Biospecimen_ID,
        -is_cluster1,
        -is_cluster2,
        -is_cluster3
      ),
      type = 'response'
    )
  
  test_predict <- as.data.frame(test_predict)
  test_predict <- test_predict %>%
    dplyr::mutate('observed_labels' = test_data[[responses[i]]],
                  'group' = 'test') %>%
    dplyr::rename('probabilities' = 'test_predict')
  
  pdf(file = file.path(
    gbm_plots_dir,
    paste0('gbm_roc_curves_', responses[i], '.pdf')
  ),
  onefile = TRUE)
  
  roc(
    train_predict$observed_labels,
    train_predict$probabilities,
    plot = TRUE,
    legacy.axes = FALSE,
    percent = TRUE,
    col = "salmon",
    lwd = 2,
    print.auc = TRUE
  )
  
  plot.roc(
    test_predict$observed_labels,
    test_predict$probabilities,
    percent = TRUE,
    col = "goldenrod",
    lwd = 2,
    print.auc = TRUE,
    add = TRUE,
    print.auc.y = 40
  )
  
  legend(
    "bottomright",
    legend = c("train", "test"),
    text.col = c("salmon", "goldenrod")
  )
  
  dev.off()
  
  var_influence <-
    relative.influence(gbm_fit, 500, scale. = FALSE, sort. = TRUE)
  var_influence <-
    data.table('variable' = names(var_influence),
               'influence_score' = var_influence)
  fwrite(var_influence,
         file = file.path(gbm_output_dir, paste0(
           'influence_score_', names(formulas)[i], '.txt'
         )),
         sep = '\t')
  var_influence_sub <- var_influence %>%
    slice_head(n = 20)
  
  bar_plot <-
    ggplot(var_influence_sub,
           aes(
             y = reorder(
               var_influence_sub$variable,
               var_influence_sub$influence_score
             ),
             x = var_influence_sub$influence_score
           )) +
    geom_bar(stat = "identity", fill = colorRampPalette(brewer.pal(3, 'Blues'))(nrow(var_influence_sub))) +
    labs(
      title = paste0("Top 20 Features by Influence Score (", names(formulas)[i], ")"),
      y = "",
      x = "Influence Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot as a PDF file in the "plots" directory
  ggsave(
    file.path(
      gbm_plots_dir,
      paste0("gbm_top_20_features_", names(formulas)[i], ".pdf")
    ),
    plot = bar_plot,
    width = 15,
    height = 6
  )
}
