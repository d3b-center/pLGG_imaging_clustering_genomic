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
dir.create(output_dir,
           showWarnings = F,
           recursive = T)

# plots dirs
plot_dir <-
  file.path(analysis_dir, "plots", "imaging_cluster_prediction")
dir.create(plot_dir, showWarnings = F, recursive = T)

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
  dplyr::mutate(race = ifelse(race == 'Reported Unknown', ethnicity, race)) %>%
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
ssgsea_scores$REACTOME_REPRODUCTION <- NULL

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
  'age_at_diagnosis_days',
  'race'
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

get_best_result <- function(caret_fit) {
  best <-
    which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result <- caret_fit$results[best,]
  rownames(best_result) <- NULL
  best_result
}


for (i in 1:length(formulas)) {
  foi <- formulas[[i]]
  ctr_cv <-
    caret::trainControl(method = "cv",
                        allowParallel = TRUE,
                        number = 10)
  
  # Elastic net log reg
  fit_glmnet <- caret::train(
    form = as.formula(foi),
    data = as.data.frame(train_data),
    method = "glmnet",
    trControl = ctr_cv,
    tuneLength = 50,
    search = 'random'
  )
  
  ### write out accuracy metrics and tuning plots
  best_result <- get_best_result(fit_glmnet)
  fwrite(best_result,
         file = file.path(
           output_dir,
           paste0('elnet_estimated_test_err_cv_', names(formulas)[i], '.txt')
         ),
         sep = '\t')
  
  p <- ggplot(data = fit_glmnet) +
    theme(legend.position = "none") +
    theme_classic()
  ggsave(
    filename = file.path(plot_dir, paste0(
      'elnet_metrics_', names(formulas)[i], '.pdf'
    )),
    plot = p,
    width = 15
  )
  
  ### Predict on training data
  train_predict <- fit_glmnet %>% predict(train_data)
  confusion_mat <-
    confusionMatrix(data = train_predict,
                    reference = train_data[[responses[i]]],
                    mode = "prec_recall")
  confusion_mat_overall <- confusion_mat$overall
  write.table(
    confusion_mat_overall,
    file = file.path(
      output_dir,
      paste0("elnet_CM_overall_train_", names(formulas)[i], ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
  
  confusion_mat_class <- confusion_mat$byClass
  write.table(
    confusion_mat_class,
    file = file.path(
      output_dir,
      paste0("elnet_CM_class_train_", names(formulas)[i], ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
  
  # Make predictions on the test data
  test_predict <- fit_glmnet %>% predict(test_data)
  confusion_mat <-
    confusionMatrix(data = test_predict,
                    reference = test_data[[responses[i]]],
                    mode = "prec_recall")
  confusion_mat_overall <- confusion_mat$overall
  write.table(
    confusion_mat_overall,
    file = file.path(
      output_dir,
      paste0("elnet_CM_overall_test_", names(formulas)[i], ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
  
  confusion_mat_class <- confusion_mat$byClass
  write.table(
    confusion_mat_class,
    file = file.path(
      output_dir,
      paste0("elnet_CM_class_test_", names(formulas)[i], ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
  
  # Important features
  
  # Print the random forest model summary
  glm_features <- coef(fit_glmnet$finalModel, fit_glmnet$bestTune$lambda)
  glm_features <- data.table('pathways' = glm_features@Dimnames[[1]][glm_features@i],
                             'coefficients' = glm_features@x[2:length(glm_features@x)]) %>%
    dplyr::arrange(desc(coefficients))
  
  write.table(
    glm_features,
    file = file.path(
      output_dir,
      paste0('elnet_full_feature_table_', names(formulas)[i], '.txt')
    ),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
  )
  
  glm_features_up <- glm_features %>% 
    dplyr::filter(coefficients > 0) %>%
    dplyr::slice_head(n = 10)
  
  glm_features_down <- glm_features %>% 
    dplyr::filter(coefficients < 0) %>%
    dplyr::slice_tail(n = 10)
  
  glm_features_filt <- rbind(glm_features_up, glm_features_down)
  
  bar_plot <-
    ggplot(glm_features_filt, aes(y = reorder(pathways, coefficients), x = coefficients)) +
    geom_bar(stat = "identity", fill = colorRampPalette(rev(brewer.pal(3, 'Blues')))(nrow(glm_features_filt))) +
    labs(title = "Top Features by GLM Coefficient",
         y = "Features",
         x = "GLM Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 18)) + 
    labs(title = '')
  
  # Save the plot as a PDF file in the "plots" directory
  ggsave(
    file.path(
      plot_dir,
      paste0("elnet_top_features_", names(formulas)[i], ".pdf")
    ),
    plot = bar_plot,
    width = 15,
    height = 9
  )
  
  # Training and testing ROC_AUC
  train_predict <-
    predict(
      fit_glmnet,
      newdata = dplyr::select(
        train_data,
        -cluster_assignment,
        -Kids_First_Biospecimen_ID,
        -is_cluster1,
        -is_cluster2,
        -is_cluster3
      ),
      type = 'prob'
    )
  
  train_predict <- as.data.frame(train_predict) %>%
    dplyr::select(-rest)
  train_predict <- train_predict %>%
    dplyr::mutate('observed_labels' = train_data[[responses[i]]],
                  'group' = 'train')
  
  test_predict <-
    predict(
      fit_glmnet,
      newdata = dplyr::select(
        test_data,
        -cluster_assignment,
        -Kids_First_Biospecimen_ID,
        -is_cluster1,
        -is_cluster2,
        -is_cluster3
      ),
      type = 'prob'
    )
  
  test_predict <- as.data.frame(test_predict) %>%
    dplyr::select(-rest)
  test_predict <- test_predict %>%
    dplyr::mutate('observed_labels' = test_data[[responses[i]]],
                  'group' = 'test')
  
  colnames(train_predict) <-
    c('probabilities', 'observed_labels', 'group')
  colnames(test_predict) <-
    c('probabilities', 'observed_labels', 'group')
  
  pdf(file = file.path(
    plot_dir,
    paste0('glm_roc_curves_', responses[i], '.pdf')
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
  
  sink(file = file.path(output_dir, paste0('ci_auc_', responses[i], '.txt')))
  
  cat('Training: AUC Confidence Intervals \n')
  ci.auc(roc(
    train_predict$observed_labels,
    train_predict$probabilities,
    method = 'b'
  ))
  
  cat('\n Testing: AUC Confidence Intervals \n')
  
  ci.auc(roc(
    test_predict$observed_labels,
    test_predict$probabilities,
    method = 'b'
  ))
  
  sink()
}
