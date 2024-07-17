# function for doing chisq test + downstream plots
corr_plots <- function(lgg_cluster_file, prefix){
  
  # subset file for comparison at the level of cohort participant id
  cohort_level <- lgg_cluster_file %>% 
    dplyr::select(cohort_participant_id, molecular_subtype, cluster_assigned) %>% 
    unique()
  
  # chi-square test of independence across molecular_subtype and imaging_cluster
  capture.output(table(cohort_level$molecular_subtype, cohort_level$cluster_assigned), file = file.path(output_dir, paste0(prefix, "_chisq_test.txt")))
  capture.output(chisq.test(x = cohort_level$molecular_subtype, y = cohort_level$cluster_assigned), file = file.path(output_dir, paste0(prefix, "_chisq_test.txt")), append = T)

  # generate balloon plot of rows with at least 3 samples in a group
  dat <- cohort_level %>%
    filter(!is.na(molecular_subtype)) %>%
    group_by(molecular_subtype, cluster_assigned)  %>%
    summarise(n = n()) %>%
    # mutate(nmax = max(n)) %>%
    # filter(nmax >= 3) %>%
    # ungroup() %>%
    # dplyr::select(-c(nmax)) %>%
    spread(key = cluster_assigned, value = n, fill = 0) %>%
    column_to_rownames("molecular_subtype")
  dat <- dat[grep("To be classified", rownames(dat), invert = T),] # remove To be classified
  pdf(file = file.path(plots_dir, paste0(prefix, "_balloonplot.pdf")))
  balloonplot(x = as.table(as.matrix(t(dat))),
              main = "Imaging clusters vs Molecular subtypes",
              xlab = "", ylab = "",
              label = T, 
              label.size = 0.5, 
              rowmar = 4,
              show.margins = FALSE)
  dev.off()
  
  # generate corrplot of rows with at least 5 samples in a group to show pearson residuals
  chisq <- chisq.test(dat)
  pdf(file = file.path(plots_dir, paste0(prefix, "_corrplot.pdf")))
  corrplot(chisq$residuals, is.cor = FALSE, tl.srt = 360, tl.offset = 1, mar = c(1, 2, 1, 1), cl.align.text = "l")
  dev.off()
}
