# Functions for generating count and enrichment heatmaps
#
# Ryan Corbett
#
# 2024

# Function to create frequency count and enrichment heatmap for two selected variables. Significance of enrichment is calculated using hypergeometric tests

plot_enr <- function(df, var1, var2,
                     var1_names, var2_names,
                     padjust = FALSE) {
  
  enr <- matrix(0, nrow(unique(df[,var1])),
                        nrow(unique(df[,var2])),
                dimnames = list(var1_names,
                                var2_names))
  pval <- enr
  ct <- enr
  
  # loop through cancer groups to calculate enrichment 
  for (i in 1:nrow(enr)){
    no_var1 <- sum(unlist(df[,var1]) == rownames(enr)[i])
    for (j in 1:ncol(enr)){
      no_var2 <- sum(unlist(df[,var2]) == colnames(enr)[j] & !is.na(unlist(df[,var2])))
      no_var1_var2 <- sum(unlist(df[,var1]) == rownames(enr)[i] & unlist(df[,var2]) == colnames(enr)[j])
      ct[i,j] <- no_var1_var2
      enr[i,j] <- (no_var1_var2/no_var2)/(no_var1/nrow(df))
      pval[i,j] <- phyper(no_var1_var2, no_var1, nrow(df) - no_var1, no_var2, lower.tail = F)
    }
  }
  
  if (padjust == TRUE) {
    
    fdr <- t(apply(pval, 1, function(x) p.adjust(x, "fdr")))
    sig_mat <- ifelse(fdr < 0.05 & enr > 1 & ct > 1, "*", "")
    
  } else {
    
    sig_mat <- ifelse(pval < 0.05 & enr > 1 & ct > 1, "*", "")
    
  }
  
  fill_mat <- matrix(glue::glue("{round(enr, 1)}{sig_mat}"), 
                     nrow(enr), ncol(enr))
  
  ct_enr_mat <- matrix(glue::glue("{ct}\n({fill_mat})"),
                       nrow(ct), ncol(ct))
  
  col_fun = colorRamp2(c(0, ceiling(max(enr))), c("white", "orangered"))
  
  region_ht <- Heatmap(enr,
                       name = "Odds ratio",
                       cluster_rows = F,
                       cluster_columns = F,
                       rect_gp = gpar(col = "black", lwd = 2),
                       col = col_fun,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%s", ct_enr_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                       },
                       column_names_side = "top",
                       column_names_rot = 25,
                       column_names_centered = FALSE)
  
  return(region_ht)
  
}
