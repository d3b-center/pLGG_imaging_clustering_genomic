suppressPackageStartupMessages({
  library(tidyverse)
  library(GSAR)
  library(org.Hs.eg.db)
  library(EGSEA)
  library(DGCA)
})

######################### prepare for GSNCA tests-filter out lowly expressed genes 
filter_low_expr_df <- function(expr_df){
  # filter lowly expressed genes by DGSA
  expr_df_filtered <- filterGenes(expr_df, 
                                  filterTypes = c("central", "dispersion"),
                                  filterDispersionType = "cv", 
                                  filterDispersionPercentile = 0.2,
                                  sequential= TRUE)
  
  # also additionally filter on sd to avoid error message on GSNCA step
  expr_df_filtered_sd <- apply(expr_df_filtered, 1, sd, na.rm = TRUE)
  expr_df_filtered_sd_filter <- which(expr_df_filtered_sd > 0.015) %>% as.data.frame() %>% rownames()
  expr_df_filtered <- expr_df_filtered[(row.names(expr_df_filtered) %in% expr_df_filtered_sd_filter),]
  
  return(expr_df_filtered)
}

####### get the function of plotMST2.pathway and modify for return of hub genes
return_hub_gene <- function(object, group, name=NULL, cor.method="pearson", min.sd=1e-3){
  nv <- ncol(object)
  object <- object[,c(which(group == 1), which(group == 2))]
  nv1 <- sum(group == 1)
  if(length(rownames(object)) < nrow(object)) 
    gnames <- as.character(c(1:nrow(object))) else gnames <- rownames(object)
  
  objt <- aperm(object, c(2,1))
  group1 <- objt[1:nv1,]
  group2 <- objt[(nv1+1):nv,]
  
  cormat1 <- abs(cor(group1, method=cor.method))
  cormat2 <- abs(cor(group2, method=cor.method))
  e1 <- eigen(cormat1)
  e2 <- eigen(cormat2)
  p1 <- matrix(abs(e1$vectors[,1]))
  p2 <- matrix(abs(e2$vectors[,1]))
  p1 <- p1 * norm(p1)
  p2 <- p2 * norm(p2)
  colnames(p1) <- "class1"
  colnames(p2) <- "class2"
  rownames(p1) <- rownames(p2) <- gnames
  major1.val <- max(p1)
  major2.val <- max(p2)
  major1.ind <- which.max(p1)
  major2.ind <- which.max(p2)
  MST2.group1 <- findMST2(object[,c(1:nv1)], cor.method, min.sd, TRUE)
  MST2.group2 <- findMST2(object[,c((nv1+1):nv)], cor.method, min.sd, TRUE)
  
  return(gnames[major1.ind])
}

#### function to create barplots for top pathways -----------------------------
pathway_barplots <- function(dat, title){
  # calculate log score
  dat <- dat %>% 
    mutate(log_score = (-1)*log10(pvalue)) %>%
    arrange(log_score, descending = TRUE)
  
  dat[,"pathway_description"] <- factor(dat[,"pathway_description"], levels = unique(dat[,"pathway_description"]))
  
  p <- ggplot(dat, aes(x = pathway_description, 
                       y = log_score,
                       fill = log_score)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 P-Value") + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
    ggtitle(title) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    guides(fill = "none")
  
  return(p)
}


######################### Run GSNCA and output plots and text files --------------
gsnca_analysis_plot <- function(tpm_matrix, cluster_anno, pathway_df, comparison, output_file_dir, output_plot_dir, top_bar=20, top_net=5){
  # get the list of pathways 
  # pathway_list <- pathway_df %>%
  #   pull(term) %>% 
  #   unique() 
  
  # only look at pathways that has 15-500 members
  # also only look at pathways that have genes overlapping with the tpm matrix
  pathway_list <- pathway_df %>% 
    group_by(term) %>% 
    dplyr::mutate(n=n()) %>% 
    filter(n >= 15 & n <= 500,
           gene %in% rownames(tpm_matrix)) %>% 
    pull(term) %>% 
    unique()

  # define matrix to store results
  gsnca_results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(gsnca_results) <- c("pathway_description", "pvalue", "hub_gene")
  
  # results file name
  fname <- file.path(output_file_dir, paste0(comparison, "_GSNCA_anaysis.tsv"))
  if(!file.exists(fname)){
    for(i in 1:length(pathway_list)){
      # iterate through path list 
      pathway_of_interest <- pathway_list[i]
      print(i)
      
      # find genes in pathway of interest
      genes_in_pathway <- pathway_df %>% 
        filter(term == pathway_of_interest) %>% 
        pull(gene) %>% unique()
      
      # filter to genes in target pathway
      tpm_matrix_per_pathway <- tpm_matrix[row.names(tpm_matrix) %in% genes_in_pathway,,drop = F] 
      
      # run GSNCA test on filtered 
      if(nrow(tpm_matrix_per_pathway) > 2){
        result_pval<-GSNCAtest(object = as.matrix(tpm_matrix_per_pathway), 
                               # since the matrix is selected by order of row, the group will match
                               group = cluster_anno$group, 
                               nperm = 500, 
                               cor.method = "spearman", 
                               check.sd = TRUE, 
                               min.sd = 1e-3, 
                               max.skip = 10)
      } else {
        result_pval <- 1
      }
      
      
      # store the pathway name and results in the results table
      gsnca_results[i,1] <- pathway_of_interest
      gsnca_results[i,2] <- result_pval
      
      # output the hubgene for each pathway analysis
      if(nrow(tpm_matrix_per_pathway) > 2){
        gsnca_results[i,3] <- return_hub_gene(object=as.matrix(tpm_matrix_per_pathway),
                                              # since the matrix is selected by order of row, the group will match
                                              group=cluster_anno$group,
                                              cor.method="spearman")
      } else {
        gsnca_results[i,3] <- NA
      }
    }
    
    # write out results
    gsnca_results <- gsnca_results %>% 
      dplyr::filter(!is.na(pathway_of_interest)) %>%
      arrange(pvalue, descending = FALSE) %>%
      # filter to pval < 0.05
      dplyr::filter(pvalue < 0.05) %>% 
      dplyr::mutate(comparison = comparison)
    gsnca_results %>% 
      readr::write_tsv(fname)
  } else {
    # read results file instead of re-running as it takes a long time
    gsnca_results <- read_tsv(fname)
  }
  
  # we plot out top n for networks
  gsnca_top_net <- gsnca_results %>%
    slice_head(n = top_net) %>%
    as.data.frame()
  
  pdf(file = file.path(output_plot_dir, paste0(comparison, "_GSNCA_plots.pdf")), width = 15, height = 10)
  
  #### For pathways with top 5 pval, we plot out the network plots for them
  for(k in 1:nrow(gsnca_top_net)){
    # gather genes in the pathway of interest 
    pathway_network <- gsnca_top_net[k,1]
    
    # find genes in pathway of interest
    genes_in_network_pathway <- pathway_df %>% 
      filter(term == pathway_network) %>% 
      pull(gene) %>% 
      unique()
    
    # filter to genes in target pathway
    tpm_matrix_per_network_pathway <- tpm_matrix[row.names(tpm_matrix) %in% genes_in_network_pathway,]
    
    plotMST2.pathway(object = as.matrix(tpm_matrix_per_network_pathway),
                     # since the matrix is selected by order of row, the group will match
                     group = cluster_anno$group,
                     cor.method = "spearman",
                     legend.size = 0.9,
                     label.size = 1.2,
                     name = paste0(pathway_network, " ", comparison)) 
  }
  
  #### For pathways with top n pval, we plot out as bar plot
  gsnca_top_bar <- gsnca_results %>%
    slice_head(n = top_bar) %>%
    as.data.frame()
  barplot <- pathway_barplots(gsnca_top_bar, title = paste0("GSNCA: TOP ", top_bar, " ", comparison))
  print(barplot)
  
  dev.off()
}



