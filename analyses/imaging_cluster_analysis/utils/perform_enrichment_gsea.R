# identify differentially expressed genes and perform pathway enrichment for input deseq2 results

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(ggplot2)
  library(ggpubr)
})

perform_enrichment_gsea <- function(diffexpr_res, pathways, prefix, plots_dir, results_dir){
  
  # output directory
  dir.create(plots_dir, recursive = T, showWarnings = F)
  dir.create(results_dir, recursive = T, showWarnings = F)

  # perform enrichment using clusterProfiler
  ranks <- diffexpr_res %>%
    arrange(desc(log2FC)) %>%
    pull(log2FC)
  names(ranks) <- diffexpr_res$gene
  path_res <- clusterProfiler::GSEA(
    geneList = ranks,
    minGSSize = 1,
    maxGSSize = length(ranks)-1,
    TERM2GENE = pathways
  )
  write_tsv(x = path_res %>% as.data.frame(), file = file.path(results_dir, paste0(prefix, "_gsea.tsv")))
  
  # barplot of top 10 significant pathways (padj < 0.05)
  # top 10 upregulated
  top10_up <- path_res %>% 
    as.data.frame() %>%
    filter(p.adjust < 0.05, 
           NES > 0) %>%
    mutate(direction = "Up") %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
  
  # top 10 downregulated
  top10_down <- path_res %>% 
    as.data.frame() %>% 
    filter(p.adjust < 0.05,
           NES < 0) %>%
    mutate(direction = "Down") %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
  
  # combine both pathways
  top10_output <- rbind(top10_up, top10_down)
  top10_output <- top10_output %>% 
    mutate(ID = gsub("_", " ", ID))
  
  # generate plots
  if(nrow(top10_output) > 0){
    # barplot
    pdf(file = file.path(plots_dir, paste0(prefix, "_gsea_barplot.pdf")), width = 8, height = 8)
    top10_output$ID <- factor(top10_output$ID, levels = unique(top10_output$ID))
    top10_output$direction <- factor(top10_output$direction, levels = c("Up", "Down"))
    
    p <- ggplot(top10_output, aes(ID, y = (-1)*log10(p.adjust), fill = direction)) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw() +
      xlab("") + 
      ylab("-log10 Adj. P-Value") + 
      scale_fill_manual(name = "Direction", values = c("Down" = "forest green", "Up" = "red")) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
      ggpubr::theme_pubr(base_size = 8) + 
      ggtitle("Top 10 Upregulated/Downregulated Pathways")
    print(p)
    dev.off()
  }
}
