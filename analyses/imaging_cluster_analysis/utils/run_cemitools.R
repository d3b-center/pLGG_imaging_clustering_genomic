suppressPackageStartupMessages({
  library(sva)
  library(CEMiTool)
})

# run CEMItools
run_cemitools_functions <- function(expr_df, annot_df, 
                                    cor_method, network_type, tom_type, n = 100, 
                                    output_dir, plots_dir, gmt_file){

  # if no values are provided, then use default values
  if(is.null(cor_method)){
    cor_method = "pearson"
  }
  if(is.null(network_type)){
    network_type = "signed"
  }
  if(is.null(tom_type)){
    tom_type = "signed"
  }
  
  # run cemitool
  cem <- cemitool(expr = expr_df, 
                  annot = annot_df, 
                  filter = T, 
                  cor_function = "bicor", 
                  cor_method = cor_method,
                  network_type = network_type,
                  tom_type = tom_type,
                  sample_name_column = "Kids_First_Biospecimen_ID",
                  class_column = "cluster_assigned",
                  merge_similar = T,
                  apply_vst = T, 
                  verbose = F)
  
  hubs <- get_hubs(cem, n, method = "kME")
  summary <- mod_summary(cem)
  
  # save hubs output
  saveRDS(hubs, file = file.path(output_dir, "hubs.rds"))
  
  # generate heatmap of gene set enrichment analysis
  cem <- mod_gsea(cem)
  cem <- plot_gsea(cem)
  
  # plot gene expression within each module
  cem <- plot_profile(cem)
  
  # read GMT file - reactome file
  # gmt_in <- read_gmt(gmt_file)
  gmt_in <- gmt_file
  
  # perform over representation analysis
  cem <- mod_ora(cem, gmt_in)
  
  # plot ora results
  cem <- plot_ora(cem)
  
  # read interactions
  int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
  int_df <- read.delim(int_fname)
  
  # plot interactions
  interactions_data(cem) <- int_df # add interactions
  cem <- plot_interactions(cem) # generate plot
  
  # diagnostic_report
  CEMiTool::diagnostic_report(cem = cem, directory = output_dir, force = T)
  
  # generate report
  CEMiTool::generate_report(cem = cem, directory = output_dir, force = T)
  
  # write all files
  write_files(cem = cem, directory = output_dir, force = T)
  
  # save plots
  save_plots(cem = cem, "all", directory = plots_dir, force = T)
}
