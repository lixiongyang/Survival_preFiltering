
###################################################################################### #
# Goal of this script: compute the optimal thresholds for miRNA dataset
###################################################################################### #

if(!file.exists(paste0("data_fit/", cancer, "/EN/pred_df_pre_filtering_miRNA.RData"))){
  
  optimize_flt_miRNA <- optimize_preFiltering(miRNA_data, y_cox, clinical_data_srv, "EN", k_fold = K_folds, 
                                              flds = NULL, IQR_thrs_miRNA, thrs_p_val_miRNA)
  
  save(optimize_flt_miRNA,
       file =  paste0("data_fit/", cancer, "/EN/pred_df_pre_filtering_miRNA.RData"))
}

# load the data and plot the results
load(file =  paste0("data_fit/", cancer, "/EN/pred_df_pre_filtering_miRNA.RData"))

print(plot_preFiltering(optimize_flt_miRNA$C_ary, optimize_flt_miRNA$n_genes_ary))
