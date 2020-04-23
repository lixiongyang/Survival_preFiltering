###################################################################################### #
# Goal of this script: study the stability of the prognostic indices by computing 
# several PIs for each patient with boostrapping method.
###################################################################################### #


# Compute the PIs by cross-validation and bootstrap ---

# start index of the method
ind <- 1

for(method in methods_vect_mult){
  
  if(!file.exists(paste0("data_fit/", cancer, "/", method, "/PI_bootstrap.RData"))){
    
    print(paste0("*** Start learning for method: ", method, " ***"))
    
    PI_test_list <- list()
    PI_test_list_flt <- list()
    
    for(k in 1:K_folds){
      
      print(paste0("*** Start learning for fold: ", k, " ***"))
      
      # build training and testing dataset
      id_test <- flds[[k]]
      
      gene_data_train <- gene_data[-id_test,]
      gene_data_test <- gene_data[id_test,]
      y_cox_train <- y_cox[-id_test]
      y_cox_test <- y_cox[id_test]
      
      # standardize the data
      mean_train <- apply(gene_data_train, 2, mean)
      sd_train <- apply(gene_data_train, 2, sd)
      gene_data_train_std <- sweep(gene_data_train, 2, mean_train, FUN="-")
      gene_data_train_std <- sweep(gene_data_train_std,2, sd_train, FUN="/")
      
      gene_data_test_std <- sweep(gene_data_test, 2, mean_train, FUN="-")
      gene_data_test_std <- sweep(gene_data_test_std,2, sd_train, FUN="/")
      
      id_NA_gene_train <- apply(gene_data_train_std, 2, function(x) sum(is.na(x)) >= 1)
      gene_data_train_std <- gene_data_train_std[, !id_NA_gene_train]
      gene_data_test_std <- gene_data_test_std[, !id_NA_gene_train]
      gene_data_train <- gene_data_train[,!id_NA_gene_train]
      
      # compute p-values and IQR
      p_val_cox_train <- p_val_univCox_func(gene_data_train_std, y_cox_train)
      p_val_cox_train_BH <- p.adjust(p_val_cox_train, method = "BH")
      id_NA <- is.na(p_val_cox_train_BH)
      gene_data_train_std <- gene_data_train_std[,!id_NA]
      gene_data_train <- gene_data_train[,!id_NA]
      gene_data_test_std <- gene_data_test_std[, !id_NA]
      p_val_cox_train_BH <-  p_val_cox_train_BH[!id_NA]
      
      IQR_vect_train <- apply(gene_data_train, 2, IQR)
      
      # pre-filter the data
      id_non_flt_genes <- p_val_cox_train_BH <= opt_thrs_sup[ind] & IQR_vect_train >= opt_thrs_unsup[ind]
      gene_data_train_std_flt <- gene_data_train_std[, id_non_flt_genes]
      gene_data_test_std_flt <- gene_data_test_std[, id_non_flt_genes]
      
      PI_test_df <- matrix(ncol = n_PIs, nrow = nrow(gene_data_test_std))
      row.names(PI_test_df) <- row.names(gene_data_test_std)
      
      PI_test_df_flt <- matrix(ncol = n_PIs, nrow = nrow(gene_data_test_std_flt))
      row.names(PI_test_df_flt) <- row.names(gene_data_test_std_flt)
      
      for(i in 1:n_PIs){
        
        # bootstrap training datasets
        id_bootstrap <- sample(1:nrow(gene_data_train_std_flt), replace = T)
        gene_data_train_bootstrap_flt <- gene_data_train_std_flt[id_bootstrap, ]
        gene_data_train_bootstrap <- gene_data_train_std[id_bootstrap, ]
        y_cox_train_bootstrap <- y_cox_train[id_bootstrap]
        
        # learn models
        res <- try(cv_fit_train_bootstrap <- learn_models(gene_data_train_bootstrap, 
                                                          y_cox_train_bootstrap, method, 1)[[1]])
        
        if(inherits(res, "try-error")){next}
        res <- try(cv_fit_train_bootstrap_flt <- learn_models(gene_data_train_bootstrap_flt, 
                                                              y_cox_train_bootstrap, method, 1)[[1]])
        if(inherits(res, "try-error")){next}
        
        # compute the PIs on the testing dataset
        names_genes_selected <- genes_selected_func(cv_fit_train_bootstrap, "lambda.min")
        if(length(names_genes_selected) >= 1){
          # compute the C-index
          beta <- as.numeric(coef(cv_fit_train_bootstrap, "lambda.min"))
          names(beta) <- colnames(gene_data_train_std) 
          
          PI_test_df[,i] <-  gene_data_test_std %*% beta 
          
        }
        
        names_genes_selected_flt <- genes_selected_func(cv_fit_train_bootstrap_flt, "lambda.min")
        if(length(names_genes_selected_flt) >= 1){
          # compute the C-index
          beta_flt <- as.numeric(coef(cv_fit_train_bootstrap_flt, "lambda.min"))
          names(beta_flt) <- colnames(gene_data_train_std_flt) 
          
          PI_test_df_flt[,i] <-  gene_data_test_std_flt %*% beta_flt
          
        }
        
        if(i%%5 == 0){
          print(paste0("Number of models learned with boostrap: ",i))
        }
      }
      
      PI_test_list[[k]] <- PI_test_df
      PI_test_list_flt[[k]] <- PI_test_df_flt
      
    }
    
    save(PI_test_list_flt, PI_test_list, file = paste0("data_fit/", cancer, "/", method, "/PI_bootstrap.RData"))
  }
  
  ind <- ind + 1
}

