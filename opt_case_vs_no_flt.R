###################################################################################### #
# Goal of this script: Compare the C-indices and the IBS obtained without
# pre-filtering and in the optimal case.
###################################################################################### #

if(!file.exists(paste0("data_fit/", cancer, "/best_vs_no_flt.RData"))){
  
  # data frames that will contain the results
  C_df_flt = C_df <- matrix(nrow = n_rep*K_folds, ncol = length(methods_vect))
  colnames(C_df_flt) = colnames(C_df) = methods_vect
  IBS_df_flt = IBS_df <- matrix(nrow = n_rep*K_folds, ncol = length(methods_vect))
  colnames(IBS_df_flt) = colnames(IBS_df) = methods_vect
  
  ind <- 1
  
  for(i in 1:n_rep){
    
    print(paste0("*** Start learning for repetition number: ", i, " ***"))
    
    # create the folds
    flds_tmp <- createFolds(1:nrow(gene_data), k = K_folds, list = TRUE, returnTrain = FALSE)
    
    for(k in 1:K_folds){
      
      print(paste0("Start learning for fold: ", k, " ***"))
      
      # build training and testing dataset
      id_test <- flds_tmp[[k]]
      
      # training and testing dataset
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
      gene_data_train <- gene_data_train[, !id_NA_gene_train]
      
      # compute the IQR and univariate p-value
      p_val_cox_train <- p_val_univCox_func(gene_data_train_std, y_cox_train)
      p_val_cox_train_BH <- p.adjust(p_val_cox_train, method = "BH")
      id_NA <- is.na(p_val_cox_train_BH)
      gene_data_train_std <- gene_data_train_std[,!id_NA]
      gene_data_train <- gene_data_train[,!id_NA]
      gene_data_test_std <- gene_data_test_std[, !id_NA]
      p_val_cox_train_BH <-  p_val_cox_train_BH[!id_NA]
      
      IQR_vect_train <- apply(gene_data_train, 2, IQR)
      
      # learn models for each method
      for(j in 1:length(methods_vect)){
        
        method <- methods_vect[j]
        print(paste0("Start learning for: ", method))
        
        # pre-filter the data
        id_non_flt_genes <- p_val_cox_train_BH <= opt_thrs_sup[j] & IQR_vect_train >= opt_thrs_unsup[j]
        gene_data_train_flt_std <- gene_data_train_std[, id_non_flt_genes]
        gene_data_test_flt_std <- gene_data_test_std[, id_non_flt_genes]
        
        dim(gene_data_train_flt_std)
        dim(gene_data_test_flt_std)
        dim(gene_data_train_std)
        dim(gene_data_test_std)
        
        if(method != "univCox"){
          
          # with pre-filtering
          res <- try(fit_flt <- learn_models(gene_data_train_flt_std, y_cox_train, method, 1)[[1]])
          if(!inherits(res, "try-error")){
            if(length(genes_selected_func(fit_flt, "lambda.min")) >= 1){
              beta_flt <- as.numeric(coef(fit_flt, "lambda.min"))
              names(beta_flt) <- colnames(gene_data_train_flt_std) 
              PI_test_flt <-  gene_data_test_flt_std %*% beta_flt
              PI_train_flt <- gene_data_train_flt_std %*% beta_flt
              C_df_flt[ind,method] <- concordance.index(PI_test_flt, surv.time = clinical_data_srv[id_test,"time"], 
                                                        surv.event = clinical_data_srv[id_test,"status"])$c.index 
              IBS_df_flt[ind,method] <- predErr(y_cox_train, y_cox_test, PI_train_flt, PI_test_flt, times = seq(0, 3, by=1/365.25), 
                                                type = "brier")$ierror
            }
            
          }
          
          
          # without pre-filtering
          res <- try(fit <- learn_models(gene_data_train_std, y_cox_train, method, 1)[[1]])
          if(!inherits(res, "try-error")){
            if(length(genes_selected_func(fit, "lambda.min")) >= 1){
              beta <- as.numeric(coef(fit, "lambda.min"))
              names(beta) <- colnames(gene_data_train_std) 
              PI_test <-  gene_data_test_std %*% beta
              PI_train <- gene_data_train_std %*% beta
              C_df[ind,method] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                    surv.event = clinical_data_srv[id_test,"status"])$c.index 
              IBS_df[ind,method] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                            type = "brier")$ierror
            }
            
          }
          
          
        }else{ # univariate Cox
          
          # with pre-filtering
          res <- try(cv_p_val_univCox_flt <- p_val_univCox_func(gene_data_train_flt_std, y_cox_train))
          
          if(all(is.na(cv_p_val_univCox_flt))){
            C_df_flt[ind, method] <- 0.5
            IBS_df_flt[ind, method] <- 0.25
          }else if(is.null(cv_p_val_univCox_flt)){
            C_df_flt[ind, method] <- 0.5
            IBS_df_flt[ind, method] <- 0.25
          }else if(!inherits(res, "try-error")){
            # gene with minimum p-value
            gene_min_p_val_flt <- names(cv_p_val_univCox_flt)[which.min(cv_p_val_univCox_flt)]
            
            res <- try(fit_flt <- coxph(y_cox_train~gene_data_train_flt_std[, gene_min_p_val_flt]))
            if(inherits(res, "try-error")){next}
            
            # compute the concordance on the testing data with this gene
            PI_test_flt <- gene_data_test_flt_std[, gene_min_p_val_flt] * fit_flt$coefficients
            PI_train_flt <- gene_data_train_flt_std[, gene_min_p_val_flt] * fit_flt$coefficients
            
            C_df_flt[ind, method] <- concordance.index(PI_test_flt, surv.time = clinical_data_srv[id_test,"time"], 
                                                       surv.event = clinical_data_srv[id_test,"status"])$c.index 
            IBS_df_flt[ind, method] <- predErr(y_cox_train, y_cox_test, PI_train_flt, PI_test_flt, times = seq(0, 3, by=1/365.25), #times,
                                               type = "brier")$ierror
          }else{
            C_df[ind, method] <- NA
            IBS_df[ind, method] <- NA
          }
          
          # without pre-filtering
          res <- try(cv_p_val_univCox <- p_val_univCox_func(gene_data_train_std, y_cox_train))
          
          if(all(is.na(cv_p_val_univCox))){
            C_df[ind, method] <- 0.5
            IBS_df[ind, method] <- 0.25
          }else if(is.null(cv_p_val_univCox)){
            C_df[ind, method] <- 0.5
            IBS_df[ind, method] <- 0.25
          }else if(!inherits(res, "try-error")){
            # gene with minimum p-value
            gene_min_p_val <- names(cv_p_val_univCox)[which.min(cv_p_val_univCox)]
            
            res <- try(fit <- coxph(y_cox_train~gene_data_train_std[, gene_min_p_val]))
            if(inherits(res, "try-error")){next}
            
            # compute the concordance on the testing data with this gene
            PI_test <- gene_data_test_std[, gene_min_p_val] * fit$coefficients
            PI_train <- gene_data_train_std[, gene_min_p_val] * fit$coefficients
            
            C_df[ind, method] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                   surv.event = clinical_data_srv[id_test,"status"])$c.index 
            IBS_df[ind, method] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), #times,
                                           type = "brier")$ierror
          }else{
            C_df[ind, method] <- NA
            IBS_df[ind, method] <- NA
          }
          
        }
        
      }
      
      ind <- ind + 1
    }
  }
  
  
  save(C_df, C_df_flt, IBS_df, IBS_df_flt,
       file = paste0("data_fit/", cancer, "/best_vs_no_flt.RData"))
}






