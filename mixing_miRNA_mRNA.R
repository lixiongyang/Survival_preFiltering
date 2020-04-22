###################################################################################### #
# Goal of this script: Compute the C-index without pre-filtering and in the optimal
# case for mRNA, miRNA and mixed datasets.
# The mRNA and miRNA datasets are integrated with the mRNA and the miRNA selected
# by the elastic net independantly.
###################################################################################### #


if(!file.exists(paste0("data_fit/", cancer, "/mixOmics.RData"))){
  
  # matrix with the C-index and IBS - miRNA + mRNA
  C_df_flt <- matrix(nrow = n_rep, ncol = K_folds)
  C_df <-  matrix(nrow = n_rep, ncol = K_folds)
  
  IBS_df_flt <- matrix(nrow = n_rep, ncol = K_folds)
  IBS_df <-  matrix(nrow = n_rep, ncol = K_folds)
  
  
  # matrix with the C-index and IBS - mRNA
  C_df_flt_mRNA <- matrix(nrow = n_rep, ncol = K_folds)
  C_df_mRNA <-  matrix(nrow = n_rep, ncol = K_folds)
  
  IBS_df_flt_mRNA <- matrix(nrow = n_rep, ncol = K_folds)
  IBS_df_mRNA <-  matrix(nrow = n_rep, ncol = K_folds)
  
  
  # matrix with the C-index and IBS - miRNA
  C_df_flt_miRNA <- matrix(nrow = n_rep, ncol = K_folds)
  C_df_miRNA <-  matrix(nrow = n_rep, ncol = K_folds)
  
  IBS_df_flt_miRNA <- matrix(nrow = n_rep, ncol = K_folds)
  IBS_df_miRNA <-  matrix(nrow = n_rep, ncol = K_folds)
  
  
  for(i in 1: n_rep){
    
    print(paste0("*** Start learning for repetition number: ", i, " ***"))
    
    flds <- createFolds(1:nrow(mRNA_data), k = K_folds, list = TRUE, returnTrain = FALSE)
    
    for(k in 1:K_folds){
      
      
      print(paste0("*** Start learning for fold: ", k, " ***"))
      
      # build training and testing dataset
      id_test <- flds[[k]]
      
      mRNA_data_train <- mRNA_data[-id_test,]
      mRNA_data_test <- mRNA_data[id_test,]
      miRNA_data_train <- miRNA_data[-id_test,]
      miRNA_data_test <- miRNA_data[id_test,]
      y_cox_train <- y_cox[-id_test]
      y_cox_test <- y_cox[id_test]
      
      # pre-filter the data
      mRNA_data_train_flt <- pre_filter_func(mRNA_data_train, y_cox_train, thrs_sup, thrs_unsup)
      mRNA_data_test_flt <- mRNA_data_test[, colnames(mRNA_data_train_flt)]
      miRNA_data_train_flt <- pre_filter_func(miRNA_data_train, y_cox_train, thrs_sup_miRNA, thrs_unsup_miRNA)
      miRNA_data_test_flt <- miRNA_data_test[, colnames(miRNA_data_train_flt)]
      
      # standardize the data ---
      # mRNA - with filtering
      mean_train <- apply(mRNA_data_train_flt, 2, mean)
      sd_train <- apply(mRNA_data_train_flt, 2, sd)
      mRNA_data_train_flt_std <- sweep(mRNA_data_train_flt, 2, mean_train, FUN="-")
      mRNA_data_train_flt_std <- sweep(mRNA_data_train_flt_std,2, sd_train, FUN="/")
      mRNA_data_test_flt_std <- sweep(mRNA_data_test_flt, 2, mean_train, FUN="-")
      mRNA_data_test_flt_std <- sweep(mRNA_data_test_flt_std,2, sd_train, FUN="/")
      
      # miRNA - with filtering
      mean_train <- apply(miRNA_data_train_flt, 2, mean)
      sd_train <- apply(miRNA_data_train_flt, 2, sd)
      miRNA_data_train_flt_std <- sweep(miRNA_data_train_flt, 2, mean_train, FUN="-")
      miRNA_data_train_flt_std <- sweep(miRNA_data_train_flt_std,2, sd_train, FUN="/")
      miRNA_data_test_flt_std <- sweep(miRNA_data_test_flt, 2, mean_train, FUN="-")
      miRNA_data_test_flt_std <- sweep(miRNA_data_test_flt_std,2, sd_train, FUN="/")
      
      # mRNA - without filtering
      mean_train <- apply(mRNA_data_train, 2, mean)
      sd_train <- apply(mRNA_data_train, 2, sd)
      mRNA_data_train_std <- sweep(mRNA_data_train, 2, mean_train, FUN="-")
      mRNA_data_train_std <- sweep(mRNA_data_train_std,2, sd_train, FUN="/")
      mRNA_data_test_std <- sweep(mRNA_data_test, 2, mean_train, FUN="-")
      mRNA_data_test_std <- sweep(mRNA_data_test_std,2, sd_train, FUN="/")
      
      id_NA_gene_train <- apply(mRNA_data_train_std, 2, function(x) sum(is.na(x)) >= 1)
      mRNA_data_train_std <- mRNA_data_train_std[, !id_NA_gene_train]
      mRNA_data_test_std <- mRNA_data_test_std[, !id_NA_gene_train]
      
      # miRNA - without filtering
      mean_train <- apply(miRNA_data_train, 2, mean)
      sd_train <- apply(miRNA_data_train, 2, sd)
      miRNA_data_train_std <- sweep(miRNA_data_train, 2, mean_train, FUN="-")
      miRNA_data_train_std <- sweep(miRNA_data_train_std,2, sd_train, FUN="/")
      miRNA_data_test_std <- sweep(miRNA_data_test, 2, mean_train, FUN="-")
      miRNA_data_test_std <- sweep(miRNA_data_test_std,2, sd_train, FUN="/")
      
      dim(mRNA_data_train_flt_std)
      dim(mRNA_data_test_flt_std)
      dim(mRNA_data_train_std)
      dim(mRNA_data_test_std)
      
      dim(miRNA_data_train_flt_std)
      dim(miRNA_data_test_flt_std)
      dim(miRNA_data_train_std)
      dim(miRNA_data_test_std)
      
      # select genes and compute the C-index on signle omics ---
      # mRNA with filtering
      res <- try(fit_mRNA_flt <- learn_models(mRNA_data_train_flt_std, y_cox_train, "EN", 1)[[1]])
      if(inherits(res, "try-error")){next}
      genes_mRNA_flt <- genes_selected_func(fit_mRNA_flt, "lambda.min")
      if(length(genes_mRNA_flt) >= 1){
        beta <- as.numeric(coef(fit_mRNA_flt, "lambda.min"))
        names(beta) <- colnames(mRNA_data_train_flt_std) 
        PI_test <-  mRNA_data_test_flt_std %*% beta
        PI_train <- mRNA_data_train_flt_std %*% beta
        C_df_flt_mRNA[i, k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                 surv.event = clinical_data_srv[id_test,"status"])$c.index 
        IBS_df_flt_mRNA[i, k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                         type = "brier")$ierror
      }
      
      # miRNA with filtering
      res <- try(fit_miRNA_flt <- learn_models(miRNA_data_train_flt_std, y_cox_train, "EN", 1)[[1]])
      if(inherits(res, "try-error")){next}
      genes_miRNA_flt <- genes_selected_func(fit_miRNA_flt, "lambda.min")
      if(length(genes_miRNA_flt) >= 1){
        beta <- as.numeric(coef(fit_miRNA_flt, "lambda.min"))
        names(beta) <- colnames(miRNA_data_train_flt_std) 
        PI_test <-  miRNA_data_test_flt_std %*% beta
        PI_train <- miRNA_data_train_flt_std %*% beta
        C_df_flt_miRNA[i,k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                 surv.event = clinical_data_srv[id_test,"status"])$c.index 
        IBS_df_flt_miRNA[i,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                         type = "brier")$ierror
      }
      
      # mRNA without filtering
      res <- try(fit_mRNA <- learn_models(mRNA_data_train_std, y_cox_train, "EN", 1)[[1]])
      if(inherits(res, "try-error")){next}
      genes_mRNA <- genes_selected_func(fit_mRNA, "lambda.min")
      if(length(genes_mRNA) >= 1){
        beta <- as.numeric(coef(fit_mRNA, "lambda.min"))
        names(beta) <- colnames(mRNA_data_train_std) 
        PI_test <-  mRNA_data_test_std %*% beta
        PI_train <- mRNA_data_train_std %*% beta
        C_df_mRNA[i,k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                            surv.event = clinical_data_srv[id_test,"status"])$c.index 
        IBS_df_mRNA[i,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                    type = "brier")$ierror
      }
      
      # miRNA without filtering
      res <- try(fit_miRNA <- learn_models(miRNA_data_train_std, y_cox_train, "EN", 1)[[1]])
      if(inherits(res, "try-error")){next}
      genes_miRNA <- genes_selected_func(fit_miRNA, "lambda.min")
      if(length(genes_miRNA) >= 1){
        beta <- as.numeric(coef(fit_miRNA, "lambda.min"))
        names(beta) <- colnames(miRNA_data_train_std) 
        PI_test <-  miRNA_data_test_std %*% beta
        PI_train <- miRNA_data_train_std %*% beta
        C_df_miRNA[i,k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                             surv.event = clinical_data_srv[id_test,"status"])$c.index 
        IBS_df_miRNA[i,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                     type = "brier")$ierror
      }
      
      # learn models and compute the C-index on multi-omics ---
      # with pre-filtering
      gene_data_train_mix_flt <- cbind(mRNA_data_train_flt_std[, genes_mRNA_flt], miRNA_data_train_flt_std[, genes_miRNA_flt])
      gene_data_test_mix_flt <- cbind(mRNA_data_test_flt_std[, genes_mRNA_flt], miRNA_data_test_flt_std[, genes_miRNA_flt])
      
      if(ncol(gene_data_train_mix_flt) <= 1 ){next}
      
      res <- try(fit_mix_flt <- learn_models(gene_data_train_mix_flt, y_cox_train, "ridge", 1)[[1]])
      if(!inherits(res, "try-error")){
        if(length(genes_selected_func(fit_mix_flt, "lambda.min")) >= 1){
          beta <- as.numeric(coef(fit_mix_flt, "lambda.min"))
          names(beta) <- colnames( gene_data_train_mix_flt) 
          PI_test_flt <-  gene_data_test_mix_flt %*% beta
          PI_train_flt <- gene_data_train_mix_flt %*% beta
          C_df_flt[i,k] <- concordance.index(PI_test_flt, surv.time = OS_tcga_df[id_test,"time"], 
                                             surv.event = OS_tcga_df[id_test,"status"])$c.index 
          IBS_df_flt[i,k] <- predErr(y_cox_train, y_cox_test, PI_train_flt, PI_test_flt, times = seq(0, 3, by=1/365.25), 
                                     type = "brier")$ierror
        }
        
      }
      
      # without pre-filtering
      gene_data_train_mix <- cbind(mRNA_data_train_std[, genes_mRNA], miRNA_data_train_std[, genes_miRNA])
      gene_data_test_mix <- cbind(mRNA_data_test_std[, genes_mRNA], miRNA_data_test_std[, genes_miRNA])
      
      if(ncol(gene_data_train_mix) <= 1 ){next}
      
      res <- try(fit_mix <- learn_models(gene_data_train_mix, y_cox_train, "ridge", 1)[[1]])
      if(!inherits(res, "try-error")){
        if(length(genes_selected_func(fit_mix, "lambda.min")) >= 1){
          beta <- as.numeric(coef(fit_mix, "lambda.min"))
          names(beta) <- colnames( gene_data_train_mix) 
          PI_test <-  gene_data_test_mix %*% beta
          PI_train <-  gene_data_train_mix %*% beta
          C_df[i,k] <- concordance.index(PI_test, surv.time = OS_tcga_df[id_test,"time"], 
                                         surv.event = OS_tcga_df[id_test,"status"])$c.index 
          IBS_df[i,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), 
                                 type = "brier")$ierror
        }
      }
    }
  }
  
  
  save(C_df, C_df_flt, IBS_df, IBS_df_flt,
       C_df_mRNA, C_df_flt_mRNA, IBS_df_mRNA, IBS_df_flt_mRNA,
       C_df_miRNA, C_df_flt_miRNA, IBS_df_miRNA, IBS_df_flt_miRNA,
       file = paste0("data_fit/", cancer, "/mixOmics.RData"))
}



