
# @description: compute C-indices by K-fold cross-validation for each pair of
#               thresholds (supervised, unsupervised)
#
# @parameters: 
#   - data: genetic data (mRNA or miRNA dataset)
#   - y_cox: Surv object containing follow-up times and status
#   - method: "EN" for elastic net, "lasso" for lasso, "ridge" for ridge and "AEN"
#             for adaptive elastic net
#   - k_fold: number of folds for the cross-validation
#   - flds: folds constructed externally (if needed)
#   - IQR_thr: thresholds on the IQR (unsupervised)
#   - thrs_p_val: thresholds on the p-values of the univariate Cox (supervised)
#
# @return: all the arrays returned have length(thrs_p_val) rows, length(IQR_thrs) columns,
#          and a third dimension corresponding to the number of folds.
#   - C_ary: array containing the C-indices
#   - IBS_ary: array containing the Integrated Brier Score (IBS)
#   - p_val_ary: p-value of the logrank test between patients with low PIs, 
#                (lower than the median) and patients with high PIs.
#   - n_genes_ary: number of genes retained after the pre-filtering
optimize_preFiltering <- function(data, y_cox, method, k_fold = 5, flds = NULL, 
                                  IQR_thrs = 0:4, 
                                  thrs_p_val = c(0.01, 0.05, 0.1, 0.5, 1)){
  
  
  # create folds
  if(is.null(flds)){
    flds <- createFolds(1:nrow(data), k = k_fold, list = TRUE, returnTrain = FALSE)
  }else{
    k_fold <- length(flds)
  }
  
  
  # array that will contain the results
  C_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds)),
                      dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds)))
  IBS_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds)),
                        dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds)))
  p_val_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds)),
                          dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds)))
  n_genes_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds)),
                            dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds)))
  
  # loop on each fold
  for(k in 1:K_folds){
    
    print(paste0("*** Start learning for fold: ", k, " ***"))
    
    # build training and testing dataset
    id_test <- flds[[k]]
    
    data_train <- data[-id_test,]
    data_test <- data[id_test,]
    y_cox_train <- y_cox[-id_test]
    y_cox_test <- y_cox[id_test]
    
    # standardize the data
    mean_train <- apply(data_train, 2, mean)
    sd_train <- apply(data_train, 2, sd)
    data_train_std <- sweep(data_train, 2, mean_train, FUN="-")
    data_train_std <- sweep(data_train_std,2, sd_train, FUN="/")
    
    data_test_std <- sweep(data_test, 2, mean_train, FUN="-")
    data_test_std <- sweep(data_test_std,2, sd_train, FUN="/")
    
    # remove features with potential NA values after standardization
    id_NA_train <- apply(data_train_std, 2, function(x) sum(is.na(x)) >= 1)
    data_train_std <- data_train_std[, !id_NA_train]
    data_test_std <- data_test_std[, !id_NA_train]
    data_train <- data_train[,!id_NA_train]
    
    # compute p-values with benjamini-Hoshberg correction on the training dataset for the pre-filtering
    p_val_cox_train <- p_val_univCox_func(data_train_std, y_cox_train)
    p_val_cox_train_BH <- p.adjust(p_val_cox_train, method = "BH")
    
    # remove potential NA values computed with univariate Cox
    id_NA <- is.na(p_val_cox_train_BH)
    data_train_std <- data_train_std[,!id_NA]
    data_train <- data_train[,!id_NA]
    data_test_std <- data_test_std[, !id_NA]
    p_val_cox_train_BH <-  p_val_cox_train_BH[!id_NA]
    
    # compute IQR on the training dataset for the pre-filtering
    IQR_vect_train <- apply(data_train, 2, IQR)
    
    # loop on each univariate cox p-values threshold
    for(i in 1:length(thrs_p_val)){
      
      print(paste0("*** Start learning for p_val_univ_Cox: ", thrs_p_val[i], "***"))
      
      # loop on each IQR threshold
      for(j in 1:length(IQR_thrs)){
        
        print(paste0("*** Start learning for IQR_thrs: ", IQR_thrs[j], "***"))
        
        # at least two features have to be selected by the pre-filtering step
        if(sum(p_val_cox_train_BH <= thrs_p_val[i] & 
               IQR_vect_train >= IQR_thrs[j]) >= 2){
          
          # number of genes remaining after pre-filtering
          n_genes_ary_final[i,j,k] <- sum(p_val_cox_train_BH <= thrs_p_val[i] & 
                                           IQR_vect_train >= IQR_thrs[j])
          
          # pre-filter the data and learn the models
          data_train_tmp <- data_train_std[, p_val_cox_train_BH <= thrs_p_val[i] & 
                                                       IQR_vect_train >= IQR_thrs[j] ]
          data_test_tmp <- data_test_std[, p_val_cox_train_BH <= thrs_p_val[i] & 
                                                     IQR_vect_train >= IQR_thrs[j] ]
          
          # univariate cox
          if(method == "univ_cox"){
            cv_p_val_univCox <-p_val_univCox_func(data_train_tmp, y_cox_train)
            
            if(all(is.na(cv_p_val_univCox))){
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }else if(is.null(cv_p_val_univCox)){
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }else{
              # gene with minimum p-value on the training dataset
              min_p_val <- names(cv_p_val_univCox)[which.min(cv_p_val_univCox)]
              
              # compute the concordance on the testing data with this gene
              PI_test <- data_test_tmp[, min_p_val]
              PI_train <- data_train_tmp[, min_p_val]
              
              C_ary_final[i,j,k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                     surv.event = clinical_data_srv[id_test,"status"])$c.index 
              IBS_ary_final[i,j,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), #times,
                                             type = "brier")$ierror
              p_val_ary_final[i,j,k] <- logrank_test_func(y_cox_test, PI_test, 0.5, 0.5)
            }
            
            # multivariate cox
          }else{
            
            # fit a model
            cv_fit_tmp <- learn_models(data_train_tmp, y_cox_train, method, 1)[[1]]
            
            # features selected
            names_genes_selected <- genes_selected_func(cv_fit_tmp, "lambda.min")
            
            # at least one features have to be selected for computing prognostic indices
            if(length(names_genes_selected) >= 1){
              
              # compute the C-index
              beta <- as.numeric(coef(cv_fit_tmp, "lambda.min"))
              names(beta) <- colnames(data_train_tmp) 
              
              PI_test <-  data_test_tmp %*% beta 
              PI_train <- data_train_tmp %*% beta 
              
              C_ary_final[i,j,k] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                     surv.event = clinical_data_srv[id_test,"status"])$c.index 
              IBS_ary_final[i,j,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), #times,
                                             type = "brier")$ierror
              p_val_ary_final[i,j,k] <- logrank_test_func(y_cox_test, PI_test, 0.5, 0.5)
            }else{
              
              warning("Less than 2 features selected by the pre-filtering step. 
                      Default values: C-index (0.5), IBS (0.25), p-val logrank test (1)")
              
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }
          }
        } else{
          print("Less than two genes selected after pre-filtering.")}
      }
    } 
  }
  
  return(list(C_ary = C_ary_final, IBS_ary = IBS_ary_final, p_val_ary = p_val_ary_final,
              n_genes_ary = n_genes_ary_final))
}

# @description: plot the median C-indices for each pair of thresholds as in
#               figure 2A
#
# @parameters: 
#   - ary: array of the C-indices computed with the function 'optimize_pre_filtering'
#   - n_genes_ary: number of genes retained after the pre-filtering (output of the 
#                 function 'optimize_pre_filtering')
#
# @return: plot a grid with the median C-indices for each pair of thresholds as in
#               figure 2A
plot_preFiltering <- function(ary, n_genes_ary){
  
  # compute the medians
  med_ary <- apply(ary, c(1,2), median)
  med_n_genes_df <- apply(n_genes_ary, c(1,2), median)
  
  # data frame for the plot
  med_ary <- melt(med_ary)
  med_ary[,c(1,2)] <- apply(med_ary[,c(1,2)], 2, as.character)
  med_n_genes_df <- melt(med_n_genes_df)
  med_n_genes_df[,c(1,2)] <- apply(med_n_genes_df[,c(1,2)], 2, as.character)
  
  # thresholds for the optimal pre-filtering
  higher_C <- best_thrs(ary)
  id_p_val_best <- which(thrs_p_val %in% higher_C[1]) 
  id_IQR_best <- which(IQR_thrs %in% higher_C[2]) 
  
  # plot
  ggplot() + 
    geom_tile(data = med_ary, aes(x=Var1, y=Var2, fill=value))  + xlab("Threshold (p-value)") + 
    ylab("Threshold (IQR)") + theme_Publication_legend_right() +
    scale_fill_gradient(high = "#132B43", low = "#56B1F7", name = "C-index") +
    geom_text(data = med_n_genes_df, aes(x=Var1, y=Var2, label=value), color = "white", size = 4) +
    geom_rect(size=1., fill=NA, colour="grey",
              aes(xmin=6 - 0.5, xmax=6 + 0.5, ymin=1 - 0.5, ymax=1+ 0.5)) +
    geom_rect(size=1., fill=NA, colour="red",
              aes(xmin=id_p_val_best - 0.5, xmax=id_p_val_best + 0.5, 
                  ymin=id_IQR_best - 0.5, ymax=id_IQR_best+ 0.5))
}

# @description: compute the optimal thresholds
#
# @parameters: 
#   - ary: array of the C-indices computed with the function 'optimize_pre_filtering'
#   
#
# @return: 
#   - thrs_p_val_opt: optimal supervised threshold
#   - thrs_IQR_opt: optimal unsupervised threshold
#   - obj: highest median C-index
best_thrs <- function(ary){
  
  # compute the medians
  med_ary <- apply(ary, c(1,2), median)
  
  # melted data frame
  med_ary <- melt(med_ary)
  
  id_opt <- which.max(med_ary$value)
  
  return(c(thrs_p_val_opt = med_ary[id_opt,1],
              thrs_IQR_opt = med_ary[id_opt,2],
              obj = med_ary[id_opt,3]))
  
}




