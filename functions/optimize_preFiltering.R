
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
optimize_preFiltering <- function(data, y_cox, method, k_fold = 5, flds_tmp = NULL, 
                                  IQR_thrs = 0:4, 
                                  thrs_p_val = c(0.01, 0.05, 0.1, 0.5, 1)){
  
  
  # create folds
  if(is.null(flds_tmp)){
    flds_tmp <- createFolds(1:nrow(data), k = k_fold, list = TRUE, returnTrain = FALSE)
  }else{
    k_fold <- length(flds_tmp)
  }
  
  
  # array that will contain the results
  C_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds_tmp)),
                      dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds_tmp)))
  IBS_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds_tmp)),
                        dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds_tmp)))
  p_val_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds_tmp)),
                          dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds_tmp)))
  n_genes_ary_final <- array(dim=c(length(thrs_p_val), length(IQR_thrs),  length(flds_tmp)),
                            dimnames = list(thrs_p_val, IQR_thrs, 1:length(flds_tmp)))
  
  # loop on each fold
  for(k in 1:k_fold){
    
    print(paste0("Start learning for fold: ", k, " ***"))
    
    # build training and testing dataset
    id_test <- flds_tmp[[k]]
    
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
    
    id_NA_gene_train <- apply(data_train_std, 2, function(x) sum(is.na(x)) >= 1)
    data_train_std <- data_train_std[, !id_NA_gene_train]
    data_test_std <- data_test_std[, !id_NA_gene_train]
    data_train <- data_train[,!id_NA_gene_train]
    
    # compute p-values and IQR
    p_val_cox_train <- p_val_univCox_func(data_train_std, y_cox_train)
    p_val_cox_train_BH <- p.adjust(p_val_cox_train, method = "BH")
    id_NA <- is.na(p_val_cox_train_BH)
    data_train_std <- data_train_std[,!id_NA]
    data_train <- data_train[,!id_NA]
    data_test_std <- data_test_std[, !id_NA]
    p_val_cox_train_BH <-  p_val_cox_train_BH[!id_NA]
    
    IQR_vect_train <- apply(data_train, 2, IQR)
    
    for(i in 1:length(thrs_p_val)){
      
      print(paste0("Start learning for p_val_univCox: ", thrs_p_val[i], "***"))
      
      for(j in 1:length(IQR_thrs)){
        
        print(paste0("Start learning for IQR_thrs: ", IQR_thrs[j], "***"))
        
        if(sum(p_val_cox_train_BH <= thrs_p_val[i] & 
               IQR_vect_train >= IQR_thrs[j]) >= 2){
          
          # number of genes
          n_genes_ary_final[i,j,k] <- sum(p_val_cox_train_BH <= thrs_p_val[i] & 
                                           IQR_vect_train >= IQR_thrs[j])
          
          # pre-filter the data and learn the models
          data_train_tmp <- data_train_std[, p_val_cox_train_BH <= thrs_p_val[i] & 
                                                       IQR_vect_train >= IQR_thrs[j] ]
          data_test_tmp <- data_test_std[, p_val_cox_train_BH <= thrs_p_val[i] & 
                                                     IQR_vect_train >= IQR_thrs[j] ]
          if(method == "univCox"){
            cv_p_val_univCox <- p_val_univCox_func(data_train_tmp, y_cox_train)
            
            if(all(is.na(cv_p_val_univCox))){
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }else if(is.null(cv_p_val_univCox)){
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }else{
              
              # gene with minimum p-value
              gene_min_p_val <- names(cv_p_val_univCox)[which.min(cv_p_val_univCox)]
              
              res <- try(fit <- coxph(y_cox_train~data_train_tmp[, gene_min_p_val]))
              if(inherits(res, "try-error")){next}
              
              # compute the concordance on the testing data with this gene
              PI_test <- data_test_tmp[, gene_min_p_val] * fit$coefficients
              PI_train <- data_train_tmp[, gene_min_p_val] * fit$coefficients
              
              C_ary_final[i,j,k] <- max(concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                         surv.event = clinical_data_srv[id_test,"status"])$c.index,
                                       concordance.index(-PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                         surv.event = clinical_data_srv[id_test,"status"])$c.index) 
              IBS_ary_final[i,j,k] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), #times,
                                             type = "brier")$ierror
              p_val_ary_final[i,j,k] <- logrank_test_func(y_cox_test, PI_test, 0.5, 0.5)
            }
            
          }else{
            res <- try(cv_fit_tmp <- learn_models(data_train_tmp, y_cox_train, method, 1)[[1]])
            
            if(!inherits(res, "try-error")){
              
              names_genes_selected <- genes_selected_func(cv_fit_tmp, "lambda.min")
              
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
                C_ary_final[i,j,k] <- 0.5
                IBS_ary_final[i,j,k] <- 0.25
                p_val_ary_final[i,j,k] <- 1
              }
            }else{
              C_ary_final[i,j,k] <- 0.5
              IBS_ary_final[i,j,k] <- 0.25
              p_val_ary_final[i,j,k] <- 1
            }
            
          }
        }else{
          print("Less than two genes selected after pre-filtering.")
        }
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
  id_p_val_best <- which(row.names(ary) %in% higher_C[1]) 
  id_IQR_best <- which(colnames(ary) %in% higher_C[2]) 
  
  # plot
  ggplot() + 
    geom_tile(data = med_ary, aes(x=Var1, y=Var2, fill=value))  + xlab("Threshold (p-value)") + 
    ylab("Threshold (IQR)") + theme_Publication_legend_right() +
    scale_fill_gradient(high = "#132B43", low = "#56B1F7", name = "C-index") +
    geom_text(data = med_n_genes_df, aes(x=Var1, y=Var2, label=value), color = "white", size = 4) +
    geom_rect(size=1., fill=NA, colour="grey",
              aes(xmin= dim(ary)[1] - 0.5, xmax=dim(ary)[1] + 0.5, ymin=1 - 0.5, ymax=1+ 0.5)) +
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

# @description: pre-filter the data with both supervised and unsupervised
# filtering
#
# @parameters: 
#   - gene_data: data to be filtered
#   - y_cox: Surv object containing follow-up times and status
#   - thrs_sup: threshold for the supervised (p-value of the univariate
#     Cox model) pre-filtering 
#   - thrs_unsup: threshold for the unsupervised (IQR) pre-filtering
#   
# @return: 
#   - gene_data_flt: data after the pre-filtering
pre_filter_func <- function(gene_data, y_cox, thrs_sup, thrs_unsup){
  
  # standardize the data
  gene_data_std <- scale(gene_data)
  
  id_NA_gene <- apply(gene_data_std, 2, function(x) sum(is.na(x)) >= 1)
  gene_data_std <- gene_data_std[, !id_NA_gene]
  gene_data <- gene_data[,!id_NA_gene]
  
  # compute p-values and IQR
  p_val_cox <- p_val_univCox_func(gene_data_std, y_cox)
  p_val_cox_BH <- p.adjust(p_val_cox, method = "BH")
  id_NA <- is.na(p_val_cox_BH)
  gene_data_std <- gene_data_std[,!id_NA]
  gene_data <- gene_data[,!id_NA]
  p_val_cox_BH <-  p_val_cox_BH[!id_NA]
  
  IQR_vect <- apply(gene_data, 2, IQR)
  
  # pre-filter the data
  id_non_flt_genes <- p_val_cox_BH <= thrs_sup & IQR_vect >= thrs_unsup
  gene_data_flt <- gene_data[, id_non_flt_genes]
  
  return(gene_data_flt)
}


