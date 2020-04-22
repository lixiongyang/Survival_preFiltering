

# Multivariate Cox models -------------------------------------------------

# @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - y_cox: Surv object containing follow-up times and status
#   - method: "EN" for elastic net, "lasso" for lasso, "ridge" for ridge and "AEN"
#             for adaptive elastic net
#   - n_rep: number of models to learn
#
# @return:
#   - list of length 'n_rep' containing 'cv.glmnet' object
learn_models <- function(gene_data, y_cox, method, n_rep){
  
  cv_fits <- list()
  
  # Elastic Net
  if(method == "EN"){
    for(i in 1:n_rep){
      cv_fits[[i]] <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 0.3, nfolds = 5, grouped = T, standardize = F)
      if(i%%5 == 0){
        print(paste0("Number of models learned: ",i))
      }
    }
    
    # lasso
  }else if(method == "lasso"){
    for(i in 1:n_rep){
      cv_fits[[i]] <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 1, nfolds = 5, grouped = T, standardize = F)
      if(i%%5 == 0){
        print(paste0("Number of models learned: ",i))
      }
    }
    
    # ridge
  }else if(method == "ridge"){
    for(i in 1:n_rep){
      cv_fits[[i]] <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 0, nfolds = 5, grouped = T, standardize = F)
      if(i%%5 == 0){
        print(paste0("Number of models learned: ",i))
      }
    }
    
  }else if(method == "coxph"){
    for(i in 1:n_rep){
      cv_fits[[i]] <- coxph(y_cox~gene_data)
    }
    if(i%%5 == 0){
      print(paste0("Number of models learned: ",i))
    }
    
  # Adaptive Elastic Net
  }else if(method == "AEN"){
    for(i in 1:n_rep){
      
      # penalty terms
      fit_ridge <- cv.glmnet(gene_data, y_cox, family = "cox",
                             alpha = 0, nfolds = 5, grouped = T, standardize = F)
      betas_ridge <- as.numeric(coef(fit_ridge, "lambda.min"))
      
      w_AEN <- abs(1 / betas_ridge)
      
      # apply the method to the new data set
      cv_fits[[i]] <- cv.glmnet(gene_data, y_cox, family = "cox", 
                                nfolds = 5, standardize = F, alpha = 0.3, 
                                grouped = T, penalty.factor = w_AEN) 
      if(i%%5 == 0){
        print(paste0("Number of models learned: ",i))
      }
    }
  }
  
  print("Models learned")
  return(cv_fits)
}

# @description: return a vector containing the names of the genes selected
#
# @parameters:
#   - fit: 'cv.glmnet' object
#   - lambda: weight of the penalty
#
# @return:
#   - vector containing the names of the genes selected
genes_selected_func <- function(fit, lambda){
  
  if(class(fit) == "cv.glmnet"){
    # coefficients and active coefficients
    coefs <- coef(fit, s = lambda)
    names_active_coefs <- row.names(coefs)[coefs@i+1]
    return(names_active_coefs)
  }else{
    return(NULL)
  }
  
  
  
}


# Univariate Cox selection ------------------------------------------------

# @description: univariate Cox for one feature
#
# @parameters:
#   - x: feature (one gene)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_cox_func <- function(x, y_cox){
  
  res <- try(fit <- coxph(y_cox ~ x))
  
  if(!inherits(res, "try-error")){
    
    test <- summary(fit)
    return(test$coefficients[5])
    
  }else{
    return(NA)
  }
  
}


# @description: univariate Cox for an entire matrix
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_univCox_func <- function(gene_data, y_cox){
  
  # compute the p-values
  p_val_univCox_vect <- apply(gene_data, 2, function(x) p_val_cox_func(x, y_cox)) 
  
  print("Model learned")
  
  return(p_val_univCox_vect)
}

# @description: compute the p-value of the log-rank test with a 'survdiff' 
#               object as input
#
# @parameters:
#   - x: 'survdiff' object
#
# @return:
#   - the p-value of logrank test
pvalue.survdiff <- function(x, log_p=FALSE, ...){
  pchisq(x$chisq, length(x$n) - 1, lower.tail=FALSE, log.p=log_p)
}

# @description: compute the p-value of the log-rank test with a risk
#               vector
#
# @parameters:
#   - y_cox: Surv object containing follow-up times and status
#   - pred_vect: risk vector used to classify the patients
#   - prob_low: quantile level for low group patients according to pred_vect 
#   - prob_high: quantile level for high group patients according to pred_vect 
#
# @return:
#   - the p-value of logrank test
logrank_test_func <- function(y_cox, pred_vect, prob_low, prob_high){
  
  # quantiles
  q_low <- quantile(pred_vect, probs = prob_low)
  q_high <- quantile(pred_vect, probs = prob_high)
  
  # id of patients in the low and high group
  id_low <- pred_vect <= q_low
  id_high <- pred_vect > q_high
  
  y_cox_tmp <- c(y_cox[id_low], y_cox[id_high])
  group <- c(rep("low", sum(id_low)), rep("high", sum(id_high)))
  
  surv_obj <- survdiff(y_cox_tmp ~ group)
  return(pvalue.survdiff(surv_obj))
}


# Other functions ---------------------------------------------------------

# @description: compute the inices of the n lowest number in a vector
#               used to compute the indices of the n lowest p-values
#
# @parameters: 
#   - x: a numeric vector
#   - n: the number of indices with lowest values in x to return
#
# @return:
#   - indices of the n lowest number in the vector x
which.minn <- function(x,n=1){
  if (n==1)
    which.min(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=FALSE)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}

