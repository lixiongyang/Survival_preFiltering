########################################################################################
# Goal of this script: study the characteristics of cancers in TCGA with more
# than 500 patients
# Cancer name, Number of patients, of genes, censoring rate, median concordance (ridge with all genes),
# number of genes with IQR >= 2, survival at 3 years, survival at 5 years
######################################################################################


# loop to compute all the characteristics
if(!file.exists(file = "data_fit/C_all_cancers.RData")){
  
  # data frame with all the characteristics
  clinical_feat <- c("n patients", "p genes", "censoring rate",
                     "survival - 3 years", "survival - 5 years",
                     "Female/Male/NA", "Age at diagnosis (q25/q75/NA)")
  data_ch <- data.frame(matrix(ncol = length(clinical_feat), nrow = length(cancers_vect)))
  row.names(data_ch) <- cancers_vect
  colnames(data_ch) <- clinical_feat
  
  # data frames with the C-indices and the IBS
  C_df_all_can <- matrix(ncol = length(cancers_vect), nrow = K_folds*n_rep)
  colnames(C_df_all_can) <- cancers_vect
  
  IBS_df_all_can <- matrix(ncol = length(cancers_vect), nrow = K_folds*n_rep)
  colnames(IBS_df_all_can) <- cancers_vect
  
  for(cancer in cancers_vect){
    
    print(paste0("*** Cancer under charaterization: ", cancer, " ***"))
    
    source(file = "load_data.R")  
    
    # gender
    id_male <- clinical_data$patient.gender == "male"
    id_female <- clinical_data$patient.gender == "female"
    id_gender_NA <- is.na(clinical_data$patient.gender)
    data_ch[cancer, "Female/Male/NA"] <- paste(sum(id_female), sum(id_male), sum(id_gender_NA), sep = " / ")  
    print(paste("Female/Male/NA:", paste(sum(id_female), sum(id_male), sum(id_gender_NA), sep = " / ")))
    
    # age
    age <- as.numeric(as.character(clinical_data$patient.age_at_initial_pathologic_diagnosis))
    data_ch[cancer, "Age at diagnosis (q25/q75/NA)"] <- paste0(median(age, na.rm = T), " (",
                                                               quantile(age, probs = 0.25, na.rm = T), "/", 
                                                               quantile(age, probs = 0.75, na.rm = T), "/",
                                                               sum(is.na(age)), ")")
    print(paste("Age at diagnosis:", median(age, na.rm = T)))
    
    # remove genes with NA values
    id_NA_gene <- apply(gene_data, 2, function(x) sum(which(is.na(x))))
    id_NA_gene <- id_NA_gene > 0
    gene_data <- gene_data[, !id_NA_gene]
    print(paste0(sum(id_NA_gene), " gene(s) removed due to NA values"))
    
    # remove genes with constant values
    id_cst_rm <- apply(gene_data, 2, function(x) var(x) == 0)
    gene_data <- gene_data[, !id_cst_rm]
    print(paste0(sum(id_cst_rm), " gene(s) removed due to constant values"))
    
    # log2 transformation
    if(min(gene_data) < 0) print("Negative values in genetics dataset")
    gene_data <- log2(gene_data + 1)
    gene_data[1:3, 1:3]
    
    # surv obect vector build with clinical data
    y_cox <- Surv(clinical_data_srv$time, clinical_data_srv$status)
    
    # interquantile range
    IQR_vect <- apply(gene_data, 2, IQR)
    hist(gene_data, breaks = 50, xlab = "Gene expression", main = paste0(cancer, ": gene expression"))
    hist(IQR_vect, breaks = 50, xlab = "IQR", main = paste0(cancer, ": IQR"))
    
    # dimension and censoring rate
    print(paste0("n patients: ", nrow(gene_data), ", p genes: ", ncol(gene_data)))
    print(paste0("Censoring rate: ", signif(1 - sum(clinical_data_srv$status) / nrow(gene_data),3)))
    data_ch[cancer, "n patients"] <- nrow(gene_data)
    data_ch[cancer, "p genes"] <- ncol(gene_data)
    data_ch[cancer, "censoring rate"] <- signif(1 - sum(clinical_data_srv$status) / nrow(gene_data),3)
    
    # --- survival at 3 and 5 years --- 
    km_fit <- survfit(y_cox~1)
    km_fit_summary <- summary(km_fit, times = c(3,5))
    data_ch[cancer, "survival - 3 years"] <- paste(signif(km_fit_summary$surv[1],3), "+/-",
                                                   signif(km_fit_summary$std.err[1],3))
    data_ch[cancer, "survival - 5 years"] <- paste(signif(km_fit_summary$surv[2],3), "+/-",
                                                   signif(km_fit_summary$std.err[2],3))
    print(paste("Survival at 3 years:", paste(signif(km_fit_summary$surv[1],3), "+/-",
                                              signif(km_fit_summary$std.err[1],3))))
    print(paste("Survival at 5 years:", paste(signif(km_fit_summary$surv[2],3), "+/-",
                                              signif(km_fit_summary$std.err[2],3))))
    
    # --- compute the concordance ---
    if(!file.exists(paste0("data_fit/", cancer, "/ridge/C_all_genes_ridge.RData"))){
      
      ind <- 1
      
      for(i in 1: n_rep){
        
        print(paste0("*** Start learning for repetition number: ", i, " ***"))
        
        flds <- createFolds(1:nrow(gene_data), k = K_folds, list = TRUE, returnTrain = FALSE)
        
        for(k in 1:K_folds){
          
          print(paste0("*** Start learning for fold: ", k, " ***"))
          
          # build training and testing dataset
          id_test <- flds[[k]]
          
          # standardize the data
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
          
          # learn a model
          res <- try(fit_ridge <- learn_models(gene_data_train_std, y_cox_train, "ridge", 1)[[1]])
          
          if(!inherits(res, "try-error")){
            # compute the concordance
            beta <- as.numeric(coef(fit_ridge, "lambda.min"))
            names(beta) <- colnames(gene_data_train_std) 
            
            PI_test <-  gene_data_test_std %*% beta 
            PI_train <- gene_data_train_std %*% beta 
            
            C_df_all_can[ind, cancer] <- concordance.index(PI_test, surv.time = clinical_data_srv[id_test,"time"], 
                                                           surv.event = clinical_data_srv[id_test,"status"])$c.index 
            IBS_df_all_can[ind, cancer] <- predErr(y_cox_train, y_cox_test, PI_train, PI_test, times = seq(0, 3, by=1/365.25), #times,
                                                   type = "brier")$ierror
          }
          
          ind <- ind + 1
          
        }
      }
      
      # save the C- indices and the IBS
      C_vect_ridge <- C_df_all_can[, cancer]
      IBS_vect_ridge <- IBS_df_all_can[, cancer]
      
      save(C_vect_ridge, IBS_vect_ridge, 
           file = paste0("data_fit/", cancer, "/C_all_genes_ridge.RData"))
      
    }else{
      load(file = paste0("data_fit/", cancer, "/C_all_genes_ridge.RData"))
      
      C_df_all_can[ind, cancer] <- C_vect_ridge
      IBS_df_all_can[ind, cancer] <- IBS_vect_ridge
    }
    
  }
  
  # save the characteristics of the data as csv file
  write.csv(data_ch, file = "cancer_characteristics.csv")
  print(paste("csv file with characteristics of the different cancers saved"))
  
  # save the c-indices of all cancers
  save(C_df_all_can, IBS_df_all_can, data_ch, file = "data_fit/C_IBS_ch_all_cancers.RData")
  
}else{
  load(file = "data_fit/C_IBS_ch_all_cancers.RData")
}

C_df_ch <- melt(C_df_all_can)
C_df_ch <- C_df_ch[, -1]
C_df_ch <- cbind(C_df_ch, rep(NA, nrow(C_df_ch)))
colnames(C_df_ch) <- c("cancer", "C", "color")
C_df_ch[C_df_ch$cancer %in% c("KIRC", "BRCA", "LUSC", "LGG"), 'color'] <- "lightblue"


# Boxplot of the concordance for all cancers
C_all_cancers_ggplot <- 
  ggplot(data = C_df_ch, aes(x = reorder(cancer, C, FUN = function(x) median(x, na.rm = T)), 
                                                   y=C, fill=color)) +
  geom_boxplot() + theme_Publication() + scale_fill_Publication() +
  geom_hline(yintercept=0.5, col="red", linetype="dashed", size = 0.5) +
  scale_fill_manual(values="lightblue") +
  theme(axis.text.x = element_text(size = 10)) + 
  ylab("C-index") +
  xlab("Cancer") + scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) 
print(C_all_cancers_ggplot)

# median C-indices for all cancers
apply(C_df_all_can, 2, function(x) median(x, na.rm = T))





