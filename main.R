########################################################################## #
#
# This script sources each script allowing to create each figure of the
# paper. 
#
########################################################################## #


########################################################################## #
#
# This section allows for loading the packages needed, sourcing the 
# functions, and choose if the graphics have to be saved or not.
#
########################################################################## #


# packages needed (to be downloaded from CRAN)
require(survival)
require(survminer)
require(survcomp)
require(glmnet)
require(ggplot2)
require(ggpubr)
require(caret)
library(survAUC)
library(reshape2)
require(scales)
require(zipfR)

# packages needed (to be downloaded from bioconductor)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("survminer")
# BiocManager::install("survcomp")
library(survminer)
library(survcomp)

# load functions
invisible(lapply(list.files(path = "functions/", pattern = "[.]R$", 
                            recursive = TRUE, full.names = T), source))

########################################################################## #
#
# --- Characteristics of the different cancers ----
# 
# This section allows for studying the cancer characteristics.
# 
# It will create the file 'cancer_characteristics.csv' and return a 
# graphics with the concordances for all cancers (Fig. 1).
#
# !!! 'n_models' will be learned for each cancers, and it takes time !!!
#
# Data related to each cancer have to be dowloaded and saved in 
# "data_cancer/'cancer-name'".
# - gene epxression data: 'cancer-name'.uncv2.mRNAseq_RSEM_all.txt
#       e.g. KIRC.uncv2.mRNAseq_RSEM_all.txt
# - clinical data: 'cancer-name'.clin.merged.txt
#       e.g. KIRC.clin.merged.txt
########################################################################## #

# choose if you want to study all cancer characteristics
study_all_can <- F

# Choose the cancers to study
cancers_vect <- "KIRC"

if(study_all_can){
  # Choose the cancers to study
  cancers_vect <- c("BRCA", "COADREAD", "HNSC", "KIRC", "LGG", "LUAD", 
                    "LUSC", "PRAD", "STES", "THCA", "UCEC")
}

# number of repetition and number of folds
K_folds <- 5
n_rep <- 10

# Compute all the characteristics for each cancers
source(file = "cancers_characteristics.R")



########################################################################## #
#
# --- Choose cancer, IQR and p-value thresholds, and load the 
#     data (Fig. 1) ----
#
# This section allows for:
#   - loading the data (mRNA-seq and clinical) for a desired cancer
#   - choosing the different IQR and p-value thresholds to study 
#   - choosing the number of folds for the cross-validation and create them
#   - choosing the methods to study
#
########################################################################## #

# choose the cancer and the IQR threshold
cancer <- "KIRC"
IQR_thrs <- c(0, 2, 4)
thrs_p_val <- c(0.01, 0.5, 1)

# load the data of the desired cancer
source(file = "load_data.R")

# create folds
K_folds <- 5
if(!file.exists(paste0("data_fit/", cancer, "/folds.RData"))){
  set.seed(1234)
  flds <- createFolds(1:nrow(gene_data), k = K_folds, list = TRUE, returnTrain = FALSE)
  save(flds, file = paste0("data_fit/", cancer, "/folds.RData"))
}else{
  load(file = paste0("data_fit/", cancer, "/folds.RData"))
}

# choose the methods to study
methods_vect <- c("ridge", "lasso", "EN", "AEN", "univCox")


########################################################################## #
#
# --- Determine the optimal thresholds (Fig. 2A) ----
#
# This section allows for fixing the optimal thresholds by K-fold cross
# validation
#
########################################################################## #

# learn the models by K-fold cross-validation for each pair of thresholds
for(method in methods_vect){
  
  # test is the C-indices have already been learned and savec
  if(!file.exists(paste0("data_fit/", cancer, "/", method, "/optimize_flt_mRNA.RData"))){
    
    print(paste0("*** Start learning for ", method, " ***"))
    
    # call the function to compute the C-indices by cross-validation
    optimize_flt_mRNA <- optimize_preFiltering(gene_data, y_cox, method, k_fold = 5, 
                                               flds = flds, IQR_thrs, thrs_p_val)
    
    # save the C-indices
    save(optimize_flt_mRNA,
         file =  paste0("data_fit/", cancer, "/", method, "/optimize_flt_mRNA.RData"))
  }
}

# plot the result for the method of your choice
method <- "ridge"
load(file = paste0("data_fit/", cancer, "/", method, "/optimize_flt_mRNA.RData"))

plot_preFiltering(optimize_flt_mRNA$C_ary, optimize_flt_mRNA$n_genes_ary)


########################################################################## #
#
# --- Compare the C-index and the IBS obtained without pre-filtering
#     and in the optimal case (Fig. 2B) ----
#
# This section allows for estimating if the pre-filtering increases the
# prediction accuracy of the data.
#
########################################################################## #

# choose the number of repetitions of the K-fold cross-validation
n_rep <- 10

# choose the optimal supervised (p-value univariate Cox) and unsupervised 
# (IQR) thresholds according to the previous step for each method 
opt_thrs_sup <- c(0.5, 0.5, 0.5, 0.5, 0.01)
opt_thrs_unsup <- c(2, 4, 3, 4, 5)
names(opt_thrs_sup) = names(opt_thrs_unsup) <- methods_vect

# learn the models
source(file = "optimal_case_vs_no_filtering.R")

# plot the results
source(file = "plot/optimal_case_vs_no_filtering_plot.R")


########################################################################## #
#
# --- Study the stability of the prognostic indices (Fig. 3) ----
#
# This section allows for estimating variance of the prognostic index for 
# all patients
#
########################################################################## #

# choose the number of PIs to estimate for each patient 
#   (number of boostrap repetitions)
# number of PIs to estimate per patient
n_rep <- 50

# choose the methods to study (only multivariate)
methods_vect <- c("ridge", "lasso", "EN", "AEN")

# learn the models
source(file = "stability_PIs.R")

# plot the results
source(file = "plot/stability_PIs_plot.R")


########################################################################## #
#
# --- Integration of mRNA data together with miRNA data (Fig. 4) ----
#
# This section allows for mRNA data together with miRNA data with and 
# without pre-filtering
#
########################################################################## #

# load miRNA dataset together with mRNA dataset
source(file = "load_data_miRNA.R")

# choose the optimal thresholds for miRNA dataset and the elastic net
IQR_thrs_miRNA <- c(0, 1, 2)
thrs_p_val_miRNA <- c(0.01, 0.5, 1)
source(file = "opt_thresholds_miRNA.R")

# best thresholds for the elastic net for mRNA dataset
thrs_sup <- 0.5
thrs_unsup <- 3

# best thresholds for the elastic net for miRNA dataset
thrs_sup_miRNA <- 0.5
thrs_unsup_miRNA <- 1

# choose the number of folds and the number of repetitions
K_folds <- 5
n_rep <- 10

# learn the models
source(file = "mixing_miRNA_mRNA.R")

# plot the results
source(file = "plot/mixing_miRNA_mRNA_plot.R")