# Optional - to install all packages at once

#from CRAN
install.packages('survival')
install.packages('survminer')
install.packages('glmnet')
install.packages('ggplot2')
install.packages('ggthemes')
install.packages('ggpubr')
install.packages('caret')
install.packages('survAUC')
install.packages('reshape2')
install.packages('scales')
install.packages('zipfR')

#from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("survcomp")
