###################################################################################### #
# Goal of this script: show the results with boxplots comparing the C-indices 
# obtained without pre-filtering and in the optimal case using the output of the 
# scirpt 'optimal_case_vs_no_filtering.R".
###################################################################################### #

# load the data learned in 'optimal_case_cv_no_filtering.R'
load(file = paste0("data_fit/", cancer, "/best_vs_no_flt.RData"))

# data frames for ggplot
C_best_df <- rbind(cbind(melt(C_df_flt), flt = "flt"),
                   cbind(melt(C_df), flt = "no_flt"))
C_best_df$flt <- factor(C_best_df$flt, levels = c("no_flt", "flt"), ordered = TRUE)

IBS_best_df <- rbind(cbind(melt(IBS_df_flt), flt = "flt"),
                     cbind(melt(IBS_df), flt = "no_flt"))
IBS_best_df$flt <- factor(IBS_best_df$flt, levels = c("no_flt", "flt"), ordered = TRUE)

# p-value of the one-sided wilcoxon test
p_val_C <- rep(NA, length(methods_vect))
p_val_IBS <- rep(NA, length(methods_vect))
for(i in 1:length(methods_vect)){
  
  # C-index
  p_C_tmp <- wilcox.test(C_df_flt[,i], C_df[,i], alternative = "greater", paired = T)$p.value
  if(p_C_tmp < 1.e-2){
    p_val_C[i] <- sprintf("p = %1.2e", p_C_tmp)
  }else{
    p_val_C[i] <- sprintf("p = %1.2f", p_C_tmp)
  }
  
  
  # IBS
  p_IBS_tmp <- wilcox.test(IBS_df_flt[,i], IBS_df[,i], alternative = "less", paired = T)$p.value
  if(p_IBS_tmp < 1.e-2){
    p_val_IBS[i] <- sprintf("p = %1.2e", p_IBS_tmp)
  }else{
    p_val_IBS[i] <- sprintf("p = %1.2f", p_IBS_tmp)
  }
}


# median number of genes
n_genes_bestC <- rep(NA, length(methods_vect))
ind <- 1

for(method in methods_vect){
  
  # load the data
  load(file = paste0("data_fit/", cancer, "/", method, "/optimize_flt_mRNA.RData"))
  
  # best median C-index
  higher_C <- best_thrs(optimize_flt_mRNA$C_ary)
  id_p_val_best <- which(thrs_p_val %in% higher_C[1]) 
  id_IQR_best <- which(IQR_thrs %in% higher_C[2]) 
  
  med_n_genes_df <- apply(optimize_flt_mRNA$n_genes_ary, c(1,2), median)
  n_genes_bestC[ind] <- med_n_genes_df[id_p_val_best, id_IQR_best]
  
  ind <- ind + 1
}

# ggplot
IBS_best_ggplot <- ggplot(data = IBS_best_df, aes(x = Var2, y = value, fill = flt)) +
  geom_boxplot(position = position_dodge(0.75)) + 
  xlab("Method") + ylab("IBS") + theme_Publication() + 
  geom_hline(yintercept=0.25, col="red", linetype="dashed", size = 0.5) +
  discrete_scale("fill","Publication",manual_pal(values = c("grey","red"))) +
  annotate("text", x = 1:length(methods_vect), y = 0, label = n_genes_bestC, col = "red", size = 5) +
  annotate("text", x = 1:length(methods_vect), y = 0.4, label = p_val_IBS, size = 4.5)
print(IBS_best_ggplot)

C_best_ggplot <- ggplot(data = C_best_df, aes(x = Var2, y = value, fill = flt)) +
  geom_boxplot(position = position_dodge(0.75)) + 
  xlab("Method") + ylab("C-index") + theme_Publication() + ylim(0.38,NA) +
  geom_hline(yintercept=0.5, col="red", linetype="dashed", size = 0.5) +
  discrete_scale("fill","Publication",manual_pal(values = c("grey","red"))) +
  annotate("text", x = 1:length(methods_vect), y = 0.38, label = n_genes_bestC, col = "red", size = 5) +
  annotate("text", x = 1:length(methods_vect), y = 0.9, label = p_val_C, size = 4.5)
print(C_best_ggplot)


# values for table 2 and 3 ---
apply(C_df_flt, 2, function(x) median(x, na.rm = T))

apply(C_df_flt-C_df, 2, function(x) median(x, na.rm = T)) 