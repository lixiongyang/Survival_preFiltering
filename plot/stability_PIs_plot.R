###################################################################################### #
# Goal of this script: show the results of the stability of the prognostic indices
# (median and IQR of the PIs of each patient) unsing the results of the script
# 'stability_PIs.R'.
###################################################################################### #

IQR_bootstrap_list = IQR_bootstrap_list_flt = med_bootstrap_list = med_bootstrap_list_flt <- list()
ind <- 1

for(method in methods_vect){
  
  # load the PIs
  load(file = paste0("data_fit/", cancer, "/", method, "/PI_bootstrap.RData"))
  PI_test_df <- do.call(rbind, PI_test_list)
  PI_test_df_flt <- do.call(rbind, PI_test_list_flt)
  
  # compute the IQRs and the medians for each patient
  IQR_bootstrap_list[[ind]] <- apply(PI_test_df, 1, IQR)
  IQR_bootstrap_list_flt[[ind]] <- apply(PI_test_df_flt, 1, IQR)
  
  med_bootstrap_list[[ind]] <- apply(PI_test_df, 1, median)
  med_bootstrap_list_flt[[ind]] <- apply(PI_test_df_flt, 1, median)
  
  ind <- ind + 1
}

names(IQR_bootstrap_list) = names(IQR_bootstrap_list_flt) =
  names(med_bootstrap_list) = names(med_bootstrap_list_flt) <- methods_vect

# data frames with the IQR and the medians
IQR_bootstrap_df <- rbind(cbind(stack(IQR_bootstrap_list), flt = "No pre-filtering"),
                          cbind(stack(IQR_bootstrap_list_flt), flt = "Pre-filtering"))
med_bootstrap_df <- rbind(cbind(stack(med_bootstrap_list), flt = "No pre-filtering"),
                          cbind(stack(med_bootstrap_list_flt), flt = "Pre-filtering"))

dim(IQR_bootstrap_df)

# boxplot of the median PI for each patient
med_bootstrap_ggplot <- ggplot(data = med_bootstrap_df, aes(x = ind, y = values)) +
  geom_boxplot(aes(fill = flt), position = position_dodge(0.75)) + ylim(-5,5) +
  xlab("Method") + ylab("median (PI)") +
  theme_Publication() + 
  discrete_scale("fill","Publication",manual_pal(values = c("grey","red")))
print(med_bootstrap_ggplot)

# boxplot of the variability of PI (IQR) for each patient
IQR_bootstrap_ggplot <- ggplot(data = IQR_bootstrap_df, aes(x = ind, y = values)) +
  geom_boxplot(aes(fill = flt), position = position_dodge(0.75)) + ylim(0,4) +
  xlab("Method") + ylab("IQR (PI)") +
  theme_Publication() + 
  discrete_scale("fill","Publication",manual_pal(values = c("grey","red")))
print(IQR_bootstrap_ggplot)