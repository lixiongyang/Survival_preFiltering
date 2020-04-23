###################################################################################### #
# Goal of this script: show the results of 'mixing_miRNA_mRNA.R' script as in
# figure 4. 
###################################################################################### #

# load the files
load(file = paste0("data_fit/", cancer, "/mixOmics.RData"))


# IBS ---
IBS_mixOmics_df <- as.data.frame(rbind(cbind(as.vector(IBS_df), "Without pre-filtering", "mix"),
                                     cbind(as.vector(IBS_df_mRNA), "Without pre-filtering", "mRNA"),
                                     cbind(as.vector(IBS_df_miRNA), "Without pre-filtering", "miRNA"),
                                     cbind(as.vector(IBS_df_flt), "With pre-filtering", "mix_flt"),
                                     cbind(as.vector(IBS_df_flt_mRNA), "With pre-filtering", "mRNA_flt"),
                                     cbind(as.vector(IBS_df_flt_miRNA), "With pre-filtering", "miRNA_flt")),
                               stringsAsFactors = F)

colnames(IBS_mixOmics_df) <- c("IBS", "flt", "data_type")
str(IBS_mixOmics_df)
IBS_mixOmics_df$IBS <-as.numeric(IBS_mixOmics_df$IBS)
str(IBS_mixOmics_df)
IBS_mixOmics_df$data_type <- factor(IBS_mixOmics_df$data_type, levels = c("mix", "mRNA", "miRNA",
                                                                      "mix_flt", "mRNA_flt", "miRNA_flt"))

# comparison to compute the p-values of the one-sided wilcoxon signed rank test
my_comparison <- list(c("mix_flt", "mRNA_flt"), c("mix", "mRNA"), c("mix_flt", "mix"))

IBS_mixOmics_ggplot <- ggplot(data = IBS_mixOmics_df, aes(x = data_type, y = IBS, fill = data_type)) +
  geom_boxplot(position = position_dodge(0.75)) + 
  xlab("") + ylab("IBS") + theme_Publication() +
  stat_compare_means(comparisons = my_comparison, paired = T, method.args = list(alternative = "greater"))+
  discrete_scale("fill","Publication", manual_pal(values = c("dimgray", "darkgray", "gainsboro", "darkred", "red", "#ffcccb"))) +
  scale_x_discrete(breaks=c("mRNA_flt", "mRNA"),
                   labels=c("With pre-filtering", "Without pre-filtering")) +
  annotate("text", x = 1:6, y = 0.4, label = c("Mix", "mRNA", "miRNA", "Mix", "mRNA", "miRNA"), 
           col = c("dimgray", "dimgray", "dimgray", "red", "red", "red"), size = 5)
print(IBS_mixOmics_ggplot)


# C-index ---
C_mixOmics_df <- as.data.frame(rbind(cbind(as.vector(C_df), "Without pre-filtering", "mix"),
                                     cbind(as.vector(C_df_mRNA), "Without pre-filtering", "mRNA"),
                                     cbind(as.vector(C_df_miRNA), "Without pre-filtering", "miRNA"),
                                     cbind(as.vector(C_df_flt), "With pre-filtering", "mix_flt"),
                                     cbind(as.vector(C_df_flt_mRNA), "With pre-filtering", "mRNA_flt"),
                                     cbind(as.vector(C_df_flt_miRNA), "With pre-filtering", "miRNA_flt")),
                               stringsAsFactors = F)

colnames(C_mixOmics_df) <- c("C", "flt", "data_type")
str(C_mixOmics_df)
C_mixOmics_df$C <-as.numeric(C_mixOmics_df$C)
str(C_mixOmics_df)
C_mixOmics_df$data_type <- factor(C_mixOmics_df$data_type, levels = c("mix", "mRNA", "miRNA",
                                                                      "mix_flt", "mRNA_flt", "miRNA_flt"))

# comparison to compute the p-values of the one-sided wilcoxon signed rank test
my_comparison <- list(c("mix_flt", "mRNA_flt"), c("mix", "mRNA"), c("mix_flt", "mix"))

C_mixOmics_ggplot <- ggplot(data = C_mixOmics_df, aes(x = data_type, y = C, fill = data_type)) +
  geom_boxplot(position = position_dodge(0.75)) + 
  xlab("") + ylab("C-index") + theme_Publication() +
  stat_compare_means(comparisons = my_comparison, paired = T, method.args = list(alternative = "greater"))+
  discrete_scale("fill","Publication", manual_pal(values = c("dimgray", "darkgray", "gainsboro", "darkred", "red", "#ffcccb"))) +
  scale_x_discrete(breaks=c("mRNA_flt", "mRNA"),
                   labels=c("With pre-filtering", "Without pre-filtering")) +
  annotate("text", x = 1:6, y = 0.4, label = c("Mix", "mRNA", "miRNA", "Mix", "mRNA", "miRNA"), 
           col = c("dimgray", "dimgray", "dimgray", "red", "red", "red"), size = 5)
print(C_mixOmics_ggplot)