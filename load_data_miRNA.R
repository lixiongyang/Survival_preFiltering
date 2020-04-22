###################################################################################### #
# Goal of this script: load mRNA together with miRNA and clinical data
###################################################################################### #

#  --- Working directories ---
wd_cancer <- paste0("data_cancer/", cancer, "/")
wd_fit <- paste0("data_fit/", cancer, "/")

# clinical data ---
clinical_data <- t(read.csv(file = paste0(wd_cancer, cancer, ".clin.merged.txt"),
                            stringsAsFactors = F, sep = "\t", header = F))

colnames(clinical_data) <- clinical_data[1,]
clinical_data <- clinical_data[-1,]
row.names(clinical_data) <- clinical_data[, "patient.bcr_patient_barcode"]
clinical_data <- as.data.frame(clinical_data)

# extract information usefull for survival (srv) analysis (time, status)
clinical_data_srv <- clinical_data_prep(clinical_data)


# mRNA data ---
mRNA_data_all <- t(read.csv(file = paste0(wd_cancer, cancer, ".uncv2.mRNAseq_RSEM_all.txt"),
                            stringsAsFactors = F, sep = "\t", header =F))
mRNA_data_all[1:5, 1:5]

# colnames, rownames, NA values
colnames(mRNA_data_all) <- mRNA_data_all[1,]
mRNA_data_all <- mRNA_data_all[-1,]
row.names(mRNA_data_all) <-  tolower(mRNA_data_all[,1])
mRNA_data_all <- mRNA_data_all[,-1]
mRNA_data_all <- as.matrix(mRNA_data_all)
class(mRNA_data_all) <- "numeric"

# keep only primary solid tumor samples 
# (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
freq_sample_types <- as.data.frame(table(substr(row.names(mRNA_data_all), 14, 15)))
freq_sample_types[,2] <- paste0("(", freq_sample_types[,2], ")") 
print(paste("Sample types:", paste(paste(freq_sample_types[,1], freq_sample_types[,2]), collapse = ", ")))

id_prim_tumor <- as.numeric(substr(row.names(mRNA_data_all), 14, 15)) == 1
mRNA_data <- mRNA_data_all[id_prim_tumor, ]
row.names(mRNA_data) <- substr(row.names(mRNA_data),1,12)
dim(mRNA_data)

id_normal <- as.numeric(substr(row.names(mRNA_data_all), 14, 15)) == 11
print(paste(sum(id_normal), "normal patients"))
mRNA_data_normal <- mRNA_data_all[id_normal, ]
dim(mRNA_data_normal)

# remove genes with NA values
id_NA_gene <- apply(mRNA_data, 2, function(x) sum(which(is.na(x))))
id_NA_gene <- id_NA_gene > 0

id_NA_gene_normal <- apply(mRNA_data_normal, 2, function(x) sum(which(is.na(x))))
id_NA_gene_normal <- id_NA_gene_normal > 0

mRNA_data <- mRNA_data[, !id_NA_gene & !id_NA_gene_normal]
mRNA_data_normal <- mRNA_data_normal[, !id_NA_gene & !id_NA_gene_normal]
print(paste0(sum(id_NA_gene & id_NA_gene_normal), " gene(s) removed due to NA values"))

# remove genes with constant values
id_cst_rm <- apply(mRNA_data, 2, function(x) var(x) == 0)
mRNA_data <- mRNA_data[, !id_cst_rm]
mRNA_data_normal <- mRNA_data_normal[, !id_cst_rm]
print(paste0(sum(id_cst_rm), " gene(s) removed due to constant values"))

# log2 transformation
if(min(mRNA_data) < 0) print("Negative values in genetics dataset")
mRNA_data <- log2(mRNA_data + 1)
mRNA_data[1:3, 1:3]

if(min(mRNA_data_normal) < 0) print("Negative values in genetics dataset")
mRNA_data_normal <- log2(mRNA_data_normal + 1)
mRNA_data_normal[1:3, 1:3]

# miRNA data ---

miRNA_data <- t(read.csv(file = paste0(wd_cancer, cancer, ".miRseq_RPKM.txt"),
                         stringsAsFactors = F, sep = "\t", header =F))
miRNA_data[1:5, 1:5]

# colnames, rownames, NA values
colnames(miRNA_data) <- miRNA_data[1,]
miRNA_data <- miRNA_data[-1,]
row.names(miRNA_data) <-  tolower(miRNA_data[,1])
miRNA_data <- miRNA_data[,-1]
miRNA_data <- as.matrix(miRNA_data)
class(miRNA_data) <- "numeric"

miRNA_data[1:5, 1:5]

# keep only primary solid tumor samples 
# (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
freq_sample_types <- as.data.frame(table(substr(row.names(miRNA_data), 14, 15)))
freq_sample_types[,2] <- paste0("(", freq_sample_types[,2], ")") 
print(paste("Sample types:", paste(paste(freq_sample_types[,1], freq_sample_types[,2]), collapse = ", ")))

id_prim_tumor <- as.numeric(substr(row.names(miRNA_data), 14, 15)) == 1
miRNA_data <- miRNA_data[id_prim_tumor, ]
row.names(miRNA_data) <- substr(row.names(miRNA_data),1,12)
dim(miRNA_data)

# remove genes with NA values
id_NA_gene <- apply(miRNA_data, 2, function(x) sum(which(is.na(x))))
id_NA_gene <- id_NA_gene > 0
miRNA_data <- miRNA_data[, !id_NA_gene]
print(paste0(sum(id_NA_gene), " gene(s) removed due to NA values"))

# remove genes with constant values
id_cst_rm <- apply(miRNA_data, 2, function(x) var(x) == 0)
miRNA_data <- miRNA_data[, !id_cst_rm]
print(paste0(sum(id_cst_rm), " gene(s) removed due to constant values"))

# log2 transformation
if(min(miRNA_data) < 0) print("Negative values in genetics dataset")
miRNA_data <- log2(miRNA_data + 1)
miRNA_data[1:3, 1:3]

# Patients overlap in clinical and genetic datas ---
intersect_patients <- intersect(row.names(mRNA_data), 
                                row.names(clinical_data_srv))
intersect_patients <- intersect(intersect_patients, row.names(miRNA_data))

print(paste("Number of patients in miRNA, mRNA and clinical data:", length(intersect_patients)))

mRNA_data <- mRNA_data[intersect_patients,]
miRNA_data <- miRNA_data[intersect_patients,]
clinical_data_srv <- clinical_data_srv[intersect_patients,]
clinical_data <- clinical_data[intersect_patients, ]

# surv obect vector build with clinical data
y_cox <- Surv(clinical_data_srv$time, clinical_data_srv$status)
names(y_cox) <- row.names(clinical_data_srv)

# change name convention
colnames(mRNA_data) <- sub('.*\\|', "", colnames(mRNA_data))

# interquantile range
IQR_vect_mRNA <- apply(mRNA_data, 2, IQR)
IQR_vect_miRNA <- apply(miRNA_data, 2, IQR)

# dimension and censoring rate
print(paste0("n patients: ", nrow(mRNA_data), ", p genes: ", ncol(mRNA_data)))
print(paste0("Censoring rate: ", signif(1 - sum(clinical_data_srv$status) / nrow(mRNA_data),3)))


dim(clinical_data_srv)
dim(mRNA_data)
dim(miRNA_data)