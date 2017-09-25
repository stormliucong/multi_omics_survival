# read luad data.
##########################################
# read clinical file.
##########################################
rm(list=ls())
setwd(getwd())
source("surTools.r")
clinical_data <- read.table("gdac.broadinstitute.org_LUAD.Clinical_Pick_Tier1.Level_4.2016012800.0.0/LUAD.clin.merged.picked.txt",sep = "\t",header = T,row.names = 1)
patient_id <- colnames(clinical_data)
length(patient_id) # 522
time <- rep(-1,length(patient_id))
censor <- rep(-1,length(patient_id))
# when vital status = 1 (i.e. death), check days_to_death first.
death_idx <- which(clinical_data["vital_status",]==1) 
tmp <- clinical_data["days_to_death",death_idx]
time[death_idx] <- as.numeric(as.matrix(tmp))
no_follow_up_idx <- c(which(is.na(time)),which(time==0))
time[no_follow_up_idx]
clinical_data[,no_follow_up_idx]
censor[death_idx] <- rep(1,length(death_idx))
# when vital status = 0 (i.e. not death), check days_to_last_followup.
tmp <- clinical_data["days_to_last_followup",-death_idx]
time[-death_idx] <- as.numeric(as.matrix(tmp))
censor[-death_idx] <- rep(0,length(censor)-length(death_idx))
# remove NA and 0 (i.e. no follow_up data available).
no_follow_up_idx <- c(which(is.na(time)),which(time==0))
time[no_follow_up_idx]
time <- time[-no_follow_up_idx]
censor <- censor[-no_follow_up_idx]
length(time)
patient_id <- patient_id[-no_follow_up_idx]
survival_data <- data.frame(patient_id,time,status=censor)
save(survival_data,file="survival_data.rda")

##########################################
# read mRNA file. log2(RSEM)
# TBD: How to deal w/ NA entries?
##########################################
rm(list=ls())
setwd(getwd())
source("surTools.r")
load("survival_data.rda")
setwd(getwd())
mRNA_data <- read.table("gdac.broadinstitute.org_LUAD.mRNAseq_Preprocess.Level_3.2016012800.0.0/LUAD.uncv2.mRNAseq_RSEM_normalized_log2.txt",header = T,row.names = 1,sep = "\t")
dim(mRNA_data)
patient_id <- colnames(mRNA_data)
# only select tumor data.
normal_tumor_indicator <- unlist(strsplit(patient_id,split = "\\."))[seq(4,length(patient_id)*4,by=4)]
tumor_idx <- which(normal_tumor_indicator=="01")
patient_id <- patient_id[tumor_idx]
length(patient_id)
mRNA_matrix <- mRNA_data[,patient_id]
gene_name <- rownames(mRNA_data)
gene_symbol <- unlist(strsplit(gene_name,split="\\|"))[seq(1,length(gene_name)*2,by=2)]
gene_id <- unlist(strsplit(gene_name,split="\\|"))[seq(1,length(gene_name)*2,by=1)]
draw_ge_survival_plot(gene_name = "EGFR|1956",survival_data = survival_data,mRNA_matrix = mRNA_matrix)
mg_mRNA <- cox_marginal_correlation(matrix = mRNA_matrix,survival_data = survival_data)
save(mRNA_matrix,file = "mRNA_matrix.rda")
save(mg_mRNA,file = "mg_mRNA.rda")
##########################################
# read CNV files.
##########################################
rm(list=ls())
source("surTools.r")
load("survival_data.rda")
# read CNV data.
CNV_data <- read.table("gdac.broadinstitute.org_LUAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt",head=T,sep="\t")
# map gene_id data
cnv_matrix <- CNV_data[,-c(1:3)] # no NA entry.
rownames(cnv_matrix) <- paste(CNV_data[,1],CNV_data[,2],sep=";")
draw_cnv_survival_plot(gene_name = "EGFR;1956",survival_data = survival_data,cnv_matrix = cnv_matrix)
mg_cnv <- cox_marginal_correlation(matrix = cnv_matrix,survival_data = survival_data)
save(cnv_matrix,file = "cnv_matrix.rda")
save(mg_cnv,file = "mg_cnv.rda")
##########################################
# read methylation files.
##########################################
rm(list=ls())
source("surTools.r")
load("survival_data.rda")
# read methyl data.
me_data <- read.table("gdac.broadinstitute.org_LUAD.Methylation_Preprocess.Level_3.2016012800.0.0/LUAD.meth.by_mean.data.txt",sep="\t",header = T)
me_matrix <- me_data[-1,-1]
me_matrix <- apply(me_matrix, 2, as.numeric)
rownames(me_matrix) <- me_data[-1,1]
patient_id <- colnames(me_matrix)
# only select tumor data.
normal_tumor_indicator <- unlist(strsplit(patient_id,split = "\\."))[seq(4,length(patient_id)*4,by=4)]
tumor_idx <- which(normal_tumor_indicator=="01")
patient_id <- patient_id[tumor_idx]
length(patient_id)
me_matrix <- me_matrix[,patient_id]
# impute NA in me_matrix.
k <- which(is.na(me_matrix), arr.ind=TRUE)
me_matrix[k] <- rowMeans(as.matrix(me_matrix), na.rm=TRUE)[k[,1]]
draw_me_survival_plot(gene_name = "A2LD1",survival_data = survival_data,me_matrix = me_matrix)
mg_me <- cox_marginal_correlation(matrix = me_matrix,survival_data = survival_data)
save(mg_me,file="mg_me.rda")
save(me_matrix,file="me_matrix.rda")
##########################################
# read miRNA files.
# TBD: How to link miRNA to mRNA?
##########################################
rm(list=ls())
source("surTools.r")
load("survival_data.rda")
# read methyl data.
mi_data <- read.table("gdac.broadinstitute.org_LUAD.miRseq_Preprocess.Level_3.2016012800.0.0/LUAD.miRseq_RPKM_log2.txt",header = T,sep="\t")
mi_matrix <- mi_data[,-1]
rownames(mi_matrix) <- mi_data[,1]
patient_id <- colnames(mi_matrix)
# only select tumor data.
normal_tumor_indicator <- unlist(strsplit(patient_id,split = "\\."))[seq(4,length(patient_id)*4,by=4)]
tumor_idx <- which(normal_tumor_indicator=="01")
patient_id <- patient_id[tumor_idx]
length(patient_id)
mi_matrix <- mi_matrix[,patient_id]
# impute NA in me_matrix.
which(is.na(mi_matrix))
k <- which(is.na(mi_matrix), arr.ind=TRUE)
mi_matrix[k] <- rowMeans(as.matrix(mi_matrix), na.rm=TRUE)[k[,1]]
which(is.na(mi_matrix))
draw_mi_survival_plot(gene_name = "hsa-mir-548v",survival_data = survival_data,mi_matrix = mi_matrix)
mg_mi <- cox_marginal_correlation(matrix = mi_matrix,survival_data = survival_data)
save(mg_mi,file="mg_mi.rda")
save(mi_matrix,file="mi_matrix.rda")



