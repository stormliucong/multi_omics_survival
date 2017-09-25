# contain all function may need.
library(glmnet)
library(survival)
library(SIS)
library(ggplot2)
patient_id_convert <- function(x){
  tmp <- unlist(strsplit(x,"\\."))
  tmp[1] <- tolower(tmp[1])
  id <- paste(tmp[1],tmp[2],tmp[3],sep = ".")
  return(id)
}

cox_marginal_correlation <- function(matrix,survival_data){
  matrix <- t(matrix)
  # convert patient id.
  rownames(matrix) <- unlist(lapply(rownames(matrix),patient_id_convert))
  # X each row is a sample with rownames as TCGA sample name.
  ref_idx <- match(rownames(matrix),survival_data$patient_id)
  survival_data_map <- survival_data[ref_idx,]
  NA_idx <- which(is.na(survival_data_map$patient_id)) 
  survival_data_map <- survival_data_map[-NA_idx,]
  matrix <- matrix[-NA_idx,]
  
  y <- Surv(survival_data_map$time,survival_data_map$status)
  # for coxph model fitting.
  mg <- rep(NA,dim(matrix)[2])
  p <- rep(NA,dim(matrix)[2])
  for(i in 1:dim(matrix)[2]){
    fit <- coxph(y ~ matrix[,i])
    summary <- summary(fit)
    coef <- summary$coefficients[1]
    pval <- summary$coefficients[5]
    mg[i] <- coef
    p[i] <- pval
  }
  q <- p.adjust(p,method = "fdr")
  marginal_correlation <- data.frame(gene=colnames(matrix),mg=mg,p=p,q=q)
  return(marginal_correlation)
}

draw_cnv_survival_plot <- function(gene_name,survival_data,cnv_matrix){
  cnv <- cnv_matrix[gene_name,]
  cnv <- t(cnv)
  rownames(cnv) <- unlist(lapply(rownames(cnv),patient_id_convert))
  ref_idx <- match(rownames(cnv),survival_data$patient_id)
  survival_data_map <- survival_data[ref_idx,]
  NA_idx <- which(is.na(survival_data_map$patient_id)) 
  survival_data_map <- survival_data_map[-NA_idx,]
  cnv <- cnv[-NA_idx,]
  dat <- data.frame(cnv = cnv,time=survival_data_map$time,status=as.factor(survival_data_map$status))
  p <- ggplot(dat,aes(x = cnv,y = time))
  p <- p + geom_point(aes(colour = as.factor(status)))
  p
  return(p)
}
draw_ge_survival_plot <- function(gene_name,survival_data,mRNA_matrix){
  ge <- mRNA_matrix[gene_name,]
  ge <- t(ge)
  rownames(ge) <- unlist(lapply(rownames(ge),patient_id_convert))
  ref_idx <- match(rownames(ge),survival_data$patient_id)
  survival_data_map <- survival_data[ref_idx,]
  NA_idx <- which(is.na(survival_data_map$patient_id)) 
  survival_data_map <- survival_data_map[-NA_idx,]
  ge <- ge[-NA_idx,]
  dat <- data.frame(ge = ge,time=survival_data_map$time,status=as.factor(survival_data_map$status))
  p <- ggplot(dat,aes(x = ge,y = time))
  p <- p + geom_point(aes(colour = as.factor(status)))
  p
  return(p)
}
draw_me_survival_plot <- function(gene_name,survival_data,me_matrix){
  colnames(me_matrix) <- unlist(lapply(colnames(me_matrix),patient_id_convert))
  ref_idx <- match(colnames(me_matrix),survival_data$patient_id)
  survival_data_map <- survival_data[ref_idx,]
  NA_idx <- which(is.na(survival_data_map$patient_id)) 
  survival_data_map <- survival_data_map[-NA_idx,]
  me_matrix <- me_matrix[,-NA_idx]
  me <- me_matrix[gene_name,]
  dat <- data.frame(me = as.numeric(me),time=as.numeric(survival_data_map$time),status=as.factor(survival_data_map$status))
  p <- ggplot(dat,aes(x = me,y = time))
  p <- p + geom_point(aes(colour = as.factor(status)))
  p
  return(p)
}
draw_mi_survival_plot <- function(gene_name,survival_data,mi_matrix){
  mi <- mi_matrix[gene_name,]
  mi <- t(mi)
  rownames(mi) <- unlist(lapply(rownames(mi),patient_id_convert))
  ref_idx <- match(rownames(mi),survival_data$patient_id)
  survival_data_map <- survival_data[ref_idx,]
  NA_idx <- which(is.na(survival_data_map$patient_id)) 
  survival_data_map <- survival_data_map[-NA_idx,]
  mi <- mi[-NA_idx,]
  dat <- data.frame(mi = mi,time=survival_data_map$time,status=as.factor(survival_data_map$status))
  p <- ggplot(dat,aes(x = mi,y = time))
  p <- p + geom_point(aes(colour = as.factor(status)))
  p
  return(p)
}

