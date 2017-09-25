rm(list=ls())
source("surTools.r")
load("cnv_matrix.rda")
load("survival_data.rda")
load("mRNA_matrix.rda")
load("me_matrix.rda")
load("mg_mRNA.rda")
load("mg_cnv.rda")
load("mg_me.rda")

mg_mRNA[mg_mRNA$gene == "FLJ40504|284085",]
mg_cnv[mg_cnv$gene == "EGFR;1956",]
mg_me[mg_me$gene == "KIAA1841",]

mg_mRNA[mg_mRNA$gene == "PRDM16|63976",]
mg_cnv[mg_cnv$gene == "PRDM16;63976",]
mg_me[mg_me$gene == "PRDM16",]

draw_me_survival_plot(gene_name = "KIAA1841",survival_data = survival_data,me_matrix = me_matrix)
draw_cnv_survival_plot(gene_name = "PRDM16;63976",survival_data = survival_data,cnv_matrix = cnv_matrix)
draw_ge_survival_plot(gene_name = "FLJ40504|284085",survival_data = survival_data,mRNA_matrix = mRNA_matrix)


x <-matrix[,"PRDM16|63976"]
survival_data_map
y <- Surv(survival_data_map$time,survival_data_map$status)
fit <- coxph(y ~ x)
