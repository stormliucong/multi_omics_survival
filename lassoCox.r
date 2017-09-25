rm(list=ls())
library(glmnet)
library(survival)

riskScoreCal <- function(x,beta){
  return(sum(x%*%beta))
}
load("simu_data_cor05_X1_b1.rda")
simu_data <- simu_data_cor05_X1_b1
X = simu_data$X
y = Surv(simu_data$survival$y,simu_data$survival$delta)
beta <- simu_data$beta
b = 0
boostrap = 1
pval <- rep(NA,boostrap)
tp <- rep(NA,boostrap)
tn <- rep(NA,boostrap)
fp <- rep(NA,boostrap)
fn <- rep(NA,boostrap)
while(b < boostrap){
  b <- b + 1
  cat(b)
  train_idx <- sample(1:dim(X)[1],replace = F,size = round(dim(X)/5*4))
  cv.fit <- cv.glmnet(x = X[train_idx,],y = y[train_idx,],family = "cox")
  beta <- coef(cv.fit, s = "lambda.1se")
  p <- sum((beta!=0))
  n <- sum(beta==0)
  tp <- sum(which(beta!=0) %in% c(1:20))
  fp <- p - tp
  tn <- sum(which(beta==0) %in% c(21:10000))
  fn <- n - tn
  tp[b] <- tp
  tn[b] <- tn
  fp[b] <- fp
  fn[b] <- fn
  # test result 
  risk_score <- apply(X[-train_idx,],1,riskScoreCal,beta=beta)
  low_risk_group <- c((1:dim(X)[1])[-train_idx])[which(risk_score < median(risk_score))]
  high_risk_group <- c((1:dim(X)[1])[-train_idx])[which(risk_score >= median(risk_score))]
  group <- rep(NA,dim(X)[1]-length(train_idx))
  ytest <- Surv(c(y[low_risk_group,1],y[high_risk_group,1]),c(y[low_risk_group,2],y[high_risk_group,2]))
  group <- c(rep("low",length(low_risk_group)),rep("high",length(high_risk_group)))
  low_fit <- survfit(y[low_risk_group] ~ 1)
  high_fit <- survfit(y[high_risk_group] ~ 1)
  group.fit <- survdiff(ytest~group)
  pval[b] <- 1 - pchisq(group.fit$chisq, length(group.fit$n) - 1)
  plot(low_fit,conf.int = "none",mark.time = T, col = 'blue', xlab ='Time', ylab = 'Survival Probability')
  lines(high_fit, conf.int = "none",mark.time = T,col = 'red')
  legend("topright",c('high risk', 'low risk'), col = c('red','blue'), lty = 1)
}


