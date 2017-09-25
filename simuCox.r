rm(list=ls())
library(stats)
library(survsim)
n = 500
p = 20000
beta = rep(0,p)
beta[1:10] = runif(n = 10,-1,-0.1)
beta[11:20] = runif(n = 10,0.1,1)
generateTtime <- function(x,beta,lambda,v){
  # lambda is the scale.
  # v is the shape.
  # of baseline hazard ratio.
  u <- runif(n=1,0,1)
  T <- (log(u)/(-exp(sum(x*beta)*lambda)))^(1/v)
  return(T)
}
maxCor <- function(x,ref){
  A = apply(ref,2,cor,y=y)
  return(max(A))
}

# A Weibull distribustion with the shape parameter of 5 
# and the scale parameter of 2 is used for the baseline hazard function, 
# and a uniform U (2, 10) is used for simulating the censoring times.
shape = 5
scale = 2
v = shape
# lambda = scale^shape
lambda = 5
C_time = round(runif(n,2,10)*12)

# generate X matrix.
# we first generate a 500 × 10000 dataset X from a uniform
# U(−1.5,1.5) distribution.
# generate X1:20
X = matrix(runif(n = n*p,-1.5,1.5),nrow = n,ncol = p)
X.qr = qr(X)
Q <- qr.Q(X.qr)
alpha <- Q[,1:20]
gamma <- Q[,21:n]
# generate X21:n w/ correlation of X1:20.
N = floor((p-20)/(n-20))
for(i in 1:N){
  # the correlation is controled by eigenvalue of T'T.
  #T <- matrix(0,nrow = 20,ncol = n-20) w/o correlation.
  T = matrix(runif(n = 20*(n-20),-1.5,1.5),nrow = 20,ncol = n-20)
  C = gamma + alpha %*% T
  X[,c(21+(i-1)*480):(20+i*480)] = C
}
T = matrix(runif(n = 20*(n-20),-1.5,1.5),nrow = 20,ncol = n-20)
C = gamma + alpha %*% T
X[,c((21+N*480):dim(X)[2])] = C[,1:length(dim(X)[2] -(21+N*480))]
dim(X)
#max_cor = apply(C, 2, FUN = maxCor, ref = X)

T_time = apply(X,1,generateTtime,beta=beta,lambda=lambda,v=v)
T_time = round((T_time+1)*12)
T_time
tmp <- cbind(T_time,C_time)
y <- apply(tmp,1,min)
delta <- y==T_time
survival <- data.frame(C_time,T_time,y,delta)
head(survival)
sum(delta)/n
simu_data_cor05_X1_b1 <- list(survival=survival,X=X,lambda=lambda,v=v,beta=beta)
save(simu_data_cor05_X1_b1,file = "simu_data_cor05_X1_b1.rda")
