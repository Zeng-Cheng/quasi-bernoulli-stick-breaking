rm(list=ls())

setwd("~/git/quasi-bernoulli-stick-breaking/")
library(LaplacesDemon)


K = 10  # original number of clusters
N = 1000    # sample size

p = 2


mu_truth = cbind(c(-4, 0, 5),c(1, 2, 3))


n_truth = round(N* c(0.3, 0.3, 0.4))
y = matrix(0, N,p)

for(i in 1:N){
  y[i,] = mu_truth[sum(i >cumsum(n_truth))+1,] + rnorm(p)
}

plot(y[,1],y[,2])


source("qbsb_gibbs_multivariate_gau.R")

# initialize the parameters
mu = matrix(rnorm(K*p),ncol = p)
I_p = diag(1,p)

SigmasInv = array(0,dim=c(p,p,K))

for(k in c(1:K)){
  SigmasInv[,,k] = rwishart(p,I_p)
}


for (i in 1:10000){

  C = updateC(mu, SigmasInv, w)
  C_label = colSums((t(C) * c(1:K)))
  n_C <- colSums(C)
  
  mu = updateMu(C,SigmasInv)
  SigmasInv= updateSigmasInv(C, mu)
  
  b <- updateB(C, eps = eps, p_b = p_b)
  beta <- updateBeta(C, b, eps = eps)
  w <- updateW(b, beta)
  
  # MH step to re-order components
  idx = MH_order_idx(C, eps = eps, p_b = p_b)
  mu = mu[idx,]
  SigmasInv = SigmasInv[,,idx]
  C = C[, idx]
  # b <- updateB(C, eps = eps, p_b = p_b)
  # beta <- updateBeta(C, b, eps = eps)
  # w <- updateW(b, beta)
  
  if (i %% 1 == 0) {
    print(n_C)
  }
}
