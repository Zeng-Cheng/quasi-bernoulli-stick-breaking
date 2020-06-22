rm(list=ls())

setwd("~/git/quasi-bernoulli-stick-breaking/")
library(LaplacesDemon)


K = 20  # original number of clusters
N = 1000    # sample size

y = rt(N,df=3)

hist(y,breaks = 100)




source("qbsb_gibbs.R")

# initialize the parameters
mu = rep(0,K)

sigma2 = 1 / rgamma(K, 2, 1)
hist(y, breaks = 100)

w = rdirichlet(1, alpha = rep(1, K))
# b = rep(1,K)

eps <- 1 / N ^ 1.01
p_b = 0.9


n_C_trace=numeric()
sigma2_trace=numeric()

for (i in 1:10000){
  C <- updateC(mu, sigma2, w)
  # C_label = colSums((t(C) * c(1:K)))
  n_C <- colSums(C)
  # mu <- updateMu(C, sigma2, mu_prior_mean = (max(y) - min(y)) / 2, mu_prior_sigma2 = (max(y) - min(y)) ^ 2)
  beta_sigma2 <- updateBeta_sigma2(sigma2, sigma2_prior_par1 = 2, beta_sigma2_prior_par1 = 0.2, beta_sigma2_prior_par2 = 10 / (max(y) - min(y)) ^ 2)
  sigma2 <- updateSigma2(C, mu, sigma2_prior_par1 = 2, sigma_prior_par2 = beta_sigma2)
  # b <- updateB(C, eps = eps, p_b = p_b)
  # beta <- updateBeta(C, b, eps = eps)
  # w <- updateW(b, beta)
  
  # MH step to re-order components
  idx = MH_order_idx(C, eps = eps, p_b = p_b)
  mu = mu[idx]
  sigma2 = sigma2[idx]
  C = C[, idx]
  b <- updateB(C, eps = eps, p_b = p_b)
  beta <- updateBeta(C, b, eps = eps)
  w <- updateW(b, beta)
  
  if (i %% 100 == 0) {
    print(n_C)
  }
  
  if(i>1000){
    if(i%%10==0){
      n_C_trace= c(n_C_trace,sum(n_C>0))
      sigma2_trace=rbind(sigma2_trace,sigma2)
    }
  }
}


barplot(table(n_C_trace))


par(mfrow=c(1,2))
ts.plot(n_C_trace)
acf(n_C_trace,lag.max = 40)
