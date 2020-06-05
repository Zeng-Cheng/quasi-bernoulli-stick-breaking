rm(list=ls())

rtbeta <- function(n, alpha, beta,a=0, b=1)
{
  stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
  x <- runif(n)
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  y <-  x*Fb
  y[y<1E-16]=1E-16
  return(qbeta(y, alpha, beta))
}

# rtbeta(K, par1_beta, par2_beta,0,eps)



setwd("~/git/sparse_dirichlet/")
library(LaplacesDemon)

K.true = 3
K = 30  # original number of clusters
N = 1000    # sample size
y <- rnormm(N, p = c(0.3, 0.3, 0.4), mu = c(-3, 0, 4), sigma = c(1, 1, 1))

# initialize the parameters
mu = rnorm(K)
sigma2 = 1/rgamma(K,2,1)
hist(y,breaks = 100)

w = rdirichlet(1, alpha = rep(1,K))
s = 0.1

b = rep(1,K)
beta = rbeta(K,1,s)
beta



################################


updateC <- function(mu,sigma2,w){
  diff2 = (outer(mu,y,"-"))**2
  diff2_sigma2 = diff2/sigma2
  loglik_normal = - diff2_sigma2/2 - log(sigma2)/2
  
  loglik = t(loglik_normal + log(c(w)))
  gumbel = -log(-log(runif(N*K,min = 0,max=1)))
  C<- t(apply(loglik + gumbel, 1,function(x){x==max(x)}))
  return(C)
}


updateMu <- function(C,sigma2){
  mu_prior_mean = 0
  mu_prior_sigma2 = 1
  
  m_var = 1/(   colSums(C)/sigma2 + 1/mu_prior_sigma2)
  m_part2 = c(y%*%C/sigma2 +  mu_prior_mean/mu_prior_sigma2)
  m_mean = m_var*m_part2
  
  mu = rnorm(K, m_mean, sqrt(m_var))
  
  return(mu)
}


updateSigma2<- function(C,mu){
  
  C_label = colSums((t(C)*c(1:K)))
  
  par1_ig = colSums(C)/2 + 2
  par2_ig =  c(((y - mu[C_label])**2)%*%C) / 2 + 1
  
  sigma2 = 1/rgamma(K, par1_ig, rate = par2_ig)
  
  if(any(is.na(sigma2))){
    sigma2 = par2_ig/par1_ig
  }
  
  return(sigma2)
}

updateBeta<- function(C, b, alpha_beta = 1 ,eps   = 1E-3){
  
  n_C <- colSums(C)
  par1_beta = (N - cumsum(n_C)) + alpha_beta 
  par2_beta = n_C + 1
  
  beta_1 = rbeta(K, par1_beta, par2_beta)
  beta_2 = rtbeta(K, par1_beta, par2_beta,0,eps)/eps
  
  beta= beta_1*(b==1) + beta_2*(b!=1)
  
  # beta[K]=1E-8
  if(any(is.na(beta))){
    # print("beta err")
    beta = updateBeta(C,b, alpha_beta,eps)
  }
  return(beta)
}


updateB<- function(C,    eps=1E-3, p_b = 1E-1,alpha_beta = 1) {
  
  n_C <- colSums(C)
  m_C = N - cumsum(n_C)
  
  choice1 = log(p_b)
  choice2 = log(1-p_b) + pbeta(eps, alpha_beta+ m_C, n_C+1, log.p =T) - alpha_beta * log(eps) 
  
  gumbel = -log(-log(runif(K*2,min = 0,max=1)))
  
  prob_choice = cbind(choice1,choice2)  + gumbel
  b<- colSums((apply(prob_choice,1,function(x)x==max(x))) * c(1,eps))
  return(b)
}


updateW<- function(b,beta){
  
  v = 1- b*beta
  w = v * (cumprod(c(1, 1-v))[1:K])

  return(w)
  
}






eps = 1/N**(1.5)

C<- updateC(mu,sigma2,w)
C_label = colSums((t(C)*c(1:K)))
n_C <- colSums(C)
mu<- updateMu(C,sigma2)
sigma2<- updateSigma2(C,mu)
b<- updateB(C, eps=eps)
beta<- updateBeta(C,b,eps = eps)
w<- updateW(b,beta) 





partition_prob <- function(n_C, eps=1E-3, p_b = 1E-1,alpha_beta = 1 ){
  logsumexp<- function(x) log(sum(exp(x - max(x)))) + max(x)
  
  
  m_C = N - cumsum(n_C)
  part1 = sum(lbeta(m_C+alpha_beta, n_C+1))
  choice1 = log(p_b)
  choice2 = log(1-p_b) + pbeta(eps, alpha_beta+ m_C, n_C+1, log.p =T) - alpha_beta * log(eps) 
  part2 = sum( apply(cbind(choice1,choice2),1,logsumexp))
  
  part1+part2
}



MH_order_idx<- function(C, eps=1E-3, p_b = 1E-1,alpha_beta = 1 ){
  
  n_C <- colSums(C)
  
  idx_prop = order(n_C, decreasing = T)
  n_C_prop = n_C[idx_prop]
  
  if (log(runif(1)< (partition_prob(n_C) - partition_prob(n_C_prop)))){
    idx = idx_prop
  }else{
    idx= c(1:K)
  }
  
  idx
}


for (i in c(1:5000)){
  
  C<- updateC(mu,sigma2,w)
  C_label = colSums((t(C)*c(1:K)))
  n_C <- colSums(C)
  mu<- updateMu(C,sigma2)
  sigma2<- updateSigma2(C,mu)
  b<- updateB(C, eps=eps)
  beta<- updateBeta(C,b,eps = eps)
  w<- updateW(b,beta) 
  
  # MH step to re-order components
  {
  idx = MH_order_idx(C,eps,p_b,alpha_beta)
  mu = mu[idx]
  sigma2 = sigma2[idx]
  C = C[,idx]
  b<- updateB(C, eps=eps)
  beta<- updateBeta(C,b,eps = eps)
  w<- updateW(b,beta) 
  }
  
  if(i %% 100==0){
    print(n_C)
    print(sum(b))
  }
}
