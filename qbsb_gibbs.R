rm(list=ls())

rtbeta <- function(n, alpha, beta, a=0, b=1)
{
  stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
  x <- runif(n)
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  y <- (1-x)*Fa + x*Fb
  return(qbeta(y, alpha, beta))
}



rtbeta(100,1,100,0, 1E-8)


setwd("~/git/sparse_dirichlet/")
library(LaplacesDemon)

K.true = 3
K = 10  # original number of clusters
N = 3000    # sample size
y <- rnormm(N, p = c(0.3, 0.3, 0.4), mu = c(-3, 0, 4), sigma = c(1, 1, 1))


w = rdirichlet(1, alpha = rep(1,K))
s = 0.1

u = rep(1,K)
beta = rbeta(K,1,s)
beta

mu = rnorm(K)
sigma2 = 1/rgamma(K,2,1)

hist(y,breaks = 100)

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

updateBeta<- function(C,u,  eps = 1E-12){
  
  n_C <- colSums(C)
  alpha_beta = 1
  par1_beta = n_C + 1
  par2_beta = (N - cumsum(n_C)) + alpha_beta
  
  beta_1 = rbeta(K, par1_beta, par2_beta)
  beta_2 = rtbeta(K, par1_beta, par2_beta,0,eps)/eps
  
  beta= beta_1*(u==1) + beta_2*(u!=1)
  
  # beta[K]=1E-8
  if(any(is.na(beta))){
    # print("beta err")
    beta = updateBeta(C,u)
  }
  return(beta)
}


updateU<- function(C, beta,   eps = 1E-12) {
  alpha_beta = 1
  p_u = 1E-3
  n_C <- colSums(C)
  n_gtK = N - cumsum(n_C)
  
  choice1 = log(p_u) 
  choice2 = log(1-p_u) + pbeta(eps, n_C+1, alpha_beta+ n_gtK,log.p =T) - log(eps) 
  
  gumbel = -log(-log(runif(K*2,min = 0,max=1)))
  
  prob_choice = cbind(choice1,choice2)  + gumbel
  u<- colSums((apply(prob_choice,1,function(x)x==max(x))) * c(1,eps))
  return(u)
}

updateW<- function(u,beta){
  
  ubeta = u*beta
  
  w = ubeta * (cumprod(c(1, 1-ubeta))[1:K])

  return(w)
  
}

C<- updateC(mu,sigma2,w)
C_label = colSums((t(C)*c(1:K)))
n_C <- colSums(C)
mu<- updateMu(C,sigma2)
sigma2<- updateSigma2(C,mu)

beta<- updateBeta(C,u)
u<- updateU(C,beta)
w<- updateW(u,beta) 



for (i in c(1:5000)){
  C<- updateC(mu,sigma2,w)
  C_label = colSums((t(C)*c(1:K)))
  n_C <- colSums(C)
  mu<- updateMu(C,sigma2)
  sigma2<- updateSigma2(C,mu)
  
  beta<- updateBeta(C,u)
  u<- updateU(C,beta)
  w<- updateW(u,beta) 
  if(i %% 100==0){
    print(n_C)
    print(sum(u))
  }
}

sum(u)

plot(w)
