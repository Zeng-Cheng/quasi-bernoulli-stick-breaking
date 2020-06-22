rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
  stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
  x <- runif(n)
  Fa <- pbeta(a, alpha, beta)
  Fb <- pbeta(b, alpha, beta)
  y <- (1 - x) * Fa + x * Fb
  y[y < 1E-16] = 1E-16
  return(qbeta(y, alpha, beta))
}


w = rdirichlet(1, alpha = rep(1, K))
# b = rep(1,K)

eps <- 1 / (N ^ 1.01)
p_b = 0.9

################################




updateC <- function(mu, SigmasInv, w) {
  
  loglik_normal = matrix(0,N,K)
  
  for(k in c(1:K)){
    diff = y - rep(mu[k,],each=N)
    logdet1 = logdet(SigmasInv[,,k])
    loglik_normal[,k] = -diag(diff%*%SigmasInv[,,k]%*%t(diff))/2 + logdet1/2
  }
  
  loglik = t( t(loglik_normal) + log(c(w))) # R is adding by columns, so I used tranpose's twice 
  gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
  C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
  return(C)
}


updateMu <- function(C, SigmasInv, mu_prior_mean = 0, mu_prior_sigma2 = 1) {
  
  ySum = t(C)%*%y
  n_C = colSums(C)
  
  for(k in c(1:K)){
    m_part2 = ySum[k,]%*%SigmasInv[,,k]+ mu_prior_mean / mu_prior_sigma2
    
    m_var = solve(SigmasInv[,,k]*n_C[k] + diag(1,p)/mu_prior_sigma2)
    m_mean = m_var%*%t(m_part2)
    
    mu[k,] = t(chol(m_var)) %*%rnorm(p) + m_mean
  }
  
  
  return(mu)
}


updateSigmasInv <- function(C, mu) {
  
  C_label = colSums((t(C) * c(1:K)))
  
  diff = (y- mu[C_label,])
  
  for(k in c(1:K)){
    pick = (C_label==k)
    
    if(sum(pick)>1){
      par2 = solve( t(diff[pick,])%*% (diff[pick,])+ diag(p))
    }
    
    if(sum(pick)==1){
      vec = matrix(diff[pick,],1,p)
      par2 = solve( t(vec)%*% (vec)+ diag(p))
    }
    if(sum(pick)==0){
      par2 =  diag(p)
    }
    
    par2 = (par2+t(par2))/2
    par1 = sum(C_label==k) + p
    
    SigmasInv[,,k] = rwishart(par1,par2)
  }
  
  return(SigmasInv)
}

# updateBeta_sigma2 <- function(sigma2, sigma2_prior_par1 = 2, beta_sigma2_prior_par1 = 0.2, beta_sigma2_prior_par2 = 10) {
#   
#   par1_gamma <- K * sigma2_prior_par1 + beta_sigma2_prior_par1
#   par2_gamma <- beta_sigma2_prior_par2 + sum(1 / sigma2)
#   
#   beta_sigma2 <- rgamma(1, par1_gamma, rate = par2_gamma)
#   
#   if(any(is.na(beta_sigma2))) {
#     beta_sigma2 = par1_gamma / par2_gamma
#   }
#   return(beta_sigma2)
# }

updateBeta <- function(C, b, alpha_beta = 1, eps = 1E-3) {
  # N = dim(C)[1]
  # K = dim(C)[2]
  n_C <- colSums(C)
  par1_beta = (N - cumsum(n_C)) + alpha_beta
  par2_beta = n_C + 1
  
  beta_1 = rbeta(K, par1_beta, par2_beta)
  beta_2 = rtbeta(K, par1_beta, par2_beta, 0, eps) / eps #truncated Beta distribution
  
  beta = beta_1 * (b == 1) + beta_2 * (b != 1)
  
  # beta[K]=1E-8
  if (any(is.na(beta))) {
    # print("beta err")
    beta = updateBeta(C, b, alpha_beta, eps)
  }
  return(beta)
}


updateB <- function(C, eps = 1E-3, p_b = 1E-1, alpha_beta = 1) {
  # K = dim(C)[2]
  # N = dim(C)[1]
  
  n_C <- colSums(C)
  m_C = N - cumsum(n_C)
  
  choice1 = log(p_b)
  choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps)
  
  gumbel = -log(-log(runif(K * 2, min = 0, max = 1)))
  
  prob_choice = cbind(choice1, choice2) + gumbel
  b <- colSums((apply(prob_choice, 1, function(x) x == max(x))) * c(1, eps))
  return(b)
}

updateW <- function(b, beta) {
  v = 1 - b * beta
  # K = length(ubeta)
  w = v * (cumprod(c(1, 1 - v))[1:K])
  return(w)
}


partition_prob <- function(n_C, eps = 1E-3, p_b = 1E-1, alpha_beta = 1) {
  logsumexp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  
  m_C = N - cumsum(n_C)
  part1 = sum(lbeta(m_C + alpha_beta, n_C + 1))
  choice1 = log(p_b)
  choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps) 
  part2 = sum(apply(cbind(choice1, choice2), 1, logsumexp))
  return(part1 + part2)
}


MH_order_idx <- function(C, eps = 1E-3, p_b = 1E-1, alpha_beta = 1) {
  n_C <- colSums(C)
  idx_prop = order(n_C, decreasing = T)
  n_C_prop = n_C[idx_prop]
  
  if (log(runif(1)) < (partition_prob(n_C_prop, eps, p_b, alpha_beta) - partition_prob(n_C, eps, p_b, alpha_beta))) {
    idx = idx_prop
  }
  else {
    idx = c(1:K)
  }
  return(idx)
}

