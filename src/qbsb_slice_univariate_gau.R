require("LaplacesDemon")

runQBSB_slice_efficient_UvGaussian <- function(
    y, xi, steps = 5E4, burnins = 1E4,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    mu_prior_mean = (max(y) + min(y)) / 2, mu_prior_sigma2 = (max(y) - min(y)) ^ 2, sigma2_prior_par1 = 2,
    beta_sigma2_prior_par1 = 0.2, beta_sigma2_prior_par2 = 10 / (max(y) - min(y)) ^ 2
    ) {

    stopifnot(eps > 0 & (alpha_beta + d_beta) > 0 & d_beta < 1 & d_beta >= 0 & p_b >= 0 & p_b <= 1)
    stopifnot(p_b == 1 | d_beta == 0)
    
    y = c(y)
    N = length(y) 
    rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
        stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
        x <- runif(n)
        Fa <- pbeta(a, alpha, beta)
        Fb <- pbeta(b, alpha, beta)
        y <- (1 - x) * Fa + x * Fb
        y[y < 1E-16] = 1E-16
        return(qbeta(y, alpha, beta))
    }
 
    updateC <- function(u, mu, sigma2, w) {
        diff2 = (outer(mu, y, "-")) ** 2
        diff2_sigma2 = diff2 / sigma2
        loglik_normal = -diff2_sigma2 / 2 - log(sigma2) / 2

        loglik = t(loglik_normal + log(outer(xi[1:K], u, ">")) + log(c(w)) - log(xi[1:K]))
        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
        return(C)
    }
    
    updateMu <- function(C, sigma2) {
        
        m_var = 1 / (colSums(C) / sigma2 + 1 / mu_prior_sigma2)
        m_part2 = c(y %*% C / sigma2 + mu_prior_mean / mu_prior_sigma2)
        m_mean = m_var * m_part2
        
        mu = rnorm(K, m_mean, sqrt(m_var))

        return(mu)
    }
       
    updateSigma2 <- function(C, mu, sigma2_prior_par2 = 1) {
        
        C_label = colSums((t(C) * c(1:K)))       
        par1_ig = colSums(C) / 2 + sigma2_prior_par1
        par2_ig = c(((y - mu[C_label]) ** 2) %*% C) / 2 + sigma2_prior_par2
        
        sigma2 = 1 / rgamma(K, par1_ig, rate = par2_ig)       
        if (any(is.na(sigma2))) {
            sigma2 = par2_ig / par1_ig
        }
        
        return(sigma2)
    }
    
    updateBeta_sigma2 <- function(sigma2) {
        
        par1_gamma <- K * sigma2_prior_par1 + beta_sigma2_prior_par1
        par2_gamma <- beta_sigma2_prior_par2 + sum(1 / sigma2)
        
        beta_sigma2 <- rgamma(1, par1_gamma, rate = par2_gamma)
        
        if (any(is.na(beta_sigma2))) {
            beta_sigma2 = par1_gamma / par2_gamma
        }
        return(beta_sigma2)
    }
    
    updateBeta <- function(C, b) {

        n_C <- colSums(C)
        par1_beta = (N - cumsum(n_C)) + alpha_beta + d_beta * c(1:K)
        par2_beta = n_C + 1 - d_beta
        
        beta_1 = rbeta(K, par1_beta, par2_beta)
        if (p_b == 1) {
            beta = beta_1
            if (any(is.na(beta))) {
                # print("beta err")
                beta = updateBeta(C, b, alpha_beta, eps, d_beta)
            }
            return(beta)
        }
        beta_2 = rtbeta(K, par1_beta, par2_beta, 0, eps) / eps # truncated Beta distribution
        
        beta = beta_1 * (b == 1) + beta_2 * (b != 1)
        # beta[K]=1E-8
        if (any(is.na(beta))) {
        # print("beta err")
            beta = updateBeta(C, b, alpha_beta, eps, d_beta)
        }
        return(beta)
    }
    
    
    updateB <- function(C) {
        
        n_C <- colSums(C)
        m_C = N - cumsum(n_C)
        
        choice1 = log(p_b)
        choice2 = log(1 - p_b) + pbeta(eps, alpha_beta + m_C, n_C + 1, log.p = T) - alpha_beta * log(eps)
        
        gumbel = -log(-log(runif(K * 2, min = 0, max = 1)))
        
        prob_choice = cbind(choice1, choice2) + gumbel
        b <- colSums((apply(prob_choice, 1, function(x) x == max(x))) * c(1, eps))
        return(b)
    }

    updateU <- function(C) {

        C_label = colSums((t(C) * c(1:K)))

        u = runif(N, max = xi[C_label])
        #if(1 - sum(w) > min(u)) print("ERROR! Sum W greater than min U")
        return(u)
    }
    
    updateW <- function(b, beta) {
        v = 1 - b * beta
        w = v * (cumprod(c(1, 1 - v))[1:K])
        return(w)
    }

    updateK <- function(b, beta, u, mu, sigma2, sigma2_prior_par2 = 1) {
        k = length(xi[xi > min(u)])
        if(k == length(xi)) stop("xi is too short")
        if(k <= K) {
            b = b[1:k]
            beta = beta[1:k]
            mu = mu[1:k]
            sigma2 = sigma2[1:k]
        }
        else {
            b = c(b, eps + (1 - eps) * rbinom(k - K, 1, p_b))
            beta = c(beta, rbeta(k - K, alpha_beta + d_beta * c((K + 1):k), 1 - d_beta))
            mu = c(mu, rnorm(k - K, mu_prior_mean, sqrt(mu_prior_sigma2)))
            sigma2 = c(sigma2, 1 / rgamma(k - K, sigma2_prior_par1, rate = sigma2_prior_par2))
        }
        return(list(k = k, b = b, beta = beta, mu = mu, sigma2 = sigma2))
    }
    
    
    
    # initialize the parameters
    K = 20 # will be changed for slice sampler
    mu = rnorm(K, 0, 1)
    sigma2 = 1 / rgamma(K, 2, 1)
    beta_sigma2 = 1
    w = rdirichlet(1, alpha = rep(1, K))
    b = rep(1, K)
    u = rep(0.01, N)

    # record trace
    trace_nC = list()
    starttime = proc.time()
    
    for (i in 1:steps){
        
        C <- updateC(u, mu, sigma2, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        mu <- updateMu(C, sigma2)
        sigma2 <- updateSigma2(C, mu, sigma2_prior_par2 = beta_sigma2)
        beta_sigma2 <- updateBeta_sigma2(sigma2)    

           
        if(p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)

        u <- updateU(C)
        uK_res = updateK(b, beta, u, mu, sigma2, sigma2_prior_par2 = beta_sigma2)
        K <- uK_res$k
        b <- uK_res$b
        beta = uK_res$beta
        mu = uK_res$mu
        sigma2 = uK_res$sigma2
        w <- updateW(b, beta)

            
        if (i > burnins) {
            idx = i - burnins
            trace_nC[[idx]] = n_C
        }

        if(i == 5000) print((proc.time() - starttime)[3] / 5000)
        if(i %% 5000 == 0) print(n_C)
        
    }   
    
    M = length(trace_nC)
    nT <- sapply(1:M, function(i) sum(trace_nC[[i]] != 0)) # number of clusters
    return(nT)
}
