require("LaplacesDemon")
require("msm")

runQBSBUvLaplace <- function(
    y, steps = 5E4, burnins = 3E4, K = 40,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    mu_prior_mean = (max(y) + min(y)) / 2, mu_prior_sigma2 = (max(y) - min(y)) ^ 2,
    lambda_prior_par1 = 2, lambda_prior_par2 = 1
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
 
    updateC <- function(mu, lambda, w) {
        diff2 = abs(outer(mu, y, "-"))
        diff2_lambda = diff2 / lambda
        loglik_normal = -diff2_lambda - log(lambda)

        loglik = t(loglik_normal + log(c(w)))
        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
        return(C)
    }
    
    updateMu <- function(C, lambda) {
        
        for(k in 1:K) {
            x <- sort(y[C[, k]])
            nk <- sum(C[, k])
            lower <- c(-Inf, x)
            upper <- c(x, Inf)
            components_mean = (nk - 2 * c(0:nk)) * mu_prior_sigma2 / lambda[k] + mu_prior_mean
            components_truncate <- rtnorm(nk + 1, components_mean, sqrt(mu_prior_sigma2), lower, upper)

            w_part1 = (components_mean ^ 2 - mu_prior_mean ^ 2) / 2 / mu_prior_sigma2
            w_part2 = log(pnorm(upper, components_mean, sqrt(mu_prior_sigma2)) - pnorm(lower, components_mean, sqrt(mu_prior_sigma2)))
            w_part3 <- (2 * c(0, cumsum(x)) - sum(x)) / lambda[k]
            weights <- w_part1 + w_part2 + w_part3

            gumbel = -log(-log(runif(nk + 1, min = 0, max = 1)))

            mu[k] <- components_truncate[which.max(weights + gumbel)]
        }

        return(mu)
    }
       
    updateLambda <- function(C, mu) {
        
        C_label = colSums((t(C) * c(1:K)))       
        par1_ig = colSums(C) + lambda_prior_par1
        par2_ig = c(abs(y - mu[C_label]) %*% C) + lambda_prior_par2
        
        lambda = 1 / rgamma(K, par1_ig, rate = par2_ig)       
        if (any(is.na(lambda))) {
            lambda = par2_ig / par1_ig
        }
        
        return(lambda)
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
    
    updateW <- function(b, beta) {
        v = 1 - b * beta
        w = v * (cumprod(c(1, 1 - v))[1:K])
        return(w)
    }
    

    # initialize the parameters
    mu = rep(0, K)
    lambda = 1 / rgamma(K, 2, 1)
    w = rdirichlet(1, alpha = rep(1, K))
    b = rep(1, K)

    # record trace
    trace_mu = list()
    trace_lambda = list() 
    trace_nC = list()
    trace_label = list()
    trace_b = list()
    
    for (i in 1:steps){
        
        C <- updateC(mu, lambda, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        mu <- updateMu(C, lambda)
        lambda <- updateLambda(C, mu)
        
        if (p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)
        
        if (i > burnins) {
            idx = i - burnins
            trace_mu[[idx]] = mu
            trace_b[[idx]] = b
            trace_lambda[[idx]] = lambda
            trace_nC[[idx]] = n_C
            trace_label[[idx]] = C_label
        }

        if(i %% 500 == 0) print(n_C)
    }
    return(list(mu = trace_mu, b = trace_b, lambda = trace_lambda, nC = trace_nC, label = trace_label))  
}
