require("LaplacesDemon")
require("msm")

runQBSB_slice_efficient_UvLaplace <- function(
    y, xi, steps = 5E4, burnins = 2E4,
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
 
    updateC <- function(u, mu, lambda, w) {
        diff2 = abs(outer(mu, y, "-"))
        diff2_lambda = diff2 / lambda
        loglik_normal = -diff2_lambda - log(lambda)

        loglik = t(
            loglik_normal +
            log(outer(xi[1:K], u, ">")) +
            log(c(w)) -
            log(xi[1:K])
            )
        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) {x == max(x)}))
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
    
    updateK <- function(b, beta, u, mu, lambda) {
        k = length(xi[xi > min(u)])
        if(k == length(xi)) stop("xi is too short")
        if(k <= K) {
            b = b[1:k]
            beta = beta[1:k]
            mu = mu[1:k]
            lambda = lambda[1:k]
        }
        else {
            b = c(b, eps + (1 - eps) * rbinom(k - K, 1, p_b))
            beta = c(beta, rbeta(k - K, alpha_beta + d_beta * c((K + 1):k), 1 - d_beta))
            mu = c(mu, rnorm(k - K, mu_prior_mean, sqrt(mu_prior_sigma2)))
            lambda = c(lambda, 1 / rgamma(k - K, lambda_prior_par1, rate = lambda_prior_par2))
        }
        return(list(k = k, b = b, beta = beta, mu = mu, lambda = lambda))
    }
    

    
    # initialize the parameters
    K = 20
    mu = rep(0, K)
    lambda = 1 / rgamma(K, 2, 1)
    w = rdirichlet(1, alpha = rep(1, K))
    b = rep(1, K)
    u = rep(0.01, N)

    # record trace
    trace_nC = list()
    starttime = proc.time()
    for (i in 1:steps){
        
        C <- updateC(u, mu, lambda, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        mu <- updateMu(C, lambda)
        lambda <- updateLambda(C, mu)

        if(p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)

        u <- updateU(C)
        uK_res = updateK(b, beta, u, mu, lambda)
        K <- uK_res$k
        b <- uK_res$b
        beta = uK_res$beta
        mu = uK_res$mu
        lambda = uK_res$lambda
        w <- updateW(b, beta)

                 
        if (i > burnins) {
            idx = i - burnins
            trace_nC[[idx]] = n_C
        }

        if(i == 1000) print((proc.time() - starttime)[3] / 1000)
        if(i %% 1000 == 0) print(n_C)
    }
    M = length(trace_nC)
    nT <- sapply(1:M, function(i) sum(trace_nC[[i]] != 0)) # number of clusters
    return(nT)
}
