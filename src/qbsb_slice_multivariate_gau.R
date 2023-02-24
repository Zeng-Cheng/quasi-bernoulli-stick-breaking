require("LaplacesDemon")

runQBSB_slice_efficient_MvGaussian <- function(
    y, xi, steps = 5E4, burnins = 2E4,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    mu_prior_mean = colMeans(y), mu_prior_SigmasInv = solve(cov(y)),
    SigmasInv_prior_nu = ncol(y), SigmasInv_prior_VInv = SigmasInv_prior_nu * cov(y)
    ) {

    stopifnot(eps > 0 & (alpha_beta + d_beta) > 0 & d_beta < 1 & d_beta >= 0 & p_b >= 0 & p_b <= 1)
    stopifnot(p_b == 1 | d_beta == 0)

    N = nrow(y)
    p = ncol(y)
    rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
        stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
        x <- runif(n)
        Fa <- pbeta(a, alpha, beta)
        Fb <- pbeta(b, alpha, beta)
        y <- (1 - x) * Fa + x * Fb
        y[y < 1E-16] = 1E-16
        return(qbeta(y, alpha, beta))
    }
    
    updateC <- function(u, mu, SigmasInv, w) {
        
        loglik_normal = matrix(0, N, K)
        
        for(k in c(1:K)) {
            diff = y - rep(mu[k, ], each = N)
            logdet1 = logdet(SigmasInv[, , k])
            loglik_normal[, k] = -diag(diff %*% SigmasInv[, , k] %*% t(diff)) / 2 + logdet1 / 2
        }
        
        loglik = t(t(loglik_normal) + log(outer(xi[1:K], u, ">")) + log(c(w)) - log(xi[1:K])) # R is adding by columns, so I used tranpose's twice 
        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
        return(C)
    }
    
    
    updateMu <- function(C, SigmasInv) {
        
        ySum = t(C) %*% y # for each cluster
        n_C = colSums(C)
        
        for(k in c(1:K)) {
            m_part2 = ySum[k, ] %*% SigmasInv[, , k] + mu_prior_mean %*% mu_prior_SigmasInv
            
            m_var = solve(SigmasInv[, , k] * n_C[k] + mu_prior_SigmasInv)
            m_mean = m_var %*% t(m_part2)
            
            mu[k, ] = t(chol(m_var)) %*% rnorm(p) + m_mean
        }
        return(mu)
    }
    
    
    updateSigmasInv <- function(C, mu) {
        
        C_label = colSums((t(C) * c(1:K)))
        
        diff = (y - mu[C_label, ])
        
        for(k in c(1:K)) {
            pick = (C_label == k)
            
            if (sum(pick) > 1) {
                par2 = solve(t(diff[pick, ]) %*% (diff[pick, ]) + SigmasInv_prior_VInv)
            }
            
            if (sum(pick) == 1) {
                vec = matrix(diff[pick, ], 1, p)
                par2 = solve(t(vec) %*% (vec)  + SigmasInv_prior_VInv)
            }
            if (sum(pick) == 0) {
                par2 = solve(SigmasInv_prior_VInv)
            }
            
            par2 = (par2 + t(par2)) / 2
            par1 = sum(C_label == k) + SigmasInv_prior_nu           
            SigmasInv[, , k] = rwishart(par1, par2)
        }        
        return(SigmasInv)
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
        beta_2 = rtbeta(K, par1_beta, par2_beta, 0, eps) / eps #truncated Beta distribution
        
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

    updateK <- function(b, beta, u, mu, SigmasInv) {
        k = length(xi[xi > min(u)])
        if(k == length(xi)) stop("xi is too short")
        if(k <= K) {
            b = b[1:k]
            beta = beta[1:k]
            mu = mu[1:k, ]
            SInew = SigmasInv[, , 1:k]
        }
        else {
            b = c(b, eps + (1 - eps) * rbinom(k - K, 1, p_b))
            beta = c(beta, rbeta(k - K, alpha_beta + d_beta * c((K + 1):k), 1 - d_beta))

            m_part2 = mu_prior_mean %*% mu_prior_SigmasInv            
            m_var = solve(mu_prior_SigmasInv)
            m_mean = m_var %*% t(m_part2)
            for(kk in (K + 1):k) {
                mu = rbind(mu, t(t(chol(m_var)) %*% rnorm(p) + m_mean))
            }

            par2 = solve(SigmasInv_prior_VInv)
            par2 = (par2 + t(par2)) / 2
            par1 = SigmasInv_prior_nu
            
            SInew = array(0, dim = c(p, p, k))
            SInew[, , 1:K] = SigmasInv
            for(kk in (K + 1):k) {
                SInew[, , kk] = rwishart(par1, par2)
            }
        }
        return(list(k = k, b = b, beta = beta, mu = mu, SigmasInv = SInew))
    }
    
    
       
    # initialize the parameters
    K = 20 # will be changed for slice sampler

    mu = matrix(rnorm(K * p), ncol = p)
    I_p = diag(1, p)   
    SigmasInv = array(0, dim = c(p, p, K))    
    for(k in c(1:K)) {
        SigmasInv[, , k] = rwishart(p, I_p)
    }
   
    w = rdirichlet(1, alpha = rep(1, K))
    b = rep(1, K)
    u = rep(0.01, N)
    
    # record trace    
    trace_nC = list()
    starttime = proc.time()
    
    for (i in 1:steps){
        
        C = updateC(u, mu, SigmasInv, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        mu = updateMu(C, SigmasInv)
        SigmasInv = updateSigmasInv(C, mu)

        if(p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)

        u <- updateU(C)
        uK_res = updateK(b, beta, u, mu, SigmasInv)
        K <- uK_res$k
        b <- uK_res$b
        beta = uK_res$beta
        mu = uK_res$mu
        SigmasInv = uK_res$SigmasInv
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