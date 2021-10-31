
require("LaplacesDemon") # rdirichlet
require("rstiefel") # rbing.matrix.gibbs
require("truncnorm") # rtruncnorm

runQBSBProbitLowrank <- function(
    A, v, d, steps = 5000, burnins = 1000, K = 30,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    D_prior_mu = 0, D_prior_sigma2 = 0.1, tau = 10
    ) {
    
    # v is the original dimension of the matrix form by each row A[s, ]

    N = nrow(A)
    p = ncol(A) # number of elements of the lower triangle matrix of v * v

    rtbeta <- function(n, alpha, beta, a = 0, b = 1) {
            stopifnot(n > 0 & all(beta > 0) & all(alpha > 0))
            x <- runif(n)
            Fa <- pbeta(a, alpha, beta)
            Fb <- pbeta(b, alpha, beta)
            y <- (1 - x) * Fa + x * Fb
            y[y < 1E-16] = 1E-16
            return(qbeta(y, alpha, beta))
    }


    updateC <- function(M, D, w) {
        
        # loglik_bernoulli = matrix(0, N, K)   
        log1mat = matrix(0, p, K)
        log2mat = matrix(0, p, K)

        for(k in c(1:K)) {
            
            MDtri <- M[, , k][lower.tri(M[, , k])] + D[lower.tri(D)]
            log1mat[, k] <- pnorm(MDtri, log.p = TRUE) # log(p)
            log2mat[, k] <- pnorm(MDtri, lower.tail = FALSE, log.p = TRUE) # log(1 - p)
            
            # loglik_bernoulli[, k] = colSums(t(A) * log1 + (1 - t(A)) * log2)
        }

        loglik_bernoulli = A %*% log1mat + (1-A) %*% log2mat
        loglik = t(t(loglik_bernoulli) + log(c(w)))

        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) { x == max(x) }))
        return(C)
    }

    updateQ <- function(C, Z, D, Lambda) {

        n_C <- colSums(C)

        original_idx = c(1:d)
        for(k in 1:K) {
            sorted_idx = order(Lambda[k, ], decreasing = T)
            rev_idx = match(original_idx, sorted_idx)
            Q[, , k] <- rbing.matrix.gibbs(
                apply(Z[, , C[, k]], c(1, 2), sum) / 2 - D * n_C[k] / 2,
                diag(Lambda[k, sorted_idx]),
                Q[, sorted_idx, k])
            Q[, , k]<- Q[, rev_idx, k]
        }
        return(Q)
    }


    updateLambda <- function(C, Z, Q, D) {
        n_C <- colSums(C)
        for(k in 1:K) {
            m_var <- 1 / (n_C[k] / 2 + 1 / tau ^ 2)
            m_mean = m_var * diag(t(Q[, , k]) %*% (apply(Z[, , C[, k]], c(1, 2), sum) - D * n_C[k]) %*% Q[, , k]) / 2
            Lambda[k, ] <- rtruncnorm(d, mean = m_mean, sd = sqrt(m_var), a = 0)
        }
        return(Lambda)
    }


    updateD <- function(C, Z, M) {

        C_label = colSums((t(C) * c(1:K)))
        diff = Z - M[, , C_label]
        sumdiff <- apply(diff, c(1, 2), sum)
        D = matrix(0, v, v)

        # D_prior_mu_tri <- D_prior_mu[lower.tri(D_prior_mu)]
        # D_prior_sigma2_tri <- D_prior_sigma2[lower.tri(D_prior_sigma2)]

        d_var_tri <- 1 / (N + 1 / D_prior_sigma2)
        d_var_diag <- 1 / (N / 2 + 1 / D_prior_sigma2)

        d_mean_tri <- d_var_tri * (sumdiff[lower.tri(sumdiff)] + D_prior_mu / D_prior_sigma2)
        d_mean_diag <- d_var_diag * (diag(sumdiff) / 2 + D_prior_mu / D_prior_sigma2)

        D[lower.tri(D)] <- rnorm(p, d_mean_tri, sqrt(d_var_tri))
        D <- D + t(D)
        diag(D) <- rnorm(v, d_mean_diag, sqrt(d_var_diag))

        return(D)
    }


    updateZ <- function(C, M, D) {

        C_label = colSums((t(C) * c(1:K)))
        Mtri <- c(sapply(C_label, function(k) M[, , k][lower.tri(M[, , k])] + D[lower.tri(D)])) # (M + D)tri for each Z
        Mdiag <- c(sapply(C_label, function(k) diag(M[, , k] + D)))
        as <- c(t(A))
        upper <- rep(0, N * p)
        lower <- rep(0, N * p)
        upper[as] <- Inf
        lower[!as] <- -Inf

        # sum(as==1) + sum(as==0)
        # sum(is.na(rtruncnorm(100, mean = Mtri[1:100], sd = 1, a = lower[1:100], b = upper[1:100])))


        Ztri <- matrix(rtruncnorm(p * N, mean = Mtri, sd = 1, a = lower, b = upper), p, N)
        Zdiag <- matrix(rnorm(v * N, mean = Mdiag, sd = sqrt(2)), v, N)
        for(s in 1:N) {
            Z[, , s] = 0

            Z[, , s][lower.tri(Z[, , s])] <- Ztri[, s]

            Z[, , s] = Z[, , s] + t(Z[, , s])

            diag(Z[, , s]) <- Zdiag[, s]
        }
        return(Z)
    }


    updateM <- function(Q, Lambda) {
        for(k in 1:K)
            M[, , k] <- Q[, , k] %*% diag(Lambda[k, ]) %*% t(Q[, , k])
        
        return(M)
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

    updateW <- function(b, beta) {
        v = 1 - b * beta
        w = v * (cumprod(c(1, 1 - v))[1:K])
        return(w)
    }
        

    ### initialize the parameters

    C_label <- kmeans(A, centers = K, iter.max = 100, nstart = 3)$cluster
    C <- sapply(1:K, function(k) C_label == k)

    w = c(rdirichlet(1, alpha = rep(1, K)))
    b = rep(1, K)

    Z = array(0, dim = c(v, v, N))
    Lambda <- matrix(rnorm(K * d), K, d)

    Q <- array(0, dim = c(v, d, K))
    for(k in 1:K)
        Q[, , k] <- rustiefel(v, d)

    M <- array(0, dim = c(v, v, K))

    A_smooth = A
    A_smooth = A_smooth * 0.99 + 0.005
    D = matrix(0, v, v)
    D[lower.tri(D)] = qnorm(colMeans(A_smooth))
    D = D + t(D)

    # record trace    

    trace_M = list()
    trace_nC = list()
    trace_label = list()
    trace_D = list()

    for (i in 1:steps) { 

        M <- updateM(Q, Lambda)
        Z <- updateZ(C, M, D)
        
        Q <- updateQ(C, Z, D, Lambda)
        Lambda <- updateLambda(C, Z, Q, D)
        D <- updateD(C, Z, M)

        if(i < 50) {
            n_C <- colSums(C)
            idx = order(n_C, decreasing = T)
            # re-order components to find a good start point for MCMC
            M = M[, , idx]
            Lambda = Lambda[idx, ]
            Q = Q[, , idx]
            C = C[, idx]
        }
         
        if (p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)

        C = updateC(M, D, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)

        if (i > burnins) {
            idx = i - burnins
            trace_M[[idx]] = t(sapply(1:K, function(k) M[, , k][lower.tri(M[, , k])]))
            trace_nC[[idx]] = n_C
            trace_label[[idx]] = C_label
            trace_D[[idx]] = D
        }
        
        if(i %% 5 == 0) {print(n_C); print(i)}       
    }
    
    return(list(M = trace_M, nC = trace_nC, label = trace_label, D = trace_D))
}