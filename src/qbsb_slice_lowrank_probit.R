require("LaplacesDemon") # rdirichlet
require("rstiefel") # rbing.matrix.gibbs, rustiefel
require("truncnorm") # rtruncnorm

# This function will run slice Gibbs sampling for the lowrank Probit model 

runQBSB_slice_efficient_ProbitLowrank <- function(
    A, v, xi, d, steps = 30000, burnins = 10000,
    p_b = 0.9, d_beta = 0, alpha_beta = 1, eps = 1e-3,
    mu_prior_mu = 0, mu_prior_sigma2 = 100, tau = 7
    ) {

    # v is the original dimension of the matrix form by each row A[s, ] 
    # tau is the standard deviation of the normal prior on lambda

    stopifnot(eps > 0 & (alpha_beta + d_beta) > 0)
    stopifnot(d_beta < 1 & d_beta >= 0 & p_b >= 0 & p_b <= 1)
    stopifnot(p_b == 1 | d_beta == 0)

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

    updateC <- function(u, M, mu, w) {

        # get K * N matrix of log likelihood of Bernoulli
        log1mat = matrix(0, p, K)
        log2mat = matrix(0, p, K)

        for(k in c(1:K)) {            
            Mmutri <- M[, , k][lower.tri(M[, , k])] + mu[lower.tri(mu)]
            log1mat[, k] <- pnorm(Mmutri, log.p = TRUE) # log(p)
            log2mat[, k] <- pnorm(
                Mmutri, lower.tail = FALSE, log.p = TRUE) # log(1 - p)
        }

        loglik_bernoulli = A %*% log1mat + (1-A) %*% log2mat

        loglik = t(
            t(loglik_bernoulli) +
            log(c(w)) +
            log(outer(xi[1:K], u, ">")) -
            log(xi[1:K])
            )

        gumbel = -log(-log(runif(N * K, min = 0, max = 1)))
        C <- t(apply(loglik + gumbel, 1, function(x) {x == max(x)}))
        return(C)
    }
    

    updateQ <- function(C, Z, mu, Lambda) {

        n_C <- colSums(C)
        original_idx = c(1:d)

        for(k in 1:K) {
            # rbing.matrix.gibbs needs ordered lambda
            sorted_idx = order(Lambda[k, ], decreasing = T)
            rev_idx = match(original_idx, sorted_idx)

            Q[, , k] <- rbing.matrix.gibbs(
                apply(Z[, , C[, k]], c(1, 2), sum) / 2 - mu * n_C[k] / 2,
                diag(Lambda[k, sorted_idx]),
                Q[, sorted_idx, k])
            Q[, , k]<- Q[, rev_idx, k]
        }
        return(Q)
    }


    updateLambda <- function(C, Z, Q, mu) {
        n_C <- colSums(C)

        for(k in 1:K) {
            m_var <- 1 / (n_C[k] / 2 + 1 / tau ^ 2)
            m_mean = m_var * diag(
                t(Q[, , k]) %*% (apply(
                    Z[, , C[, k]], c(1, 2), sum) - mu * n_C[k]) %*% Q[, , k]) / 2
            Lambda[k, ] <- rnorm(d, mean = m_mean, sd = sqrt(m_var))
        }
        return(Lambda)
    }


    updateMu <- function(C, Z, M) {

        # This function use a single number for all elements of mu

        C_label = colSums((t(C) * c(1:K)))
        diff = Z - M[, , C_label]
        sumdiff <- apply(diff, c(1, 2), sum)
        
        mu_var <- 1 / (N * (p + v / 2) + 1 / mu_prior_sigma2)
        mu_mean <- mu_var * (sum(sumdiff[lower.tri(sumdiff)]) + sum(diag(sumdiff) / 2) + mu_prior_mu / mu_prior_sigma2)

        mu = matrix(rnorm(1, mu_mean, sqrt(mu_var)), v, v)

        return(mu)
    }


    updateZ <- function(C, M, mu) {

        C_label = colSums((t(C) * c(1:K)))
        Mmutri <- c(sapply(1:N, function(s) {
            k = C_label[s]
            M[, , k][lower.tri(M[, , k])] + mu[lower.tri(mu)]
        })) # (M + mu)tri for each Z

        Mmudiag <- c(sapply(1:N, function(s) diag(M[, , C_label[s]] + mu)))

        as <- c(t(A))
        upper <- rep(0, N * p)
        lower <- rep(0, N * p)
        upper[as] <- Inf
        lower[!as] <- -Inf

        Ztri <- matrix(rtruncnorm(p * N, mean = Mmutri, a = lower, b = upper), p, N)
        Zdiag <- matrix(rnorm(v * N, mean = Mmudiag, sd = sqrt(2)), v, N)
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

    updateK <- function(b, beta, u, Q, Lambda) {
        k = length(xi[xi > min(u)])
        if(k == length(xi)) stop("xi is too short")
        if(k <= K) {
            b = b[1:k]
            beta = beta[1:k]
            Qnew = Q[, , 1:k]
            Lambda = Lambda[1:k, ]
        }
        else {
            b = c(b, eps + (1 - eps) * rbinom(k - K, 1, p_b))
            beta = c(beta, rbeta(k - K, alpha_beta + d_beta * c((K + 1):k), 1 - d_beta))

            Qnew = array(0, dim = c(v, d, k))
            Qnew[, , 1:K] = Q
            for(kk in (K+1):k) {
                Qnew[, , kk] = rustiefel(v, d)
            }

            Lambda = rbind(
                Lambda,
                matrix(rnorm((k - K) * d, sd = tau), k - K, d)
                )
        }
        return(list(k = k, b = b, beta = beta, Q = Qnew, Lambda = Lambda))
    }

    ### initialize the parameters

    K = 30
    C_label <- kmeans(A, centers = 10, iter.max = 100, nstart = 3)$cluster
    C <- sapply(1:K, function(k) {C_label == k})

    w = c(rdirichlet(1, alpha = rep(1, K)))
    b = rep(1, K)
    u = rep(0.01, N)

    Z = array(0, dim = c(v, v, N))
    Lambda <- matrix(rnorm(K * d, sd = tau), K, d)

    Q <- array(0, dim = c(v, d, K))
    for(k in 1:K)
        Q[, , k] <- rustiefel(v, d)

    M <- array(0, dim = c(v, v, K))
    M <- updateM(Q, Lambda)
    mu = matrix(qnorm(mean(A)), v, v)
        
    # record trace    

    trace_M = list()
    trace_nC = list()
    trace_label = list()
    trace_mu = list()
    starttime = proc.time()

    for (i in 1:steps) {
        
        Z <- updateZ(C, M, mu)

        Q <- updateQ(C, Z, mu, Lambda)
        Lambda <- updateLambda(C, Z, Q, mu)
        M <- updateM(Q, Lambda)
        mu <- updateMu(C, Z, M)

        if(i < 50) {
            n_C <- colSums(C)
            idx = order(n_C, decreasing = T)
            # re-order components
            M = M[, , idx]
            Lambda = Lambda[idx, ]
            Q = Q[, , idx]
            C = C[, idx]
        }

        if (p_b < 1) b <- updateB(C)
        beta <- updateBeta(C, b)
        w <- updateW(b, beta)

        u <- updateU(C)
        uK_res = updateK(b, beta, u, Q, Lambda)
        K <- uK_res$k
        b <- uK_res$b
        beta = uK_res$beta
        Q = uK_res$Q
        Lambda = uK_res$Lambda 
        M <- array(0, dim = c(v, v, K))
        M <- updateM(Q, Lambda)
        w <- updateW(b, beta)

        C <- updateC(u, M, mu, w)
        C_label = colSums((t(C) * c(1:K)))
        n_C <- colSums(C)
        
        if (i > burnins) {
            idx = i - burnins
            trace_M[[idx]] = t(sapply(1:K, function(k) M[, , k][lower.tri(M[, , k])]))
            trace_nC[[idx]] = n_C
            trace_label[[idx]] = C_label
            trace_mu[[idx]] = mu[1, 1]
        }

        if(i %% 10 == 0) {
            print(n_C);
            # print(round(mu[1, 1], 3))
            print(paste("iter", i))
        }
        if(i == 1000) print((proc.time() - starttime)[3] / 1000)
    }
    return(list(M = trace_M, nC = trace_nC, label = trace_label, mu = trace_mu))
}