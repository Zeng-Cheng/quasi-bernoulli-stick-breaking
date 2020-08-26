p = 2
N = 200

mu_truth = cbind(c(-4, 0, 5), c(1, 2, 3))
n_truth = round(N * c(0.3, 0.3, 0.4))

y = matrix(0, N, p)

for(i in 1:N){
    y[i,] = mu_truth[sum(i > cumsum(n_truth)) + 1, ] + rnorm(p)
}

source("qbsb_gibbs_multivariate_gau.R")

res = runQBSBMvGaussian(y, steps = 2E4, burnins = 1E4, K = 20, p_b = 0.9, eps = 1 / nrow(y) ^ 1.1)