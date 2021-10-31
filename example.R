p = 2 # dimension of the data
N = 250 # sample size

mu_truth = cbind(c(-4, 0, 5), c(1, 2, 3)) # true mean of the three components
n_truth = round(N * c(0.3, 0.3, 0.4)) # true number of samples from the three components

y = matrix(0, N, p)

# generate the data from the three-component normal mixtrue distribution
for(i in 1:N){
    y[i,] = mu_truth[sum(i > cumsum(n_truth)) + 1, ] + rnorm(p)
}

source("qbsb_gibbs_multivariate_gau.R")
res = runQBSBMvGaussian(y, steps = 5E4, burnins = 3E4, K = 40, p_b = 0.9, eps = 1 / nrow(y) ^ 2.1)