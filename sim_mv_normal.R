### create data from the 3-component mixture of Gaussian ###
set.seed(507)
dir.create("simdata")
for(N in c(250, 1000, 2500)) {
    mu_truth = cbind(c(-4, 0, 5), c(1, 2, 3))
    n_truth = round(N * c(0.3, 0.3, 0.4))
    y = matrix(0, N, 2)

    for(i in 1:N){
        y[i, ] = mu_truth[sum(i > cumsum(n_truth)) + 1, ] + rnorm(2)
    }

    write.table(y, paste("simdata/data_mv_normal_N_", N, ".txt", sep = ""),
        row.names = FALSE, col.names = FALSE)
}

source("src/qbsb_slice_multivariate_gau.R")
kappa = 0.5
xi = (1 - kappa) * kappa ^ (c(0:250))
steps = 5E4
burnins = 2E4
dir.create("simres")

set.seed(15)
# Quasi-Bernoulli mixture model
res = matrix(0, steps - burnins, 3)
i = 1
for(N in c(250, 1000, 2500)) {
    y = as.matrix(read.table(paste("simdata/data_mv_normal_N_", N, ".txt", sep = "")))
    res[, i] = runQBSB_slice_efficient_MvGaussian(
        y, xi, steps, burnins, p_b = 0.9, eps = 1 / nrow(y) ^ 2.1)    
    i = i + 1
}
save(res, file = "simres/mv_normal_qb.Rda")


set.seed(15)
# Dirichlet process mixture model
res = matrix(0, steps - burnins, 3)
alpha = c(0.67, 0.63, 0.61)
i = 1
for(N in c(250, 1000, 2500)) {
    y = as.matrix(read.table(paste("simdata/data_mv_normal_N_", N, ".txt", sep = "")))
    res[, i] = runQBSB_slice_efficient_MvGaussian(
        y, xi, steps, burnins, p_b=1, alpha_beta=alpha[i])
    i = i + 1
}
save(res, file = "simres/mv_normal_dp.Rda")


set.seed(15)
# Pitman--Yor process mixture model
res = matrix(0, steps - burnins, 3)
alpha = c(0.38, 0.29, 0.29)
d = c(0.10, 0.11, 0.10)
i = 1
for(N in c(250, 1000, 2500)) {
    y = as.matrix(read.table(paste("simdata/data_mv_normal_N_", N, ".txt", sep = "")))
    res[, i] = runQBSB_slice_efficient_MvGaussian(
        y, xi, steps, burnins, p_b = 1, alpha_beta = alpha[i], d_beta = d[i])
    i = i + 1
}
save(res, file = "simres/mv_normal_py.Rda")



#######################################
#######################################
### run the following code in Julia ###
#######################################
#######################################

cd("quasi-bernoulli-stick-breaking") # the working directory

using BayesianMixtures
B = BayesianMixtures

n_total = 5E4
n_burn  = 2E4
t_max = 50

using DelimitedFiles
for N in [250, 1000, 2500]
    y_r = readdlm("simdata\\data_mv_normal_N_" * string(N) * ".txt")
    y = [y_r[i,:] for i = 1:N]
    options = B.options("MVN", "MFM", y, n_total, use_hyperprior = false, n_burn = n_burn, t_max = t_max)
    res = B.run_sampler(options)
    writedlm("simres\\mv_normal_mfm_Sim_N_" * string(N) * ".txt", res.t)
end

#########################
### end code in Julia ###
#########################





#############################
#############################
#### plot posterior of T ####
library(ggplot2)
library(Hmisc)
library(gridExtra)

allT = 1:10

### input results of posterior distribution of T for all models

postT <- c()
thinning = 50

load("simres/mv_normal_dp.Rda")
for(i in 1:3) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}

load("simres/mv_normal_py.Rda")
for(i in 1:3) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}

load("simres/mv_normal_qb.Rda")
for(i in 1:3) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}

## load MFM results

for(N in c(250, 1000, 2500)) {
    nT <- unlist(read.table(
        paste("simres/mv_normal_mfm_Sim_N_", N, ".txt", sep = "")
        ))[-(1:20000)]
    M = length(nT) # length of chain
    thinned_idx = seq(50, M, by = 50)
    nT <- nT[thinned_idx]
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}


width = 8
height = 4.5
units = "in"

postT_mv_normal = data.frame(
    Model = factor(rep(c("DP Mixture", "PY Mixture", "QB Mixture", "MFM"),
        each = length(allT) * 3),
        levels = c("DP Mixture", "PY Mixture", "QB Mixture", "MFM")),
    Pr = postT,
    t = allT,
    n = factor(
        rep(c("n=250", "n=1000", "n=2500"),
        each = length(allT)),
        levels = c("n=250", "n=1000", "n=2500")) 
)


ggplot(postT_mv_normal, mapping = aes(x = t, y = Pr, color = Model)) + 
geom_point() + geom_line() + facet_grid(row = vars(Model), col = vars(n)) + theme_bw() +
theme(legend.position = "none") + theme(strip.text = element_text(size = 12)) + 
scale_x_continuous(breaks = allT) + scale_y_continuous(limits = c(0, 1)) + 
ylab("Pr(T=t|y)") + xlab("T: Number of Clusters") +
scale_color_manual(values = c("#F8766D", "#ded709", "#6E9CF8", "#53B64C"))

ggsave("postT_mv_normal.png", width = width, height = height, units = units)