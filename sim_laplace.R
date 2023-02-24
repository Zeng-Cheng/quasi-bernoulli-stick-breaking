### create data from the three-component mixture of Laplace distribution ###
set.seed(101)
dir.create("simdata")
for(N in c(50, 200, 500, 1000)) {
    mu_truth = c(-10, 0, 10)
    lambda_truth = c(1, 1.5, 0.5)
    n_truth = round(N * c(0.35, 0.3, 0.35))

    y = rep(0, N)

    for(i in 1:N){
        ci <- sum(i > cumsum(n_truth)) + 1
        y[i] = mu_truth[ci] + rexp(1, 1 / lambda_truth[ci]) - rexp(1, 1 / lambda_truth[ci])
    }

    write.table(y, paste("simdata/data_laplace_N_", N, ".txt", sep = ""),
        row.names = FALSE, col.names = FALSE)
}



source("src/qbsb_slice_univariate_lap.R")
steps = 5E4
burnins = 2E4
kappa = 0.4
xi = (1 - kappa) * kappa ^ (c(0:200))
dir.create("simres")


set.seed(101)
# Quasi-Bernoulli mixture model
res = matrix(0, steps - burnins, 4)
i = 1
for(N in c(50, 200, 500, 1000)) {
    y = c(unlist(read.table(paste("simdata/data_laplace_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvLaplace(
        y, xi, steps, burnins, p_b = 0.9, eps = 1 / length(y) ^ 2.1)
    i = i + 1
}
save(res, file = "simres/laplace_qb.Rda")


set.seed(101)
# Dirichlet process mixture model
res = matrix(0, steps - burnins, 4)
alpha = c(0.71, 0.68, 0.65, 0.63)
i = 1
for(N in c(50, 200, 500, 1000)) {
    y = c(unlist(read.table(paste("simdata/data_laplace_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvLaplace(
        y, xi, steps, burnins, p_b = 1, alpha_beta = alpha[i])
    i = i + 1
}
save(res, file = "simres/laplace_dp.Rda")


set.seed(101)
# Pitman--Yor process mixture model
res = matrix(0, steps - burnins, 4)
alpha = c(0.48, 0.39, 0.35, 0.29)
d = c(0.09, 0.10, 0.10, 0.11)
i = 1
for(N in c(50, 200, 500, 1000)) {
    y = c(unlist(read.table(paste("simdata/data_laplace_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvLaplace(
        y, xim steps, burnins, p_b = 1, alpha_beta = alpha[i], d_beta = d[i])
    i = i + 1
}
save(res, file = "simres/laplace_py.Rda")



#############################
#############################
#### plot posterior of T ####

library(ggplot2)
library(Hmisc)
library(gridExtra)

allT = 1:15

### input results of posterior distribution of T for all models

postT <- c()
thinning = 50

load("simres/laplace_qb.Rda")
for(i in 1:4) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}

load("simres/laplace_dp.Rda")
for(i in 1:4) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}


load("simres/laplace_py.Rda")
for(i in 1:4) {    
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    nT <- res[thinned_idx, i] # number of clusters
    probt <- sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    postT <- c(postT, probt)
}

postT_laplace = data.frame(
    Model = factor(rep(
        c("QB Mixture", "DP Mixture", "PY Mixture"),
        each = length(allT) * 4
        )),
    Pr = postT,
    t = allT,
    n = factor(
        rep(c("n=50", "n=200", "n=500", "n=1000"),
        each = length(allT)),
        levels = c("n=50", "n=200", "n=500", "n=1000")) 
)


ggplot(postT_laplace, mapping = aes(x = t, y = Pr, color = Model)) + 
geom_point() + geom_line() + facet_grid(row = vars(Model), col = vars(n)) + theme_bw() +
theme(legend.position = "none") + 
theme(strip.text = element_text(size = 10), axis.text.x = element_text(size = 6)) + 
scale_x_continuous(breaks = allT) + scale_y_continuous(limits = c(0, 1)) + 
ylab("Pr(T=t|y)") + xlab("T: Number of Clusters") +
scale_color_manual(values = c("#F8766D", "#ded709", "#6E9CF8"))

ggsave("postT_laplace.png", width = 8, height = 4, units = "in")