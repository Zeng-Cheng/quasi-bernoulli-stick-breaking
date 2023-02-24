source("src/qbsb_slice_univariate_gau.R")
sample_sizes = c(100, 250, 1000, 2500)

### create data from the 3-component mixture of Gaussian ###
set.seed(1000)
dir.create("simdata")
for (N in sample_sizes) {
    mu_truth = c(-2, 0, 2)
    sigma_truth = c(1, 1, 1)
    n_truth = round(N * c(0.3, 0.4, 0.3))
    y = rep(0, N)

    for (i in 1:N){
        y[i] = (
            mu_truth[sum(i > cumsum(n_truth)) + 1] +
            rnorm(1) * sigma_truth[sum(i > cumsum(n_truth)) + 1]
        )
    }

    write.table(y, paste("simdata/data_uv_normal_N_", N, ".txt", sep = ""),
        row.names = FALSE, col.names = FALSE)
}


kappa = 0.3
xi = (1 - kappa) * kappa ^ (c(0:200))
steps = 5E4
burnins = 2E4
dir.create("simres")


set.seed(100)
# Quasi-Bernoulli mixture model case 1
res = matrix(0, steps - burnins, 4)
i = 1
for(N in sample_sizes) {
    y = c(unlist(read.table(paste(
        "simdata/data_uv_normal_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b = 0.9, eps = 1 / N ^ 2.1)
    i = i + 1
}
save(res, file = "simres/uv_normal_qb_power2.Rda")


set.seed(100)
# Quasi-Bernoulli mixture model case 2
res = matrix(0, steps - burnins, 4)
i = 1
for(N in sample_sizes) {
    y = c(unlist(read.table(paste(
        "simdata/data_uv_normal_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b = 0.9, eps = 1 / N ^ 3.1)
    i = i + 1
}
save(res, file = "simres/uv_normal_qb_power3.Rda")


set.seed(1000)
# Dirichlet process mixture model log rate
res <- matrix(0, steps - burnins, 4)
alpha <- 4 / log(sample_sizes)
i <- 1
for(N in sample_sizes) {
    y <- c(unlist(read.table(
        paste("simdata/data_uv_normal_N_", N, ".txt", sep = "")
        )))
    res[, i] <- runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b = 1, alpha_beta = alpha[i])
    i <- i + 1
}
save(res, file = "simres/uv_normal_dp_log.Rda")


set.seed(1000)
# Dirichlet process mixture model reciprocal rate
res = matrix(0, steps - burnins, 4)
alpha = 20 / sample_sizes
i = 1
for(N in sample_sizes) {
    y = c(unlist(read.table(
        paste("simdata/data_uv_normal_N_", N, ".txt", sep = "")
        )))
    res[, i] = runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b=1, alpha_beta=alpha[i])
    i = i + 1
}
save(res, file = "simres/uv_normal_dp_power1.Rda")


set.seed(1000)
# Dirichlet process mixture model exponential rate
res = matrix(0, steps - burnins, 4)
alpha = exp(-sample_sizes / 10)
i = 1
for(N in sample_sizes) {
    y = c(unlist(read.table(paste("simdata/data_uv_normal_N_", N, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b=1, alpha_beta=alpha[i])
    i = i + 1
}
save(res, file = "simres/uv_normal_dp_exp.Rda")


#############################
#############################
#### plot posterior of T ####

library(ggplot2)
library(Hmisc)
library(gridExtra)

### plot posterior distribution of T

allT <- 1:10

### input results of posterior distribution of T for all models
postT <- c()
thinning = 50

calc_post_T <- function(res) {
    M = dim(res)[1] # length of chain
    thinned_idx = seq(thinning, M, by = thinning)
    probt <- sapply(data.frame(res[, 1:4]), function(x) {
        nT <- x[thinned_idx] # number of clusters
        sapply(allT, function(T) sum(nT == T)) / length(thinned_idx)
    })
    return(probt)
}

load("simres/uv_normal_qb_power2.Rda")
postT <- c(postT, calc_post_T(res))

load("simres/uv_normal_qb_power3.Rda")
postT <- c(postT, calc_post_T(res))

load("simres/uv_normal_dp_exp.Rda")
postT <- c(postT, calc_post_T(res))

load("simres/uv_normal_dp_log.Rda")
postT <- c(postT, calc_post_T(res))

load("simres/uv_normal_dp_power1.Rda")
postT <- c(postT, calc_post_T(res))

width = 11
height = 7
units = "in"

label <- c(
    bquote("QBM with "*epsilon[1]*"(n)"),
    bquote("QBM with "*epsilon[2]*"(n)"),
    bquote("DPM with "*alpha[1]*"(n)"),
    bquote("DPM with "*alpha[2]*"(n)"),
    bquote("DPM with "*alpha[3]*"(n)")
    )

levels <- label[c(3, 4, 5, 1, 2)]

sslabels <- c(
    bquote("n"*"="*"100"),
    bquote("n"*"="*"250"),
    bquote("n"*"="*"1000"),
    bquote("n"*"="*"2500")
)

postT_uv_normal = data.frame(
    Model = factor(
        rep(label, each = length(allT) * 4),
        levels = levels
        ),
    Pr = postT,
    t = allT,
    n = factor(
        rep(sslabels, each = length(allT)),
        levels = sslabels
        )
)


ggplot(postT_uv_normal, mapping = aes(x = t, y = Pr, color = Model)) +
geom_point() + geom_line() + 
facet_grid(row = vars(Model), col = vars(n), labeller = label_parsed) +
theme_bw() +
theme(legend.position = "none") + theme(strip.text = element_text(size = 11)) +
scale_x_continuous(breaks = allT) + scale_y_continuous(limits = c(0, 1)) +
ylab("Pr(T=t|y)") + xlab("T: Number of Clusters") +
scale_color_manual(values = c("#a31dd4", "#d51c0f", "#cb26a8", "#6E9CF8", "#0fdcf7"))

ggsave("postT_uv_normal_ano.png", width = width, height = height, units = units)
