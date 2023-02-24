### create 5 batches data with sample size N ###
set.seed(10)
dir.create("simdata")
N = 1000
for (k in 1:5) {
    mu_truth = c(-4, 0, 5)
    n_truth = round(N * c(0.3, 0.3, 0.4))
    y = rep(0, N)

    for (i in 1:N){
        y[i] = mu_truth[sum(i > cumsum(n_truth)) + 1] + rnorm(1)
    }

    write.table(y, paste("simdata/data_uv_normal_1000_batch_", k, ".txt", sep = ""),
        row.names = FALSE, col.names = FALSE)
}


source("qbsb_slice_univariate_gau.R")
kappa = 0.5
xi = (1 - kappa) * kappa ^ (c(0:200))
steps = 25000
burnins = 5000
dir.create("simres")


set.seed(10)
# Quasi-Bernoulli mixture model
res = matrix(0, steps - burnins, 5)
for(i in 1:5) {
    y = c(unlist(read.table(paste("simdata/data_uv_normal_1000_batch_", i, ".txt", sep = ""))))
    res[, i] = runQBSB_slice_efficient_UvGaussian(
        y, xi, steps, burnins, p_b = 0.9, eps = 1 / length(y) ^ 2.1)
}
save(res, file = "simres/uv_normal_qb_1000_5batch.Rda")


#######################################
#######################################
### run the following code in Julia ###
#######################################
#######################################

cd("quasi-bernoulli-stick-breaking") # the working directory

using BayesianMixtures
B = BayesianMixtures

n_total = 25000
n_burn  = 5000
t_max = 50

using DelimitedFiles
for i in [1, 2, 3, 4, 5]
    y = vec(readdlm("simdata\\data_uv_normal_1000_batch_" * string(i) * ".txt"))
    options = B.options("Normal", "MFM", y, n_total, n_burn = n_burn, t_max = t_max)
    res = B.run_sampler(options)
    writedlm("simres\\uv_normal_mfm_1000_batch_" * string(i) * ".txt", res.t)
end


#########################
### end code in Julia ###
#########################




#############################
#############################
#### analysis of results ####

library(ggplot2)
library(coda)

#### input results
load("simres/uv_normal_qb_1000_5batch.Rda") # load QB results

### load MFM results
res_mfm <- c()
for(i in 1:5) {
     res_mfm <- cbind(res_mfm, unlist(read.table(paste(
        "simres/uv_normal_mfm_1000_batch_", i, ".txt", sep = "")))[-(1:5000)])    
}


#############################
### calculate average ESS ###

thinned_idx = seq(50, nrow(res), by = 50)
mean(effectiveSize(res[thinned_idx, ]) / nrow(res[thinned_idx, ]))
mean(effectiveSize(res_mfm[thinned_idx, ]) / nrow(res_mfm[thinned_idx, ]))


##############################
### plot traces of a chain ###

width = 5
height = 3
units = "in"

trace_uv_normal = data.frame(
    t = c(res_mfm[, 3], res[, 3]),
    Iteration = c(1:nrow(res_mfm)),
    Model = factor(rep(c("MFM", "QB Mixture"), each = nrow(res_mfm)))
)

ggplot(data = trace_uv_normal, aes(x = Iteration, y = t)) +
geom_line(aes(color = Model)) + facet_grid(rows = vars(Model)) + theme_bw() +
scale_color_manual(values = c("#53B64C", "#6E9CF8")) +
theme(legend.position = "none") + theme(strip.text = element_text(size = 12)) +
ylab("Number of Clusters")

ggsave("trace_uv_normal.png", width = width, height = height, units = units)


############################
### plot acfs of a chain ###

lag.max = 40
thinned_idx = seq(50, nrow(res), by = 50)
acf_uv_normal = data.frame(
    ACF = c(        
        as.numeric(acf(res_mfm[thinned_idx, 3], lag.max = lag.max, plot = FALSE)[[1]]),
        as.numeric(acf(res[thinned_idx, 3], lag.max = lag.max, plot = FALSE)[[1]])
    ),
    Lag = c(0:lag.max),
    Model = factor(rep(c("MFM", "QB Mixture"), each = lag.max + 1))
)

ggplot(data = acf_uv_normal, aes(x = Lag, y = ACF)) +
geom_bar(stat = "identity", aes(fill = Model), width = 0.5) +
facet_grid(rows = vars(Model)) + theme_bw() +
theme(legend.position = "none") + theme(strip.text = element_text(size = 12)) +
scale_fill_manual(values=c("#53B64C", "#6E9CF8"))

ggsave("acf_uv_normal.png", width = width, height = height, units = units)