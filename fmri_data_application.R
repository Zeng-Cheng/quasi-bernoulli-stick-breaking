### FMRI data application ###

y <- as.matrix(read.table("data/cor_data_HCP_PTN1200_recon2_nodets_3T_MSMALL_d50_ts2.txt"))

A = abs(y) > 0.5
N = nrow(A)
v = 50 # dimension of the network matrix
p = ncol(A) # number of elements in the lower triangle matrix

source("src/qbsb_slice_lowrank_probit.R")
steps = 30000
burnins = 10000
kappa = 0.6
xi = (1 - kappa) * kappa ^ (c(0:200))
dir.create("qbsb_hcp")


set.seed(101)
# Quasi-Bernoulli mixture model
res_qb <- runQBSB_slice_efficient_ProbitLowrank(
    A, v, xi, d = 2, steps, burnins, p_b = 0.9, eps = 1 / N ^ 2.1)
save(res_qb, file = "qbsb_hcp/hcp_res_qb.Rda")


set.seed(101)
# Dirichlet process mixture model
res_dp <- runQBSB_slice_efficient_ProbitLowrank(
    A, v, xi, d = 2, steps, burnins, p_b = 1, alpha_beta = 0.63)
save(res_dp, file = "qbsb_hcp/hcp_res_dp.Rda")


set.seed(101)
# Pitman--Yor process mixture model
res_py <- runQBSB_slice_efficient_ProbitLowrank(
    A, v, xi, d = 2, steps, burnins, p_b = 1, d_beta = 0.11, alpha_beta = 0.29)
save(res_py, file = "qbsb_hcp/hcp_res_py.Rda")





###############################
###############################
###############################
### analysis of the results ###

load("qbsb_hcp/hcp_res_qb.Rda")
load("qbsb_hcp/hcp_res_dp.Rda")
load("qbsb_hcp/hcp_res_py.Rda")

v = 50
lc = length(res_qb$nC) # length of chain
N = length(res_qb$label[[1]]) # sample size

########################################################
### print the proportion of subjects in each cluster ###

# print components weights of QB model 

round(rowMeans(sapply(1:lc, function(i) {
    idx = order(res_qb$nC[[i]], decreasing = TRUE)
    res_qb$nC[[i]][idx]
})) / N, 3)

# print components weights of PY model

round(rowMeans(sapply(1:lc, function(i) {
    idx = order(res_py$nC[[i]], decreasing = TRUE)
    res_py$nC[[i]][idx]
})) / N, 3)

# print components weights of DP model

round(rowMeans(sapply(1:lc, function(i) {
    idx = order(res_dp$nC[[i]], decreasing = TRUE)
    res_dp$nC[[i]][idx]
})) / N, 3)


###########################
### plot posterior of T ###

allT = 1:12

### posterior distribution of T
nT <- sapply(1:lc, function(i) sum(res_qb$nC[[i]] != 0)) #number of clusters
probt_QB <- sapply(allT, function(i) sum(nT == i) / lc)

nT <- sapply(1:lc, function(i) sum(res_dp$nC[[i]] != 0)) #number of clusters
probt_DP <- sapply(allT, function(i) sum(nT == i) / lc)

nT <- sapply(1:lc, function(i) sum(res_py$nC[[i]] != 0)) #number of clusters
probt_PY <- sapply(allT, function(i) sum(nT == i) / lc)

postT_fmri <- data.frame(
        Model = factor(rep(c("QB Mixture", "DP Mixture", "PY Mixture"), each = length(allT))),
        Pr = c(probt_QB, probt_DP, probt_PY), t = allT)

ggplot(postT_fmri, mapping = aes(x = t, y = Pr, color = Model)) + theme_bw() +
geom_point(size = 3) + geom_line(size = 1) + facet_grid(col = vars(Model)) + 
theme(legend.position = "none") + theme(strip.text = element_text(size = 12)) +
scale_x_continuous(breaks = allT) + scale_y_continuous(limits = c(0, 1)) + 
ylab("Pr(T=t|y)") + xlab('T: Number of Clusters') +
scale_color_manual(values = c("#F8766D", "#ecf847", "#6E9CF8")) +
theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
theme(axis.title = element_text(size = 12))

ggsave("postT_fmri.png", width = 8, height = 3, unit = "in")




##############################################################
##############################################################
##############################################################
### plot the edge connectivity probablities for each group ###

library(ggplot2)

### load the results
load("qbsb_hcp/hcp_res_qb.Rda")
v = 50
lc = length(res_qb$nC) # length of chain


##############################################################################
### print the proportion of edge probabilities which are greater than 0.05 ###
for(k in 1:6) {
    proportion_non_zero <- mean(sapply(1:lc, function(i) {
        idx = order(res_qb$nC[[i]], decreasing = TRUE)
        corMatrix <- diag(0, v)
        corMatrix[lower.tri(corMatrix)] <- pnorm(res_qb$M[[i]][idx[k], ] + res_qb$mu[[i]])
        AA = corMatrix + t(corMatrix)
        sum(AA > 0.05) / v ^ 2
    }))
    print(proportion_non_zero)
}


##########################################################################
### plot the heat map of edge connectivity probablities for each group ###

library(reshape2)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

for(k in 1:6) {
    postMean <- rowMeans(sapply(1:lc, function(i) {
        idx = order(res_qb$nC[[i]], decreasing = TRUE)
        pnorm(res_qb$M[[i]][idx[k], ] + res_qb$mu[[i]])
    }))
    corMatrix <- diag(0, v)
    corMatrix[lower.tri(corMatrix)] <- postMean
    corMatrix <- corMatrix + t(corMatrix) + diag(0, v)
    corMatrix <- melt(corMatrix)

    ggplot(corMatrix, aes(y = Var1, x = Var2)) + geom_raster(aes(fill = value)) + 
    scale_fill_gradientn(colours = jet.colors(100), breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 16, title = NULL)) +
    labs(x = NULL, y = NULL) + scale_y_reverse(breaks = seq(0, v, length.out = 6)) + 
    scale_x_continuous(breaks = seq(0, v, length.out = 6)) +
    theme_void() + theme(legend.justification = c(0, 0.6), legend.text = element_text(size = 11)) +
    theme(axis.text.x = element_text(size = 11, vjust = 5), axis.text.y = element_text(size = 11, hjust = 1.6))
    
    ggsave(paste("eeg", k, ".png", sep = ""), width = 4.4, height = 3.64, unit = "in")
}





##############################################################
##############################################################
##############################################################
### Mixture of Factor Analyzers to compare with our method ###

library(HDclassif)
mfa <- hddc(A, model = 7)
save(mfa, file = "qbsb_hcp/hcp_res_mfa.Rda")
load("qbsb_hcp/hcp_res_mfa.Rda")

library(reshape2)
library(ggplot2)

################################
### plot the posterior means ###

jet.colors <- colorRampPalette(c(
    "#00007F", "blue", "#007FFF", "#7FFF7F",
    "yellow", "#FF7F00", "red", "#7F0000"
    ))

for(k in 1:3) {
    corMatrix <- diag(0, v)
    corMatrix[lower.tri(corMatrix)] <- mfa$mu[k, ]
    corMatrix <- corMatrix + t(corMatrix) + diag(0, v)
    corMatrix <- melt(corMatrix)

    ggplot(corMatrix, aes(y = Var1, x = Var2)) + 
    geom_raster(aes(fill = value)) + 
    scale_fill_gradientn(
        colours = jet.colors(100),
        breaks = seq(0, 1, by = 0.2),
        limits = c(0, 1)
        ) + 
    guides(fill = guide_colourbar(
        barwidth = 0.5, barheight = 16, title = NULL
        )) +
    labs(x = NULL, y = NULL) +
    scale_y_reverse(breaks = seq(0, v, length.out = 6)) + 
    scale_x_continuous(breaks = seq(0, v, length.out = 6)) +
    theme_void() + theme(
        legend.justification = c(0, 0.6),
        legend.text = element_text(size = 11)
        ) +
    theme(
        axis.text.x = element_text(size = 11, vjust = 5),
        axis.text.y = element_text(size = 11, hjust = 1.6)
        )
    
    ggsave(paste("mfa", k, ".png", sep = ""), width = 4.4, height = 3.64, unit = "in")
}