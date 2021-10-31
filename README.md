# Quasi-Bernoulli Stick-breaking Process Mixture Model

This repo uses R for Bayesian posterior computation for the quasi-Bernoulli stick-breaking (QBSB) process mixture model, which is proposed by Zeng and Duan ([arXiv: 2008.09938](https://arxiv.org/abs/2008.09938)). Since the Dirichlet process is a special case of the QBSB process, the functions in this repo can also be used for posterior computation for the Dirichlet process mixture model. In addition, using some parameters, the functions can implement MCMC for the Pitman-Yor process mixture model.

The QBSB mixture model has the following form:
$$
\begin{aligned}
    & Y_i \mid c_i,{\theta} \sim \mathcal{F}(\theta_{c_i})\text{ independently for } 			i=1,\dots,n,\\
    & c_1,\dots,c_n \mid {w} \sim \text{Categorical}({w}) \text{ i.i.d.},\\
	& \theta_1,\theta_2,\ldots \sim \mathcal{G} \text{ i.i.d.},\\
    &  w_1 = v_1, \ w_k = v_k\prod_{l=1}^{k-1} (1-v_l), \text{ for } k\geq 2, \\
    &  v_k = 1- b_k\beta_k ,\\
    &  b_k \sim p\delta_1(\cdot) + (1-p) \delta_\epsilon(\cdot),\\
    &  \beta_k \sim  \text{Beta}(\alpha, 1).\\
        
\end{aligned}
$$
The $\mathcal{G}$ is the base measure or, in other word, the prior given to the parameter $\theta$ of the component distribution. The following component distributions $\mathcal{F}$ are implemented:

- univariate normal, with base measure normal on the mean and the inverse-Gamma (hyperprior used) on the variance,
- multivariate normal, with base measure normal on the mean and the inverse-Wishart on the variance-covariance matrix,
- Laplace distribution, with base measure normal on the location parameter and inverse-Gamma on the scale parameter,
- probit model for the network data clustering, to read the paper (Zeng and Duan, arXiv: 2008.09938) for the details.



## Installation

- Install [R](https://www.r-project.org/).

- Install packages: "LaplacesDemon", "rstiefel" and "truncnorm".

- Put the R files into the working directory of R. Each file contains a function for posterior computation for the quasi-Bernoulli mixture model under each component distribution.

  

## Usage

Source the R file for using the corresponding function.

```R
source("qbsb_gibbs_univariate_gau.R")
runQBSBUvGaussian
source("qbsb_gibbs_multivariate_gau.R")
runQBSBMvGaussian
source("qbsb_gibbs_univariate_lap.R")
runQBSBUvLaplace
source("qbsb_gibbs_lowrank_probit.R")
runQBSBProbitLowrank
```

These functions return the Markov chains of the posterior sampling of the parameters.
