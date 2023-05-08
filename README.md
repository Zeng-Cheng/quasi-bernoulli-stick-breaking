# Quasi-Bernoulli Stick-breaking Process Mixture Model

This repo uses R for Bayesian posterior computation for the quasi-Bernoulli stick-breaking (QBSB) process mixture model, which is proposed by the paper *Consistent Model-based Clustering using the Quasi-Bernoulli Stick-Breaking Process* ([arXiv: 2008.09938](https://arxiv.org/abs/2008.09938)). Since the Dirichlet process is a special case of the QBSB process, the functions in this repo can also be used for posterior computation for the Dirichlet process mixture model. In addition, using some parameters, the functions can implement MCMC for the Pitman-Yor process mixture model.

The QBSB mixture model has the following form:

$$
\begin{align}
     Y_i \mid c_i,{\theta} & \sim \mathcal{F}(\theta_{c_i})\mbox{ independently for } 			i=1,\dots,n, \\
     c_1,\dots,c_n \mid {w} & \sim \mbox{Categorical}({w}) \mbox{ i.i.d.},\\
	 \theta_1,\theta_2,\ldots & \sim \mathcal{G} \mbox{ i.i.d.},\\
      w_1 & = v_1, \ w_k = v_k\prod_{l=1}^{k-1} (1-v_l), \mbox{ for } k\geq 2, \\
      v_k & = 1- b_k\beta_k ,\\
      b_k & \sim p\delta_1(\cdot) + (1-p) \delta_\epsilon(\cdot),\\
      \beta_k & \sim  \mbox{Beta}(\alpha, 1).
\end{align}
$$

The $\mathcal{G}$ is the base measure or, in other word, the prior given to the parameter $\theta$ of the component distribution. The following component distributions $\mathcal{F}$ are implemented:

- univariate normal, with base measure normal on the mean and the inverse-Gamma (hyperprior used) on the variance,
- multivariate normal, with base measure normal on the mean and the inverse-Wishart on the variance-covariance matrix,
- Laplace distribution, with base measure normal on the location parameter and the inverse-Gamma on the scale parameter,
- lowrank Probit model for the network data clustering, to read our paper for the details.



## Installation

- Install [R](https://www.r-project.org/).
- Install packages: "LaplacesDemon", "rstiefel", "msm" and "truncnorm".
- Put the R files into the working directory of R. Each file contains a function for posterior computation for the quasi-Bernoulli mixture model under each component distribution.



## Usage

Source the R file for using the corresponding function.

```R
source("src/qbsb_slice_univariate_gau.R")
runQBSB_slice_efficient_UvGaussian
source("src/qbsb_slice_multivariate_gau.R")
runQBSB_slice_efficient_MvGaussian
source("src/qbsb_slice_univariate_lap.R")
runQBSB_slice_efficient_UvLaplace
source("src/qbsb_slice_lowrank_probit.R")
runQBSB_slice_efficient_ProbitLowrank
```

These functions return the Markov chains of the posterior sampling of the parameters (essentially, the number of clusters). When using the first three functions, the data are passed to the argument `y`, the `xi` is the decreasing positive real numbers used for the slice sampler's latent variables. The remaining arguments are parameters of the models (priors) and are self-explained. Please find the details in our paper.



## Steps for replicating the results in our paper

- Follow the Section Installation.
- Install packages: "ggplot2", "reshape2", "HDclassif", "Hmisc", "gridExtra", "coda", "gtools".
- To compare with the MFM model ([Miller and Harrison, 2018](https://www.tandfonline.com/doi/full/10.1080/01621459.2016.1255636?casa_token=tlesuHLlGLcAAAAA%3AO1CempPvL1drsvhYidJuQDFn2QCpL4aJFkhAyS7iCseMCXE7mhF-IeQVOQoCVt23wsKOCXVhobT2PRE)), we produce the results from MFM by using the Julia package [BayesianMixtures](https://github.com/jwmi/BayesianMixtures.jl) which is created by J. W. Miller and M. T. Harrison for nonparametric Bayesian mixture models. Please follow the instruction of this package.
- The corresponding relations between the R files and the produced Figures are as follows:

| R files                    | Figures     |
| -------------------------- | ----------- |
| prob_adding_cluster.R      | 1 and 11    |
| sim_uv_normal.R            | 2           |
| efficiency_mcmc.R          | 3           |
| sim_laplace.R              | 4           |
| fmri_data_application.R    | 5, 6 and 10 |
| sim_mv_normal.R            | 7           |
| efficiency_mcmc_mv.R       | 8           |
| sim_uv_normal_additional.R | 9           |

