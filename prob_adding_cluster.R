library(ggplot2)
library(gtools)

# all ways of putting n different balls into m urns
balls_urns <- function(n, m) {

    if(n == 1) {
        return(diag(1, m))
    }
    else {
        X <- NULL
        for (i in 1:m) {
            Y = Recall(n - 1, m)
            Y[, i] = Y[, i] + 1
            X <- rbind(X, Y)         
        }
        return(X)
    }
}

# probability of setting a new table if there are m new people
pr_new_table <- function(nt, m, p, alpha, eps) {

    T = length(nt)
    N = sum(nt)

    allputm <- balls_urns(m, T) # all scenarios of putting m balls into T urns
    eppf <- c()

    if(m > 6) print(paste("m =", m, " number of iterations =", dim(allputm)[1]))

    for(k in 1:dim(allputm)[1]) {

        if(k %% 4000 == 0) print(k)
        allpmtt <- permutations(T, T) # for calculating EPPF
        nt1 <- nt + allputm[k, ]
        for(j in 1:nrow(allpmtt)) {
            mk = N + m - cumsum(nt1[allpmtt[j, ]])
            nk = nt1[allpmtt[j, ]]
            nume <- p + (1 - p) * pbeta(eps, mk + alpha, nk + 1) / eps ^ alpha
            denome <- nk + mk + alpha * (1 - p) * (1 - eps ^ (nk + mk))
            allpmtt[j, ] <- nume / denome
        }
        eppf[k] = prod(
            exp(lfactorial(nt1) - lfactorial(nt))) * sum(
                apply(allpmtt, 1, prod))
    }

    allpmtt <- permutations(T, T) # for calculating EPPF
    nt1 <- nt
    for(j in 1:nrow(allpmtt)) {
        mk = N - cumsum(nt1[allpmtt[j, ]])
        nk = nt1[allpmtt[j, ]]
        nume <- p + (1 - p) * pbeta(eps, mk + alpha, nk + 1) / eps ^ alpha
        denome <- nk + mk + alpha * (1 - p) * (1 - eps ^ (nk + mk))
        allpmtt[j, ] <- nume / denome
    }
    eppf_old = sum(apply(allpmtt, 1, prod))

    return(1 - sum(eppf) / eppf_old / exp(
        lgamma(N + m + alpha) - lgamma(N + alpha)))
}


nt = c(50, 50) # existing tables
p.test <- c(1, 0.9, 0.5)
m.test <- c(1:15)

pr = sapply(m.test, function(m) sapply(p.test, function(p) {
    if (p == 1) alpha = 0.69
    else alpha = 1
    pr_new_table(nt, m, p, alpha = alpha, eps = 1 / (sum(nt) + m) ^ 2.1)
    }))


label <- c(
    "Dirichlet Process",
    bquote("Quasi-Bernoulli " * tilde("p") * "=0.9"),
    bquote("Quasi-Bernoulli " * tilde("p") * "=0.5"))

pr_res <- data.frame(
    pr = c(pr),
    Model = factor(label, level = label),
    m = rep(m.test, each = length(p.test))
    )

ggplot(pr_res, mapping=aes(x=m, y=pr, color=Model)) +
geom_line(size = 0.8) + geom_point(size = 1.5) + theme_bw() +
scale_x_continuous(breaks = m.test) + scale_y_continuous(limits = c(0, 0.1)) +
xlab("m") + ylab("Probability of adding new cluster(s)") +
theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)) +
theme(
    legend.text=element_text(),
    legend.title=element_text(size=14),
    panel.grid.minor=element_blank()) +
scale_colour_discrete(labels=label)

ggsave("add_tables_prob.png", width=6, height=3, unit="in")




###########################
### another case of eps ###

nt = c(50, 50) # existing tables
p.test <- c(1, 0.9, 0.5)
m.test <- c(1:15)

pr = sapply(m.test, function(m) sapply(p.test, function(p) {
    if (p == 1) alpha = 0.69
    else alpha = 1
    pr_new_table(nt, m, p, alpha = alpha, eps = 1 / (sum(nt) + m) ^ 5)
    }))


label <- c(
    "Dirichlet Process",
    bquote("Quasi-Bernoulli " * tilde("p") * "=0.9"),
    bquote("Quasi-Bernoulli " * tilde("p") * "=0.5"))

pr_res <- data.frame(
    pr = c(pr),
    Model = factor(label, level = label),
    m = rep(m.test, each = length(p.test))
    )

ggplot(pr_res, mapping=aes(x=m, y=pr, color=Model)) +
geom_line(size = 0.8) + geom_point(size = 1.5) + theme_bw() +
scale_x_continuous(breaks = m.test) + scale_y_continuous(limits = c(0, 0.1)) +
xlab("m") + ylab("Probability of adding new cluster(s)") +
theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)) +
theme(
    legend.text=element_text(),
    legend.title=element_text(size=14),
    panel.grid.minor=element_blank()) +
scale_colour_discrete(labels=label)

ggsave("add_tables_prob2.png", width=6, height=3, unit="in")