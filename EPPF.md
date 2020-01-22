Exchangeable Partition Probability Function

Let $g_k = n_k+m_k$, and $M(z)= \max\{z_1,z_2,\dots,z_n\}$.
$$
\begin{aligned} 
L(z_1, z_2, \dots, z_n)  ={} &  \prod_{k=1}^{\infty} \bigg\{ w \text{B}(n_k+1,  g_{k+1} + 1) \\
& + (1-w)\frac{1}{\epsilon}\text{B}(n_k+1,g_{k+1}+1)\int_{0}^{\epsilon}\frac{1}{\text{B}(n_k+1,g_{k+1}+1)}x^{n_k}(1-x)^{g_{k+1}}\,dx \bigg \}  \\ 
\stackrel{(a)}{=}{} & \prod_{k=1}^{M(z)} \bigg\{ w \text{B}(n_k+1,  g_{k+1} + 1) \\
& + (1-w)\frac{1}{\epsilon}\text{B}(n_k+1,g_{k+1}+1)F_{Beta(n_k+1,g_{k+1}+1)}(\epsilon) \bigg \}\\
={} & \prod_{k=1}^{M(z)} \frac{\Gamma(n_k+1)\Gamma(g_{k+1}+1)}{\Gamma(n_k+g_{k+1}+2)}Q_k\\
={} & \prod_{c\in C_z}\Gamma(|c|+1) \prod_{k=1}^{M(z)}\frac{\Gamma(g_{k+1} + 1)}{(g_k+1)\Gamma(g_k+1)}\prod_{k=1}^{M(z)}Q_k\\
={} & \frac{1}{(n+1)\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \prod_{k=1}^{M(z)}\frac{Q_k}{g_{k+1}+1}
\end{aligned}
$$
where $(a)$ is using $B(1,1)=1$ when $n_k=m_k=0$. For simplicity, let $Q_k = w + (1-w)\frac{1}{\epsilon}F_{Beta(n_k+1,g_{k+1}+1)}(\epsilon)$.

Let $k_1<k_2<\cdots<k_t$ denote the unique values of $(z_1,z_2,\dots,z_n)$, and set $k_0=0$.

For $k_{i-1}<k<k_{i}$, we have $n_k=0$ and $g_{k+1} = g_{k_i}$. This leads to

$$
Q_k = w + (1-w)\frac{1}{\epsilon}\int_{0}^{\epsilon}(g_{k_i}+1)(1-x)^{g_{k_i}}\,dx.
$$

For $k=k_{i}$, we have $n_k=n_{k_i}$ and $g_{k+1} = g_{k_{i+1}}$. This leads to

$$
Q_k = w + (1-w)\frac{1}{\epsilon} F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon).
$$

Combining these two we have,

$$
\begin{aligned}
L  ( z_1,z_2,\dots,z_n) ={} & \frac{1}{(n+1)\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \\
& \times \prod_{i=1}^{t} \biggl[ \biggl\{ \frac{w}{g_{k_i} + 1} + (1 - w) \frac{1}{\epsilon}\int_{0}^{\epsilon} (1-x) ^ {g_{k_i}}\,dx \biggr\}^{d_i}\\
& \times \biggl\{
\frac{w}{g_{k_{i+1}} + 1} + (1-w)\frac{1}{\epsilon(g_{k_{i+1}}+1)}F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)
\biggr \}\biggr],
\end{aligned}
$$
where $d_i = k_i-k_{i-1}-1\in \mathbb{N}$. Summing out all possible values of $d_i$,

$$
\begin{aligned}
& \sum_{d\in\mathbb{N}^t}L(z_1,z_2,\dots,z_n) \\
={}                 & \frac{1}{(n+1)\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \\
& \times \prod_{i=1}^t \bigg [ \bigg\{
\frac{w}{g_{k_{i+1}} + 1} + (1-w)\frac{1}{\epsilon(g_{k_{i+1}}+1)}F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)
\bigg \} \\
& \times \sum_{{d_i=0}}^{\infty} \bigg\{ \frac{w}{g_{k_i} + 1} + (1 - w) \frac{1}{\epsilon}\int_{0}^{\epsilon} (1-x) ^ {g_{k_i}}\,dx \bigg\}^{d_i} \bigg ] \\
\stackrel{(a)}{=}{} & \frac{1}{(n+1)\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \\
& \times \prod_{i=1}^t\bigg [ \bigg\{
\frac{w}{g_{k_{i+1}} + 1} + (1-w)\frac{1}{\epsilon(g_{k_{i+1}}+1)}F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)
\bigg \}  \\
& \times \bigg\{ 1 - \frac{w}{g_{k_i} + 1} - (1 - w) \frac{1}{\epsilon}\int_{0}^{\epsilon} (1-x) ^ {g_{k_i}}\,dx \bigg\}^{-1} \bigg ]\\
={}	& \frac{1}{(n+1)\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \\
& \times \prod_{i=1}^{t} \biggl\{
\frac{g_{k_i}+1}{g_{k_i+1}+1}\cdot
\frac{w + (1-w)\frac{1}{\epsilon}F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)}{g_{k_i} + 1 - w - (1 - w)(g_{k_i} + 1) \frac{1}{\epsilon}\int_{0}^{\epsilon} (1-x) ^ {g_{k_i}}\,dx}
\biggr\}\\
={}	& \frac{1}{\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \prod_{i=1}^{t} 
\frac{w\epsilon + (1-w)F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)}{\epsilon (g_{k_i} + 1 - w) - (1 - w)(g_{k_i} + 1) \int_{0}^{\epsilon} (1-x) ^ {g_{k_i}}\,dx}\\
={}	& \frac{1}{\Gamma(n+1)} \biggl(\prod_{c\in C_z}\Gamma(|c|+1)\biggr) \prod_{i=1}^{t} 
\frac{w\epsilon + (1-w)F_{Beta(n_{k_i}+1,g_{k_{i+1}}+1)}(\epsilon)}{\epsilon g_{k_i} - (1 - w)(1 - \epsilon)[1 - (1-\epsilon)^{g_{k_i}}]}
\end{aligned}
$$
where $(a)$ is based on $\sum_{d=0}^{\infty} x^d = (1-x)^{-1}$ with $|x|<1$.