
## Boltzmann Equations

The simplified set of Boltzmann equations with two particles ($Z$ and $\chi$) and with $BR(Z \to \chi + X) = 1$ is given below:

<p float="left">
    <img src="./figures/boltzEqs.png" alt="Boltzmann Equations" width=100% height=100% />
</p>

The more general case can be generally written as:

$$\frac{d Y_i}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-\langle\sigma v\rangle_{i i} ({Y_i}^2 - {\bar{Y}_i}^2) - \langle\sigma v\rangle_{ij} (Y_i Y_j - \bar{Y}_i \bar{Y}_j) + \langle\sigma v\rangle_{j j \rightarrow i i} \left(\frac{Y_j^2}{\bar{Y}_j^2} -  \frac{Y_i^2}{\bar{Y}_i^2} \right)- \frac{\Gamma_{i\rightarrow j}}{s} \left(Y_i -\bar{Y}_i \frac{Y_j} {\bar{Y}_j} \right) \right.\\
+ \sum_j \frac{K_1 (m_j/T)}{K_2(m_j/T)} \frac{\Gamma_j}{s} \sum_b BR(j \to i+b+...) \left(Y_j - \bar{Y}_j \frac{Y_i}{\bar{Y}_i}\frac{Y_b}{\bar{Y}_b}...\right)\\
\left.  - \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_i}{s}  \sum_{j,b} BR(i \to j +b+...)  \left(Y_i - \bar{Y}_i \frac{Y_j}{\bar{Y}_j}\frac{Y_b}{\bar{Y}_b}...\right) \right]$$

Finally, defining the "decay matrix":

$$D_{ij} \equiv \sum_b \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_i}{s}   BR(i \to j +b+...)  \left(Y_i - \bar{Y}_i \frac{Y_j}{\bar{Y}_j}\frac{Y_b}{\bar{Y}_b}...\right)$$

and the "collision terms":

$$C_{ij,lm} \equiv \langle \sigma v \rangle_{ij\to lm} \left( \frac{Y_i}{\bar{Y}_i} \frac{Y_j}{\bar{Y}_j} - \frac{Y_l}{\bar{Y}_l} \frac{Y_m}{\bar{Y}_m}\right) = -C_{lm,ij}$$
and
$$
\langle \sigma v \rangle_{iX\to jX} \equiv \bar{Y}_{i} \frac{\Gamma_{i\rightarrow j}}{s}
$$

we have:
$$\frac{d Y_i}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-C_{ii,XX} \frac{1}{\bar{Y}_i \bar{Y}_i} - C_{ij,XX} \frac{1}{\bar{Y}_i \bar{Y}_j}   - C_{ii,jj} 
 - C_{iX,jX} + \sum_j D_{ji} - \sum_j D_{ij} \right]$$

and for convenience we define $Y_0 = Y_{\rm SM} = Y_{eq} = \frac{\zeta(3)}{\pi^2}g^{*}_{S}(T) T^3/s = \frac{45}{2 \pi^4} \zeta(3) $ with $\frac{d Y_0}{dx} = 0$.