
## Boltzmann Equations

The simplified set of Boltzmann equations with two particles ($Z$ and $\chi$) and with $BR(Z \to \chi + X) = 1$ is given below:

<p float="left">
    <img src="./figures/boltzEqs.png" alt="Boltzmann Equations" width=100% height=100% />
</p>

The more general case can be generally written as:

$$\frac{d Y_{i}}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-\langle\sigma v\rangle_{i i} (Y_{i}^2 - \overline{Y}_{i}^2) - \langle\sigma v\rangle_{ij} (Y_{i} Y_j - \overline{Y}_{i} \overline{Y}_j) + \langle\sigma v\rangle_{j j \rightarrow i i} \left(\frac{Y_{j}^2}{\overline{Y}_{j}^2} -  \frac{Y_{i}^2}{\overline{Y}_{i}^2} \right)\right.\\
+ \frac{\Gamma_{j\rightarrow i}}{s} \left(Y_j -\overline{Y}_j \frac{Y_{i}} {\overline{Y}_{i}} \right) - \frac{\Gamma_{i\rightarrow j}}{s} \left(Y_i -\overline{Y}_i \frac{Y_{j}} {\overline{Y}_{j}} \right) \\
+ \sum_{j} \frac{K_1 (m_j/T)}{K_2(m_j/T)} \frac{\Gamma_{j}}{s} \sum_{b} BR(j \to i+b+...) \left(Y_j - \overline{Y}_{j} \frac{Y_{i}}{\overline{Y}_{i}}\frac{Y_{b}}{\overline{Y}_{b}}...\right)\\
\left.  - \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_{i}}{s}  \sum_{j,b} BR(i \to j +b+...)  \left(Y_i - \overline{Y}_{i} \frac{Y_{j}}{\overline{Y}_{j}}\frac{Y_{b}}{\overline{Y}_{b}}...\right) \right]$$

Finally, defining the "decay matrix":

$$D_{ij} \equiv \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_{i}}{s}  \sum_{b} BR(i \to j +b+...)  \left(Y_i - \overline{Y}_{i} \frac{Y_{j}}{\overline{Y}_{j}}\frac{Y_{b}}{\overline{Y}_{b}}...\right)$$

we have:
$$\frac{d Y_{i}}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-\langle\sigma v\rangle_{i i} (Y_{i}^2 - \overline{Y}_{i}^2) - \langle\sigma v\rangle_{ij} (Y_{i} Y_j - \overline{Y}_{i} \overline{Y}_j) + \langle\sigma v\rangle_{j j \rightarrow i i} \left(\frac{Y_{j}^2}{\overline{Y}_{j}^2} -  \frac{Y_{i}^2}{\overline{Y}_{i}^2} \right)\right.\\
+ \frac{\Gamma_{j\rightarrow i}}{s} \left(Y_j -\overline{Y}_j \frac{Y_{i}} {\overline{Y}_{i}} \right) - \frac{\Gamma_{i\rightarrow j}}{s} \left(Y_i -\overline{Y}_i \frac{Y_{j}} {\overline{Y}_{j}} \right) 
\\
\left. + \sum_{j} D_{ji} - \sum_{j} D_{ij} \right]$$
