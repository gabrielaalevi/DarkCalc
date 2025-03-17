
## Boltzmann Equations

The simplified set of Boltzmann equations with two particles ($Z$ and $\chi$) and with $BR(Z \to \chi + X) = 1$ is given below:

<p float="left">
    <img src="./figures/boltzEqs.png" alt="Boltzmann Equations" width=100% height=100% />
</p>

The more general case can be generally written as:

```math
\begin{gather*}
\frac{d Y_i}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-\langle\sigma v\rangle_{i i} ({Y_i}^2 - {\bar{Y}_i}^2) - \langle\sigma v\rangle_{ij} (Y_i Y_j - \bar{Y}_i \bar{Y}_j) + \langle\sigma v\rangle_{j j \rightarrow i i} \left(\frac{Y_j^2}{\bar{Y}_j^2} -  \frac{Y_i^2}{\bar{Y}_i^2} \right)- \frac{\Gamma_{i\rightarrow j}}{s} \left(Y_i -\bar{Y}_i \frac{Y_j} {\bar{Y}_j} \right) \right.\\
+ \sum_j \frac{K_1 (m_j/T)}{K_2(m_j/T)} \frac{\Gamma_j}{s} \sum_b BR(j \to i+b+...) \left(Y_j - \bar{Y}_j \frac{Y_i}{\bar{Y}_i}\frac{Y_b}{\bar{Y}_b}...\right)\\
\left.  - \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_i}{s}  \sum_{j,b} BR(i \to j +b+...)  \left(Y_i - \bar{Y}_i \frac{Y_j}{\bar{Y}_j}\frac{Y_b}{\bar{Y}_b}...\right) \right]
\end{gather*}
```

Finally, defining the "decay matrix":

```math
D_{ij} \equiv \sum_b \frac{K_1 (m_i/T)}{K_2(m_i/T)} \frac{\Gamma_i}{s}   BR(i \to j +b+...)  \left(Y_i - \bar{Y}_i \frac{Y_j}{\bar{Y}_j}\frac{Y_b}{\bar{Y}_b}...\right)
```

and the "collision terms":

```math
C_{ij,lm} \equiv  \langle \sigma v \rangle_{ij\to lm}  \left( Y_i Y_j - \bar{Y}_i \bar{Y}_j \frac{Y_l}{\bar{Y}_l} \frac{Y_m}{\bar{Y}_m}\right) = -C_{lm,ij}
```

since:
```math
\langle \sigma v \rangle_{lm\to ij}  = \frac{\bar{Y}_i \bar{Y}_j}{\bar{Y}_l \bar{Y}_m} \langle \sigma v \rangle_{ij\to lm}
```

In addition, we have:

```math
\frac{\Gamma_{i\rightarrow j}}{s} \equiv \bar{Y}_{SM} \langle \sigma v \rangle_{i0\to j0}
```

where we have, for convenience, defined:

```math
 Y_0 = Y_{\rm SM} = \bar{Y}_{SM} = \frac{\zeta(3)}{\pi^2}g^{*}_{S}(T) T^3/s = \frac{45}{2 \pi^4} \zeta(3)
```
 with $\frac{d Y_0}{dx} = 0$.

Therefore the Boltzmann equation can be written as:
```math
\begin{gather*}
\frac{d Y_i}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-C_{ii,00}  - C_{ij,00}   - C_{ii,jj} - C_{i0,j0} + \sum_j D_{ji} - \sum_j D_{ij} \right]\\
\Rightarrow \frac{d Y_i}{dx} = \frac{1}{3H} \left|\frac{ds}{dx} \right| \left[-\sum_{a,b,c} C_{ia,bc} + \sum_j D_{ji} - \sum_j D_{ij} \right]
\end{gather*}
```

With the above notation and taking that $a,b,c = i,j$ or $0$, we have the following table of processes which changes $Y_i$ (i.e. number of times $i$ appears in the initial state $\neq$ number of times $i$ appears in the final state).

 * Processes preserving a $Z_2$ symmetry (with $\chi$ and $Z$ being odd):

| Initial Indices | Final Indices | Process                     | Cross-section | Name
| :-------------- | :-----------: | :-------------------------: | :---------: | :---------: |
| $i,0$           | $j,0$         | $\chi + SM \to Z + SM$ | $\Gamma_{\chi \to Z}$ | Conversion
| $i,i$           | $0,0$         | $\chi + \chi \to SM + SM$ | $\langle \sigma v \rangle_{\chi \chi}$ | Self-annihilation
| $i,i$           | $j,j$         | $\chi + \chi \to Z + Z$ | $\langle \sigma v \rangle_{\chi \chi \to Z Z}$ | Scattering
| $i,i$           | $i,j$         | $\chi + \chi \to Z + \chi$ | $\langle \sigma v \rangle_{\chi \chi \to Z \chi}$ | ?
| $i,j$           | $0,0$         | $\chi + Z \to SM + SM$ |  $\langle \sigma v \rangle_{\chi Z}$ | Co-annihilation
| $i,j$           | $j,j$         | $\chi + Z \to Z + Z$ | $\langle \sigma v \rangle_{\chi Z \to Z Z}$ | ?

 * Processes without a $Z_2$ symmetry:

| Initial Indices | Final Indices | Process                     | Cross-section | Name
| :-------------- | :-----------: | :-------------------------: | :---------: | :---------: |
| $i,0$           | $0,0$         | $\chi + SM \to SM + SM$ |  $\langle \sigma v \rangle_{\chi SM \to SM SM}$ | ?  
| $i,0$           | $j,j$         | $\chi + SM \to Z + Z$ |  $\langle \sigma v \rangle_{\chi SM \to Z Z}$ | ?  
| $i,i$           | $i,0$         | $\chi + \chi \to \chi + SM$ | $\langle \sigma v \rangle_{\chi \chi \to \chi SM}$ | ?  
| $i,i$           | $j,0$         | $\chi + \chi \to Z + SM$ |  $\langle \sigma v \rangle_{\chi \chi \to Z SM}$ | ?  
| $i,j$           | $j,0$         | $\chi + Z \to Z + SM$ |  $\langle \sigma v \rangle_{\chi Z \to Z SM}$ | ?  

