# The AeroCom algorithm

The goal of the AeroCom algorithm is to calculate properties at cloud top based on the AeroCom recommendation. There are two main parts of the algorithm: probabilistically determining "cloud top" and then "calculating properties" at said cloud top.

We treat model columns independently, so we loop over all columns in parallel. We then loop over all layers in serial (due to needing an accumulative product), starting at 2 (second highest) layer because the highest is assumed to have no clouds. Let's take a photonic approach from above the model top. Let's say that $p_{k}$ is the probability of a photon passing through the layer $k$. We follow the maximum-random overlap assumption. In all cases, we assume the cloudiness (or cloudy fraction) is completely opaque.

We assume the highest layer has no clouds, thus the $p_{k} = 1$ for the highest layer. Note that $p_{k}$ is initialized as 1 for all layers. We also clip the cloudy fraction $C_{i,k}$ to ensure that $C_{i,k} \in [0+\epsilon, 1-\epsilon]$, where $\epsilon = 0.001$. Starting at the second highest layer, $k+1$, we check if some "cloudy" conditions are met. These conditions are now arbitrarily defined by a cloudiness threshold of $\epsilon$ (i.e., $C_{i,k}>\epsilon$) and a non-zero threshold on the total (both liquid and ice) droplet number concentration (i.e., $cQ_{i,k} + iQ_{i,k} > 0$). If the conditions are met, we estimate the cloud-top cloud fraction using an accumulative product following the maximum-random overlap assumption.

$$c_{i} = 1 - \prod_{k=2}^{K} p_{k} = 1 - \prod_{k=2}^{K} \frac{1 - \max(C_{i,k}, C_{i,k-1})}{1-C_{i,k-1}}$$

In order to estimate cloud-top properties, we weight by the probability of "remaining cloudiness" or $p_{k-1} - p_{k}$.

| Type | Equation |
| --- | --------- |
| cloud property | $x_{i} = \sum_{k=2}^{K} X_{i,k} \Phi_{i,k} (p_{k-1} - p_{k})$ |
| cloud content | $x_{i} = \sum_{k=2}^{K} \Phi_{i,k} (p_{k-1} - p_{k})$ |
| other property | $x_{i} = \sum_{k=2}^{K} X_{i,k} (p_{k-1} - p_{k})$ |

In the above, $\Phi_{i,k}$ is the thermodynamic phase defined by the cloud droplet number concentration ratios.

$$i\Phi_{i,k} = \frac{iQ_{i,k}}{iQ_{i,k} + cQ_{i,k}}$$

$$c\Phi_{i,k} = \frac{cQ_{i,k}}{iQ_{i,k} + cQ_{i,k}}$$

The thermodynamic phase is used only for cloud properties (e.g., cloud-top cloud droplet number concentration) or cloud content (e.g., cloud liquid content). Further, $X_{i,k}$ is the three-dimensional cloud property of interest which is needed if we are converting a property from three-dimensional ($X$) to its two-dimensional counterpart ($x$). "Other" properties here include temperature and pressure which are not dependent on the thermodynamic phase.

Most studies in this topic refer a technical report by Tiedtke et al. (1979)[@tiedtke_ecmwf_1979]. Another more recent general reference that may be of interest is that of Räisänen et al. (2004)[@raisanen2004stochastic].
