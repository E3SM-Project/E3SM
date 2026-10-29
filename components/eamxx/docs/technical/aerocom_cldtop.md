# The AeroCom algorithm

The AeroCom recommendation is to compare cloud properties where a satellite
sees them: at cloud top. In EAMxx, the [expression](../user/diags/expressions.md)
`X.aerocom('top')` gives, per column, a field `X` sampled at cloud top, and
`X.aerocom('bot')` at cloud bottom. There are two parts to it: determining
probabilistically where the cloud top is, and then weighting `X` by that
probability.

## Where the cloud top is

Columns are independent. Within a column, take a photon coming down from
above the model top, and let $p_{k}$ be the probability that it passes through
layers $0,\ldots,k$, where layer $0$ is the highest. Cloudiness is completely
opaque, and cloudy layers overlap following the maximum-random assumption.

The highest layer is assumed clear, so $p_{0} = 1$. A layer $k$ is cloudy if its
cloud fraction $C_{k} > \epsilon$, with $\epsilon = 0.001$, and it holds some
condensate, $qc_{k} + qi_{k} > 0$. A layer that is not cloudy lets the photon
through, $p_{k} = p_{k-1}$. Otherwise

$$p_{k} = p_{k-1} \frac{1 - \max(C_{k}, C_{k-1})}{1 - \min(C_{k-1}, 1 - \epsilon)}$$

so that contiguous cloudy layers overlap maximally, and layers apart overlap
randomly. The probability that layer $k$ is the cloud top is the weight

$$w_{k} = p_{k-1} - p_{k},$$

which is $0$ at every layer that is not cloudy. The weights add up to the cloud
fraction seen from above, $1 - p_{K-1}$, where $K$ is the number of layers.

The cloud bottom is the same seen from below: the lowest layer is assumed
clear, and each layer overlaps with the one *below* it, from the second-lowest
up.

## Weighting a field

A field $X$ at cloud top is

$$x = \sum_{k} X_{k} \Phi_{k} w_{k},$$

where $\Phi_{k}$ is $1$, or the share of the condensate in one phase:

$$\Phi^{liq}_{k} = \frac{qc_{k}}{qc_{k} + qi_{k}}, \qquad \Phi^{ice}_{k} = 1 - \Phi^{liq}_{k}.$$

The phase is the `phase` argument, and the cloud-top cloud fraction is $X = 1$,
written with the placeholder receiver `mask`. That makes a three by two matrix:

| `phase` | a field `X` | `X = 1` |
| ------- | ----------- | ------- |
| none | `T_mid.aerocom('top')`: temperature at cloud top | `mask.aerocom('top')`: cloud fraction seen from above |
| `'liq'` | `eff_radius_qc.aerocom('top', phase='liq')`: droplet radius at liquid cloud top | `mask.aerocom('top', phase='liq')`: liquid cloud fraction seen from above |
| `'ice'` | `eff_radius_qi.aerocom('top', phase='ice')`: crystal radius at ice cloud top | `mask.aerocom('top', phase='ice')`: ice cloud fraction seen from above |

Use `phase` for a property of one phase of the cloud (e.g., `nc` or
`eff_radius_qc` with `'liq'`, `ni` or `eff_radius_qi` with `'ice'`), and no
phase for a property of the air (e.g., `T_mid`, `p_mid`). Note that $x$ is a
weighted *sum*, not a mean: it is not divided by the cloud fraction.

$X$ can be any expression. The in-cloud droplet number concentration at
liquid cloud top, per volume of air, is

```yaml
- cdnc_top := (nc*pseudo_density/dz/gravit/cldfrac_tot.where(cldfrac_tot>0)).aerocom('top', phase='liq')
```

where `pseudo_density/dz/gravit` is the air density, which turns `nc` (per
mass of air) into a number per volume, and dividing by `cldfrac_tot` turns a
grid mean into an in-cloud value. The `where(..)` keeps that division away
from layers with no cloud, where `cldfrac_tot` is zero and the division would
be by zero.

Only the layers with a nonzero weight are summed, so `X` need not be finite
anywhere else. If `X` is masked (e.g., it is the result of `where(..)`), its
invalid layers are skipped too.

The weights $\Phi_{k} w_{k}$ are available as an expression of their own,
`cldfrac_tot.aerocom_weights('top', phase='liq')`. `X.aerocom(..)` is
computed from them, so all fields sampled at the same position and with the
same phase share one computation of the weights.

Most studies in this topic refer a technical report by Tiedtke et al. (1979)[@tiedtke_ecmwf_1979]. Another more recent general reference that may be of interest is that of Räisänen et al. (2004)[@raisanen2004stochastic].
