# Field contraction diagnostics

In EAMxx, we can automatically calculate field reductions
across the horizontal columns and across the model vertical levels.
We call these horizontal and vertical reductions.
We can also automatically calculate zonal averages and histograms.
All of these are written as [expressions](expressions.md), and compose:
e.g., `X.weighted(pseudo_density).mean('lev').weighted(area).mean('col')`. If
using weighted vertical reductions, do the vertical reduction first.

## Horizontal reduction

The horizontal reduction $C$ is the weighted average across all columns,

$$
C_{\dots k} = \frac{\sum_{i} w_{i} \cdot F_{i \dots k}}{\sum_{i} w_{i}}
$$

where $F_\text{i...k}$ is the field at column $i$ and level $k$,
and $w_{i}$ is the weight at column $i$.
We note that the field $F$ can have other dimensions ($\dots$).
As in xarray, `X.mean('col')` weights every column by 1, and
`X.weighted(W).mean('col')` by `W`, where `W` is `area`, the grid's cell area,
or any field with layout `(ncol)`:

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.weighted(area).mean('col')` | Cell area | Average over the domain |
| `X.mean('col')` | 1 | Average across all columns, each counted once |
| `X.weighted(W).mean('col')` | $W_{i}$ | Average across all columns, weighted by any field `W` |

## Vertical reduction

Vertical reductions $C$ are defined as

$$
C_{\dots} = \sum_{k} w_{k} \cdot F_{\dots k}
$$

where $F_{\dots k}$ is the field at level $k$,
and $w_{k}$ is the weight at level $k$.

(the mean divides by $\sum_k w_k$). To select one, request `X.sum('lev')` or
`X.mean('lev')`, which weight every level by 1, or, as in xarray, give the
weights first with `X.weighted(W)`, where `W` is any field or expression on the
same levels:

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.mean('lev')` | 1 | Average across all levels |
| `X.sum('lev')` | 1 | Sum across all levels |
| `X.weighted(pseudo_density).mean('lev')` | $\Delta p_{k}$ | Average across all levels, weighted by $\Delta p_{k}$ |
| `X.weighted(pseudo_density/gravit).sum('lev')` | $\Delta p_{k}/g$ | Mass-weighted vertical integral |
| `X.weighted(dz).mean('lev')` | $\Delta z_{k}$ | Average across all levels, weighted by $\Delta z_{k}$ |
| `X.weighted(dz).sum('lev')` | $\Delta z_{k}$ | Vertical integral over height |
| `X.weighted(W).sum('lev')` | $W_{k}$ | Sum across all levels, weighted by any field `W` |

`pseudo_density` is $\Delta p_{k}$ in Pa, `gravit` is g in m/s$^2$, and `dz` is
$\Delta z_{k}$ in m. `W` has layout `(lev)` or that of `X` without its
components, and must not be masked. If `X` is masked, invalid levels are
skipped; see [expressions](expressions.md) for what a column with no valid
level gives.

## Zonal reduction

We currently have a utility to calculate zonal averages online.
To select the zonal average with `Y` latitude bins, request
`X.zonal_mean(bins=Y)`. All zonal averages are calculated
using the area fraction in each bin as the weight.

For 180 latitude bins, the bins are defined
as follows: [-90, -89), [-89, -88), ..., [89, 90].
For 90 latitude bins, the bins are defined as follows:
[-90, -88), [-88, -86), ..., [88, 90].
And so on...

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.zonal_mean(bins=Y)` | Area fraction | Average across the zonal direction |

## Histograms

We currently have a utility to calculate histograms online.
To select the histogram, request `X.histogram(bins=[..])` with the edges
of the bins. For example, the histogram specified by
`T_mid.histogram(bins=[250,270,290,310])` has 5 bins effectively defined by
[$-\infty$, 250), [250, 270), [270, 290), [290, 310), and [310, $\infty$).

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.histogram(bins=[V0,V1,...,VN])` | 1 or 0 | Count of field values within each range |

## Example

```yaml
%YAML 1.1
---
filename_prefix: monthly.outputs
averaging_type: average
max_snapshots_per_file: 1
fields:
  physics_pg2:
    field_names:
      # in this example, we use T_mid in units of K
      - T_col      := T_mid.weighted(area).mean('col')        # K
      - T_dpavg    := T_mid.weighted(pseudo_density).mean('lev')        # K
      - T_dpsum    := T_mid.weighted(pseudo_density/gravit).sum('lev')  # K * kg / m2
      - T_dzavg    := T_mid.weighted(dz).mean('lev')                    # K
      - T_dzsum    := T_mid.weighted(dz).sum('lev')                     # K * m
      - T_vavg     := T_mid.mean('lev')                       # K
      - T_vsum     := T_mid.sum('lev')                        # K
      - T_vsum_col := T_mid.sum('lev').weighted(area).mean('col')  # K
      - T_col_vsum := T_mid.weighted(area).mean('col').sum('lev')  # K
      - T_dz_col   := T_mid.weighted(dz).sum('lev').weighted(area).mean('col')              # K * m
      - T_dp_col   := T_mid.weighted(pseudo_density).mean('lev').weighted(area).mean('col') # K
      - T_zon180   := T_mid.zonal_mean(bins=180)              # K
      - T_zon90    := T_mid.zonal_mean(bins=90)               # K
      - T_hist     := T_mid.histogram(bins=[250,270,290,310]) # count
output_control:
  frequency: 1
  frequency_units: nmonths
```
