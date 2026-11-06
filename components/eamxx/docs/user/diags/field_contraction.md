# Field contraction diagnostics

In EAMxx, we can automatically calculate field reductions
across the horizontal columns and across the model vertical levels.
We call these horizontal and vertical reductions.
We can also automatically calculate zonal averages and histograms.
All of these are written as [expressions](expressions.md), and compose:
e.g., `X.mean('lev',weights='dp').mean('col')`. If using dp- or dz-weighted
vertical reductions, do the vertical reduction first.

## Horizontal reduction

We currently offer only one horizontal reduction $C$, and it is defined as

$$
C_{\dots k} = \sum_{i} w_{i} \cdot F_{i \dots k}
$$

where $F_\text{i...k}$ is the field at column $i$ and level $k$,
and $w_{i}$ is the weight at column $i$.
We note that the field $F$ can have other dimensions ($\dots$).
The weight $w$ is defined as the area fraction in column $i$,
that is, the area in column $i$ divided by the total area in all columns.

To select the horizontal reduction, request `X.mean('col')`.

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.mean('col')` | Area fraction | Average across all columns |

## Vertical reduction

We currently offer three vertical reductions $C$, defined as

$$
C_{\dots} = \sum_{k} w_{k} \cdot F_{\dots k}
$$

where $F_{\dots k}$ is the field at level $k$,
and $w_{k}$ is the weight at level $k$.

To select the vertical reduction, request `X.mean('lev')` or `X.sum('lev')`,
optionally with `weights='dp'` or `weights='dz'`.

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X.mean('lev',weights='dp')` | $\Delta p_{k}$ | Average across all levels, weighted by $\Delta p_{k}$ |
| `X.sum('lev',weights='dp')` | $\Delta p_{k}$ | Sum across all levels, weighted by $\Delta p_{k}$ |
| `X.mean('lev',weights='dz')` | $\Delta z_{k}$ | Average across all levels, weighted by $\Delta z_{k}$ |
| `X.sum('lev',weights='dz')` | $\Delta z_{k}$ | Sum across all levels, weighted by $\Delta z_{k}$ |
| `X.mean('lev')` | 1 | Average across all levels |
| `X.sum('lev')` | 1 | Sum across all levels |

The supported weighting options for now are

- `pseudo_density` field in EAMxx, $\Delta p_{k}$, in units of Pa;
- `dz` field in EAMxx, $\Delta z_{k}$, in units of m;
- and no weighting, which is equivalent to using a weight of 1.

In the case of `pseudo_density`, the weighting is scaled by 1/g,
where g is the gravitational acceleration, in units of m/s$^2$.

It is important to note that in the case of dp- and dz-weighted
reductions, the field provided to the reduction utility must be
at least 2 in rank, e.g., (ncol, nlev).

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
      - T_col      := T_mid.mean('col')                       # K
      - T_dpavg    := T_mid.mean('lev',weights='dp')          # K
      - T_dpsum    := T_mid.sum('lev',weights='dp')           # K * Pa * s / (m * m)
      - T_dzavg    := T_mid.mean('lev',weights='dz')          # K
      - T_dzsum    := T_mid.sum('lev',weights='dz')           # K * m
      - T_vavg     := T_mid.mean('lev')                       # K
      - T_vsum     := T_mid.sum('lev')                        # K
      - T_vsum_col := T_mid.sum('lev').mean('col')            # K
      - T_col_vsum := T_mid.mean('col').sum('lev')            # K
      - T_dz_col   := T_mid.sum('lev',weights='dz').mean('col')  # K * m
      - T_dp_col   := T_mid.mean('lev',weights='dp').mean('col') # K
      - T_zon180   := T_mid.zonal_mean(bins=180)              # K
      - T_zon90    := T_mid.zonal_mean(bins=90)               # K
      - T_hist     := T_mid.histogram(bins=[250,270,290,310]) # count
output_control:
  frequency: 1
  frequency_units: nmonths
```
