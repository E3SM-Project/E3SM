# Field contraction diagnostics

In EAMxx, we can automatically calculate field reductions
across the horizontal columns and across the model vertical levels.
We call these horizontal and vertical reductions.
We can also automatically calculate zonal averages and histograms.
We have *experimental* support for composing diagnostics; below,
the horizontal and vertical reductions can be composed
sequentially, but if using dp- or dz-weighted reductions,
one must do the vertical reduction first.

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

To select the horizontal reduction, you only need to suffix
a field name `X` with `_horiz_avg` in the output requests.

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_horiz_avg` | Area fraction | Average across all columns |

## Vertical reduction

We currently offer three vertical reductions $C$, defined as

$$
C_{\dots} = \sum_{k} w_{k} \cdot F_{\dots k}
$$

where $F_{\dots k}$ is the field at level $k$,
and $w_{k}$ is the weight at level $k$.

To select the vertical reduction, you only need to suffix
a field name `X` with `_vert_(avg|sum)_(dp|dz)_weighted` or
`_vert_(avg|sum)` in the output yaml files.

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_vert_avg_dp_weighted` | $\Delta p_{k}$ | Average across all levels, weighted by $\Delta p_{k}$ |
| `X_vert_sum_dp_weighted` | $\Delta p_{k}$ | Sum across all levels, weighted by $\Delta p_{k}$ |
| `X_vert_avg_dz_weighted` | $\Delta z_{k}$ | Average across all levels, weighted by $\Delta z_{k}$ |
| `X_vert_sum_dz_weighted` | $\Delta z_{k}$ | Sum across all levels, weighted by $\Delta z_{k}$ |
| `X_vert_avg` | 1 | Average across all levels |
| `X_vert_sum` | 1 | Sum across all levels |

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
To select the zonal average, you need to suffix
a field name `X` with `_zonal_avg` and the
number of bins `Y` as `_Y_bins`. All zonal averages are calculated
using the area fraction in each bin as the weight.

For 180 latitude bins, the bins are defined
as follows: [-90, -89), [-89, -88), ..., [89, 90).
For 90 latitude bins, the bins are defined as follows:
[-90, -88), [-88, -86), ..., [88, 90).
And so on...

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_zonal_avg_Y_bins` | Area fraction | Average across the zonal direction |

## Histograms

We currently have a utility to calculate histograms online.
To select the histogram, you need to suffix a field name `X` with
`_histogram_` followed by the values specifying the edges of the bins,
separated by `_`. For example, the histogram specified by
`T_mid_histogram_250_270_290_310` has 5 bins effectively defined by
[$-\infty$, 250), [250, 270), [270, 290), [290, 310), and [310, $\infty$).

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_histogram_V0_V1_..._VN` | 1 or 0 | Count of field values within each range |

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
      - T_mid_horiz_avg  # K
      - T_mid_vert_avg_dp_weighted  # K
      - T_mid_vert_sum_dp_weighted  # K * Pa * s / (m * m)
      - T_mid_vert_avg_dz_weighted  # K
      - T_mid_vert_sum_dz_weighted  # K * m
      - T_mid_vert_avg  # K
      - T_mid_vert_sum  # K
      - T_mid_vert_sum_horiz_sum  # K
      - T_mid_horiz_avg_vert_sum  # K
      - T_mid_vert_sum_dz_weighted_horiz_avg  # K * m
      - T_mid_vert_avg_dp_weighted_horiz_avg  # K
      - T_mid_zonal_avg_180_bins  # K
      - T_mid_zonal_avg_90_bins  # K
      - T_mid_histogram_250_270_290_310 # count
output_control:
  frequency: 1
  frequency_units: nmonths
```
