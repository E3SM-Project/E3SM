# Field contraction diagnostics

In EAMxx, we can automatically calculate field reductions
across the horizontal columns and across the model vertical levels.
We call these horizontal and vertical reductions.

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

We currently offer two vertical reductions $C$, defined as

$$
C_{\dots} = \sum_{k} w_{k} \cdot F_{\dots k}
$$

where $F_{\dots k}$ is the field at level $k$,
and $w_{k}$ is the weight at level $k$.

To select the vertical reduction, you only need to suffix
a field name `X` with `_vert_(avg|sum)_(dp|dz)_weighted`

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_vert_avg_dp_weighted` | $\Delta p_{k}$ | Average across all levels, weighted by $\Delta p_{k}$ |
| `X_vert_sum_dp_weighted` | $\Delta p_{k}$ | Sum across all levels, weighted by $\Delta p_{k}$ |
| `X_vert_avg_dz_weighted` | $\Delta z_{k}$ | Average across all levels, weighted by $\Delta z_{k}$ |
| `X_vert_sum_dz_weighted` | $\Delta z_{k}$ | Sum across all levels, weighted by $\Delta z_{k}$ |

The only supported weighting for now is that of either
`pseudo_density` field in EAMxx, $\Delta p_{k}$, in units of Pa,
or `dz` field in EAMxx, $\Delta z_{k}$, in units of m.
In the case of `pseudo_density`, the weighting is scaled by 1/g,
where g is the gravitational acceleration, in units of m/s$^2$.

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
output_control:
  frequency: 1
  frequency_units: nmonths
```
