# Field contraction diagnostics

In EAMxx, we can automatically calculate field reductions
across the horizontal columns and across the model vertical levels.
We call these horizontal and vertical reduction.

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
a field name `X` with `_vert_(avg|sum)`

| Reduction | Weight | Description |
| --------- | ------ | ----------- |
| `X_vert_avg` | $\Delta p_{k}$ | Average across all levels, weighted by $\Delta p_{k}$ |
| `X_vert_sum` | $\Delta p_{k}$ | Sum across all levels, weighted by $\Delta p_{k}$ |

The only supported weighting for now is that of
`pseudo_density` field in EAMxx, $\Delta p_{k}$.

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
    - T_mid_horiz_avg
    - T_mid_vert_avg
    - T_mid_vert_sum
output_control:
  frequency: 1
  frequency_units: nmonths
```
