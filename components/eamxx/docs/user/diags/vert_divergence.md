# vertical divergence diagnostics

We can automatically calculate the vertical divergence of a given field
with respect to the vertical coordinate (pressure or height).
In other words, the vertical divergence is the derivative of the field
along a vertical coordinate.
To automatically get the divergence of a field X, we can request
`dX_dp_vert_divergence` or `X_dz_vert_divergence`
in the output yaml files.

## Brief description

For simplicity, the field `X` can only have two dimensions
(column and level). This covers most fields of interest, but
notably does not cover higher-dimensioned fields.
We also assume that the highest level (top of the atmosphere)
has a zero derivative.
The input field `X` must be defined on the midpoints, and its
difference is defined as the difference between two consecutive
levels, pointing downwards. That is, for a pressure coordinate,

$$
\nabla_{\text{vert}} X = \frac{X_{i, k} - X_{i, k-1}}{\Delta p_{k}}
$$

for all $k>0$ (at $k=0$, the derivative is zero).

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
      - dT_mid_dp_vert_divergence  # K / Pa
      - dT_mid_dz_vert_divergence  # K / m
output_control:
  frequency: 1
  frequency_units: nmonths
```
