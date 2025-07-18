# Vertical derivative diagnostics

We can automatically calculate the vertical derivative of a given field
with respect to a vertical coordinate (pressure or height).
In other words, the vertical derivative is the derivative of the field
along a vertical coordinate.
To automatically get the derivative of a field X, we can request
`dX_dp_vert_derivative` or `dX_dz_vert_derivative`
in the output yaml files.

*WARNING: support for the dz-based derivative is experimental!*

## Brief description of algorithm

For simplicity, the field `X` can only have two dimensions
(column and level). This covers most fields of interest, but
notably does not cover higher-dimensioned fields in radiation and MAM.
We calculate the derivative by weighting two finite differences.
The input field `X` must be defined on the midpoints, and its
difference is defined as the difference between two consecutive
levels, pointing downwards. That is, for a pressure coordinate,

$$
\nabla_{p} X = w  \frac{x}{a} + (1-w)\frac{y}{b}
$$

where

- $x = X_{i, k} - X_{i, k-1}$ and $a = p_{i,k} - p_{i,k-1}$
- $y = X_{i, k+1} - X_{i, k}$ and $b = p_{i,k+1} - p_{i,k}$
- $w = b / (a + b)$

and $x = 0$, $w=0$ for $k=0$ (highest level)
and $y = 0$, $w=1$ for $k=K$ (lowest level).

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
      - dT_mid_dp_vert_derivative  # K / Pa
      - dT_mid_dz_vert_derivative  # K / m
output_control:
  frequency: 1
  frequency_units: nmonths
```
