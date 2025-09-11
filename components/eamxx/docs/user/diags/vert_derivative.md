# Vertical derivative diagnostics

We can automatically calculate the vertical derivative of a given field
with respect to a vertical coordinate (pressure or height).
In other words, the vertical derivative is the derivative of the field
along a vertical coordinate.
To automatically get the derivative of a field X, we can request
`X_pvert_derivative` or `X_zvert_derivative`
in the output yaml files.

*WARNING: support for the dz-based derivative is experimental!*

## Brief description of algorithm

For simplicity, the field `X` can only have two dimensions
(column and level). This covers most fields of interest, but
notably does not cover higher-dimensioned fields in radiation and MAM.
We calculate the derivative by finding the difference in the
input field `X`, which must be provided on midpoints, by interpolating
its values at interfaces weighted by the pressure differential.
That is, for a pressure-based derivative

$$
\nabla_{p} X = \frac{\Delta X}{\Delta p}
$$

and for height-based derivative

$$
\nabla_{p} X = \frac{\Delta X}{\Delta z}
$$

where

$$
\Delta X = \tilde{X}_{i,k+1} - \tilde{X}_{i,k}
$$

$$
\tilde{X}_{i,k+1} = \frac{\left(
    X_{i,k} \times \Delta p_{i,k+1} +
    X_{i,k+1} \times \Delta p_{i,k}
  \right)}
  {\left(
    \Delta p_{i,k+1} + \Delta p_{i,k}
  \right)}
$$

Here, $\Delta p$ is the "pseudo density" pressure differential in Pa;
and $\Delta z$ is the height of the layer in m

Note that the definition $\tilde{X}_{i,k+1}$ assumes that all
fields are weighted by the pressure differential, which is a
prognostic parameter of the dynamics core solver. Also note that
when $\Delta p_{i,k+1} \gg \Delta p_{i,k}$,
$\tilde{X}_{i,k+1}$ is closer to $X_{i,k}$ than $X_{i,k+1}$.

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
      - T_mid_pvert_derivative  # K / Pa
      - T_mid_zvert_derivative  # K / m
output_control:
  frequency: 1
  frequency_units: nmonths
```
