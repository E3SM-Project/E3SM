# Vertical derivative diagnostics

We can calculate the vertical derivative of a field with respect to a vertical
coordinate, pressure or height, by requesting the
[expression](expressions.md) `X.differentiate('p_mid')` or
`X.differentiate('z_mid')` in the output yaml files. The name and the numbers
are xarray's: the result is what `X.differentiate('p_mid')` gives in xarray on
one column whose level coordinate is that column's `p_mid`.

## Algorithm

For simplicity, the field `X` can only have two dimensions (column and level).
This covers most fields of interest, but notably does not cover
higher-dimensioned fields in radiation and MAM. `X` and the coordinate `x`
(`p_mid` or `z_mid`) are both on midpoints.

In each column, the derivative is numpy's `gradient` of `X` with respect to
the values of `x`, which is what xarray's `differentiate` computes. At an
interior level $k$, with $h_1 = x_k - x_{k-1}$ and $h_2 = x_{k+1} - x_k$, it is
second-order accurate on uneven spacing:

$$
\left.\frac{\partial X}{\partial x}\right|_k =
  \frac{-h_2}{h_1 (h_1+h_2)} X_{k-1}
+ \frac{h_2 - h_1}{h_1 h_2} X_k
+ \frac{h_1}{h_2 (h_1+h_2)} X_{k+1}
$$

and at the top and bottom levels it is one-sided and first-order, as with
xarray's default `edge_order=1`:

$$
\left.\frac{\partial X}{\partial x}\right|_0 = \frac{X_1 - X_0}{x_1 - x_0},
\qquad
\left.\frac{\partial X}{\partial x}\right|_{n-1} = \frac{X_{n-1} - X_{n-2}}{x_{n-1} - x_{n-2}}.
$$

The derivative is with respect to the coordinate as it stands. Level indices
increase downward, so `p_mid` increases with $k$ and `z_mid` decreases: in a
troposphere where temperature falls with height, `T_mid.differentiate('p_mid')`
is positive and `T_mid.differentiate('z_mid')` negative.

If `X` is masked (e.g., `X.where(..)`), a level of the result is valid only if
every level its difference reads is valid ($k-1$, $k$, $k+1$ inside; the two
levels used at the top and bottom), which is where xarray would not propagate
a NaN. Elsewhere it is missing.

The units of the result are those of `X` divided by those of the coordinate.

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
      - dTdp := T_mid.differentiate('p_mid')  # K / Pa
      - dTdz := T_mid.differentiate('z_mid')  # K / m
output_control:
  frequency: 1
  frequency_units: nmonths
```
