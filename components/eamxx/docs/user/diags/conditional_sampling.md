# Conditional sampling diagnostics

The conditional sampling diagnostic allows you to extract field values
where a specified condition is met, filling other locations with a
default fill value. This is useful for analyzing field behavior under
specific atmospheric conditions or at particular levels.
The goal of the diagnostics is to enable you to discard unwanted data
via *other* reductions, such as horizontal and vertical reductions
(see "Field contraction" section of the docs).

## Overview

Conditional sampling creates a new field where:

- Values from the input field are preserved
  where the condition is **true**
- Fill values (default: large number) are used
  where the condition is **false**

The general syntax is the [expression](expressions.md) `X.where(<lhs> <op> <rhs>)`.

## Supported operators

| Operator | Description |
| -------- | ----------- |
| `==` | Equal to |
| `!=` | Not equal to |
| `>` | Greater than |
| `>=` | Greater than or equal to |
| `<` | Less than |
| `<=` | Less than or equal to |

## Field-based conditional sampling

Compare the input field against another field at each grid point.

**Examples**:

- `T_mid.where(qv>0.01)`
- `p_mid.where(T_mid<=273.15)`
- `qv.where(p_mid<50000)`

## Level-based conditional sampling

Compare against the vertical level index (0-based indexing).
Use the special condition operand `lev`.

**Examples**:

- `T_mid.where(lev>5)`
- `p_mid.where(lev==0)`
- `qv.where(lev<=10)`
- `T_mid.where(qv.isel(lev=125)<qv.isel(lev=128))`

## Mask-only calculation

Compute a mask stating where the condition is met.
Put `mask` where the sampled field goes. It is a placeholder, not a field:
nothing named `mask` exists in the field manager, and it is recognized only as
the receiver of `.where(..)`, so it cannot be requested or operated on by
itself. The output will be `1` where the condition is satisfied and `0`
elsewhere. This is particularly useful when combined with horizontal or
vertical reductions to count occurrences of specific conditions.

**Examples**:

- `mask.where(qv>0.01)`
- `mask.where(T_mid<=273.15)`
- `mask.where(p_mid<50000)`
- `mask.where(lev>5)`

**A mask cannot currently be written to file as is.** The mask diagnostic's
output has `int` data type, and IO cannot write it: listing `mask.where(..)` in
`field_names` aborts the run with a narrowing-conversion error, and reducing it
does not help. The workaround is to promote it to `Real` by multiplying by 1.
E.g., the fraction of the earth where the surface pressure is > 1e5Pa:

```yaml
field_names:
  - wet_frac := (mask.where(ps>100000)*1.0).mean('col')
```

## Caveats

- For now, we only support 1D, 2D, and 3D fields.
  These include most fields of interest with any possible sub-layout
  of (ncol,ncmp,nlev)
  Adding support for higher-dimensioned fields is straightforward,
  but we will add it when requested.
- The condition and input fields must have the same layout
  (except for level-based sampling). In level-based sampling,
  the level index must be present in the field.
- We assume the user knows the standard unit of the lhs/rhs that are used
  internally in EAMxx field manager (there is no way to specify it).
- The RHS in the condition can be a scalar value or another field/diag
- `where` takes a single comparison. `and`/`or` are not supported; chain
  `where(..)` calls instead.

## Example configuration

```yaml
%YAML 1.1
---
filename_prefix: conditional_sampling_outputs
averaging_type: instant
fields:
  physics_pg2:
    field_names:
      # Field-based conditional sampling
      - T_wet   := T_mid.where(qv>0.01)
      - p_cold  := p_mid.where(T_mid<=273.15)
      - qv_high := qv.where(p_mid<50000)
      # Level-based conditional sampling
      - T_low   := T_mid.where(lev>5)
      - p_top   := p_mid.where(lev==0)
      - qv_top  := qv.where(lev<=10)
      # Mask-based conditional sampling (promoted to Real, see above)
      - wet     := mask.where(qv>0.01)*1.0
      - cold    := mask.where(T_mid<=273.15)*1.0
      - low     := mask.where(lev>5)*1.0
output_control:
  frequency: 6
  frequency_units: nhours
```
