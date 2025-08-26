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

The general syntax is: `<input>_where_<condition>_<operator>_<value>`

## Supported operators

| Operator | Aliases | Description |
| -------- | ------- | ----------- |
| `eq` | `==` | Equal to |
| `ne` | `!=` | Not equal to |
| `gt` | `>` | Greater than |
| `ge` | `>=` | Greater than or equal to |
| `lt` | `<` | Less than |
| `le` | `<=` | Less than or equal to |

## Field-based conditional sampling

Compare the input field against another field at each grid point.

**Examples**:

- `T_mid_where_qv_gt_0.01`
- `p_mid_where_T_mid_le_273.15`
- `qv_where_p_mid_lt_50000`

## Level-based conditional sampling

Compare against the vertical level index (0-based indexing).
Use the special condition field name `lev`.

**Examples**:

- `T_mid_where_lev_gt_5`
- `p_mid_where_lev_eq_0`
- `qv_where_lev_le_10`

## Count-based conditional sampling

Count the number of grid points where a condition is met.
Use the special input field name `count`. The output will be `1.0`
where the condition is satisfied and `0.0` elsewhere.
This is particularly useful when combined with horizontal or vertical
reductions to count occurrences of specific conditions.

**Examples**:

- `count_where_qv_gt_0.01`
- `count_where_T_mid_le_273.15`
- `count_where_p_mid_lt_50000`
- `count_where_lev_gt_5`

## Caveats

- For now, we only support 1D or 2D fields.
  These include most fields of interest with the following layouts
  (ncol,), (nlev,), (ncol,nlev).  
  Adding support for higher-dimensioned fields is straightforward,
  but we will add it when requested.
- The condition and input fields must have the same layout
  (except for level-based sampling). In level-based sampling,
  the level index must be present in the field.
- The condition is provided as a triplet of information string
  (`<condition>_<operator>_<value>`). We assume the user knows
  the standard unit of the condition that is used internally
  in EAMxx field manager (there is no way to specify it).

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
      - T_mid_where_qv_gt_0.01
      - p_mid_where_T_mid_le_273.15
      - qv_where_p_mid_lt_50000
      # Level-based conditional sampling  
      - T_mid_where_lev_gt_5
      - p_mid_where_lev_eq_0
      - qv_where_lev_le_10
      # Count-based conditional sampling
      - count_where_qv_gt_0.01
      - count_where_T_mid_le_273.15
      - count_where_lev_gt_5
output_control:
  frequency: 6
  frequency_units: nhours
```
