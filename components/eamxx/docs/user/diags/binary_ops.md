# Binary operations diagnostics

In EAMxx, we can perform binary arithmetic operations on pairs of fields
to create new diagnostic outputs. The binary operations diagnostic allows
you to compute element-wise arithmetic operations between two fields.

## Supported operations

The binary operations diagnostic supports four basic arithmetic operations:

| Operator | Symbol | Description |
| -------- | ------ | ----------- |
| Addition | `plus` | Element-wise addition of two fields |
| Subtraction | `minus` | Element-wise subtraction of two fields |
| Multiplication | `times` | Element-wise multiplication of two fields |
| Division | `over` | Element-wise division of two fields |

## Requirements

For two fields to be compatible for binary operations, they must satisfy:

1. **Same layout**
2. **Same data type**
3. **Same grid**
4. **Compatible units** for addition and subtraction

## Unit handling

The resulting diagnostic field will have units determined by the operation.

## Configuration

To use the binary operations diagnostic, one can request an output
using the general syntax of `<field_1>_<binary_op>_<field_2>`:

- `field_1` is the name of the first input field
- `binary_op` is the operator: `plus`, `minus`, `times`, or `over`
- `field_2` is the name of the second input field

## Example

```yaml
field_names:
  # T_mid has units K
  - T_mid_times_p_mid  # K*Pa
  - T_mid_over_p_mid   # K/Pa 
  - T_mid_plus_T_mid   # K
  - T_mid_minus_T_mid  # K
```

## Caveats

- As the name suggests, these diagnostics were written for a single
  operation connecting two fields; doing multiple binary ops in
  succession carries risks; read on.
- Strictly speaking, multiple operations will take place in order.
  For example, `qc_plus_qv_times_p_mid` will first compute `qc + qv`,
  then multiply the result by `p_mid`. (Beware, we do NOT support
  mathematical operation precedence; it is simply about the parser
  which processes the expression left-to-right.)
- In fact, `p_mid_times_qc_plus_qv` will fail because of units!
- We only support existing fields in the Field Manager, e.g.,
  grid information like lat, lon, and area are not available.
- We do not support dimensionality broadcasting, so both fields
  must have the same shape.

Contact developers on github if you have additional questions.
