# Binary operations diagnostics

In EAMxx, we can perform binary arithmetic operations on pairs of fields
to create new diagnostic outputs. The binary operations diagnostic allows
you to compute element-wise arithmetic operations between two fields.

## Supported operations

| Operation | Expression | Description |
| --------- | ---------- | ----------- |
| Addition | `a + b` | Element-wise addition of two fields |
| Subtraction | `a - b` | Element-wise subtraction of two fields |
| Multiplication | `a * b` | Element-wise multiplication of two fields |
| Division | `a / b` | Element-wise division of two fields |

Either operand can also be a numeric literal or a named physical constant
(`Rgas`, `P0`, ...), and `-a` is shorthand for `(-1)*a`.

## Requirements

For two fields to be compatible for binary operations, they must satisfy:

1. **Same layout**
2. **Same data type**
3. **Same grid**
4. **Compatible units** for addition and subtraction

## Unit handling

The resulting diagnostic field will have units determined by the operation.

## Configuration

Binary operations are written as [expressions](expressions.md), and an
expression that is written to file must be given a name with `:=`.

## Example

```yaml
field_names:
  # T_mid has units K
  - T_times_p := T_mid*p_mid     # K*Pa
  - T_over_p  := T_mid/p_mid     # K/Pa
  - T_twice   := T_mid+T_mid     # K
  - T_zero    := T_mid-T_mid     # K
  - q_flux    := (qc+qv)*p_mid   # kg/kg*Pa
```

## Caveats

- The usual operator precedence applies, so `qc+qv*p_mid` is
  `qc+(qv*p_mid)`, which fails because of units. Parenthesize when in doubt.
- We only support existing fields in the Field Manager, e.g.,
  grid information like lat, lon, and area are not available.
- We do not support implicit dimensionality broadcasting, so both fields
  must have the same shape. Use `X.broadcast_like(Y)` to make them match.

Contact developers on github if you have additional questions.
