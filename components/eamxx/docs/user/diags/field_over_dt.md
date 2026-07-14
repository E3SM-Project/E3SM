# Field-over-dt diagnostic

Divides any input field by the model timestep `dt = t_current − t_start_of_step`.
This is useful for converting accumulated differences into per-second rates.

## Configuration

Append `_over_dt` to any field name `X`:

```yaml
field_names:
  - X_over_dt
```

The output field has the same layout and grid as `X`, with units `[X_units / s]`.

## Behavior

- At the start of each timestep (`init_timestep`), the start-of-step timestamp
  is saved internally.
- At evaluation time (`compute_diagnostic`), `dt` is computed as the difference
  between the field's current timestamp and the saved start timestamp (in seconds).
- The output is `X / dt` element-wise.
- On the very first evaluation — before `init_timestep` has been called — the
  output is filled with the standard EAMxx fill value (`constants::fill_value<Real>`).

## Example: atmospheric backward tendency

```yaml
field_names:
  - T_mid_minus_T_mid_prev_over_dt
```

This composes three diagnostics:

1. `T_mid_prev` — `FieldPrevDiag`: value of `T_mid` at the start of the step
2. `T_mid_minus_T_mid_prev` — `BinaryOpsDiag`: difference
3. `T_mid_minus_T_mid_prev_over_dt` — `FieldOverDtDiag`: rate

The shorthand `T_mid_atm_backtend` expands to the same expression via the
built-in alias system (see [Built-in aliases](builtin_aliases.md)).

## Caveats

- `dt` must be strictly positive; a zero or negative timestep triggers a runtime error.
- Only fields with `Real` data type are supported.
- The suffix `_over_dt` is matched *before* binary-op patterns, so
  `A_minus_B_over_dt` parses as `FieldOverDtDiag(A_minus_B)`, not
  `BinaryOpsDiag(A_minus_B, over, dt)`.
  See [Parsing precedence](parsing_precedence.md) for details.
