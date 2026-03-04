# Previous-timestep field diagnostic

In EAMxx, we can capture the value of any field at the beginning of the
previous timestep. This is useful for comparing the state of the model
before and after a process updates a field within a single timestep.

## Configuration

To use this diagnostic, suffix a field name `X` with `_prev` in the output
requests:

```
X_prev
```

The output field has the same layout, units, and grid as the input field `X`.

## Behavior

At the start of each timestep (during `init_timestep`), the current value of
`X` is saved. When the diagnostic is evaluated (during `compute_diagnostic`),
the saved value is written to output. This means `X_prev` always reflects the
state of `X` at the beginning of the current timestep (i.e., before any
processes have modified it in the current step).

On the very first diagnostic evaluation — before `init_timestep` has been
called — the output is filled with the standard EAMxx fill value
(`constants::fill_value<Real>`).

## Example

```yaml
field_names:
  # T_mid has units K
  - T_mid        # current value of T_mid
  - T_mid_prev   # value of T_mid at the start of the current timestep
```

## Caveats

- The `_prev` suffix is stripped to identify the tracked field. If the model
  has a registered field whose name literally ends in `_prev`, the diagnostic
  parser will attempt to resolve it as a `FieldPrevDiag` rather than a direct
  field lookup. In practice this is unlikely to cause conflicts because literal
  model fields are resolved before diagnostics are considered.
- Only fields with `Real` data type are supported.

Contact developers on github if you have additional questions.
