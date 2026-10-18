# Previous-timestep field diagnostic

In EAMxx, we can capture the value of any field at the beginning of the
previous timestep. This is useful for comparing the state of the model
before and after a process updates a field within a single timestep.

## Configuration

To use this diagnostic, suffix a field name `X` with `_prev` in the output
requests:

```none
X_prev
```

The output field has the same layout, units, and grid as the input field `X`.

## Behavior

At the start of each timestep (during `init_timestep`), the current value of
`X` is saved. When the diagnostic is evaluated (during `compute_diagnostic`),
the saved value is written to output. This means `X_prev` always reflects the
state of `X` at the beginning of the current timestep (i.e., before any
processes have modified it in the current step).

### Instant streams and derived diagnostics

An instant stream evaluates diagnostics on output steps only. When `X` is
itself a computed diagnostic, `X_prev` in an instant stream is therefore only
available if the stream writes every step (`frequency_units: nsteps` with
`frequency: 1`); otherwise the request is rejected when the stream is set up. An averaged stream evaluates its diagnostics every
step, and a model field `X` is always available, so neither is restricted.

### The previous value of a previous value

`X_prev_prev` is rejected when the stream is set up, and so is any other
previous value of a previous value, e.g. `X.shift(time=1).shift(time=1)`, the
same through an alias, or `X.shift(time=1).tend()`. It would not be two steps
back: `init_timestep` runs for the inner `X_prev` first, so by the time the
outer one saves its input, that input already holds the value of the current
step's start, and the result is one step back, like `X_prev`. Any other
diagnostic is evaluated after `init_timestep`, so to look two steps back, put
an operation between the two: `(X.shift(time=1)*1.0).shift(time=1)`.

### First-step initialisation for derived diagnostics

When `X` is itself a computed diagnostic (rather than a raw model field) it
has not been evaluated yet when `init_timestep` runs at the start of the
first step. In that case the `init_timestep` call has nothing to capture (the
source has no valid timestamp), and `X_prev` is invalid on the first step: it
is written as the fill value, and so is anything computed from it, such as
`X - X_prev`. From the second step onward `X_prev` reflects the normally
captured start-of-step value.

## Example

```yaml
field_names:
  # T_mid has units K
  - T_mid        # current value of T_mid
  - T_mid_prev   # value of T_mid at the start of the current timestep
```

## Caveats

- The `_prev` suffix is stripped to identify the tracked field. If the model
  has a registered field whose name literally ends in `_prev`, it will be
  output directly as a model field (model fields are resolved before
  diagnostics are considered), not treated as a `FieldPrevDiag`.
- Only fields with `Real` data type are supported.

Contact developers on GitHub if you have additional questions.
