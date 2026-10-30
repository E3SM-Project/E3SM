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

### The first step

`init_timestep` captures `X` as it is at the start of the step, so on the
first step `X_prev` holds a value only if `X` had one before the run started.
A model field always does (its initial condition). A computed diagnostic has
one only if it was evaluated at t0, which happens when an instant stream writes
the initial snapshot (the default; see `skip_t0_output`) and requests `X` or
anything computed from it. Otherwise `X_prev` is invalid on the first step: it
is written as the fill value, and so is anything computed from it, such as
`X - X_prev`. From the second step onward `X_prev` reflects the normally
captured start-of-step value.

The t0 snapshot itself has no previous step, and there the shift of a model
field holds the field's current value. A shift of an expression built on a
shift carries that into the first step. So when an instant stream writes the
t0 snapshot, on the first step `(X.shift(time=1)*1.0).shift(time=1)` holds `X`
at t0 rather than the fill value, and with `bt := X - X_prev`,
`bt - bt.shift(time=1)` equals `bt` rather than a second difference. Both are
right from the second step on. If the first step matters, set
`skip_t0_output: true` in the instant streams.

### Restarts

No diagnostic state is saved in a restart file, and no t0 snapshot is written
after a restart. If `X` is a model field in the restart set, a restart gives
it a valid value before the first step, so `X_prev` continues as if the run had
not stopped. Otherwise (a computed diagnostic, or a model field outside the
restart set, such as one imported from the coupler during the step), the first
step after a restart is like the first step of a run without the t0 snapshot:
`X_prev` is the fill value, and an averaged stream counts one sample fewer in
that output window. A shift of an expression that contains a shift, such as
`(X.shift(time=1)*1.0).shift(time=1)`, is the fill value for one step if `X`
is a model field in the restart set, and for two steps otherwise.

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
