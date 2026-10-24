# Previous-timestep field diagnostic

In EAMxx, we can capture the value of any field at the beginning of the
previous timestep. This is useful for comparing the state of the model
before and after a process updates a field within a single timestep.

## Configuration

To use this diagnostic, request the [expression](expressions.md)
`X.shift(time=1)`:

```yaml
field_names:
  - X_prev := X.shift(time=1)
```

The output field has the same layout, units, and grid as the input field `X`.

## Behavior

At the start of each timestep (during `init_timestep`), the current value of
`X` is saved. When the diagnostic is evaluated (during `compute_diagnostic`),
the saved value is written to output. This means `X.shift(time=1)` always
reflects the state of `X` at the beginning of the current timestep (i.e.,
before any processes have modified it in the current step).

### Instant streams and derived diagnostics

An instant stream evaluates diagnostics on output steps only. When `X` is
itself a computed diagnostic, `X.shift(time=1)` in an instant stream is
therefore only available if the stream writes every step
(`frequency_units: nsteps` with `frequency: 1`); otherwise the request is
rejected when the stream is set up. An averaged stream evaluates its
diagnostics every step, and a model field `X` is always available, so neither
is restricted.

### The previous value of a previous value

`X.shift(time=1).shift(time=1)` is rejected when the stream is set up, and so
is any other previous value of a previous value, e.g. the same through an
alias, or `X.shift(time=1).tend()`. It would not be two steps back:
`init_timestep` runs for the inner shift first, so by the time the outer one
saves its input, that input already holds the value of the current step's
start, and the result is one step back, like `X.shift(time=1)`. Any other
diagnostic is evaluated after `init_timestep`, so to look two steps back, put
an operation between the two: `(X.shift(time=1)*1.0).shift(time=1)`.

### First-step initialisation for derived diagnostics

When `X` is itself a computed diagnostic (rather than a raw model field) it
has not been evaluated yet when `init_timestep` runs at the start of the
first step. In that case the `init_timestep` call has nothing to capture (the
source has no valid timestamp), and `X.shift(time=1)` is invalid on the first
step: it is written as the fill value, and so is anything computed from it,
such as `X - X.shift(time=1)`. From the second step onward `X.shift(time=1)`
reflects the normally captured start-of-step value.

### Restarts

No diagnostic state is saved in a restart file. If `X` is a model field in the
restart set, a restart gives it a valid value before the first step, so
`X.shift(time=1)` continues as if the run had not stopped. Otherwise (a
computed diagnostic, or a model field outside the restart set, such as one
imported from the coupler during the step), the first step after a restart is
like the first step of a run: `X.shift(time=1)` is the fill value, and an
averaged stream counts one sample fewer in that output window. A shift of an
expression that contains a shift, such as
`(X.shift(time=1)*1.0).shift(time=1)`, is the fill value for two steps.

## Example

```yaml
field_names:
  # T_mid has units K
  - T_mid                              # current value of T_mid
  - T_mid_prev := T_mid.shift(time=1)  # value of T_mid at the start of the current timestep
```

## Caveats

- Only `time=1` is available, since only one step of history is kept.
- Only fields with `Real` data type are supported.

Contact developers on GitHub if you have additional questions.
