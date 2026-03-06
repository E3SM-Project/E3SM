# Diagnostic parsing precedence

EAMxx resolves diagnostic field names through a series of pattern tests applied
in a fixed priority order.  Understanding this order lets you predict — and
control — how composite names are parsed.

## Priority order (highest to lowest)

1. **Plain model field** — the output layer first checks whether the name
   exists directly in the model field manager (FM).  If it does, no diagnostic
   parsing is performed at all; the field is used as-is.

2. **User-defined `:=` aliases** — if the name was declared as a `:=` alias in
   the output YAML, it is resolved to the aliased field or expression before
   any diagnostic parsing.

   Items 3–9 below apply only when the name is not found by steps 1–2, i.e.
   they describe the priority order *inside* the diagnostic parser
   (`create_diagnostic`).

3. **Built-in aliases** — recognized shorthands expand before anything else.
   E.g., `X_atm_backtend` → `X_minus_X_prev_over_dt`.
   See [Built-in aliases](builtin_aliases.md).

4. **`_over_dt` suffix** — checked before binary ops to prevent `X_over_dt`
   from being misinterpreted as `BinaryOpsDiag(X, over, dt)`.

5. **Specific named patterns** — `_at_<level>`, `_at_<pressure>`, `_at_<height>`,
   `_horiz_avg`, `_vert_avg`, `_vert_sum`, `_zonal_avg`, `_pvert_derivative`,
   `_zvert_derivative`, `_where_..._op_val` (conditional sampling), etc.

6. **Binary ops** — `A_(plus|minus|times|over)_B`.
   The first capture group is *greedy*, so the left operand extends as far as
   possible: `A_minus_B_over_C` parses as `(A_minus_B) over C`, i.e., `(A−B)/C`.
   Binary ops are checked **before** `_prev` so that `X_minus_X_prev` correctly
   resolves as `BinaryOpsDiag(X, minus, X_prev)` rather than
   `FieldPrevDiag(X_minus_X)`.

7. **`_prev` suffix** — checked after binary ops so that the right-hand operand
   of a binary op can itself be a `_prev` field.

8. **Histograms** — `X_histogram_<bin_config>`.

9. **Plain diagnostic name** — the string is looked up directly in the
   atmosphere-diagnostic factory.

## Key rules

### Rule 1: greedy left operand for binary ops

The binary-ops regex is:

```none
([A-Za-z0-9_.+\-\*\÷]+)_(plus|minus|times|over)_([A-Za-z0-9_.+\-\*\÷]+)$
```

The first group is greedy.  For a string with multiple operator tokens, the
**rightmost** operator becomes the outermost (lowest-precedence) operation:

```none
A_minus_B_over_C  →  (A_minus_B) over C  →  (A − B) / C
```

### Rule 2: `_prev` is resolved after binary ops

`binary_ops` is tested before `field_prev`, so `X_minus_X_prev` parses
as `BinaryOpsDiag(X, minus, X_prev)`.  Then `X_prev` recurses and is
resolved as `FieldPrevDiag(X)`.  A bare `X_prev` (no operator keyword)
still resolves as `FieldPrevDiag(X)` because `binary_ops` does not match
(there is no operator keyword in the string).

## Worked example: the `bt` family

```yaml
fields:
  physics_pg2:
    aliases:
      # bt1 is an intermediate — needed by bt2 and bt_prod but not written to NC
      - bt1:=f_minus_f_prev_over_dt
      # Parsing:
      #   _over_dt suffix matched first → FieldOverDtDiag( f_minus_f_prev )
      #   f_minus_f_prev → BinaryOpsDiag( f, minus, f_prev )
      #   f_prev → FieldPrevDiag( f )
      # Result: (f(t1) − f(t0)) / (t1 − t0)

      # bt2 = tendency of the tendency (also intermediate)
      - bt2:=bt1_minus_bt1_prev
      # Only one "minus" token → BinaryOpsDiag( bt1, minus, bt1_prev )
      # bt1_prev → FieldPrevDiag( bt1 )

    field_names:
      # bt_prod = product of the two tendencies — IS written to NC
      - bt_prod:=bt1_times_bt2

      # bt_osc_count = indicator that bt_prod is negative (tendency oscillating)
      - bt_osc_count:=count_where_bt_prod_lt_0
      # ConditionalSampling: input=count, condition_field=bt_prod, op=lt, value=0
```

The `aliases` section names intermediate sub-expressions that are needed
as inputs to other diagnostics but should not appear in the output file.
The `field_names` section lists what actually gets written to NetCDF.
Use `:=` in both sections to give a convenient name to a complex expression.

## Avoiding ambiguity

Use `:=` aliases to break complex expressions into named sub-expressions.
This is always clearer than relying on implicit precedence.
