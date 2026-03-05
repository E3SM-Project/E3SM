# Built-in aliases

EAMxx recognises certain diagnostic name patterns as **built-in aliases**.
A built-in alias is expanded to a canonical composable expression *before*
any other regex is tested.  The expansion recurses, so the resulting name
is itself parsed normally.

## Reference table

| User request     | Expands to               | Meaning                                               |
|------------------|--------------------------|-------------------------------------------------------|
| `X_atm_backtend` | `X_minus_X_prev_over_dt` | Atmospheric backward tendency of X: (X − X_prev) / dt |

## Example

```yaml
field_names:
  - T_mid_atm_backtend   # expands to T_mid_minus_T_mid_prev_over_dt
```

This is equivalent to requesting:

```yaml
field_names:
  - T_mid_minus_T_mid_prev_over_dt
```

which in turn chains:

- `FieldPrevDiag(T_mid)` → `T_mid_prev`
- `BinaryOpsDiag(T_mid, minus, T_mid_prev)` → `T_mid_minus_T_mid_prev`
- `FieldOverDtDiag(T_mid_minus_T_mid_prev)` → `T_mid_minus_T_mid_prev_over_dt`

## Notes

- Built-in aliases are checked before all other pattern rules (including suffix
  operators and binary ops), so they always take precedence.
- The match is greedy on the field-name part, so `T_mid_at_500hPa_atm_backtend`
  expands to `T_mid_at_500hPa_minus_T_mid_at_500hPa_prev_over_dt`.
- Additional built-in aliases may be added in `eamxx_io_utils.cpp` in the
  dedicated block at the top of `create_diagnostic()`.
