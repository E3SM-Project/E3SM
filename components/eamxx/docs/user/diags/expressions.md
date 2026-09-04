# Diagnostic expressions

As the underscore-based system is being deprecated, EAMxx can now
build diagnostics from **expressions**:

```yaml
field_names:
  - q_flux := (qc+qv)*p_mid
  - cold_T  := T_mid.where(T_mid<273.15).mean('lev')
```

Expressions are parsed by [`dexpr`](https://github.com/E3SM-Project/E3SM/tree/master/share/dexpr),
the expression parser shared across E3SM components. They exist because the
underscore syntax cannot express grouping: `qc_plus_qv_times_p_mid` is grouped
by regex greediness, not by operator precedence (see
[Binary arithmetics](binary_ops.md) and [Parsing precedence](parsing_precedence.md)).
Written as an expression, it means what it looks like.

Both syntaxes work, and they are resolved independently: a name is matched
against the underscore patterns first, and only names that match none of them
are read as an expression. Nothing you write today changes meaning.

## Expressions must be given a name

An expression is not a usable NetCDF variable name, and `field_names` entries
are written to file verbatim. So any expression you want written must be given
a name with `:=`:

```yaml
field_names:
  - q_flux := (qc+qv)*p_mid    # good
  - (qc+qv)*p_mid              # error: "An expression must be given an output name"
```

The expression itself is recorded in the variable's `alias_of` attribute, so
the provenance stays in the file.

This applies only to what is *written*. Expressions nested inside other
expressions, and expressions in the `aliases:` section (intermediates that other
diagnostics use but that are not written), need no name of their own.

## Operators

| Operator | Meaning |
| -------- | ------- |
| `a + b`, `a - b`, `a * b`, `a / b` | element-wise arithmetic |
| `-a` | negation, shorthand for `(-1)*a` |
| `>`, `>=`, `<`, `<=`, `==`, `!=` | comparison, **only** inside `where(..)` |
| `( )` | grouping |

Operands may be fields, other expressions, named physical constants (`Rgas`,
`P0`, ...), or numeric literals:

```yaml
field_names:
  - p_hPa      := p_mid/100
  - T_scaled   := T_mid*0.5
  - dry_air    := p_mid - qv*p_mid
  - p_ratio    := p_mid / P0
```

Precedence is the usual one, so `qc+qv*p_mid` is `qc+(qv*p_mid)`. Note this
differs from the underscore spelling `qc_plus_qv_times_p_mid`, which means
`(qc+qv)*p_mid`. When in doubt, parenthesize. (That pair only illustrates
grouping: `qc+(qv*p_mid)` adds a mixing ratio to a pressure, and is rejected at
run time for incompatible units, while `(qc+qv)*p_mid` is fine.)

## Functions

Functions are written as methods on the field they act on, `X.f(..)`, and are
named after their [xarray](https://docs.xarray.dev) counterparts wherever there
is an honest analogue.

| Expression | Equivalent underscore name |
| ---------- | -------------------------- |
| `X.isel(lev=10)` | `X_at_lev_10` |
| `X.isel(lev=0)` | `X_at_model_top` |
| `X.isel(lev=-1)` | `X_at_model_bot` |
| `X.interp(p_mid=500, units='hPa')` | `X_at_500hPa` |
| `X.interp(z_mid=10, reference='surface')` | `X_at_10m_above_surface` |
| `X.mean('col')` | `X_horiz_avg` |
| `X.mean('lev')` | `X_vert_avg` |
| `X.sum('lev')` | `X_vert_sum` |
| `X.mean('lev', weights='dp')` | `X_vert_avg_dp_weighted` |
| `X.where(qv>0.01)` | `X_where_qv_gt_0.01` |
| `X.derivative('p_mid')` | `X_pvert_derivative` |
| `X.histogram(bins=[0,1,2])` | `X_histogram_0_1_2` |
| `X.zonal_mean(bins=20)` | `X_zonal_avg_20_bins` |
| `X.shift(time=1)` | `X_prev` |
| `X.over_dt()` | `X_over_dt` |
| `X.tend()` | `X_atm_backtend` |

Notes:

- `isel` indexes levels the way xarray does: `lev=0` is the model top and
  `lev=-1` the bottom. Other negative indices are rejected, since they cannot
  be resolved without knowing the level count.
- `interp` takes exactly one of `p_mid` or `z_mid`, named after the coordinate
  as in xarray. `units` defaults to `Pa` and `m` respectively; `reference`
  (`'surface'` or `'sealevel'`, default `'sealevel'`) applies to `z_mid` only.
- **`shift` looks BACK in time**, following xarray's sign convention: a positive
  shift moves data forward, so `shift(time=1)` is the value from the *previous*
  step. `shift(time=-1)` would be the value from a step *ahead*, which a running
  model cannot know, and is rejected -- do not read it as "one step back" by
  analogy with `isel(lev=-1)`. Only `time=1` is available, since only one step
  of history is kept.
- `mean('col')` is always area weighted, so `weights` does not apply to it.
  `weights` for `'lev'` is `'dp'` or `'dz'`.
- `where` takes a single comparison. `and`/`or` are not supported; chain
  `where(..)` calls instead.
- **`mask` and `lev` are placeholders, not fields.** Neither exists in the field
  manager, and neither means anything on its own. `mask` is only legal as the
  receiver of `.where(..)`: `mask.where(cond)` is the 0/1 indicator of `cond`
  itself, rather than some field sampled by it. `lev` is only legal as the left
  operand of the comparison inside `where(..)`, where it means the level index.
  So `mask.where(lev>5)` is a request, while `mask`, `mask*1.0`,
  `mask.mean('lev')` and `lev>5` are not: each of those leaves `mask` or `lev`
  as an operand to be resolved on its own, and the run aborts with
  `The key 'mask' is not associated to any registered product`, since it is
  neither a field nor a diagnostic. Both names are inherited from the underscore
  syntax (`mask_where_lev_gt_5`), not from xarray.
- **`over_dt` and `tend` cannot go in an `instant` stream.** An instant stream
  takes a snapshot at `t0`, before any step has run, and `dt` is zero there, so
  the diagnostic aborts with "dt must be positive". Write them to an averaged
  stream (`averaging_type: average`, `max`, `min`) instead. This is a property
  of the diagnostic, not of the syntax: the underscore spellings `X_over_dt` and
  `X_atm_backtend` behave the same way.
- **A `mask.where(..)` result cannot currently be written to file as is.** The
  diagnostic is built, but its output has `int` data type, which IO cannot
  handle: listing `mask.where(..)` (or the underscore `mask_where_..`) in
  `field_names` fails with a narrowing conversion error, and reducing it first
  does not help. Multiplying by 1 promotes it to `Real`, and
  `mask.where(ps>100000)*1.0` writes fine. Note the `*1.0` applies to the field
  that `where` produced, which is a perfectly ordinary field; it is not an
  operation on `mask`, which per the bullet above cannot be multiplied.
- `X.tend()` is shorthand for `(X - X.shift(time=1)).over_dt()`, and inherits
  the `over_dt` restriction above.
- Histogram bin edges must be non-negative, and must be writable without an
  exponent: the diagnostic joins them with `_` and splits them back apart, so
  `1e3` is fine (it is written `1000.0`) but `1e30` is not.

Named quantities are not operations, and stay plain names in an expression:
`LiqWaterPath`, `RelativeHumidity`, `z_mid`, `dz`, `precip_liq_surf_mass_flux`,
and so on.

### Where this departs from xarray

`isel`, `interp`, `mean`, `sum`, `where` and `shift` mean what they mean in
xarray, on the same arguments. The rest have no xarray spelling we could adopt:

| | |
| --- | --- |
| `histogram` | not in xarray; the keyword follows [xhistogram](https://xhistogram.readthedocs.io) |
| `zonal_mean` | xarray would be `groupby_bins('lat', 20).mean()`, a call chain we cannot represent |
| `over_dt`, `tend` | no xarray equivalent |
| `mean('lev', weights=..)` | xarray would be `weighted(dp).mean('lev')`, again a chain |

`X.mean('lev', weights='dp')` is `(X*dp).sum('lev')/dp.sum('lev')`, but computed
in one sweep rather than as three diagnostics.

Two operand names have no xarray analogue at all, because they are not fields:
`mask` and `lev`, both inherited from the underscore syntax and both legal only
inside `where(..)`. See the note above for what they mean and where they may
appear.

## Errors

An expression that does not parse, calls a function that does not exist, or
uses an operator no diagnostic implements is rejected with a message pointing
at the offending part. Unknown function names come back with the list of
available ones.

```text
Error! Unsupported expression: no diagnostic implements the operator '**'.
 - subexpression: (T_mid**2)
```

## Worked example

```yaml
fields:
  physics_pg2:
    aliases:
      # An intermediate: used below, but not written to file. No name needed
      # for the expression itself beyond the alias.
      - warm_T := T_mid.where(T_mid>273.15)

    field_names:
      # Column-mean temperature, but only where it is above freezing
      - warm_T_colavg := warm_T.mean('lev')
      # Surface temperature where the air is moist
      - sfc_moist_T   := T_mid.where(qv>0.0).isel(lev=-1)
      # Moisture flux
      - q_flux        := (qc+qv)*p_mid
```

## Adding a function

The set of callable functions is registered by EAMxx, in
`components/eamxx/src/share/io/eamxx_dexpr_diags.cpp`, not by `dexpr` itself.
Adding one means adding a `FunctionSpec` to the registry there and a case to the
translator that maps it onto a diagnostic; nothing under `share/dexpr` changes.

Each `FunctionSpec` carries an `example` of how the call is written. That is not
just documentation: two checks use it, so a function that does not hang together
fails fast rather than at the moment a user first writes it.

- `dexpr::validate_registry`, called where the registry is built, proves the
  example parses, matches the spec that declared it (arity, keyword names,
  required keywords), and really calls that function. Declare `mean` as taking
  no positional arguments while writing `T_mid.mean('lev')` and this catches it.
- The `dexpr_every_registered_function_is_buildable` case in
  `src/share/io/tests/create_diag.cpp` runs every example through
  `create_diagnostic`, so a function registered without a translator case is
  caught too.

So adding a function means: register it with an example, translate it, and both
checks come along for free. The `dexpr` command line tool does the same for the
generic vocabulary -- `dexpr check "<expr>"` and `dexpr check-registry`.

Those two checks stop at "a diagnostic was built". Every expression on this page
is also run end to end -- computed and written to NetCDF -- by the standalone
test in
`components/eamxx/tests/multi-process/dynamics_physics/homme_shoc_cld_p3_rrtmgp_pg2`,
which carries two extra output streams for the purpose: `output_diags_ins.yaml`
and, for the ones an instant stream cannot hold, `output_diags_avg.yaml`. Keep
them in step with this page: an example added here belongs in one of those.
