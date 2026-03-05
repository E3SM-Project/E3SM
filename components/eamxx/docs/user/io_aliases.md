# Field Aliasing

This document demonstrates the field aliasing feature for EAMxx I/O operations.

## Overview

The field aliasing feature allows users to specify custom names for
variables in netcdf output files while maintaining the original
internal field names in the model. This is useful for:

- Creating shorter, more convenient variable names for output
- Maintaining compatibility with external tools that expect specific variable names
- Providing user-friendly names for complex internal field names

## Syntax

The alias syntax uses the delimiter `:=` to separate the alias name
from the internal field name:

```yaml
alias_name:=internal_field_name
```

## YAML Configuration Examples

### Basic Usage

```yaml
field_names:
  - "LWP:=LiqWaterPath"          # Alias LWP for LiqWaterPath
  - "SWP:=SolidWaterPath"        # Alias SWP for SolidWaterPath  
  - "T:=T_mid"                   # Alias T for T_mid
  - "qv"                         # Regular field name (no alias)
```

### Mixed Usage

You can mix aliased and non-aliased fields in the same configuration:

```yaml
field_names:
  - "T_mid"                        # Regular field name
  - "LWP:=LiqWaterPath"            # Aliased field
  - "p_mid"                        # Regular field name  
  - "temp:=T_mid"                  # Aliased field
  - "surf_temp:=temp_at_model_bot" # Alias of diagnostic of an alias
```

## Output Behavior

When using aliases:

1. **NetCDF Variables**: The netcdf file will contain variables
named according to the aliases

    - `LWP` instead of `LiqWaterPath`
    - `T` instead of `T_mid`
    - `RH` instead of `RelativeHumidity`

2. **Internal Processing**: The FieldManager used by IO stores
both the original and the alias fields, with the latter being an
alias of the former.

3. **Metadata**: Variable attributes (units, long_name, etc.)
are preserved from the original fields, and `alias_of`
is added to the netcdf files to document aliasing

## Intermediate-only fields (`:=:` syntax)

Sometimes a diagnostic is only needed as a building block for another
diagnostic, and you do not want it written to the NetCDF output file.
Use `alias:=:original` (note the extra colon) to declare such fields:

```yaml
field_names:
  - "MyDiag:=:T_mid_times_p_mid"   # computed but NOT written to NC
  - "MyDiagSqrd:=MyDiag_times_MyDiag"  # uses MyDiag; IS written to NC
```

`MyDiag` is created and registered in the internal field manager so that
`MyDiagSqrd` can depend on it, but it does not appear in the output file.
Only `MyDiagSqrd` is written.

Rules:

- A name declared with `:=:` must not also appear as a regular output field.
- The same alias name cannot be used for both `:=:` and `:=`.
- The expression after `:=:` is resolved with the same composable-diagnostic
  parser as any other field name.

## Caveats

Currently, a field can be requested only once in a single stream,
but we do allow requesting the same field again but with an alias

```yaml
field_names:
  - "LiqWaterPath"          # OK: a known diag field
  - "LWP:="                 # Error: empty field name
  - ":=LiqWaterPath"        # Error: empty alias name  
  - "LWP:=LiqWaterPath"     # OK: an alias
  - "liq_w_path:=Field2"    # OK: another alias, but different name
  - "LWP:=LiqWaterPath"     # Error: duplicate field LWP (even if the same as before)
  - "LWP:=T_mid"            # Error: duplicate field LWP
```
