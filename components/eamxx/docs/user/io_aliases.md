# EAMxx Field Aliasing Feature

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
  - "LWP:=LiqWaterPath"      # Alias LWP for LiqWaterPath
  - "SWP:=SolidWaterPath"    # Alias SWP for SolidWaterPath  
  - "T:=T_mid"               # Alias T for T_mid
  - "qv"                     # Regular field name (no alias)
```

### Mixed Usage

You can mix aliased and non-aliased fields in the same configuration:

```yaml
field_names:
  - "T_mid"                        # Regular field name
  - "LWP:=LiqWaterPath"            # Aliased field
  - "p_mid"                        # Regular field name  
  - "RH:=RelativeHumidity"         # Aliased field
```

## Output Behavior

When using aliases:

1. **NetCDF Variables**: The netcdf file will contain variables
named according to the aliases

   - `LWP` instead of `LiqWaterPath`
   - `T` instead of `T_mid`
   - `RH` instead of `RelativeHumidity`

2. **Internal Processing**: All internal model operations use the
original field names

   - Field validation uses `LiqWaterPath`, `T_mid`, etc.
   - Diagnostic calculations use original names
   - Memory management uses original field structures

3. **Metadata**: Variable attributes (units, long_name, etc.)
are preserved from the original fields, and `eamxx_name`
is added to the netcdf files to document aliasing

## Caveats

Currently, a field can be requested only once in a single stream,
and either the original name or the alias name counts.

```yaml
field_names:
  - "LWP:="                   # Error: empty field name
  - ":=LiqWaterPath"          # Error: empty alias name  
  - "LWP:=Field1"             # OK
  - "LWP:=Field2"             # Error: duplicate alias LWP
  - "LWP1:=LiqWaterPath"      # OK
  - "LWP2:=LiqWaterPath"      # Error: duplicate field LiqWaterPath
```
