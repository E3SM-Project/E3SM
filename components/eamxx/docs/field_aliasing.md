# EAMxx Field Aliasing Feature

This document demonstrates the field aliasing feature for EAMxx I/O operations.

## Overview

The field aliasing feature allows users to specify custom names for variables in netcdf output files while maintaining the original internal field names in the model. This is useful for:

- Creating shorter, more convenient variable names for output
- Maintaining compatibility with external tools that expect specific variable names
- Providing user-friendly names for complex internal field names

## Syntax

The alias syntax uses the delimiter `:=:` to separate the alias name from the internal field name:

```yaml
alias_name:=:internal_field_name
```

## YAML Configuration Examples

### Basic Usage

```yaml
output_streams:
  stream1:
    filename_prefix: output
    averaging_type: INSTANT
    fields:
      physics_gll:
        field_names:
          - "LWP:=:LiqWaterPath"      # Alias LWP for LiqWaterPath
          - "SWP:=:SolidWaterPath"    # Alias SWP for SolidWaterPath  
          - "T:=:T_mid"               # Alias T for T_mid
          - "qv"                      # Regular field name (no alias)
```

### Mixed Usage

You can mix aliased and non-aliased fields in the same configuration:

```yaml
field_names:
  - "T_mid"                          # Regular field name
  - "LWP:=:LiqWaterPath"            # Aliased field
  - "ps"                            # Regular field name  
  - "RH:=:RelativeHumidity"         # Aliased field
```

## Output Behavior

When using aliases:

1. **NetCDF Variables**: The netcdf file will contain variables named according to the aliases
   - `LWP` instead of `LiqWaterPath`
   - `T` instead of `T_mid`
   - `RH` instead of `RelativeHumidity`

2. **Internal Processing**: All internal model operations use the original field names
   - Field validation uses `LiqWaterPath`, `T_mid`, etc.
   - Diagnostic calculations use original names
   - Memory management uses original field structures

3. **Metadata**: Variable attributes (units, long_name, etc.) are preserved from the original fields

## Input Reading

The same alias syntax works for reading input files:

```yaml
input_streams:
  initial_conditions:
    filename: ic_file.nc
    field_names:
      - "LWP:=:LiqWaterPath"    # Read variable LWP from file into field LiqWaterPath
      - "T:=:T_mid"             # Read variable T from file into field T_mid
```

## Error Handling

The system provides comprehensive error checking:

- **Duplicate Aliases**: Each alias must be unique within a field list
- **Malformed Syntax**: Proper error messages for invalid alias specifications
- **Empty Names**: Both alias and field names must be non-empty

Example error cases:

```yaml
field_names:
  - "LWP:=:"                    # Error: empty field name
  - ":=:LiqWaterPath"          # Error: empty alias name  
  - "LWP:=:Field1"             # OK
  - "LWP:=:Field2"             # Error: duplicate alias LWP
```

## Implementation Notes

- **Backward Compatible**: Existing configurations without aliases continue to work unchanged
- **Minimal Performance Impact**: Alias processing occurs only during initialization  
- **Thread Safe**: All alias operations are performed during setup phase
- **Memory Efficient**: No data copying, only name mapping during I/O operations

## Example Use Cases

### 1. Climate Model Output
```yaml
field_names:
  - "T:=:T_mid"
  - "RH:=:RelativeHumidity"  
  - "P:=:ps"
```

### 2. Atmospheric Chemistry
```yaml
field_names:
  - "QC:=:qc"
  - "QI:=:qi"
  - "QV:=:qv"
```

### 3. Diagnostic Output
```yaml
field_names:
  - "CLDTOT:=:cldfrac_tot_for_analysis"
  - "CLDLIQ:=:cldfrac_liq"
  - "CLDICE:=:cldfrac_ice_for_analysis"
```
