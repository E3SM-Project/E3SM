# EAMxx Column Extraction Tool

This directory contains tools for extracting single columns from EAMxx NetCDF output files.

## extract_single_column.py

Extract a single atmospheric column from an EAMxx output file for use in single-column tests.

### Requirements

- Python 3.6+
- netCDF4-python (`pip install netCDF4`)
- numpy (`pip install numpy`)

### Usage

```bash
python extract_single_column.py INPUT.nc OUTPUT.nc --col-idx COL_INDEX [OPTIONS]
```

### Arguments

- `INPUT.nc`: Path to the input EAMxx NetCDF output file
- `OUTPUT.nc`: Path for the output single-column NetCDF file
- `--col-idx COL_INDEX`: Column index to extract (0-based, required)
- `--time-idx TIME_INDEX`: Time index to extract (0-based, default: last time step)
- `-v, --verbose`: Print detailed information during extraction

### Examples

```bash
# Extract column 0 from the last time step
python extract_single_column.py \
    mam4_aero_microphys_standalone_output.INSTANT.nsteps_x1.np4.2021-10-12-45000.nc \
    single_col_ic.nc \
    --col-idx 0

# Extract column 100 from time step 5
python extract_single_column.py \
    model_output.nc \
    single_col_ic.nc \
    --col-idx 100 \
    --time-idx 5

# Extract with verbose output
python extract_single_column.py \
    model_output.nc \
    single_col_ic.nc \
    --col-idx 50 \
    --verbose
```

### Output

The output file will contain:
- All variables from the input file
- `ncol` dimension set to 1
- `time` dimension set to 1 (single time snapshot)
- Provenance metadata attributes:
  - `column_extraction_source`: Original input file path
  - `column_extraction_col_idx`: Extracted column index
  - `column_extraction_time_idx`: Extracted time index

### Using Extracted Columns as Initial Conditions

The output file can be used as an initial condition file for single-column EAMxx tests by:

1. Placing the file in the appropriate data directory
2. Updating the test's `input.yaml` to reference the file:
   ```yaml
   initial_conditions:
     filename: /path/to/single_col_ic.nc
   ```
3. Setting `number_of_global_columns: 1` in the grid configuration

See `../tests/single-process/mam/aero_microphys_scm/` for an example single-column test.
