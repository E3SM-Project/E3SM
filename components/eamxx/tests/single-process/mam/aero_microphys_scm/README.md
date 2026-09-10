# MAM4 Aerosol Microphysics Single Column Model (SCM) Test

This test is a single-column version of the MAM4 aerosol microphysics standalone test.
It runs aerosol microphysical processes (nucleation, coagulation, condensation, mode merging) 
for a single atmospheric column.

## Differences from Multi-Column Test

Compared to the `aero_microphys` test (218 columns), this test:
- Uses only **1 column** (`number_of_global_columns: 1`)
- Uses a single-column initial condition file extracted from the multi-column IC
- Has the same physics configuration and input data files
- Runs much faster due to reduced computational domain

## Purpose

Single-column tests are useful for:
- Rapid development and debugging of aerosol microphysics schemes
- Detailed analysis of aerosol processes at specific locations
- Testing code changes before running expensive multi-column tests
- Creating test cases for specific atmospheric conditions

## Initial Conditions

The initial condition file (`EAMxx_tests_IC_FILE_MAM4xx_72lev_SCM`) should be a single-column
NetCDF file with `ncol=1` dimension. This can be created by extracting a column from the
multi-column IC file using the extraction tool:

```bash
cd components/eamxx/scripts/column-extraction
python extract_single_column.py \
    $SCREAM_DATA_DIR/init/screami_unit_tests_mam4xx_ne2np4L72_20240329.nc \
    $SCREAM_DATA_DIR/init/screami_unit_tests_mam4xx_ne2np4L72_scm_20240329.nc \
    --col-idx 0
```

## Running the Test

From the build directory:

```bash
cd components/eamxx
./scripts/test-all-eamxx -m <MACHINE> -t dbg -p mam4_aero_microphys_scm
```

Or run the test directly:

```bash
cd <build_dir>/tests/single-process/mam/aero_microphys_scm
ctest -R mam4_aero_microphys_scm
```

## Configuration

### Time Stepping
- Time step: 1800 seconds (30 minutes)
- Number of steps: Depends on test size (2/5/48 for small/medium/large)
- Start time: October 12, 2021

### Grid
- 1 horizontal column
- 72 vertical levels

### Physics
Same MAM4 aerosol microphysics configuration as the multi-column test:
- Condensation (mam4_do_cond: true)
- New particle nucleation (mam4_do_newnuc: true)
- Coagulation (mam4_do_coag: true)
- Mode merging/renaming (mam4_do_rename: true)

## Output

Generates NetCDF output files with format:
```
mam4_aero_microphys_scm_output.INSTANT.nsteps_x1.np<RANKS>.2021-10-12-45000.nc
```

Output includes aerosol species concentrations for all MAM4 modes and species.
