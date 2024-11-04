# Quality Control (QC) testing for MPAS-Seaice

This testing script to determine if the answers are climate-changing as determined by a 2 stage t-test was adapted from [CICE](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_testing.html#code-validation-test-non-bit-for-bit-validation).

## Setup

1. Setup and build a baseline and test case with recommended minimum resolution `RES=TL319_EC30to60E2r2`, and `COMPSET=DTESTM-JRA1p5`.
2. Turn on the ice thickness analysis member with the following namelist options in `user_nl_mpassi`:

   config_AM_thicknesses_enable = .true.
   config_AM_thicknesses_compute_on_startup = .true.
   config_AM_thicknesses_write_on_startup = .true.

3. In the run directory's `streams.seaice` file, modify the `output` stream to have daily snapshot output by setting `output_interval="00-00-01_00:00:00"` (default is "none"). Copy this modified `streams.seaice` file to the case directory in `SourceMods/src.mpassi/.`.
4. Run each case for (at least) 5 years.

## Usage

    python mpas-seaice.t-test.py $BASE $TEST

where `$BASE`, `$TEST` are the paths to the run directories of the two tests containing the `mpassi.hist.*` files. A sample batch script for Chrysalis and Anvil is provided, `job_script.cice-qc.anvil`, `job_script.cice-qc.chrysalis`.

## Output

A `qc_log.txt` file provides testing progress and results. Map *png files are also generated (currently produces blank output, support still needed).