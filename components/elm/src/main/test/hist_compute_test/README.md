# histComputeMod unit test

Checks `histComputeMod` (the per-step history accumulation behind `histFileMod`)
against `hist_reference_mod`. The reference is the algorithm ELM used before
`histComputeMod`: `hist_update_hbuf_field_1d/2d`, `hist_set_snow_field_2d` and
subgridAveMod's `p2g/c2g/l2g/t2g`, applied one field at a time and serially. It
differs from the old code in two deliberate fixes, which `histComputeMod` makes too:

- every averaging flag indexes the model array from its own lower bound (before,
  only 1d `'A'` fields did, so 1d `'I'/'X'/'M'` fields read the wrong elements when
  the pointer's lower bound had been reset to 1);
- `'X'` and `'M'` follow the `'A'` rule: only valid samples count, and a missing
  sample never overwrites the extreme found so far (before, the gridcell path
  replaced it with spval, and a buffer whose first sample was missing kept spval
  for the rest of the interval).

The test runs on a synthetic subgrid (`fake_subgrid_mod`). It is laid out like
ELM's: landunits are ordered by type across gridcells, so a gridcell's points
aren't contiguous, and proc bounds don't start at 1. It has random weights
(some zero), active flags, snow layers and urban/landunit-type scale tables.
About 180 fields cover:

- every subgrid level, 1d and 2d;
- model arrays whose lower bound was reset to 1;
- snow fields dimensioned `-nlevsno+1:0`;
- all averaging flags and scale types;
- gridcell-averaged and native output.

Every step gets new values (10% missing), and the subgrid changes every 7 steps.
After every update, and after every end-of-interval normalize, all buffers and
counters must equal the reference bit for bit. The whole test runs twice: once
with normal snow-layer order and 256-point tiles, then with `extrasnowlayers`
and 3-point tiles.

```bash
cmake -S . -B build && cmake --build build && ctest --test-dir build --output-on-failure
# nvfortran, OpenACC on host threads
source ~/.local/nvhpc-libs/env.sh
cmake -S . -B build-mc -DCMAKE_Fortran_COMPILER=nvfortran -DACC=MULTICORE && cmake --build build-mc
build-mc/test_hist_compute
# OpenMP, and the timing comparison with the old algorithm
cmake -S . -B build-omp -DOPENMP=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build-omp
OMP_NUM_THREADS=8 build-omp/test_hist_compute bench
```
