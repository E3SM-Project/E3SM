# histComputeMod / subgridAveViewMod unit test

Checks the per-step history accumulation (`histComputeMod`) and the per-gridcell
averages it shares with `subgridAveMod` (`subgridAveViewMod`) against
`hist_reference_mod`: the original serial ELM algorithm (`hist_update_hbuf_field_1d/2d`,
`hist_set_snow_field_2d` and the serial `p2g/c2g/l2g/t2g` loops), applied one field at a
time. The reference keeps the original rules for `'X'`/`'M'` (described in
`histComputeMod`), implements the new `'U'`/`'L'` (extremes of the valid samples), and
makes one deliberate fix that `histComputeMod` makes too: every flag indexes the model
array from its own lower bound (the original code did so only for 1d `'A'` fields, so 1d
`'I'/'X'/'M'` fields read the wrong elements when the pointer's lower bound had been
reset to 1).

The test runs on a synthetic subgrid (`fake_subgrid_mod`), passed to the code as a
`subgrid_ave_view` instead of ELM's global subgrid instances, with ELM's index arrays
(`topi/topf`, `coli/colf`, `pfti/pftf`) and the landunits of each gridcell. It is laid out like ELM's:
landunits are ordered by type across gridcells, so a gridcell's points aren't contiguous,
gridcells have one or two topounits, and proc bounds don't start at 1. It has random
weights (some zero), active flags, snow layers and urban/landunit-type scale tables.
About 180 fields cover:

- every subgrid level, 1d and 2d;
- model arrays whose lower bound was reset to 1;
- snow fields dimensioned `-nlevsno+1:0`;
- all time operations (`I`, `A`, `X`, `M`, `U`, `L`) and scale types;
- gridcell-averaged and native output.

Every step gets new values (10% missing), and the subgrid changes every 7 steps. After
every update, and after every end-of-interval normalize, all buffers and counters must
equal the reference bit for bit. Every step also checks `p2g_range`, `c2g_range`,
`l2g_range` and `t2g_range` against the serial loops: over all gridcells on the host, and in
chunks of gridcells in a kernel, as `subgridAveMod` calls them. The whole test runs twice:
with normal snow-layer order, then with `extrasnowlayers`.

With OpenACC the model arrays and the subgrid are put on the device, as in ELM; the test
copies new values there each step and the buffers back before comparing.

```bash
cmake -S . -B build && cmake --build build && ctest --test-dir build --output-on-failure
# nvfortran, OpenACC on host threads
cmake -S . -B build-mc -DCMAKE_Fortran_COMPILER=nvfortran -DACC=MULTICORE && cmake --build build-mc
build-mc/test_hist_compute
# nvfortran, OpenACC on a GPU (e.g. a Perlmutter login node's A100)
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DACC=GPU && cmake --build build-gpu
build-gpu/test_hist_compute
# timing of the original algorithm vs histComputeMod
build/test_hist_compute bench
```

The GPU build is compiled with `-Mnofma`: the reference runs on the host, and fused
multiply-adds on the device change the averages in the last bits (up to ~3e-15
relative). nvfortran 26.5 ignores `-gpu=nofma` at `-O0`; `-Mnofma` works at any
optimization level.

For gfortran multithreading, translate the OpenACC directives to OpenMP with SPEL's
`acc2omp` utility.
