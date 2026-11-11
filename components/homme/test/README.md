# Homme tests

Homme's tests are plain CMake/CTest tests: `ctest` runs everything, and each test runs a CMake
script (`cmake -P`) or an executable directly. There are no generated shell scripts.

## Layout

| Folder    | Contents |
|-----------|----------|
| `execs/`  | All the executables used by the tests (and by `cases/`), in one file. Each one is a build of a model (e.g., `theta-l`) with given compile-time sizes (levels, tracers, ...), and is built in `<build>/test/execs/<name>`. |
| `unit/`   | Unit tests of the C++ code (Catch2). Many of them compare the C++ implementation against the F90 one. |
| `system/` | Tests that run an executable with one or more namelists, and check the results. One folder per model, or per group of test cases. The namelists of a test are next to it. `system/README` describes the tests. |
| `cases/`  | Setups for hand-run simulations (e.g., the DCMIP cases), with job scripts and plotting scripts. These are not ctest tests. |
| `data/`   | Input files used by tests and cases (vertical coordinates, meshes). |
| `cmake/`  | `HommeTesting.cmake` (the functions used to create the tests), and the scripts that ctest runs (`cmake/testing`). Also `ctest_script.cmake`, used by `test_homme.py`. |

`test_homme.py` configures, builds and runs the tests for a supported machine (it is what the GitHub CI uses).

## Running the tests

```bash
cmake -C <homme>/cmake/machineFiles/<machine>.cmake <homme>
make -j8 test-execs    # build all the executables needed by the tests
ctest                  # or: make check
```

The tests must run inside an allocation. The MPI launcher is set by the standard CMake variables
`MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS` and `MPIEXEC_POSTFLAGS` (see the machine files),
and `USE_NUM_PROCS` caps the number of ranks of each test.

Useful commands:

```bash
ctest -N                    # list the tests
ctest -N -L theta-l_kokkos  # list the tests with a label
ctest -R swtc1 -V           # run swtc1 (and its checks), printing the commands
```

Labels: `system`, `unit`, `baseline_gen` (tests whose outputs are baselines), `baseline_cmp` (comparisons against
baselines), `cxx_vs_f90`, the model of the executable (`sweqx`, `preqx`, `preqx_kokkos`, `theta-l`, `theta-l_kokkos`,
`tool`), and the testing profile (`dev`, `short`, `nightly`) for the tests created for each profile.

A system test `X` runs in `<build>/tests/X`, which contains its inputs and its outputs (the log `X.log`, and
`movies/*.nc`). The test runs

```bash
<mpi-launcher> <exe> X.nl X.log   # Homme reads the namelists from X.nl, and writes its log in X.log
```

(`ctest -R X -V` prints the exact command). A test with N namelists is split in N tests `X_run1`, ..., `X_runN`,
which ctest runs in order. Its checks are separate tests, which use `X` as a ctest fixture (so running a check
also runs `X` first):

| Test                | Check |
|---------------------|-------|
| `X_baseline_cmp`    | The outputs match the baselines (cprnc) |
| `X_cxx_vs_f90`      | The netcdf outputs of the F90 test `X` and of its C++ counterpart `X-kokkos` are bit-for-bit identical (BFB builds only) |
| `X_ref_cmp`         | Some outputs match some reference outputs of the same test, within a tolerance |
| `X_l2err`           | The error norms printed in the log are within bounds |
| `X_omp`             | Not a check: runs the hybrid MPI+OpenMP namelists of `X` (only if `ENABLE_HORIZ_OPENMP=ON`) |

## Baselines

Baselines are the outputs of the `baseline_gen` tests. The comparisons look in `${HOMME_BASELINE_DIR}/tests/baseline/X/movies`
(by default, `HOMME_BASELINE_DIR` is the build dir).

```bash
# Generate baselines in <dir>/tests/baseline
cmake -C <machine-file> -DHOMME_BASELINE_DIR=<dir> -DHOMME_GENERATE_BASELINES=ON <homme>
make -j8 test-execs && ctest -L baseline_gen

# Compare against them
cmake -C <machine-file> -DHOMME_BASELINE_DIR=<dir> <homme>
make -j8 test-execs && ctest
```

The `baseline` target (used by CIME's `HOMME` test) runs `ctest -L baseline_gen`, and copies the outputs in
`<build>/tests/baseline`. Use `-DHOMME_SKIP_BASELINE_COMPARE=ON` to not create the comparisons against baselines.

## Adding a test

1. If no executable in `execs/CMakeLists.txt` has the right configuration, add one with `homme_add_model_exec`.
2. Add the namelist(s) in the right folder of `system/`, and declare the test in its `CMakeLists.txt`:

   ```cmake
   homme_add_system_test (my-test
     EXEC         theta-l-nlev20-native
     RANKS        16
     NAMELISTS    my-test.nl
     VCOORD_FILES ${HOMME_TEST_DATA_DIR}/vcoord/*20*
     OUTPUTS      my_output1.nc)
   ```

   Namelists are configured with `configure_file`, so they can use CMake variables set before the call
   (e.g., `ne = ${HOMME_TEST_NE}`). If `HOMME_TEST_QSIZE` is set, `${HOMME_TEST_QVARS}` expands to the
   tracers to output (e.g., `output_varnames1 = 'ps','T'${HOMME_TEST_QVARS}`). Use `homme_add_f90_cxx_system_tests`
   to also create the C++ version (`my-test-kokkos`, using the executable `<EXEC>-kokkos`) and, optionally, the
   `cxx_vs_f90` comparison. Since that comparison only looks at the netcdf files, output all the fields you want
   to compare (including the tracers).
   See `cmake/HommeTesting.cmake` for all the options.
