# Standalone EAMxx Testing

In this section we describe our testing methodology for standalone EAMxx
configurations. We use several types of tests

* **Unit tests** are individual test programs that demonstrate that a small set
  of code performs a single function or a set of related functions. We use
  a C++ unit testing framework called [Catch2](https://catch2-temp.readthedocs.io/en/latest/index.html)
  to implement unit tests.
* **Property (verification) tests** are test programs that configure code that
  demonstrates that a part of EAMxx (for example, an atmospheric physics
  parameterization or the dynamical core) is able to produce an answer that
  satisfies some physical constraint or matches a known solution under specific
  circumstances.
* **Fortran-C++ "bit-for-bit" (BFB) tests** are test programs, often implemented
  as unit tests, that demonstrate that a set of C++ code ported from Fortran
  produces bit-for-bit identical results to its Fortran counterpart, provided
  certain compiler options are enabled (such as "strict" floating-point
  arithmetic).
* **Test Suites** are named collections of tests that can be run on demand using
  the [ctest](https://cmake.org/cmake/help/latest/manual/ctest.1.html) command.

We also support a `test-all-scream` configuration that runs all of the
standalone tests for an EAMxx configuration. Note, your current machine
must be known to EAMxx before `test-all-scream` will work. A machine can
be made known to EAMxx by editing the eamxx/scripts/machines_specs.py files.
There are some instructions on what to do at the top of this file.

`test-all-scream` has a good help dump

```shell
cd $scream_repo/components/eamxx
./scripts/test-all-scream -h
```

If you are unsure of the cmake configuration for you development cycle, one
trick you can use is to run `test-all-scream` for the `dbg` test and just
copy the cmake command it prints (then ctrl-C the process).

```shell
cd $scream_repo/components/eamxx
./scripts/test-all-scream -t dbg -m $machine
# wait for a few seconds*
# Ctrl-C *
# Copy the contents of DCMAKE_COMMAND that was passed to ctest *
# Add "cmake" to beginning of contents and path to eamxx at the end. *
```

Considerations for using `test-all-scream`:

* Your machine must be known to our scripts, see above.
* If you try to run commands by-hand (outside of test-all-scream;
  cmake, make, ctest, etc), you'll need to remember to
  load the scream-env into your shell, which can be done like this:
  `cd eamxx/scripts; eval $(./scripts/scream-env-cmd $machine)`
* test-all-scream expects to be run from a compute node if you
  are on a batch machine.
* You'll need to think about your baseline situation, as many of our
  tests rely on pre-existing baselines. The -b option controls the baseline
  location and can have the following values:
  * AUTO: A common public baseline area shared by all developers
  * LOCAL: A private baseline area for the current developer in the current repo
  * $path: A specific arbitrary path
  * None: If there is no -b at all, no baseline testing will be done

## Running EAMxx's Tests with CTest

Before running the tests, generate a baseline file:

```shell
cd $RUN_ROOT_DIR
make baseline
```

The tests will run, automatically using the baseline file, which is located in
the CMake-configurable path `${SCREAM_BASELINES_DIR}`. By default, this path is
set to an invalid string. If baselines tests are enabled, we check that a valid
path has been provided.

To run all of SCREAM's tests, make sure you're in `$RUN_ROOT_DIR` and type

```shell
ctest -VV
```

This runs everything and reports results in an extra-verbose (`-VV`) manner.

You can also run subsets of the SCREAM tests. For example, to run only the
P3 regression tests (again, from the `$RUN_ROOT_DIR` directory), use

```shell
ctest -R p3_regression
```

### Grouping Tests with Labels

We can create groupings of tests by using **labels**. For example, we have a
`driver` label that runs tests for SCREAM's standalone driver. You can see a
list of available labels by typing

```shell
ctest --print-labels
```

To see which tests are associated with a given label (e.g. `driver`), use

```shell
ctest -L driver -N
```

## EAMxx Test Suites

### The `p3_regression` Suite

`p3_regression` uses a baseline file to compare any new or altered
implementations with our P3 Fortran reference implementation. If you're working
on the C++/Kokkos implementation, you can invoke any new tests to the function
`Baseline::run_and_cmp` in
`${SCREAM_SRC_DIR}/components/eamxx/p3/tests/p3_run_and_cmp.cpp`.

If the reference Fortran implementation changes enough that a new baseline file
is required, make sure to let other SCREAM team members know, in order to
minimize disruptions.
