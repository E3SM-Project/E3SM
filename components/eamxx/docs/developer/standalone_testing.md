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
standalone tests for an EAMxx configuration.

## Running EAMxx's Tests with CTest

Before running the tests, generate a baseline file:

```
cd $RUN_ROOT_DIR
make baseline
```

The tests will run, automatically using the baseline file, which is located in
the CMake-configurable path `${SCREAM_TEST_DATA_DIR}`. By default, this path is
set to `data/` within your build directory (which is `$RUN_ROOT_DIR`, in
our case).

To run all of SCREAM's tests, make sure you're in `$RUN_ROOT_DIR` and type

```
ctest -VV
```

This runs everything and reports results in an extra-verbose (`-VV`) manner.

You can also run subsets of the SCREAM tests. For example, to run only the
P3 regression tests (again, from the `$RUN_ROOT_DIR` directory), use

```
ctest -R p3_regression
```

### Grouping Tests with Labels

We can create groupings of tests by using **labels**. For example, we have a
`driver` label that runs tests for SCREAM's standalone driver. You can see a
list of available labels by typing

```
ctest --print-labels
```

To see which tests are associated with a given label (e.g. `driver`), use

```
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

