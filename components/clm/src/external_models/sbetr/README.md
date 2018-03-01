# BeTR

BeTR is a standalone reactive transport libary designed to be
integrated into land surface models such as CLM and ALM.

Jinyun Tang, jinyuntang@lbl.gov

References

Tang, J. Y., Riley, W. J., Koven, C. D., and Subin, Z. M.: CLM4-BeTR, a generic biogeochemical transport and reaction module for CLM4: model development, evaluation, and application, Geosci. Model Dev., 6, 127-140, doi:10.5194/gmd-6-127-2013, 2013.

## Building

[![Build Status](https://travis-ci.org/BeTR-biogeochemistry-modeling/sbetr.svg?branch=master)](https://travis-ci.org/BeTR-biogeochemistry-modeling/sbetr)

### Requirements

Requirements for configuring and building BeTR:

* cmake > 3.1

* fortran compiler

  * gfortran > 5.3.0. Older versions make work, but are not
    supported. At a minimum, you will have to change the compiler
    flags. The pfunit testing frame work can NOT be build with gcc <
    4.9.

  * intel - TBD

  * pgi - TBD

  * nag - TBD

* c compiler - only used for 3rd-party libraries. Any modern compiler
  should work. The following versions are tested:

  * clang - approximately 7.3

  * gcc - 5.3

* python 2.7

All third party dependancies for building and running standalone BeTR
are included in the 3rd-party directory.

### Build


BeTR uses a cmake based build system. The default build is debug. To
build using the default debug configuration:

    cd ${SBETR_ROOT}
    make config
    make all

This will do an out of source build in:

    ${SBETR_ROOT}/build/Xyz-Debug/src

where Xyz is the build configuration.

The standalone `sbetr` executable is in:

    ${SBETR_(ROOT}/build/Xyz-Debug/src/driver/sbetr

To build a release configuration of the code:

    cd ${SBETR_ROOT}
    make debug=0 config
    make debug=0 all

The following command will create local/bin/, where sbetr will be installed.

    make install

Please note, the top level make file providing `make config` etc is a
convenience for the most common use cases. You don't have to use it
and can specify all configuration and build commands manually.


### HPC machines

To run betr on cluster, one needs to load the following (take intel
compiler for example) the compiler, cmake, python (>2.7), and
mkl, the configuration command is then

    make config CC=icc CXX=icpc FC=ifort

and the install command is

    make install CC=icc CXX=icpc FC=ifort

Others should be done similarly as one run betr on a desktop or
laptop, with appropriate modifications as described above.

* yellowstone

  * intel
    ```SH
    module unload ncarbinlibs
    module load intel/16.0.2
    module load mkl/11.3.0
    module load cmake/3.3.1
    module load python/2.7.7

    make CC=icc CXX=icpc FC=ifort config all install test
    ```

  * gnu - doesn't work yet, can't link blas
    ```SH
    module swap intel gnu/5.3.0
    module load cmake/3.3.1
    module load python/2.7.7

    make config all install test
    ```

  * pgi - in progress
    ```SH
    module swap intel pgi
    ```

* edison

  * intel

        module load cmake/3.3.2
        module load intel/15.0.1.133

        make CC=icc CXX=icpc FC=ifort config all install test

* cori

  * intel

        module load cmake/3.3.2
        module load intel/16.0.0.109

        make CC=icc CXX=icpc FC=ifort config all install test

## Testing

BeTR testing inclueds [pFUnit](http://pfunit.sourceforge.net/) based
unit tests and a systems level regression test driver.

pFUnit tests are created automatical during the build. Run the unit
tests with

    make test

or by calling ctest in the build directory.

Regression tests are based on calling the standalone sbetr executable
and checking the results are within a specified epsilon of a baseline.

Regression testing will eventually be integrated into the 'make test'
command with unit tests. For now they have to be run separately.

    cd regression-tests
    make rtest




### Creating new tests

#### Regression tests

Regression tests are grouped into test suites, which are defined by
configuration files. The directory where the configuration file is
referred to as the 'suite directory'. A suite directory can contain
multiple configuratin files. Configuration files are in cfg/ini
format, and have the following information:

```INI
[default_tolerances]
#category = value type
general = 1.0e-14 absolute
concentration = 1.0e-13 relative

[mock]
# override default tolarance
concentration = 1.0e-14 absolute
timeout = 30 seconds

```

Where there is one required section: `default_tolerances`. This
contains the default tolerances for *all tests in this suite.*
Tolerances are specified by category, concentration, velocity,
general. The type of toleraces can be absolute, relative or
percent. These values can be over ridden for individual tests.

In general, tests should be kept short. Long tests are harder to debug
and they rarely provide a benefit over shorter tests. If a particular
set of conditions needs to be tested, it is better to create a special
forcing data set that exercises those conditions in a short test. The
test suite enforces a timeout limit for tests to prevent infinite
loops and other hangs from running. If a test is exceeding the default
timeout limit, then it can be increased on a case by case basis in the
test suite configuration.

All other sections are considered to define tests. The name of the
section is the test name. It is expected that a `test_name.namelist`
file will be present in the same directory as the suite configuration
file. sbetr is run from the same directory as the configuration file,
and all paths in the namelist file must be relative to this directory.
sbetr will write a `test_name.regression` file with the regression
test data. This is compared to `test_name.regression.baseline`, also
in the same directory as the suite configuration file. Keys in the
test section are used to modify the test, either by changing the
tolerance, or modifying some other functionality, e.g. we will
eventually have restart tests triggered by a setting in the test
section.

Setting tolerances is a balancing act that requires some trial and
error. Test tolerances should be set as tightly as possible to
identify small changes in behavior. But they should be loose enough to
be platform independent, so we don't get false positives when moving
to a new platform.

A requirement for testing is that the repository be self contained and
platform independent. All input data and baselines should be contained
in the repository. Inorder to be revision control friendly, the data
should be saved in the repository as plain text files. Netcdf files
should be converted to their text representation with:

    ncdump -p 9,17 ${DATA}.nc > ${DATA}.nc.cdl

These files will be automatically converted back to binary when the
test suite is run with:

    ncgen -o ${DATA}.nc ${DATA}.nc.cdl

Check the cdl file to ensure all variables have the same data type as
is expected in the code, i.e. double percision for r8. Before saving a
cdl file in the repo, verify that the round trip of cdl-nc-cdl results
in the same files. *Failure to do so will result in unreproducible
test results!*


#### Unit tests

*Document procedure for new pFUnit tests*

## Running

sbetr takes exactly one command line arguement, the name of the input
namelist file. The namelist input file specifies runtime configuration
and paths to other input data files. NOTE: paths are relative to the
directory where sbetr is executed!


    cd ${SBETR_ROOT}/example_input
    ${SBETR_ROOT}/local/bin/sbetr mock.namelist

The example is set with mock run that transport five tracers: N2, O2,
AR, CO2, CH4 and DOC


## Development

Key direcotries:

* 3rd-party - select 3rd-party sources that betr depends on.

* src - contains model code

  * src/betr - LM independent betr library

  * src/driver - standalone driver

  * src/stub\_clm, esmf\_wrf\_timemgr, src/shr - stub version of land
    model code needed to make the standalone LM interfaces compile.

  * Application - place to hold customized betr applications

* cmake - contains utilities for the configuration and build system.

* regression_tests - regressiont test input and baselines.


To configure a new bgc implementation, follow the example in
BGCReactionsMockRunType and add the new configuration name to
ApplicationsFactoryMod.
