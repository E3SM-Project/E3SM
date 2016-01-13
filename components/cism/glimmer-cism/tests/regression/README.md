========================
Build and test structure
========================

The build and test structure (BATS) is primarily intended to allow users and developers to
quickly generate a set of regression tests for use with the Land Ice Validation
and Verification (LIVV) toolkit. 

BATS is a [Python 2.7](https://www.python.org/) module that is primarily
controlled by command line options. 
BATS requires:

Python Packages 
* [python-numpy](https://pypi.python.org/pypi/numpy) 
* [python-scipy](https://pypi.python.org/pypi/scipy) 
* [python-netCDF4](https://pypi.python.org/pypi/netCDF4) 
* [python-matplotlib](https://pypi.python.org/pypi/matplotlib)

External Packages
* [NetCDF 4.3.0+](http://www.unidata.ucar.edu/software/netcdf/)
* [NCO (NetCDF Operators) 4.4.0](http://nco.sourceforge.net/)
* [HDF5 1.8.6](https://www.hdfgroup.org/HDF5/)

If you have a working install of CISM, and you installed the suggested packages
in the [CISM](http://oceans11.lanl.gov/cism/documentation.html) users manual,
you'll likely already have everything you need. If you haven't previously built
CISM on your machine, we suggest following the installation instructions as they
are laid out in the users manual first.

(Note: on some High Performance Computing platforms, a
`setup_PLATFORM.bash` script has been provided which will attempt to load the
needed modules.)

If you are having any troubles with dependencies, open an issue on the 
[LIVVkit issue tracker](https://github.com/LIVVkit/LIVVkit/issues)!

===========
   Usage
===========

If the location of your CISM source
code is stored in the environment variable `CISM`, you can go to the BATS
location by:

```sh
cd $CISM/tests/regression/
```

There you will find the main BATS run script `build_and_test.py`. BATS is
primarily controlled via options specified at the command line.  To see the full
list of options, run:

```sh
./build_and_test.py -h
```

or 

```sh
python build_and_test.py -h
```

Each LIVVkit workflow, listed on the 
[LIVVkit usage](https://github.com/LIVVkit/LIVVkit/wiki/Usage) page, shows how
BATS is used to generate data for LIVVkit, and the commands used. We suggest
looking at the workflows first to get a general idea of how BATS will be used.  

====================
   How BATS works
====================

BATS works very similar to how you would build CISM and then run one or more CISM
tests. BATS will build a version of CISM, and then either run a set of
regression tests if you are using a personal computer (PC), or setup a series of
regression tests and generate a job submission script if you are using a high
performance computer (HPC). 

That is, if CISM is located in `$CISM`  on ORNL's HPC
platform titan, for example, and you want to run the `$CISM/higher-order/dome`
test using the gnu compiler, you would typically build CISM by:

```sh
cd $CISM/builds/titan-gnu/
source titan-gnu-cmake
make -j 8
```

which would build a parallel version of CISM into the
`$CISM/builds/titan-gnu/cism_driver/`
directory. Then you would go to the dome test directory, create a link to the
CISM driver, and run the test:

```sh
cd $CISM/tests/higher-order/dome
ls -s $CISM/builds/titan-gnu/cism_driver/cism_driver cism_driver
./runDome.py
```

You would have to repeat this process for each individual test you'd like
to run. BATS, however, automates these steps. On titan, running BATS like so

```sh
cd $CISM/tests/regression/
source setup_titan.bash
./build_and_test.py -p titan -c gnu -b ./build
```

will result in BATS generating all the CMake build file into a new directory
called `build`, and building CISM into the directory `build/cism_driver/`. BATS
will then run a set of CISM's higher-order tests:

* Dome (at a variety of resolutions and processor counts)
* Circular Shelf
* Confined Shelf
* ISMIP-HOM a and c (at 20 and 80 km resolutions)
* ISMIP-HOM f
* Stream

All of the files associated with each test will be output to a directory called
`reg_test/titan-gnu/` which has a directory structure that mirrors CISM's test
directory structure: 

```sh
 reg_test/
    └── PLATFORM-COMPILER/
        ├── CMakeCache.txt
        ├── higher-order/
        │    ├── dome/
        │    │   ├── dome.RESO.pPRC.*
        │    │   └── timing/
        │    │       └── dome-t?.RESO.pPRC.*
        │    ├── ismip-hom
        │    │   ├── ismip-hom-a.RESO.pPRC.*
        │    │   ├── ismip-hom-c.RESO.pPRC.*
        │    │   └── ismip-hom-f.0100.pPRC.*
        │    ├── shelf
        │    │   ├── shelf-circular.RESO.pPRC.*
        │    │   └── shelf-confined.RESO.pPRC.*
        │    └── stream
        │        └── stream.RESO.pPRC.*
    --------------------------------------------
        ├── Jobs/
        │    ├── platform_job.small
        │    ├── platform_job.small_timing_?
        │    ├── platform_job.large
        │    └── platform_job.large_timing_?
        ├── submit_all_jobs.bash
        └── clean_timing.bash
```

where `RESO` is a four-digit number indicating the model resolution (units are
test specific), and `pPRC` is an optional three-digit number, prefixed by a
`p`, indicating the number of processors used to run the model. Everything
below the dashed line will only appear on HPC systems. `submit_all_jobs.bash`
will submit all the jobs in the `jobs/` directory and `clean_timing.bash` will
clean out any `higher-order/*/timing/` directory such that only the timing
files remain (to be used once all jobs have finished). This
`reg_test/titan-gnu/` directory is formatted to be used with LIVVkit directly. 

BATS is designed to be flexible and work with any LIVVkit usage scenario. In
order to do that, BATS provides a number of options to configure which system
you are using, when/where/how CISM is built, the destination of the output
directory, and which tests are run. For more information on these topics, see:


* Detailed discussion of [BATS options](https://github.com/LIVVkit/LIVVkit/wiki/BATS-options)
* [Out-of-source builds](https://github.com/LIVVkit/LIVVkit/wiki/CISM-out-of-source-builds)
* [The reg_test directory](https://github.com/LIVVkit/LIVVkit/wiki/BATS-reg-test) structure
* [LIVVkit usage](https://github.com/LIVVkit/LIVVkit/wiki/Usage)

=============
   Authors
=============

Joseph H. Kennedy, ORNL 
    github : jhkennedy
