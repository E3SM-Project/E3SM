======
stream
======
This directory contains Python scripts for running idealized ice stream
experiments.

This test case simulates flow over an idealized ice stream underlain by a
subglacial till with a known and specified yield stress distribution. For the
two distributions available in this test case, analytical solutions are
available from Raymond (2000) and Schoof (2006). Some additional discussion of
the Schoof test case can be found in Bueler and Brown (2009). 

For the Raymond test case, the yield stress within the ice stream is given a
uniform value below the driving stress, and outside of the ice stream it is
given a uniform value much higher than the driving stress (i.e., the yield
stress distribution is approximated by a "step" function). For the Schoof test
case, the till yield stress across the ice stream is given by a continuously
varying function. 

In both cases, the basal properties vary in the across-flow direction only and
are symmetric about the ice stream centerline.  As a result, the velocity
solutions are also uniform along flow and symmetric about the centerline.

You will need to either copy your executable into this directory, have the
executable in your path, or create a symbolic link to your executable, using:

```sh
ln -s ../../../builds/platform-compiler/cism_driver/cism_driver ./
```

To run the stream experiment, enter the following on a terminal command line:

For a serial job: 

```sh
./runStream.py 
```

For a parallel job: 

```sh
./runStream.py -n N
```

where `N` is the number of processors to use. 

Note that parallel support in the current `runStream.py` script attempts to
determine how an mpi job is executed on your machine, but if you see an error
you may need to make minor alterations to the relevant lines of the script that
make the call to MPI.

If there are problems with running in serial on a multicore machine, try:

```sh
./runStream.py -n 1
```

The default is to run the model for the Raymond configuration and compare output
to the Raymond analytic solution. This can be changed by simply editing the text
string in quotes in line 122 of the `runStream.py` script.

Depending on the horizontal and vertical resolution of the problem, you may need
to increase the maximum number of non-linear solver iterations above the default
value in order for the model to converge. If this is necessary, increase the
value of `glissade_maxiter` in the `stream.config` file.

If the number of grid cells in the across-flow direction is too small, one may
also experience problems and for this reason we recommend not decreasing this
number below the default value of 25. 

Additional information on command line options available when running the script
(e.g., for changing the horizontal and vertical resolution) can be found using:

```sh
./runStream.py --help
```

The script performs the following three steps:
1. Parses any command line options relative to the default values in `stream.config.in`
2. Creates a netCDF input file for CISM.
3. Runs CISM, and creates a netCDF output file.

The netCDF files are written to the `output` subdirectory, which is controlled
by the `-o/--output-dir` option.


PREREQUISITES:
--------------

In order to use the `confined-shelf.py` script, you will need to have Python and 
one of the following Python modules installed:

* Scientific.IO.NetCDF
* netCDF4
* pycdf

To view the results and compare model output to the analytical solutions use:

```sh
./plotStream.py
```

You can also use ncview or another utility for viewing netCDF files to view the
`*.out.nc` output file directly.


REFERENCES:
-----------

Raymond, C. F., 2000: Energy balance of ice streams. J Glaciol, 46, 665–674.

Schoof, C., 2006: A variational approach to ice stream flow. Journal of Fluid
Mechanics, 556, 227–251, doi:10.1017/S0022112006009591.

Bueler, E., and J. Brown, 2009: Shallow shelf approximation as a “sliding law”
in a thermomechanically coupled ice sheet model. J. Geophys. Res, 114, 1–21,
doi:10.1029/2008JF001179.
