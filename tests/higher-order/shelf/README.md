=====
Shelf
=====

This directory contains Python scripts for running idealized ice shelf
experiments.


=========================
Confined Shelf Experiment
=========================

Note that this test case has been slightly altered from its previous form. 
That is, the setup is now consistent with tests 3 & 4 from the more "simple"
(i.e. not Ross) EISMINT-shelf test cases. More info on these can be found at:
<http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html>

You will need to either copy your executable into this directory, have the
executable in your path, or create a symbolic link to your executable, using:

```sh
ln -s ../../../builds/platform-compiler/cism_driver/cism_driver ./
```

To run the confined shelf experiment, enter the following on a terminal 
command line:

For a serial job: 

```sh
./runShelfConfined.py 
```

For a parallel job: 

```sh
./runShelfConfined.py -n N
```

where `N` is the number of processors to use.

Note that parallel support in the current `runShelfConfined.py` script attempts to 
determine how an mpi job is executed on your machine, but if you see an error
you may need to make minor alterations to the relevant lines of code that make the call
to MPI.

If there are problems with running in serial on a multicore machine, try:

```sh
./runShelfConfined.py -n 1
```

For details of options for running the script use:

```sh
./runShelfConfined.py --help
```

Unlike in previous releases, there is only one configuration file which uses
the Payne/Price higher-order dynamical core (the Pattyn/Bocek/Johnson higher-
order core, while still present in the code, is no longer supported).


The script performs the following three steps:

1. Create a netCDF input file for Glimmer.
2. Run Glimmer, creating a netCDF output file.

The netCDF files are written into the `output` subdirectory, which is controlled
by the `-o/--output-dir` option.

PREREQUISITES:
--------------

In order to use the `runShelfConfined.py` script, you will need to have Python and 
one of the following Python modules installed:

* Scientific.IO.NetCDF
* netCDF4
* pycdf

To view the results use ncview or another utility for viewing netCDF files.


=========================
Circular Shelf Experiment
=========================

To run the circular shelf experiment, enter the following on a terminal 
command line:

For a serial job: 

```sh
./runShelfCircular.py
```

For a parallel job: 

```sh
./runShelfCircular.py -n N
```

where `N` is the number of processors to use.



For details of options for running the script use:

```sh
./runShelfCircular.py --help
```

This script has a few more options:

`-a/--alpha` specifies that the ice thickness field will
have a conical top.  The default is that the ice thickness is constant (that is,
a flat top).

`-b/--beta` specifies that a Gaussian function
will be used for beta.  The default is that there is a small square region
in the center of the domain where beta is large; Beta is zero over the rest of
the domain.  There is an abrupt step from the large value (1.0e10) to zero (0.0)
when using the default.

`-d/--dirichlet` specifies that a Dirichlet 
boundary condition (velocity = zero) will be applied in a small square at the 
center of the domain.


The script performs the following three steps:

1. Create a netCDF input file for Glimmer.
2. Run Glimmer, creating a netCDF output file.

The netCDF files are written into the `output` subdirectory, which is controlled
by the `-o/--output-dir` option.


PREREQUISITES:
--------------

In order to use the `runShelfCircular.py` script, you will need to have Python and 
one of the following Python modules installed:

* Scientific.IO.NetCDF
* netCDF4
* pycdf

To view the results use ncview or another utility for viewing netCDF files.

