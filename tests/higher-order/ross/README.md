Ross Ice Shelf Experiment
=========================

This experiment was designed to model the Ross Ice Shelf off Antarctica.  For
information about the experiment and its results see:
<http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html>

The python scripts provided here (`runRoss.py` and `plotRoss.py`, referred to
in the following as the Ross scripts) were created to run the experiment using
CISM and compare the results with measured velocities.

PREREQUISITES:
--------------

In order to use the Ross scripts, you will need to have Python and one of the
following Python modules installed:
* Scientific.IO.NetCDF
* netCDF4
* pycdf

To plot the results (which is the usual way to view the output), the Python
matplotlib library also must be installed.

You will need to either copy your executable into this directory, have the
executable in your path, or create a symbolic link to your executable, using:

```sh
ln -s ../../../builds/platform-compiler/cism_driver/cism_driver ./
```

RUNNING THE EXPERIMENT:
-----------------------

Running the experiment is a two step process.  Each step invokes a python
script by entering a command on the command line.  The working directory should
be the directory which contains the Ross scripts.  This directory should also
contain the subdirectory `data`, which contains the input data for the
experiment.

1. __Run CISM by invoking the `runRoss.py` script:__

```sh
./runRoss.py
```

This command will build the necessary input files in the `./output`
subdirectory and run the model.  

For details of how to use the script to also run CISM, see:

```sh
./runRoss.py --help
```

When CISM is run, it will be use model setting specified in `ross.config` (by
default).


2. __Plot the results by invoking the plotRoss.py script:__

```sh
./plotRoss.py
```

The `plotRoss.py` script can also be invoked using command line options to plot
filled contours, change the color map, and mask the grounded region.  Example:

```sh
./plotRoss.py --ncontours=14 --vmax=1400 --cmap=gist_ncar --mask
```

Short forms of the options are also available.  Example:

```
./plotRoss.py -n 14 -v 1400 -c gist_ncar -m
```

OUTPUT FILES:
-------------

When you run `runRoss.py`, it creates two NetCDF files in the `output`
subdirectory.  These files can be examined using a tool such as ncview. Their
contents are described below.

1. `ross.nc` is an input file created by `runRoss.py` which provides input 
information to CISM.
2. `ross.out.nc` is an output file created by CISM.

