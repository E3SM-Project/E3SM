Dome test case 
==============

This directory contains a python script for running an experiment involving an
ice sheet with an parabolic "dome" shape on a flat base.  

You will first need to either copy your executable into this directory, 
have the executable in your path, or create a symbolic link to your executable,
using:

```sh
ln -s ../../../builds/platform-compiler/cism_driver/cism_driver ./
```

To run the experiment, enter the following on a terminal command line:

For a serial job: 

```sh 
./runDome.py 
```

For a parallel job: 

```sh 
./runDome.py -n N 
```

where N is the number of processors you'd like to use. For example:

```sh 
./runDome.py -n 8 
```

will run the test case on 8 processors.

Execute: 

```sh 
./runDome.py --help 
```

for a list of all options available.

If you want to run in parallel, the configuration file and number of processors
must be specified (but can be 1). If no parameters are specified, the code will
run in serial using the `dome.config` configuration file. 

The configuration file included is set to use the higher-order Glissade dycore.


PREREQUISITES: 
-------------- 

In order to use the `runDome.py` script, you will
need to have Python and one of the following Python modules installed:
* Scientific.IO.NetCDF netCDF4 pycdf

To view the results use ncview or another utility for viewing netCDF files.

--------------------------------------------------------------------------------

This directory also includes example files for how to use time-dependent forcing
with CISM. If the `dome-forcing.config` file is passed to the `runDome.py`
script using the `-c/--config` option (or, if any passed `*.config` file has a
[CF Forcing] section) it will generate a file called `dome.forcing.nc` that
contains time-varying fields for `acab`, `artm`, `uvel`, `vvel`, and
`kinbcmask`.

Then CISM will be run using the `dome-forcing.config` file which includes a
forcing section that uses this new `*.nc` file.  This can be done with, e.g.:

```sh
./runDome.py -c dome-forcing.config
```

After the simulation completes, you can view the output in ncview or some other
tool to confirm that the forcing fields have indeed changed over the simulation
and that they are being applied properly (e.g., the `thk` field is changing as
it ought to for the time-varying forcing applied).

