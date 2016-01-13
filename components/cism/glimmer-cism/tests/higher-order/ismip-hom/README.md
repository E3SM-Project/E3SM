ISMIP-HOM Experiments
=====================

Ice Sheet Model Intercomparison Project for Higher-Order Models (ISMIP-HOM)
prescribes a set of experiments meant to test the implementation of higher-order
physics.  For more information, see
<http://homepages.ulb.ac.be/~fpattyn/ismip/>.

The python scripts provided here (`runISMIP_HOM.py` and `plotISMIP_HOM.py`,
referred to in the following as the ISMIP-HOM scripts) were created to run the
experiments using Glimmer/CISM and compare the results with results from other
models.

PREREQUISITES:
--------------
In order to use the ISMIP-HOM scripts, you will need to have Python and one of
the following Python modules installed:
 
 * Scientific.IO.NetCDF
 * netCDF4
 * pycdf

To plot the results (which is the usual way to view the output), the Python
matplotlib library also must be installed.

For comparison you need to download the ISMIP-HOM results from
<http://homepages.ulb.ac.be/~fpattyn/ismip/tc-2-95-2008-supplement.zip> and
unzip to a directory named `ismip_all`. (Make `ismip_all` a subdirectory of the
directory containing the ISMIP-HOM scripts.)

You will need to either copy your executable into this directory, have the 
executable in your path, or create a symbolic link to your executable, using:

```sh
ln -s ../../../builds/platform-compiler/cism_driver/cism_driver ./
```

RUNNING THE TESTS:
------------------
Running the tests is a three step process.  The last two steps invoke python
scripts by entering commands on the command line.  The working directory should
be the directory which contains the ISMIP-HOM scripts and the configuration
files copied in step one. Within the working directory, an `output` subdirectory
will be created that will contain all the output files.

 1. __Modify the default configuration file if desired:__

All tests use `ismip-hom.config` as a template for creating the `*.config` file
that is actually used to perform the model run.  If you wish to adjust any model
settings, you should adjust them in this file.  Most users will not have any
reason to modify settings.

 2. __Run Glimmer/CISM by invoking the `runISMIP_HOM.py` script:__

```sh
./runISMIP_HOM.py
```

Invoking the script in this way runs a default set of experiments and domain
sizes.  To specify experiments and/or domain sizes you include command line
arguments when you invoke the script.  For example,

```sh
./runISMIP_HOM.py --sizes 40 80 160
```

runs the default experiments with domain sizes 40, 80, and 160 km.

```sh
./runISMIP_HOM.py -r a c
```

runs experiments a and c using the default domain sizes.  Combine these as

```sh
./runISMIP_HOM.py -r a c --sizes 40 80 160
```

to run experiments a and c with domain sizes 40, 80, and 160 km.

CISM may not converge for every experiment with the default values of the grid
size and other parameters.  Some additional parameters can be changed using
command line arguments.  Run 

```sh
./runISMIP_HOM.py --help
```

to see all of the possible command line arguments. 

3. __Plot the results by invoking the `plotISMIP_HOM.py` script:__

```sh
./plotISMIP_HOM.py
```

This will automatically plot the newest set of tests that were run (all sizes
and experiments).

If you have multiple sets of tests (e.g., via the `-m` or `-n` option) in the
same output directory, `plotISMIP_HOM.py` will automatically plot the newest set
(see the OUTPUT FILES section below for a discussion on sets). To plot an older
set, pass any of the output (`*.out.nc`) files from the older set via the `-f`
option.  


By default,  if there is more than one domain size for an experiment within a
set, they are all plotted as small plots in one file.  You may override the
sizes, and plotting one domain size at a time results in larger plots; for example

```sh
./plotISMIP_HOM.py --sizes 40
```

creates two files, each containing one (large) plot.

The plots created are saved as `*.png` files in the `output` subdirectory of the
directory containing the `plotISMIP_HOM.py` script.  You can change the type of
file created by changing the value of the variable `plotType` in the
`plotISMIP_HOM.py` script.  If you prefer to have the plots displayed on the
screen instead of saved to file, set the variable `savePlotInFile` in
`plotISMIP_HOM.py` to `False` instead of `True`.  The variables `savePlotInFile`
and `plotType` are set near the top of the `plotISMIP_HOM.py` script. 

Output flags
------------
As with `runISMIP_HOM.py`, all of the possible command line arguments for
`plotISMIP_HOM.py` can be displayed using the command line argument `--help` (or
`-h`)

```sh
./plotISMIP_HOM.py --help
```

OUTPUT FILES:
-------------
When you run the ISMIP-HOM scripts, they create files in the `output`
subdirectory.  These files are described below.

The files will all follow this naming pattern:

```sh
ismip-hom-?[-MOD].RESO.[pPROC.].EXT
```

Where `?` is a stand-in for the experiment (e.g. `a`), `[-MOD]` is an optionally
user specified filename modifier, `RESO` is the size of the experiment, `[pPROC.]`
is the number of processors the test was run with (appears only if the `-n`
option was passed to `runISMIP_HOM.py`), and `EXT` is the file extension. A set
of experiment files consist of all files that match the pattern:

```sh
ismip-hom-?[-MOD].????.[pPROC.].EXT
```

where here the `?` is a POSIX metacharacter. That is, all experiments and sizes
are grouped into sets defined by `[-MOD]` and `[pPROC]` arguments (or lack
thereof). 

Files whose names end in `.config` are configuration files read by CISM.  They
are created by copying the base configuration files (`ismip-hom.config`) and
making changes based on the command line arguments passed to `runISMIP_HOM.py`.
Configuration files are text files and can be viewed using a text editor.

Files whose names end in `.nc` are netCDF files.  The contents of these files
can be examined using a tool such as ncview.  There are two different netCDF
files created each time CISM is run: 1. An input file created by
`runISMIP_HOM.py` which provides the ice thickness, bed topography, and
(sometimes) basal friction coefficient.  2. An output file created by CISM will
have a name ending in `.out.nc`.   

Files whose names end in `.txt` are output files written in a format used by all
models participating in the ISMIP-HOM experiment.  As their name implies these
are text files that can be viewed using a text editor.  They are written by
`plotISMIP_HOM.py` from information in the netCDF output files prior to
plotting.

Files whose names end in `.png` are Portable Network Graphics files written by
`plotISMIP_HOM.py` containing the plotted data.  This is standard file type that
can be viewed using image viewer software or a web browser.


TO DO / KNOWN ISSUES
--------------------
For ISMIP-HOM tests A and C, at the 20, 10, and 5  km wavelengths, convergence
may be very slow.

In general, the code will converge faster for fewer horizontal grid points but
the results may be less accurate (i.e. the match with the ISMIP-HOM mean may be
poor). The provided configuration files have reasonable default values for the
grid spacing in order to "pass" the ISMIP-HOM tests at *most* of the
wavelengths. In some cases these may need to be adjusted.

When multiple sizes are found by `plotISMIP_HOM.py`, and a single size is not
specified, the `f` plots will be duplicated within its plot that many times. 
