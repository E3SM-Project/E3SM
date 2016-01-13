=============
Glint example
=============

INPUT DATA:
-----------

To run this glint-example test case, you will need to download the following
climate files:

* `ncep-doe_6h_climate.64x32.nc` (precip and 2-m temp on a global grid)

* `orog.igcmgrid.64x32.nc` (surface elevation on a global grid)

You will also need a 20-km Greenland input file:

* `gland20.input.nc`

These three netCDF files are part of a tar file, `glint-example.1.0.0.tar.gz`,
that can be downloaded from the CISM website:

 <http://oceans11.lanl.gov/cism/data/glint-example.1.0.0.tar.gz>

These files are also available in the CESM inputdata repository (account
required), here:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/gland20.input.nc>

and:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/glint_example/>

To run a 5-km Greenland test case, you can download one of the files here:

<http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland>

The default is `Greenland_5km_v1.1.nc`. 

A 5-km file is also available in the CESM inputdata repo:

<https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/glc/cism/gland5.input.nc>

The 5-km test case and 20-km test case use the same climate files.


BUILDING AND RUNNING THE TEST CASES:
------------------------------------

glint_example is run using the `cism_driver` executable with two input
arguments, e.g.:

```sh
./cism_driver ice_sheet.config climate.config
```

First, copy or link the executable to this test directory:

```sh
ln -s ../../builds/mac-gnu/cism_driver/cism_driver
```

As of October 2014, the repo contains three Greenland config files:
1. `greenland_20km.config.pdd`
2. `greenland_5km.config.pdd`
3. `greenland_20km.config.smb`

These files have settings that are physically reasonable but not scientifically
validated.

Files 1 and 2 are basically identical except for the grid resolution and
timestep.  Both these files use the native Glint positive-degree-day (PDD)
scheme to force the Glide dycore, which uses the shallow-ice approximation.

Files 1 and 3 are identical except that file 3 has acab_mode = 0, which means
that Glint uses an externally computed surface mass balance (SMB) instead of
computing the SMB with a PDD scheme.  Normally the SMB would be computed and
passed in by a climate model such as CESM.  But for the standalone code, the
SMB is estimated crudely based on the input temperature and precip fields.
This allows us to test the Glint interfaces that would be called by a climate
model, without actually coupling to a climate model.

There are two glint_example config files in the repo:
1. `glint_example.config.pdd`
2. `glint_example.config.smb`

These files are set up to use the global NCEP climatological data for
downscaling to the Greenland grid. These files also set the length of the model
run (by default, we run a 10-year smoke test).  The only difference is that
file 1 should be run with the PDD scheme and file 2 with the SMB scheme.

The three standard diagnostic tests in this directory can be launched as
follows:

```sh
./cism_driver greenland_20km.config.pdd glint_example.config.pdd
./cism_driver greenland_5km.config.pdd glint_example.config.pdd 
./cism_driver greenland_20km.config.smb glint_example.config.smb
```

These are currently the only standard tests that exercise the Glint part of the
code, which is required for coupling to a climate model.

Diagnostic output is written to a file called `greenland_20km.config.smb.log` or
something similar.  NetCDF output is written to various files at a frequency
set in the Greenland config file.



