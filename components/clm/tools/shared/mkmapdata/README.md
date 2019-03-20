# Overview
The routines in this directory create a mapping dataset from
SCRIP grid files to map from one grid to another. These mapping files
are used by either CLM or `mksurfdata_map` to regrid from one resolution
to another.

The script uses ESMF and requires that ESMF be built and the path
for ESMF binary files (using the program ESMF_RegridWeightGen) 
be given as input to the script. You will probably need to two versions,
one with mpiuni (no MPI) and one with mpi. Both versions also need to be built 
with NetCDF rather than the default IO version.

Details about the various script options can be found by invoking the "help" 
option of the script:
```
   ./mkmapdata.sh -help
```

# Prerequisites
The following tasks should only have to be done once for each new machine you 
want to run the code on.

## Export the input SCRIP grid files for the resolutions you'll need. Most of 
these files are on the Subversion inputdata server at

    https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/mappingdata/grids/

Supported machines also have a copy in the CESM DIN_LOC_ROOT location for that machine.

## Obtain and build the versions of ESMF required for this script

The version needs to support ESMF_RegridWeightGen and support the
options passed to it in the mkmapdata.sh script. As such it needs
to be built with NetCDF, support UGRID format, and support the 
`--netcdf4` and `--64bit_offset` options. ESMF can be obtained from

   http://www.earthsystemmodeling.org/

The version of NetCDF used with ESMF needs to be version 4.1 or higher
and compiled with the NetCDF4 file format enabled (with HDF5 compression).
That will enable the `--netcdf4` and `--64bit_offset` options to be used.

# Running the code
The following steps provide a method to create the executable and generate the
grid map dataset.

1. cd to this directory 

2. Create map dataset(s). This is done using the `mkmapdata.sh` script. Options and 
   environment variables can be shown by passing the `--help` argument to the script.
   Basic usage involves either passing the `-r` argument with a supported resolution,
   or passing the `-f` argument with a user-defined SCRIP-formatted grid file to specify
   the target grid.

   Example for standard resolutions:
   ```
        ./mkmapdata.sh -r 10x15
   ```
   Example for non-standard resolutions where you provide an input SCRIP grid file:
   ```     
       ./mkmapdata.sh -f <SCRIP_gridfile>
   ```
   
   Alternatively, the `regridbatch.sh` can be run, which runs `mkmapdata.sh` for a bunch of 
   different resolutions. The `mknoocnmap.pl` script can also be run to create a single-point/regional
   map for an area without ocean:
   ```
      ./mknoocnmap.pl -help      # for help on this script
   ```
   
3. move (and rename if appropriate) generated map datasets
   to `$DIN_LOC_ROOT/lnd/clm/mappingdata/maps`, etc.


Important files:
```
regridbatch.sh  Script to run mkmapdata.sh for many resolutions 
mvNimport.sh    Script to copy and import mapping files in for many resolutions
mkmapdata.sh    Script to create mapping datasets for a given resolution
mknoocnmap.pl   Script to create unity mapping dataset for single-point or regional studies over land-only (no ocean).
mkunitymap.ncl  NCL script to create a unity map -- ran by above script
```
