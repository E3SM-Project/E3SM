CUBE_TO_TARGET
==============

DESCRIPTION
-----------
This code performs rigorous remapping of topography variables on a cubed-sphere 
grid to any target grid. The code is documented in:
                                                                             
  Lauritzen, Nair and Ullrich, 2010, J. Comput. Phys.                        
 
USAGE
-----
Generating topography to be used in full atmosphere simulations involves 
multiple steps:

1. Interpolate high resolution topography on cubed-sphere to target grid
2. Apply dycore-specific smoothing to interpolated topography so that atmosphere 
   can run stably
3. Recompute subgrid surface roughness

This tool performs steps 1 and 3, but the smoothing needs to be handled by the 
dycore. (NOTE: there is an option to have this code perform the smoothing, but
this option is currently *untested* and *unsupported*, and it is preferred to
use the dycore to perform the smoothing separately)

To run the code, execute:
```
    ./cube_to_target <arguments>                                            
```
See below for description of arguments. To create topography for new grids, the
code will need to be run twice: once to interpolate the input topography to the
target grid, and then again after using the dycore to apply smoothing to the
interpolated topography to recompute the subgrid surface roughness.
                                                                             
REQUIRED ARGUMENTS                                                          
------------------
  --target-grid <filename>            Target grid descriptor in SCRIP format
  --input-topography <filename>       Input USGS topography on cube sphere
  --output-topography <filename>      Output topography on target grid       
                                                                             
OPTIONAL ARGUMENTS                                                          
------------------
  --smoothed-topography <filename>    Input smoothed topography (for surface 
                                      roughness calculation). If present,    
                                      output will contain estimate of subgrid
                                      surface roughness (SGH). Note that the                     
                                      variance in elevation from the 30s to
                                      3km grid (SGH30) is also downscaled,
                                      but does not depend on the smoothing.
                                                                             
AUTHOR
------
  Peter Hjort Lauritzen (pel@ucar.edu), AMP/CGD/NESL/NCAR

INSTALLING
----------
A Makefile is provided to build the code using the `make` utility. The Makefile
will grab the following environment variables if they exist:

  FC            Fortran compiler to use (i.e., ifort, gfortran, pgif90)
  LIB_NETCDF    Location of netCDF libraries (i.e., /usr/local/lib)
  INC_NETCDF    Location of netCDF include headers (i.e., /usr/local/include)

If these environment variables do not exist, some defaults are assumed, but it
is unlikely the build will succeed on your platform, especially if you are on
some kind of high performance computing environment. So, to build, setup the
environment variables and then run `make`. As an example, on an E3SM-supported
machine, it should be sufficient to run the following commands to build:
```
    ${E3SM_ROOT}/cime/tools/configure && source .env_mach_specific.sh
    export FC=ifort
    export LIB_NETCDF=${NETCDF_DIR}/lib
    export INC_NETCDF=${NETCDF_DIR}/include
    make
```
where `${E3SM_ROOT}` is set to the location of the E3SM source on your machine.
