# rfmip-rrtmgp
This directory contains programs and support infrastructure for running
the [RTE+RRTMGP](https://github.com/RobertPincus/rte-rrtmgp) radiation parameterization for the
[RFMIP](https://www.earthsystemcog.org/projects/rfmip/) cases.

To run the codes you will need to
1. Download and build the RTE+RRTMGP libraries
2. Edit the Makefile in this directory to point to the RTE+RRTMGP installation (Makefile macro RRTMGP_DIR) as well as the location of the netCDF C and Fortran libraries and module files (macros NCHOME and NFHOME).
3. `make` in this directory. (You'll probably want to use optimizing compiler flags.)
4. Copy the absorption coefficient data files from $RRTMGP_DIR/data; obtain the input file from following links from <https://www.earthsystemcog.org/projects/rfmip/resources/>; make template output files using the script at <https://github.com/RobertPincus/RFMIP-IRF-Scripts>.

The executables are intended to run without further changes.  
