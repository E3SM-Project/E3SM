CMake_Fortran_utils
===================

CMake modules dealing with Fortran-specific issues and Fortran libraries

Currently, these modules should work with CMake version 2.8.8 and later
versions. Earlier CMake versions may work but are untested.

Below is a brief listing of modules. More detailed information on the
purpose and use of these modules can be found in comments at the top of
each file.

Find modules for specific libraries:

FindNETCDF

FindpFUnit

FindPnetcdf

Utility modules:

genf90_utils - Generate Fortran code from genf90.pl templates.

pFUnit_utils - Create executables using the pFUnit parser and driver.

Sourcelist_utils - Use source file lists defined over multiple directories.

Modules that are CESM-specific and/or incomplete:

CESM_utils - Handles a few options, and includes several other modules.

Compilers - Specify compiler-specific behavior, add build types for CESM.
