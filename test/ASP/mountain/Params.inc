PLEV  = 26
NP    = 4
PCOLS = 1
PCNST = 1
PNATS = 1
GRID_STAG =_NONSTAGGER
OMP   = _NO_OMP_THREADS
MPI   = _MPI
RESTART= _PRESTART

# Movie specifies the io method, options are
# _NETCDF   output a seperate netcdf file on the native grid for each mpi task 
# _PNETCDF  output a global netcdf file on the native grid using all tasks for IO   (slow - use _PIO)
# _PIO      output a global netcdf file on the native grid using pio to rearrange data prior to output
# _PIO_INTERP   output a global netcdf file on an interpolated grid using pio to rearrange data prior to output
#  _STUBS   no output method is linked
# see the file doc/README.io for details
MOVIE = _PIO_INTERP
#MOVIE = _PNETCDF
METIS  = _NOMETIS
TRACE  = _NOTRACE
