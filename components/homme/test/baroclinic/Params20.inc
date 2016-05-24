PLEV  = 20
NP    = 8
PCOLS = 1
PCNST = 1
PNATS = 1
GRID_STAG =_NONSTAGGER
OMP   = _OMP_THREADS
MPI   = _MPI
RESTART= _PRESTART
METIS  = _NOMETIS
TRACE  = _NOTRACE

# Movie specifies the io method, options are
# _PIO      output a global netcdf file on the native grid using pio to rearrange data prior to output
# _PIO_INTERP   output a global netcdf file on an interpolated grid using pio to rearrange data prior to output
# see the file doc/README.io for details
MOVIE = _PIO_INTERP
#MOVIE = _PIO




