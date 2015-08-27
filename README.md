# ParallelIO

A high-level Parallel I/O Library for structured grid applications

## Dependencies

PIO can use NetCDF (version 4.3.3+) and/or PnetCDF (version 1.6.0+) for I/O.
Ideally, the NetCDF version should be built with MPI, which requires that it
be linked with an MPI-enabled version of HDF5.  Optionally, NetCDF can be 
built with DAP support, which introduces a dependency on CURL.  Additionally,
HDF5, itself, introduces dependencies on LIBZ and (optionally) SZIP.

## Configuring with CMake

To configure the build, PIO requires CMake version 2.8.12+.  The typical
configuration with CMake can be done as follows:

```
CC=mpicc FC=mpif90 cmake [-DOPTION1=value1 -DOPTION2=value2 ...] /path/to/pio/source
```

where `mpicc` and `mpif90` are the appropriate MPI-enabled compiler wrappers
for your system.

The `OPTIONS` section typically should consist of pointers to the install
locations for various dependencies, assuming these dependencies are not 
located in *canonical* search locations.  

For each dependency `XXX`, one can specify the location of its 
installation path with the CMake variable `XXX_PATH`.  If the `C` and
`Fortran` libraries for the dependency are installed in different locations
(such as can be done with NetCDF), then you can specify individually
`XXX_C_PATH` and `XXX_Fortran_PATH`.  Hence, you can specify the locations
of both NetCDF-C, NetCDF-Fortran, and PnetCDF with the CMake line:

```
CC=mpicc FC=mpif90 cmake -DNetCDF_C_PATH=/path/to/netcdf-c \
                         -DNetCDF_Fortran_PATH=/path/to/netcdf-fortran \
                         -DPnetCDF_PATH=/path/to/pnetcdf \
                         /path/to/pio/source
```

This works for the dependencies: `NetCDF`, `PnetCDF`, `HDF5`, `LIBZ`, `SZIP`.

## Building

Once you have configured PIO with CMake in a build directory.  From within
the build directory, build PIO with:

```
make
```

This will build the `pioc` and `piof` libraries.

## Testing

If you desire to do testing, you may build the test executables with:

```
make tests
```

Once the tests have been built, you may run tests with:

```
ctest
```

Alternatively, you may build the test executables and then run tests immediately with:

```
make check
```

which operates similarly to the `make check` Autotools target.

**NOTE:** It is important to note that these tests are designed to run in parallel.
If you are on one of the supported supercomputing platforms (i.e., NERSC, NWSC, ALCF, 
etc.), then the `ctest` command will assume that the tests will be run in an appropriately
configured and scheduled parallel job.  This can be done by requesting an interactive
session from the login nodes and then running `ctest` from within the interactive
terminal.  Alternatively, this can be done by running the `ctest` command from a
job submission script.  It is important to understand, however, that `ctest` itself
will preface all of the test executable commands with the appropriate `mpirun`/`mpiexec`/`runjob`/etc.
Hence, you should not further preface the `ctest` command with these MPI launchers.

