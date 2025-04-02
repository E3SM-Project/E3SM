# Third-party Libraries Required by EAMxx

- First, we note that if you are running EAMxx on a supported machine,
then you can expect that the third-party libraries (TPLs) are
installed, properly configured, and recognized by the EAMxx build system.
- However, if you choose to run and/or develop EAMxx on a local machine,
workstation, or unsupported cluster, then you will likely need to build some
number of the TPLs manually.
      - **Note:** The exception to this is if you are running on a machine with
      a module system that already provides them--e.g., research institutions
      that run E3SM regularly.

## List of TPLs

- **Note:** For this section, we will assume you plan to build EAMxx on a fresh
Linux[^only-linux] system with the standard development tools available.
      - E.g., a package manager (`apt`, `yum`), `git`, `make`, `cmake`, etc.

??? Abstract "TPLs You Will (almost certainly) Need to Install or Build Manually"

    - [OpenMPI](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)
        - Only required for builds employing multi-node parallelism.
        - Required for parallel configurations of the netCDF libraries.
        - Likely to be provided by your package manager or module system.
    - [HDF5](https://www.hdfgroup.org/download-hdf5/)
        - High-performance data-wrangling library.
        - The link above provides `.deb` and `.rpm` binary packages, though
        the page also provides a link to build from source.
    - [netCDF-C](https://github.com/Unidata/netcdf-c)
        - Used for I/O of simulation data in a format containing descriptive metadata.
    - [netCDF-fortran](https://github.com/Unidata/netcdf-fortran)
    - [PnetCDF](https://parallel-netcdf.github.io)
        - Only required for builds employing multi-node parallelism.
    
    **Build Tips**
    
    - As of 2024, the author found much less resistance configuring the build
    for the netCDF libraries using Autoconf, rather than CMake.

??? List "Auxiliary TPLs Your System May Already Have Installed"

      The following list of TPLs may be pre-installed on some systems, and are
      helpful, if not required, for building the netCDF libraries.
      However, in the author's experience they did need to be built on a fresh
      installation of Red Hat Enterprise Linux 8.

    **Required**

    - [LAPACK](https://www.netlib.org/lapack/) or [BLAS](https://www.netlib.org/blas/)
        - Alternatively, [OpenBLAS](http://www.openmathlib.org/OpenBLAS/),
        which is unlikely to be pre-installed.

    **Likely Required**

    - [Autoconf](https://www.gnu.org/software/autoconf/)
    - [Libtool](https://www.gnu.org/software/libtool/)

    **Maybe Required**[^computers-amirite]

    - [help2man](https://www.gnu.org/software/help2man/)

## Help & Hints

Building third-party libraries can be a pain.
To help here are some reasonable example scripts to start from that should
help getting the ***TPLs You Will Need to Build*** compiled and working.

??? Example "Basic Configuration Script Examples"

    !!! Warning "Disclaimer"

        We provide these scripts as demonstrative examples, only.
        They will almost certainly not build your libraries successfully if they
        are copy/paste-ed as they appear here.
        That being said, every hardware/software environment is different and
        will display different and sometimes strange behaviors.
        So, the best advice we can provide is to begin with something
        resembling what is below and respond to the compiler messages
        by modifying: environment variables, `PATH` or `[LD_]LIBRARY_PATH`
        entries, configuration arguments, and compiler arguments as required.
        As long as you do ***all of that*** perfectly, you should be just fine.
        :sunglasses:

    ### Assumptions

    - Source code is downloaded to a directory we assign to the variable `${TPL_ROOT}`.
    - We build the libraries in a separate ***build*** and ***install***
    directory in the root directory of the source code.
        - E.g., `${TPL_ROOT}/<lib-name>/{build/,install/}`
        - The bash scripts below are placed in the `build/` directory.
            - Hence, note that when executing the configure script from the
            `build/` directory, the final `../configure [...]` command in the
            script (that calls autoconf) is running the `configure` command in
            the source directory for that TPL--e.g., `${TPL_ROOT}/<lib-name>/`
    - `libtool` is installed and on your `PATH`.
        - I.e., `which libtool` returns the location of the executable.

    ??? Abstract "`configure-hdf5.sh`"

        ```{.shell .copy}
        #!/bin/bash

        # add MPI C/C++ compilers to the proper environment variables
        mpi_path="${TPL_ROOT}/openmpi/install/bin/"
        export CC="${mpi_path}/mpicc"
        export CXX="${mpi_path}/mpicxx"

        # configure HDF5 build
        ../configure \
          --prefix=${TPL_ROOT}/hdf5/install \
          --enable-shared \
          --enable-hl \
          --enable-parallel \
          --with-zlib=/usr/include,/usr/lib
        ```

    ??? Abstract "`configure-netcdf-c.sh`"

        ```{.shell .copy}
        #!/bin/bash

        # add MPI C/C++ to the proper environment variables
        mpi_path="${TPL_ROOT}/openmpi/install/bin"
        export CC="${mpi_path}/mpicc"
        export CXX="${mpi_path}/mpicxx"

        # historically, netcdf can be a little picky to build, so it can help to
        # make ABSOLUTELY CERTAIN that the build system knows where HDF5 is
        # via the "sledgehammer" approach
        hdf5_path="${TPL_ROOT}/hdf5/install"
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${hdf5_path}/lib"
        export LIBRARY_PATH=${LIBRARY_PATH}:${hdf5_path}/lib
        export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${hdf5_path}/lib/pkgconfig
        export CPPFLAGS="-I${hdf5_path}/include"

        ../configure --prefix=${TPL_ROOT}/netcdf-c/install
        ```

    ??? Abstract "`configure-netcdf-fortran.sh`"

        ```{.shell .copy}
        #! /bin/bash

        # add MPI Fortran compiler to the proper environment variables
        mpi_path="${TPL_ROOT}/openmpi/install/bin"
        export FC="${mpi_path}/mpif90"

        # local variables, to be used below
        export NetCDF_C_ROOT=${TPL_ROOT}/netcdf-c/install
        export NETCDF_C_LIBRARY=${NetCDF_C_ROOT}/lib
        export NETCDF_C_INCLUDE_DIR=${NetCDF_C_ROOT}/include

        # add netcdf-C bits and pieces to the proper paths
        export LD_LIBRARY_PATH="${NETCDF_C_LIBRARY}:${LD_LIBRARY_PATH}"
        export LIBRARY_PATH="${NETCDF_C_LIBRARY}:${LIBRARY_PATH}"
        export PATH=${NetCDF_C_ROOT}/bin:${PATH}

        # add netcdf-C include and library dir to the requisite compiler flags
        export CPPFLAGS="-I${NETCDF_C_INCLUDE_DIR}"
        export LDFLAGS="-L${NETCDF_C_LIBRARY}"

        # while we're at it, tell it where to find hdf5, too
        export HDF5_ROOT=${TPL_ROOT}/hdf5/install

        ../configure --prefix=${TPL_ROOT}/netcdf-fortran/install
        ```

    ??? Abstract "`configure-pnetcdf.sh`"

        ```{.shell .copy}
        #!/bin/bash

        mpi_path="${TPL_ROOT}/openmpi/install/bin"
        export CC="${mpi_path}/mpicc"
        export CXX="${mpi_path}/mpicxx"
        export FC="${mpi_path}/mpif90"

        ../configure \
          --prefix=${TPL_ROOT}/pnetcdf/install \
          --enable-shared \
          MPICC=mpicc \
          MPICXX=mpicxx \
          MPIF77=mpifort \
          MPIF90=mpifort
        ```

<!-- ======================================================================= -->

[^only-linux]: We only support build instructions for Linux, though some individuals have previously built EAMxx on MacOS with success.
We would not recommend using EAMxx on Windows, other than via WSL.
However, if you do end up building and running on MacOS or Windows,
please let the developers know or submit a PR!
<!-- markdownlint-disable-next-line MD053 -->
[^computers-amirite]: Computers are weird. While it's by no means required from
a functionality standpoint, a recently-configured RHEL8 system seemed to be
convinced it needed help2man for some reason.
