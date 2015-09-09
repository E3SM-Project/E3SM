#! /bin/csh -f 

# directory in which glc is built
set glc_dir=$EXEROOT/glc

# directory in which glc obj files are built
set glc_obj_dir=$OBJROOT/glc/obj

# directory in which glimmer-cism library is created
set cism_libdir=$glc_dir/lib

# directory in which we can find source mods
set sourcemod_dir=$CASEROOT/SourceMods/src.cism

cd $glc_obj_dir

# ----------------------------------------------------------------------
# Create Filepath
# ----------------------------------------------------------------------
# The following just gives the filepath for the cesm-specific code:
# the glimmer-cism stuff is picked up by the cmake-based build
cat >! Filepath << EOF
$sourcemod_dir
$CODEROOT/glc/cism/drivers/cpl
$CODEROOT/glc/cism/source_glc
$CODEROOT/glc/cism/mpi
EOF

# ----------------------------------------------------------------------
# Set options to cmake
# ----------------------------------------------------------------------
# Note that some other generic CMAKE options are set in the Makefile
set cmake_opts=""
set cmake_opts="$cmake_opts -D CISM_COUPLED=ON"
set cmake_opts="$cmake_opts -D CISM_USE_MPI_WITH_SLAP=ON"
# CISM_USE_GPTL_INSTRUMENTATION is unnecessary (and possibly harmful)
# when built inside CESM; for CESM we instead use -DCCSMCOUPLED, which
# also gives us timing instrumentation
set cmake_opts="$cmake_opts -D CISM_USE_GPTL_INSTRUMENTATION=OFF"
set cmake_opts="$cmake_opts -D CISM_BINARY_DIR=$glc_dir"
set cmake_opts="$cmake_opts -D CMAKE_Fortran_MODULE_DIRECTORY=$glc_obj_dir"
set cmake_opts="$cmake_opts -D GLIMMER_NETCDF_DIR="\$"(NETCDF_PATH)"
set cmake_opts="$cmake_opts -D CISM_MPI_INC_DIR="\$"(INC_MPI)"
set cmake_opts="$cmake_opts -D GLIMMER_SOURCEMOD_DIR=$sourcemod_dir/glimmer-cism"
if ($CISM_USE_TRILINOS == 'TRUE') then
    set cmake_opts="$cmake_opts -D NO_TRILINOS=OFF"
    set cmake_opts="$cmake_opts -D CISM_MPI_MODE=ON"
    set cmake_opts="$cmake_opts -D CISM_SERIAL_MODE=OFF"
    set cmake_opts="$cmake_opts -D GLIMMER_TRILINOS_DIR="\$"(TRILINOS_PATH)"
else
    set cmake_opts="$cmake_opts -D NO_TRILINOS=ON"
    set cmake_opts="$cmake_opts -D CISM_MPI_MODE=OFF"
    set cmake_opts="$cmake_opts -D CISM_SERIAL_MODE=ON"
endif

# ----------------------------------------------------------------------
# Set mkDepends to append libglimmercismfortran.a to the end of each
# .o dependency line.
#
# Rationale: Some of the source files in the cesm-specific code depend
# on files included in this library. Ideally, we would be able to
# determine the actual dependencies, but that's not easy with the
# current tools and the fact that we build the glimmer-cism code using
# a different build system than the cesm-specific code. So for now, we
# just rebuild all the cesm-specific code whenever anything in the
# libglimmercismfortran.a library changes.
#
# WJS (3-6-13): I thought we would just need to include these options
# in the call to make the complib target. But for some reason that I
# can't determine, mkDepends is called when we make $glc_dir/Makefile,
# so we also need to include these options there.
# ----------------------------------------------------------------------

set mkdepends_opts="-d $cism_libdir/libglimmercismfortran.a"

# ----------------------------------------------------------------------
# create the glimmer-cism makefile by running cmake (done via a rule
# in the system-level makefile)
# ----------------------------------------------------------------------
$GMAKE $glc_dir/Makefile MODEL=cism USER_CMAKE_OPTS="$cmake_opts" USER_MKDEPENDS_OPTS="$mkdepends_opts" GLC_DIR=$glc_dir -f $CASETOOLS/Makefile || exit 1

# ----------------------------------------------------------------------
# create the glimmer-cism library (or libraries), using the makefile
# created by cmake
# ----------------------------------------------------------------------
pushd $glc_dir
$GMAKE -j $GMAKE_J || exit 2
popd

# ----------------------------------------------------------------------
# create the cesm-specific portion of the glc library using cesm's makefile
# ----------------------------------------------------------------------
$GMAKE complib -j $GMAKE_J MODEL=cism COMPLIB=$LIBROOT/libglc.a USER_MKDEPENDS_OPTS="$mkdepends_opts" GLC_DIR=$glc_dir -f $CASETOOLS/Makefile || exit 6

exit 0

