# This file is intended to be included by build_model function. Any changes
# to CMAKE variables intended to impact the CMakeLists.txt file that called
# build_model must use PARENT_SCOPE.

# Add INCROOT to path for Depends and Include
set(MINCROOT "")
if (INCROOT)
  list(APPEND CPP_DIRS ${INCROOT})
  set(MINCROOT ${INCROOT})
endif ()

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

# Determine whether any C++ code will be included in the build;
# currently, C++ code is included if and only if we're linking to the
# trilinos library or the Albany library.
set(USE_CXX FALSE)
if (USE_TRILINOS OR USE_ALBANY OR USE_KOKKOS)
  set(USE_CXX TRUE)
endif()

if (NOT MOD_SUFFIX)
  set(MOD_SUFFIX "mod")
endif()

# Look for -crm samxx in the CAM_CONFIG_OPTS CIME variable
# If it's found, then enable USE_SAMXX
string(FIND "${CAM_CONFIG_OPTS}" "-crm samxx" HAS_SAMXX)
if (NOT HAS_SAMXX EQUAL -1)
  # The following is for the SAMXX code:
  set(USE_SAMXX TRUE)
endif()

string(FIND "${CAM_CONFIG_OPTS}" "-rrtmgpxx" HAS_RRTMGPXX)
if (NOT HAS_RRTMGPXX EQUAL -1)
  # The following is for the RRTMGPXX code:
  set(USE_RRTMGPXX TRUE)
endif()

# If samxx or rrtmgpxx is being used, then YAKL must be used as well
if (USE_SAMXX OR USE_RRTMGPXX)
    set(USE_YAKL TRUE)
else()
    set(USE_YAKL FALSE)
endif()

# If YAKL is being used, then we need to enable USE_CXX
if (${USE_YAKL})
  set(USE_CXX TRUE)
endif()

#===============================================================================
# set CPP options (must use this before any flags or cflags settings)
#===============================================================================
include(${CASEROOT}/Macros.cmake)

set(CPPDEFS "${CPPDEFS} ${USER_CPPDEFS} -D${OS}")

# SLIBS comes from Macros, so this append must come after Macros are included
if (USE_FMS)
  set(SLIBS "${SLIBS} -lfms")
endif()

if (DEBUG)
  # e3sm still has components that cannot build with -DDEBUG
  if (CIME_MODEL STREQUAL "cesm")
    set(CPPDEFS "${CPPDEFS} -DDEBUG")
  endif()
else()
  set(CPPDEFS "${CPPDEFS} -DNDEBUG")
endif()

if (USE_ESMF_LIB)
  set(CPPDEFS "${CPPDEFS} -DUSE_ESMF_LIB")
endif()

if (COMP_INTERFACE STREQUAL "nuopc")
  set(CPPDEFS "${CPPDEFS} -DNUOPC_INTERFACE")
else()
  set(CPPDEFS "${CPPDEFS} -DMCT_INTERFACE")
endif()

if (COMPARE_TO_NUOPC)
  set(CPPDEFS "${CPPDEFS} -DCOMPARE_TO_NUOPC")
endif()

if (MPILIB STREQUAL mpi-serial)
  set(CPPDEFS "${CPPDEFS} -DNO_MPI2")
else()
  set(CPPDEFS "${CPPDEFS} -DHAVE_MPI")
endif()

if (PIO_VERSION STREQUAL "1")
  set(CPPDEFS "${CPPDEFS} -DPIO1")
else()
  set(USE_CXX TRUE)
endif()

if (USE_CXX AND NOT SUPPORTS_CXX)
  message(FATAL_ERROR "Fatal attempt to include C++ code on a compiler/machine combo that has not been set up to support C++")
endif()

# Not clear how to escape commas for libraries with their own configure
# script, and they don't need this defined anyway, so leave this out of
# FPPDEFS.
if (HAS_F2008_CONTIGUOUS)
  if (CPRE)
    set(CONTIGUOUS_FLAG "${CPRE}USE_CONTIGUOUS=contiguous,")
  else()
    set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS=contiguous,")
  endif()
else()
  if (CPRE)
    set(CONTIGUOUS_FLAG "${CPRE}USE_CONTIGUOUS=")
  else()
    set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS=")
  endif()
endif()

# JGF TODO: replace with findnetcdf
if (NETCDF_C_PATH)
  if (NOT NETCDF_FORTRAN_PATH)
    message(FATAL_ERROR "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
  endif()
  set(NETCDF_SEPARATE TRUE)
  if (NOT INC_NETCDF_C)
    set(INC_NETCDF_C ${NETCDF_C_PATH}/include)
  endif()
  if (NOT INC_NETCDF_FORTRAN)
    set(INC_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/include)
  endif()
  if (NOT LIB_NETCDF_C)
    if (EXISTS ${NETCDF_C_PATH}/lib)
      set(LIB_NETCDF_C ${NETCDF_C_PATH}/lib)
    elseif (EXISTS ${NETCDF_C_PATH}/lib64)
      set(LIB_NETCDF_C ${NETCDF_C_PATH}/lib64)
    else()
      message(FATAL_ERROR "NETCDF_C_PATH does not contain a lib or lib64 directory")
    endif()
  endif()
  if (NOT LIB_NETCDF_FORTRAN)
    if(EXISTS ${NETCDF_FORTRAN_PATH}/lib)
      set(LIB_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/lib)
    elseif(EXISTS ${NETCDF_FORTRAN_PATH}/lib64)
      set(LIB_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/lib64)
    else()
      message(FATAL_ERROR "NETCDF_FORTRAN_PATH does not contain a lib or lib64 directory")
    endif()
  endif()
elseif (NETCDF_FORTRAN_PATH)
  message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
elseif (NETCDF_PATH)
  set(NETCDF_SEPARATE FALSE)
  if (NOT INC_NETCDF)
    set(INC_NETCDF ${NETCDF_PATH}/include)
  endif()
  if (NOT LIB_NETCDF)
    if (EXISTS ${NETCDF_PATH}/lib)
      set(LIB_NETCDF ${NETCDF_PATH}/lib)
    elseif(EXISTS ${NETCDF_PATH}/lib64)
      set(LIB_NETCDF ${NETCDF_PATH}/lib64)
    else()
      message(FATAL_ERROR "NETCDF_PATH does not contain a lib or lib64 directory")
    endif()
  endif()
else()
  message(FATAL_ERROR "NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml")
endif()

if (MPILIB STREQUAL mpi-serial)
  if (PNETCDF_PATH)
    unset(PNETCDF_PATH)
  endif()
else()
  if (PNETCDF_PATH)
    if (NOT INC_PNETCDF)
      set(INC_PNETCDF ${PNETCDF_PATH}/include)
    endif()
    if (NOT LIB_PNETCDF)
      set(LIB_PNETCDF ${PNETCDF_PATH}/lib)
    endif()
  endif()
endif()

# Set PETSc info if it is being used
if (USE_PETSC)
  if (PETSC_PATH)
    if (NOT INC_PETSC)
      set(INC_PETSC ${PETSC_PATH}/include)
    endif()
    if (NOT LIB_PETSC)
      set(LIB_PETSC ${PETSC_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "PETSC_PATH must be defined when USE_PETSC is TRUE")
  endif()

  # Get the "PETSC_LIB" list an env var
  set(PETSC_DIR ${PETSC_PATH})
  find_package(PETSc)
  set(PETSC_LIB ${PETSC_LIBRARIES})
endif()

if (USE_TRILINOS)
  if (TRILINOS_PATH)
    if (NOT INC_TRILINOS)
      set(INC_TRILINOS ${TRILINOS_PATH}/include)
    endif()
    if (NOT LIB_TRILINOS)
      set(LIB_TRILINOS ${TRILINOS_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "TRILINOS_PATH must be defined when USE_TRILINOS is TRUE")
  endif()

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TRILINOS_PATH} PARENT_SCOPE)
  find_package(Trilinos)
endif()

if (USE_ALBANY)
  if (ALBANY_PATH)
    if (NOT INC_ALBANY)
      set(INC_ALBANY ${ALBANY_PATH}/include)
    endif()
    if (NOT LIB_ALBANY)
      set(LIB_ALBANY ${ALBANY_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "ALBANY_PATH must be defined when USE_ALBANY is TRUE")
  endif()

  # get the "ALBANY_LINK_LIBS" list as an env var
  file(READ ${ALBANY_PATH}/export_albany.in ALBANY_OUTPUT)
  string(REPLACE "ALBANY_LINK_LIBS=" "" ALBANY_LINK_LIBS "${ALBANY_OUTPUT}")
endif()

if (USE_KOKKOS)
  # LB 01/23
  # CMake's find_package, when used with multiple PATHS and PATH_SUFFIXES,
  # follows the following rule when looking in paths:
  #  1. look in all the path suffixes of the first PATHS entry, in the order provided
  #  2. look in the first PATH provided
  #  3. repeat 1-2 with the following entry of PATHS
  #  4. look in cmake/system default paths (unless told not to).
  # So the following cmd will fist look in the KOKKOS_PATH folder and subfolders,
  # if KOKKOS_PATH is non-empty. Then will proceed to look in the lib, lib/cmake,
  # and lib64/cmake subfolders of the INSTALL_SHAREDPATH. If all of these fail,
  # it will look in INSTALL_SHAREDPATH.

  if (KOKKOS_PATH)
    set (PATHS ${KOKKOS_PATH} ${INSTALL_SHAREDPATH})
  elseif(DEFINED ENV{KOKKOS_PATH})
    set (PATHS $ENV{KOKKOS_PATH} ${INSTALL_SHAREDPATH})
  else()
    set (PATHS ${INSTALL_SHAREDPATH})
  endif()
  find_package(Kokkos REQUIRED
               PATHS ${PATHS}
               PATH_SUFFIXES lib lib/cmake lib64/cmake
               NO_DEFAULT_PATH)
endif()

# JGF: No one seems to be using this
# if (USE_MOAB)
#   if (MOAB_PATH)
#     set(CPPDEFS "${CPPDEFS} -DHAVE_MOAB")
#     if (NOT INC_MOAB)
#       set(INC_MOAB ${MOAB_PATH}/include)
#     endif()
#     if (NOT LIB_MOAB)
#       set(LIB_MOAB ${MOAB_PATH}/lib)
#     endif()
#   else()
#     message(FATAL_ERROR "MOAB_PATH must be defined when USE_MOAB is TRUE")
#   endif()

#   # # get the "IMESH_LIBS" list as an env var
#   #include $(LIB_MOAB)/iMesh-Defs.inc
# endif()

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)
string(FIND "${CPPDEFS}" "-DLINUX" HAS_DLINUX)
string(FIND "${CPPDEFS}" "DBG" HAS_DBG)
string(FIND "${CPPDEFS}" "Darwin" HAS_DARWIN)
if (NOT HAS_DLINUX EQUAL -1 AND HAS_DBG EQUAL -1 AND HAS_DARWIN EQUAL -1)
  set(CPPDEFS "${CPPDEFS} -DHAVE_SLASHPROC")
endif()

# Atleast on Titan+cray mpi, MPI_Irsends() are buggy, causing hangs during I/O
# Force PIO to use MPI_Isends instead of the default, MPI_Irsends
if (PIO_VERSION STREQUAL 2)
  set(EXTRA_PIO_CPPDEFS "-DUSE_MPI_ISEND_FOR_FC")
else()
  set(EXTRA_PIO_CPPDEFS "-D_NO_MPI_RSEND")
endif()

if (LIB_PNETCDF)
  set(CPPDEFS "${CPPDEFS} -D_PNETCDF")
  set(SLIBS "${SLIBS} -L${LIB_PNETCDF} -lpnetcdf")
endif()

# Set esmf.mk location with ESMF_LIBDIR having precedent over ESMFMKFILE
set(CIME_ESMFMKFILE "undefined_ESMFMKFILE")
if (ESMFMKFILE)
  set(CIME_ESMFMKFILE ${ESMFMKFILE})
endif()
if (ESMF_LIBDIR)
  set(CIME_ESMFMKFILE ${ESMF_LIBDIR}/esmf.mk)
endif()

# For compiling and linking with external ESMF.
# If linking to external ESMF library then include esmf.mk
# ESMF_F90COMPILEPATHS
# ESMF_F90LINKPATHS
# ESMF_F90LINKRPATHS
# ESMF_F90ESMFLINKLIBS
if (USE_ESMF_LIB)
  # include(${CIME_ESMFMKFILE}) # JGF SKIPPING FOR NOW
  # Will need something like 'make -f esmf.mk  -p 2> /dev/null | grep ESMF_F90COMPILEPATHS'
  #set(CPPDEFS "${CPPDEFS} -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR} -DESMF_VERSION_MINOR=${ESMF_VERSION_MINOR}")
  #set(FFLAGS "${FFLAGS} ${ESMF_F90COMPILEPATHS}")
  #set(SLIBS "${SLIBS} ${ESMF_F90LINKPATHS} ${ESMF_F90LINKRPATHS} ${ESMF_F90ESMFLINKLIBS}")
  message(FATAL_ERROR "ESMF not supported in CMake yet")
endif()

#===============================================================================
# User-specified INCLDIR
#===============================================================================

set(INCLDIR ".")
if (USER_INCLDIR)
  list(APPEND INCLDIR "${USER_INCLDIR}")
endif()

#===============================================================================
# Set compilers
#===============================================================================

if (MPILIB STREQUAL "mpi-serial")
  set(CC ${SCC})
  set(FC ${SFC})
  set(CXX ${SCXX})
  set(MPIFC ${SFC})
  set(MPICC ${SCC})
  set(MPICXX ${SCXX})
else()
  set(CC ${MPICC})
  set(FC ${MPIFC})
  set(CXX ${MPICXX})
  if (MPI_PATH)
    set(INC_MPI ${MPI_PATH}/include)
    set(LIB_MPI ${MPI_PATH}/lib)
  endif()
endif()
set(CSM_SHR_INCLUDE ${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include)

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
list(APPEND INCLDIR "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include")

if (NOT NETCDF_SEPARATE)
  list(APPEND INCLDIR "${INC_NETCDF}")
else()
  list(APPEND INCLDIR "${INC_NETCDF_C}" "${INC_NETCDF_FORTRAN}")
endif()

foreach(ITEM MOD_NETCDF INC_MPI INC_PNETCDF INC_PETSC INC_TRILINOS INC_ALBANY) # INC_MOAB)
  if (${ITEM})
    list(APPEND INCLDIR "${${ITEM}}")
  endif()
endforeach()

if (NOT MCT_LIBDIR)
  set(MCT_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (PIO_LIBDIR)
  if (PIO_VERSION STREQUAL ${PIO_VERSION_MAJOR})
    list(APPEND INCLDIR "${PIO_INCDIR}")
    set(SLIBS "${SLIBS} -L${PIO_LIBDIR}")
  else()
    # If PIO_VERSION_MAJOR doesnt match, build from source
    unset(PIO_LIBDIR)
  endif()
endif()
if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (NOT GPTL_LIBDIR)
  set(GPTL_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (NOT GLC_DIR)
  set(GLC_DIR "${EXEROOT}/glc")
endif()

if (NOT CISM_LIBDIR)
  set(CISM_LIBDIR "${GLC_DIR}/lib")
endif()

if (NOT GLCROOT)
  # Backwards compatibility
  set(GLCROOT "${CIMEROOT}/../components/cism")
endif()

list(APPEND INCLDIR "${INSTALL_SHAREDPATH}/include")

string(FIND "${CAM_CONFIG_OPTS}" "-cosp" HAS_COSP)
if (NOT HAS_COSP EQUAL -1)
  # The following is for the COSP simulator code:
  set(USE_COSP TRUE)
endif()

# System libraries (netcdf, mpi, pnetcdf, esmf, trilinos, etc.)
if (NOT SLIBS)
  if (NOT NETCDF_SEPARATE)
    set(SLIBS "-L${LIB_NETCDF} -lnetcdff -lnetcdf")
  else()
    set(SLIBS "-L${LIB_NETCDF_FORTRAN} -L${LIB_NETCDF_C} -lnetcdff -lnetcdf")
  endif()
endif()

if (LAPACK_LIBDIR)
  set(SLIBS "${SLIBS} -L${LAPACK_LIBDIR} -llapack -lblas")
endif()

if (LIB_MPI)
  if (NOT MPI_LIB_NAME)
    set(SLIBS "${SLIBS} -L${LIB_MPI} -lmpi")
  else()
    set(SLIBS "${SLIBS} -L${LIB_MPI} -l${MPI_LIB_NAME}")
  endif()
endif()

# Add PETSc libraries
if (USE_PETSC)
  set(SLIBS "${SLIBS} ${PETSC_LIB}")
endif()

# Add trilinos libraries; too be safe, we include all libraries included in the trilinos build,
# as well as all necessary third-party libraries
if (USE_TRILINOS)
  set(SLIBS "${SLIBS} -L${LIB_TRILINOS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARIES}")
endif()

# Add Albany libraries.  These are defined in the ALBANY_LINK_LIBS env var that was included above
if (USE_ALBANY)
  set(SLIBS "${SLIBS} ${ALBANY_LINK_LIBS}")
endif()

# Add MOAB libraries.  These are defined in the MOAB_LINK_LIBS env var that was included above
# if (USE_MOAB)
#   set(SLIBS "${SLIBS} ${IMESH_LIBS}")
# endif()

# Add libraries and flags that we need on the link line when C++ code is included
if (USE_CXX)
  if (CXX_LIBS)
    set(SLIBS "${SLIBS} ${CXX_LIBS}")
  endif()

  if (CXX_LDFLAGS)
    set(LDFLAGS "${LDFLAGS} ${CXX_LDFLAGS}")
  endif()
endif()

# Decide whether to use a C++ or Fortran linker, based on whether we
# are using any C++ code and the compiler-dependent CXX_LINKER variable
if (USE_CXX AND CXX_LINKER STREQUAL "CXX")
  set(LD "CXX")
else()
  set(LD "Fortran")
  # Remove arch flag if it exists, it break fortran linking
  string(REGEX REPLACE "-arch[^ ]+" "" LDFLAGS "${LDFLAGS}")
endif()

if (NOT IO_LIB_SRCROOT)
  if (PIO_VERSION STREQUAL 2)
    # This is a pio2 library
    set(PIOLIB "${PIO_LIBDIR}/libpiof.a ${PIO_LIBDIR}/libpioc.a")
    set(PIOLIBNAME "-lpiof -lpioc")
    set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio2")
  else()
    # This is a pio1 library
    set(PIOLIB "${PIO_LIBDIR}/libpio.a")
    set(PIOLIBNAME "-lpio")
    if (NOT EXISTS "${CIMEROOT}/src/externals/pio1/pio")
      set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio1")
    else()
      set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio1/pio")
    endif()
  endif()
else()
  set(IO_LIB_SRC_DIR "IO_LIB_v${PIO_VERSION}_SRCDIR")
  set(PIO_SRC_DIR "${IO_LIB_SRCROOT}/${IO_LIB_SRC_DIR}")
endif()

set(MCTLIBS "${MCT_LIBDIR}/libmct.a ${MCT_LIBDIR}/libmpeu.a")

set(GPTLLIB "${GPTL_LIBDIR}/libgptl.a")

#------------------------------------------------------------------------------
# Set key cmake vars
#------------------------------------------------------------------------------
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}" PARENT_SCOPE)
