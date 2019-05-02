#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/scripts/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
#===============================================================================

cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0057 NEW)

project(E3SM CXX Fortran)

# Load dependency search path.
set(DIRS ".")
file(STRINGS Filepath FILEPATH_DIRS)
list(APPEND DIRS ${FILEPATH_DIRS})

set(CPP_DIRS ${DIRS})
# Add INCROOT to path for Depends and Include
set(MINCROOT "")
if (INCROOT)
  list(APPEND CPP_DIRS ${INCROOT})
  set(MINCROOT ${INCROOT})
endif ()

# # Expand any tildes in directory names. Change spaces to colons. JGF
# VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
# VPATH := $(subst $(space),:,$(VPATH))

# RM    := rm
# CP    := cp

# exec_se: $(EXEC_SE)  Depends
# complib: $(COMPLIB)  Depends

# Determine whether to compile threaded or not
# Set the THREADDIR for the shared build
# based on the threaded build status
if (SMP OR SMP_PRESENT)
  set(THREADDIR "threads")
  set(compile_threaded TRUE)
else()
  set(THREADDIR "nothreads")
  set(compile_threaded FALSE)
endif()

# set the debug directory based on the debug status
if (DEBUG)
  set(DEBUGDIR "debug")
else()
  set(DEBUGDIR "nodebug")
endif()

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

#===============================================================================
# set CPP options (must use this before any flags or cflags settings)
#===============================================================================

set(CPPDEFS ${USER_CPPDEFS} "-D${OS}")

include(${CASEROOT}/Macros)

# SLIBS comes from Macros, so this append must come after Macros are included
if (USE_FMS)
  list(APPEND SLIBS "-lfms")
endif()

if (DEBUG)
  if (CIME_MODEL STREQUAL "cesm")
    list(APPEND CPPDEFS "-DDEBUG")
  endif()
else()
  list(APPEND CPPDEFS "-DNDEBUG")
endif()

if (USE_ESMF_LIB)
  list(APPEND CPPDEFS "-DUSE_ESMF_LIB")
endif()

if (COMP_INTERFACE STREQUAL "nuopc")
  list(APPEND CPPDEFS "-DNUOPC_INTERFACE")
else()
  list(APPEND CPPDEFS "-DMCT_INTERFACE")
endif()

if (COMPARE_TO_NUOPC)
  list(APPEND CPPDEFS "-DCOMPARE_TO_NUOPC")
endif()

if (MPILIB STREQUAL mpi-serial)
  list(APPEND CPPDEFS "-DNO_MPI2")
else()
  list(APPEND CPPDEFS "-DHAVE_MPI")
endif()

if (PIO_VERSION STREQUAL "1")
  list(APPEND CPPDEFS "-DPIO1")
else()
  set(USE_CXX TRUE)
endif()

if (NOT SHAREDPATH)
  set(SHAREDPATH ${COMPILER}/${MPILIB}/${DEBUGDIR}/${THREADDIR}/${COMP_INTERFACE})
  set(INSTALL_SHAREDPATH ${EXEROOT}/${SHAREDPATH})
endif()

# Decide whether to use a C++ or Fortran linker, based on whether we
# are using any C++ code and the compiler-dependent CXX_LINKER variable
if (USE_CXX AND CXX_LINKER STREQUAL "CXX")
  if (MPILIB STREQUAL mpi-serial)
    set(LD ${SCXX})
  else()
    set(LD ${MPICXX})
  endif()
else()
  if (MPILIB STREQUAL mpi-serial)
    set(LD ${SFC})
  else()
    set(LD ${MPIFC})
  endif()
endif()

if (USE_CXX AND NOT SUPPORTS_CXX)
  message(FATAL_ERROR "Fatal attempt to include C++ code on a compiler/machine combo that has not been set up to support C++")
endif()

# Not clear how to escape commas for libraries with their own configure
# script, and they don't need this defined anyway, so leave this out of
# FPPDEFS. JGF
if (HAS_F2008_CONTIGUOUS)
  set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS=contiguous,")
else()
  set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS")
endif()

if (CPRE)
  set(CONTIGUOUS_FLAG ${CPRE})
endif()

if (NOT AR)
  set(AR "ar")
endif()
if (NOT ARFLAGS)
  set(ARFLAGS "-r")
endif()

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
    set(LIB_NETCDF_C ${NETCDF_C_PATH}/lib)
  endif()
  if (NOT LIB_NETCDF_FORTRAN)
    set(LIB_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/lib)
  endif()
elseif (NETCDF_FORTRAN_PATH)
  message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
elseif (NETCDF_PATH)
  set(NETCDF_SEPARATE FALSE)
  if (NOT INC_NETCDF)
    set(INC_NETCDF ${NETCDF_PATH}/include)
  endif()
  if (NOT LIB_NETCDF)
    set(LIB_NETCDF ${NETCDF_PATH}/lib)
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

  # # Get the "PETSC_LIB" list an env var  JGF
  #include $(PETSC_PATH)/lib/petsc/conf/variables
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

  # # get a bunch of variables related to this trilinos installation;
  # these variables begin with "Trilinos_" JGF
  #include $(INC_TRILINOS)/Makefile.export.Trilinos
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

  # # get the "ALBANY_LINK_LIBS" list as an env var JGF
  #include $(ALBANY_PATH)/export_albany.in
endif()

if (USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/Kokkos)
  list(APPEND SLIBS ${KOKKOS_LIBS})
  list(APPEND CXXFLAGS ${KOKKOS_CXXFLAGS})
  list(APPEND CXX_LDFLAGS ${KOKKOS_LDFLAGS})
endif()

if (USE_MOAB)
  if (MOAB_PATH)
    list(APPEND CPPDEFS "-DHAVE_MOAB")
    if (NOT INC_MOAB)
      set(INC_MOAB ${MOAB_PATH}/include)
    endif()
    if (NOT LIB_MOAB)
      set(LIB_MOAB ${MOAB_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "MOAB_PATH must be defined when USE_MOAB is TRUE")
  endif()

  # # get the "IMESH_LIBS" list as an env var JGF
  #include $(LIB_MOAB)/iMesh-Defs.inc
endif()

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)
string(FIND CPPDEFS "-DLINUX" HAS_DLINUX)
string(FIND CPPDEFS "DBG" HAS_DBG)
string(FIND CPPDEFS "Darwin" HAS_DARWIN)
if (NOT HAS_DLINUX EQUAL -1 AND HAS_DBG EQUAL -1 AND HASDARWIN EQUAL -1)
  list(APPEND CPPDEFS "-DHAVE_SLASHPROC")
endif()

# Atleast on Titan+cray mpi, MPI_Irsends() are buggy, causing hangs during I/O
# Force PIO to use MPI_Isends instead of the default, MPI_Irsends
if (PIO_VERSION STREQUAL 2)
  set(EXTRA_PIO_CPPDEFS "-DUSE_MPI_ISEND_FOR_FC")
else()
  set(EXTRA_PIO_CPPDEFS "-D_NO_MPI_RSEND")
endif

if (LIB_PNETCDF)
  list(APPEND CPPDEFS "-D_PNETCDF")
  list(APPEND SLIBS "-L${LIB_PNETCDF}" "-lpnetcdf")
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
  # include(${CIME_ESMFMKFILE} # JGF
  list(APPEND CPPDEFS "-DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR}" "-DESMF_VERSION_MINOR=${ESMF_VERSION_MINOR}")
  list(APPEND FFLAGS ${ESMF_F90COMPILEPATHS})
  list(APPEND SLIBS ${ESMF_F90LINKPATHS} ${ESMF_F90LINKRPATHS} ${ESMF_F90ESMFLINKLIBS})
endif()

# JGF
# ifdef CPRE
#   FPPDEFS := $(subst $(comma),\\$(comma),$(CPPDEFS))
#   FPPDEFS := $(patsubst -D%,$(CPRE)%,$(FPPDEFS))
#   EXTRA_PIO_FPPDEFS := $(subst $(comma),\\$(comma),$(EXTRA_PIO_CPPDEFS))
#   EXTRA_PIO_FPPDEFS := $(patsubst -D%,$(CPRE)%,$(EXTRA_PIO_FPPDEFS))
# else
#   FPPDEFS := $(CPPDEFS)
#   EXTRA_PIO_FPPDEFS := $(EXTRA_PIO_CPPDEFS)
# endif

#===============================================================================
# Set config args for pio and mct to blank and then enable serial
#===============================================================================
if (NOT CONFIG_ARGS)
  set(CONFIG_ARGS "")
endif()
string(FIND MODEL "pio" HAS_PIO)
if (NOT HAS_PIO EQUAL -1)
  list(APPEND CONFIG_ARGS "--enable-timing")
  if (DEBUG)
    list(APPEND CONFIG_ARGS "--enable-debug")
  endif()
endif()

#===============================================================================
# User-specified INCLDIR
#===============================================================================

set(INCLDIR "-I.")
if (USER_INCLDIR)
  list(APPEND INCLDIR ${USER_INCLDIR})
endif()

#===============================================================================
# MPI-serial library (part of MCT)
#===============================================================================

if (MPILIB STREQUAL mpi-serial)
  set(CC ${SCC})
  set(FC ${SFC})
  set(CXX ${SCXX})
  set(MPIFC ${SFC})
  set(MPICC ${SCC})
  set(MPICXX ${SCXX})
  list(APPEND CONFIG_ARGS "MCT_PATH=${SHAREDLIBROOT}/${SHAREDPATH}/mct/mpi-serial")
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
# VPATH+=$(CSM_SHR_INCLUDE) # JGF

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================

