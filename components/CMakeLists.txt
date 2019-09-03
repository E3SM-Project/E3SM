#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/scripts/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
# We will prefer space-separated strings over lists
#
#===============================================================================

# bmpersch
# - Read modern cmake docs, use modern features
# - Use find_package for trilinos and other TPLS

cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0057 NEW)

project(E3SM C CXX Fortran)

# Set module path
set(CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

# Load dependency search path.
file(STRINGS ${CMAKE_BINARY_DIR}/Filepath FILEPATH_DIRS)
set(CPP_DIRS ${FILEPATH_DIRS})
list(APPEND CPP_DIRS ".")

# Add INCROOT to path for Depends and Include
set(MINCROOT "")
if (INCROOT)
  list(APPEND CPP_DIRS ${INCROOT})
  set(MINCROOT ${INCROOT})
endif ()

# Determine whether to compile threaded or not
# Set the THREADDIR for the shared build
# based on the threaded build status
if (SMP_PRESENT)
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

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TRILINOS_PATH})
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
  include(${INSTALL_SHAREDPATH}/kokkos_generated_settings.cmake)
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  string (REPLACE ";" " " KOKKOS_LDFLAGS_STR "${KOKKOS_LDFLAGS}")
  set(SLIBS "${SLIBS} ${KOKKOS_LIBS_LIST}")
  set(CXXFLAGS "${CXXFLAGS} ${KOKKOS_CXXFLAGS_STR}")
  set(CXX_LDFLAGS "${CXX_LDFLAGS} ${KOKKOS_LDFLAGS_STR}")
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

if (CPRE)
  foreach(ITEM ${CPPDEFS})
    string(REPLACE "-D" "" NOD ${ITEM})
    set(FPPDEFS "${FPPDEFS} ${CPRE}${NOD}")
  endforeach()
  foreach(ITEM ${EXTRA_PIO_CPPDEFS})
    string(REPLACE "-D" "" NOD ${ITEM})
    set(EXTRA_PIO_FPPDEFS "${EXTRA_PIO_FPPDEFS} ${CPRE}${NOD}")
  endforeach()
else()
  set(FPPDEFS "${CPPDEFS}")
  set(EXTRA_PIO_FPPDEFS "${EXTRA_PIO_CPPDEFS}")
endif()

#===============================================================================
# Set config args for pio and mct to blank and then enable serial
#===============================================================================
if (NOT CONFIG_ARGS)
  set(CONFIG_ARGS "")
endif()
string(FIND "${MODEL}" "pio" HAS_PIO)
if (NOT HAS_PIO EQUAL -1)
  set(CONFIG_ARGS "${CONFIG_ARGS} --enable-timing")
  if (DEBUG)
    set(CONFIG_ARGS "${CONFIG_ARGS} --enable-debug")
  endif()
endif()

#===============================================================================
# User-specified INCLDIR
#===============================================================================

set(INCLDIR "-I.")
if (USER_INCLDIR)
  set(INCLDIR "${INCLDIR} ${USER_INCLDIR}")
endif()

#===============================================================================
# MPI-serial library (part of MCT)
#===============================================================================

if (MPILIB STREQUAL "mpi-serial")
  set(CC ${SCC})
  set(FC ${SFC})
  set(CXX ${SCXX})
  set(MPIFC ${SFC})
  set(MPICC ${SCC})
  set(MPICXX ${SCXX})
  set(CONFIG_ARGS "${CONFIG_ARGS} MCT_PATH=${SHAREDLIBROOT}/${SHAREDPATH}/mct/mpi-serial")
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
# VPATH+=$(CSM_SHR_INCLUDE) # JGF - not sure if I need to do anything here

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/include -I${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include")

if (NOT NETCDF_SEPARATE)
  set(INCLDIR "${INCLDIR} -I${INC_NETCDF}")
else()
  set(INCLDIR "${INCLDIR} -I${INC_NETCDF_C} -I${INC_NETCDF_FORTRAN}")
endif()

foreach(ITEM MOD_NETCDF INC_MPI INC_PNETCDF INC_PETSC INC_TRILINOS INC_ALBANY) # INC_MOAB)
  if (${ITEM})
    set(INCLDIR "${INCLDIR} -I${${ITEM}}")
  endif()
endforeach()

if (MODEL STREQUAL "driver")
  set(INCLDIR "${INCLDIR} -I${EXEROOT}/atm/obj -I${EXEROOT}/ice/obj -I${EXEROOT}/ocn/obj -I${EXEROOT}/glc/obj -I${EXEROOT}/rof/obj -I${EXEROOT}/wav/obj -I${EXEROOT}/esp/obj -I${EXEROOT}/iac/obj")
endif()

if (NOT MCT_LIBDIR)
  set(MCT_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (PIO_LIBDIR)
  if (PIO_VERSION STREQUAL ${PIO_VERSION_MAJOR})
    set(INCLDIR "${INCLDIR} -I${PIO_INCDIR}")
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

set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/include")

#
# Use the MCT dir for the cache for all configure calls because it is the first one
#
set(CFLAGS "${CFLAGS}")
set(CXXFLAGS "${CXXFLAGS}")
set(CONFIG_ARGS "CC=\"${CC}\" FC=\"${FC}\" MPICC=\"${MPICC}\" MPIFC=\"${MPIFC}\" FCFLAGS=\"${FFLAGS} ${FREEFLAGS} ${INCLDIR}\" CPPDEFS=\"${CPPDEFS}\" CFLAGS=\"${CFLAGS} -I.. ${INCLDIR}\" LDFLAGS=\"${LDFLAGS}\"")

if (NOT NETCDF_SEPARATE)
  set(CONFIG_ARGS "${CONFIG_ARGS} NETCDF_PATH=${NETCDF_PATH}")
else()
  set(CONFIG_ARGS "${CONFIG_ARGS} NETCDF_PATH=${NETCDF_C_PATH}")
endif()

if (COMPILER STREQUAL "nag")
  set(CONFIG_ARGS "LIBS=${SLIBS}")
endif()

string(FIND "${CAM_CONFIG_OPTS}" "-cosp" HAS_COSP)
if (NOT HAS_COSP EQUAL -1)
  # The following is for the COSP simulator code:
  get_filename_component(COSP_LIBDIR ${EXEROOT}/atm/obj/cam/src/physics/cosp ABSOLUTE)
endif()


if (MODEL STREQUAL "cam")
  # These RRTMG files take an extraordinarily long time to compile with optimization.
  # Until mods are made to read the data from files, just remove optimization from
  # their compilation.
  set(NOOPT_FILES "cam/src/physics/rrtmg/ext/rrtmg_lw/rrtmg_lw_k_g.f90;cam/src/physics/rrtmg/ext/rrtmg_sw/rrtmg_sw_k_g.f90")

  if (COSP_LIBDIR)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cam/src/physics/cosp/Cosp.cmake)
  endif()

endif()

# System libraries (netcdf, mpi, pnetcdf, esmf, trilinos, etc.)
if (NOT SLIBS)
  if (NOT NETCDF_SEPARATE)
    set(SLIBS "${SLIBS} -L${LIB_NETCDF} -lnetcdff -lnetcdf")
  else()
    set(SLIBS "${SLIBS} -L${LIB_NETCDF_FORTRAN} -L${LIB_NETCDF_C} -lnetcdff -lnetcdf")
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
# We need to do these additions after CONFIG_ARGS is set, because they can sometimes break configure for mct, etc.,
# if they are added to LDFLAGS in CONFIG_ARGS.
if (USE_CXX)
  if (CXX_LIBS)
    set(SLIBS "${SLIBS} ${CXX_LIBS}")
  endif()

  if (CXX_LDFLAGS)
    set(LDFLAGS "${LDFLAGS} ${SLIBS} ${CXX_LDFLAGS}")
  endif()
endif()

# Remove arch flag if it exists
string(REGEX REPLACE "-arch[^ ]+" "" F90_LDFLAGS "${LDFLAGS}")

#------------------------------------------------------------------------------
# Drive configure scripts for support libraries (mct)
#------------------------------------------------------------------------------

# if (NOT EXISTS ${SHAREDLIBROOT}/${SHAREDPATH}/mct/Makefile.conf)
#   execute_process(COMMAND ${CIMEROOT}/src/externals/mct/configure ${CONFIG_ARGS} --srcdir ${CIMEROOT}/src/externals/mct)
# endif()

# if (NOT EXISTS ${SHAREDLIBROOT}/${SHAREDPATH}/mct/mpi-serial/Makefile.conf)
#   execute_process(COMMAND ${CIMEROOT}/src/externals/mct/mpi-serial/configure ${CONFIG_ARGS} --srcdir ${CIMEROOT}/src/externals/mct/mpi-serial)
# endif()

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

set(ULIBS "-L${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/lib -lcsm_share -L${INSTALL_SHAREDPATH}/lib ${PIOLIBNAME} -lgptl -lmct -lmpeu")

#------------------------------------------------------------------------------
# Drive cmake script for cism and pio
#------------------------------------------------------------------------------

set(CMAKE_Fortran_FLAGS "${FFLAGS} ${FPPDEFS} ${EXTRA_PIO_FPPDEFS} ${INCLDIR}")
set(CMAKE_C_FLAGS "${CFLAGS} ${CPPDEFS} ${EXTRA_PIO_CPPDEFS} ${INCLDIR}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS} ${CPPDEFS} ${EXTRA_PIO_CPPDEFS} ${INCLDIR}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(GPTL_PATH "${INSTALL_SHAREDPATH}")
set(PIO_ENABLE_TESTS FALSE)
set(PIO_USE_MALLOC TRUE)
set(USER_CMAKE_MODULE_PATH "${CIMEROOT}/src/CMake;${CIMEROOT}/src/externals/pio2/cmake")

if (NOT PNETCDF_PATH)
  set(WITH_PNETCDF FALSE)
endif()

if (MPILIB STREQUAL "mpi-serial")
  set(PIO_USE_MPISERIAL TRUE)
  set(MPISERIAL_PATH ${INSTALL_SHAREDPATH})
endif()

set(CMAKE_C_COMPILER ${CC})
set(CMAKE_CXX_COMPILER ${CXX})
set(CMAKE_Fortran_COMPILER ${FC})
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}")

# JGF: We only want to do CAM for now
# add_subdirectory(${GLCROOT}/source_cism)
# add_subdirectory(${PIO_SRC_DIR})

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------

# JGF, I think we can rely on CMake's built-in dependency analysus

# Get src files
# JGF: Cmake does not have a VPATH concept, so we need relative/absolute paths to source files
# Note: Using absolute paths seems to wreck CMake's ability to do a dep analysis on fortran
# sources and compile them in the right order.
#
# One additional subtley is that, when mkSrcfiles found multiple files with the same basename,
# only the one found first gets compiled.

#execute_process(COMMAND ${CASETOOLS}/mkSrcfiles)
#file(STRINGS ${CMAKE_BINARY_DIR}/Srcfiles SOURCES)

set(SRCROOT_REL "${CIMEROOT}/..")
set(BASENAME_SET)
file(TO_CMAKE_PATH ${SRCROOT_REL} SRCROOT_ABS)
foreach(DIRSEARCH ${FILEPATH_DIRS})
  file(GLOB MATCHES RELATIVE "${SRCROOT_ABS}/components" "${DIRSEARCH}/*.[Ffc]" "${DIRSEARCH}/*.[Ff]90" "${DIRSEARCH}/*.cpp" "${DIRSEARCH}/*.F90.in")
  if (MATCHES)
    foreach (MATCH IN LISTS MATCHES)
      get_filename_component(BASENAME ${MATCH} NAME)
      list(FIND BASENAME_SET ${BASENAME} BASENAME_WAS_FOUND)
      if (BASENAME_WAS_FOUND EQUAL -1)
        list(APPEND SOURCES ${MATCH})
        list(APPEND BASENAME_SET ${BASENAME})
      else()
        message(WARNING "Skipping repeated base filename ${BASENAME} for ${MATCH}")
      endif()
    endforeach()
  endif()
endforeach()

foreach(SOURCE_FILE IN LISTS SOURCES)
  get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
  if (SOURCE_EXT STREQUAL ".F90.in")
    string(REPLACE ".in" "" SOURCE_NO_IN ${SOURCE_FILE})
    list(APPEND GEN_F90_SOURCES ${SOURCE_NO_IN})
    list(APPEND SOURCES ${SOURCE_NO_IN})
    list(REMOVE_ITEM SOURCES ${SOURCE_FILE})
  endif()
endforeach()

foreach(ITEM IN LISTS CPP_DIRS)
  if (EXISTS ${ITEM})
    set(INCLDIR "${INCLDIR} -I${ITEM}")
  endif()
endforeach()

#-------------------------------------------------------------------------------
# create list of component libraries - hard-wired for current ccsm components
#-------------------------------------------------------------------------------

if (CIME_MODEL STREQUAL "cesm")
  if (COMP_LND STREQUAL "clm")
    set(USE_SHARED_CLM TRUE)
  else()
    set(USE_SHARED_CLM FALSE)
  endif()
else()
  set(USE_SHARED_CLM FALSE)
endif()

if (NOT USE_SHARED_CLM)
  set(LNDOBJDIR "${EXEROOT}/lnd/obj")
  set(LNDLIBDIR "${LIBROOT}")
  if (COMP_LND STREQUAL "clm")
    set(LNDLIB "libclm.a")
  else()
    set(LNDLIB "liblnd.a")
  endif()
  set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
else()
  set(LNDLIB "libclm.a")
  set(LNDOBJDIR "${SHAREDLIBROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/clm/obj")
  set(LNDLIBDIR "${EXEROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/lib")
  set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
  if (MODEL STREQUAL "clm")
    set(INCLUDE_DIR "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
  endif()
endif()

if (NOT ULIBDEP)
  if (LIBROOT)
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libatm.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libice.a")
    set(ULIBDEP "${ULIBDEP} ${LNDLIBDIR}/${LNDLIB}")
    set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libocn.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/librof.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libglc.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libwav.a")
    if (NOT COMP_INTERFACE STREQUAL "nuopc")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libiac.a")
    endif()
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libesp.a")
  endif()
endif()

if (COMP_GLC STREQUAL "cism")
  set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismfortran.a")
  if (CISM_USE_TRILINOS)
    set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismcpp.a")
  endif()
endif()
if (OCN_SUBMODEL STREQUAL "moby")
  set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libmoby.a")
endif()

# JGF, this kind of thing should be done via target_link_library in cmake
if (NOT CLIBS)
  if (ULIBDEP)
    # For each occurrence of something like /path/to/foo/libbar.a in ULIBDEP,
    # CLIBS will contain -L/path/to/foo -lbar
    foreach (ITEM ${ULIBDEP})
      get_filename_component(ULIB_DIR ${ITEM} DIRECTORY)
      get_filename_component(ULIB_LIB ${ITEM} NAME)
      set(CLIBS "${CLIBS} -L${ULIB_DIR} -l${ULIB_LIB}")
    endforeach()
  endif()
endif()

# libcsm_share.a is in ULIBDEP, but -lcsm_share is in ULIBS rather than CLIBS,
# so this needs to be added after creating CLIBS above
set(CSMSHARELIB "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/lib/libcsm_share.a")
set(ULIBDEP "${ULIBDEP} ${CSMSHARELIB}")

# Things may have been added to flags and we need to remove PIO stuff
set(CMAKE_Fortran_FLAGS "${FFLAGS} ${INCLDIR}")
set(CMAKE_C_FLAGS "${CFLAGS} ${INCLDIR}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS} ${INCLDIR}")

# do the necessary genf90s
add_custom_target(genf90
  DEPENDS ${CIMEROOT}/src/externals/genf90/genf90.pl)
foreach (SRC_FILE IN LISTS GEN_F90_SOURCES)
  get_filename_component(BASENAME ${SRC_FILE} NAME)
  add_custom_command (
    OUTPUT ${CMAKE_BINARY_DIR}/${BASENAME}
    COMMAND ${CIMEROOT}/src/externals/genf90/genf90.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in > ${CMAKE_BINARY_DIR}/${BASENAME}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in genf90)
  list(REMOVE_ITEM SOURCES ${SRC_FILE})
  list(APPEND SOURCES ${CMAKE_BINARY_DIR}/${BASENAME})
endforeach ()

# Flags are slightly different for different fortran extensions
foreach (SOURCE_FILE IN LISTS SOURCES)
  # Cosp manages its own flags
  if (NOT SOURCE_FILE IN_LIST COSP_SOURCES)
    get_filename_component(SOURCE_EXT ${SOURCE_FILE} EXT)
    if (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".f")
      set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FIXEDFLAGS}")
    elseif(SOURCE_EXT STREQUAL ".f90")
      set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FREEFLAGS}")
    elseif(SOURCE_EXT STREQUAL ".F90")
      set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FREEFLAGS} ${CONTIGUOUS_FLAG}")
    endif()
  endif()
endforeach()

# Load machine/compiler specific settings
set(COMPILER_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${COMPILER}.cmake)
set(MACHINE_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.cmake)
set(PLATFORM_SPECIFIC_DEPENDS ${CASEROOT}/Depends.${MACH}.${COMPILER}.cmake)
set(TRY_TO_LOAD ${COMPILER_SPECIFIC_DEPENDS} ${MACHINE_SPECIFIC_DEPENDS} ${PLATFORM_SPECIFIC_DEPENDS})
foreach(ITEM IN LISTS TRY_TO_LOAD)
  if (EXISTS ${ITEM})
    include(${ITEM})
  endif()
endforeach()

# Disable optimizations on some files that would take too long to compile, expect these to all be fortran files
foreach (SOURCE_FILE IN LISTS NOOPT_FILES)
  set_property(SOURCE ${SOURCE_FILE} APPEND_STRING PROPERTY COMPILE_FLAGS " ${FFLAGS_NOOPT}")
endforeach()

#-------------------------------------------------------------------------------
# build rules:
#-------------------------------------------------------------------------------

if (MPILIB STREQUAL "mpi-serial")
  set(MPISERIAL "${INSTALL_SHAREDPATH}/lib/lib-mpi-serial.a")
  set(MLIBS "${MLIBS} ${MPISERIAL}")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)

if (COMPLIB)
  add_library(${COMPLIB})
  target_sources(${COMPLIB} PRIVATE ${SOURCES})

  # Subtle: In order for fortran dependency scanning to work, our CPPFPP/DEFS must be registered
  # as COMPILE_DEFINITIONS, not simple added via CMAKE_Fortran_Flags. Also, CPPDEFS *must*
  # be provided as a list, not a whitespace-separated string; otherwise, things get wonky.
  separate_arguments(CPPDEFS_LIST UNIX_COMMAND "${CPPDEFS}")
  target_compile_definitions(${COMPLIB} PRIVATE ${CPPDEFS_LIST})
  add_dependencies(${COMPLIB} genf90)
endif()

if (EXEC_SE)
  add_executable(${EXEC_SE})
  target_link_library(${EXEC_SE} ${COMPLIB} ${ULIBDEP} ${CSMSHARELIB} ${MCTLIBS} ${PIOLIB} ${GPTLLIB} ${CLIBS} ${ULIBS} ${SLIBS} ${MLIBS} ${F90_LDFLAGS})
endif()
