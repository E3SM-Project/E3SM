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

# Look for -crm pam in the CAM_CONFIG_OPTS CIME variable
# If it's found, then enable USE_PAM
string(FIND "${CAM_CONFIG_OPTS}" "-crm pam" HAS_PAM)
if (NOT HAS_PAM EQUAL -1)
  # The following is for the PAM code:
  set(USE_PAM TRUE)
endif()

string(FIND "${CAM_CONFIG_OPTS}" "-rrtmgpxx" HAS_RRTMGPXX)
if (NOT HAS_RRTMGPXX EQUAL -1)
  # The following is for the RRTMGPXX code:
  set(USE_RRTMGPXX TRUE)
endif()

# If samxx or rrtmgpxx is being used, then YAKL must be used as well
if (USE_SAMXX OR USE_RRTMGPXX OR USE_PAM)
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

if (NOT DEBUG)
  set(CPPDEFS "${CPPDEFS} -DNDEBUG")
endif()

if (USE_ESMF_LIB)
  set(CPPDEFS "${CPPDEFS} -DUSE_ESMF_LIB")
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

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)
string(FIND "${CPPDEFS}" "-DLINUX" HAS_DLINUX)
string(FIND "${CPPDEFS}" "DBG" HAS_DBG)
string(FIND "${CPPDEFS}" "Darwin" HAS_DARWIN)
if (NOT HAS_DLINUX EQUAL -1 AND HAS_DBG EQUAL -1 AND HAS_DARWIN EQUAL -1)
  set(CPPDEFS "${CPPDEFS} -DHAVE_SLASHPROC")
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
endif()

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
list(APPEND INCLDIR "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include")

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

#------------------------------------------------------------------------------
# Set key cmake vars
#------------------------------------------------------------------------------
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}" PARENT_SCOPE)
