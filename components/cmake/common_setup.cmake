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
endif()
# The code below is what we actually want but it's currently broken.
# Once fixes are in place, uncomment the line below and remove the 3
# lines above.
# set(CPPDEFS "${CPPDEFS} -DPIO${PIO_VERSION}")

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
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
list(APPEND INCLDIR "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include")

string(FIND "${CAM_CONFIG_OPTS}" "-cosp" HAS_COSP)
if (NOT HAS_COSP EQUAL -1)
  # The following is for the COSP simulator code:
  set(USE_COSP TRUE)
endif()
