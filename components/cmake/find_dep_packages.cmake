# This file is for finding pacakges needed by E3SM. It should be included
# from the main CMakeLists.txt file.
#
# Finding the correct packages will likely depend on the ${Package}_ROOT
# environment variable being set by config_machines.xml for your machine. Note
# that is environment var is case sensitive.

# Machine env vars should follow the following pattern:
#   <env name="Package_ROOT">$SHELL{if [ -z "$Package_ROOT" ]; then echo /default/install/location; else echo "$Package_ROOT"; fi}</env>
#
# This will allow users to easily specify a different location for all their cases by
# simply setting ${Package}_ROOT in their shell.

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

if (USE_ALBANY)
  find_package(Albany REQUIRED)
endif()

# Albany depends on Trilinos
if (USE_ALBANY OR USE_TRILINOS)
  find_package(Trilinos REQUIRED)
endif()
