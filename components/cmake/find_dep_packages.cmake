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

# Kokkos' find_package needs to come before other find_packages
# that may define Kokkos targets so we can avoid duplicate target
# errors.
if (USE_KOKKOS)

  if (NOT DEFINED ENV{Kokkos_ROOT})
    # We should use the kokkos from Trilinos if we are using Trilinos
    if (USE_TRILINOS)
      set (ENV{Kokkos_ROOT} ${Trilinos_ROOT})
    else()
      # Kokkos will be built in the sharedlibs if Kokkos_ROOT is unset.
      set(ENV{Kokkos_ROOT} ${INSTALL_SHAREDPATH})
    endif()
  endif()

  find_package(Kokkos REQUIRED)
endif()

# Albany depends on Trilinos
if (USE_ALBANY OR USE_TRILINOS)
  find_package(Trilinos REQUIRED)
endif()

if (USE_MOAB)
  find_package(MOAB REQUIRED)
endif()

if (USE_ALBANY)
  find_package(Albany REQUIRED)
endif()

if (USE_PETSC)
  find_package(PETSc REQUIRED)
endif()

find_package(PIO REQUIRED)
find_package(MCT REQUIRED)
find_package(CsmShare REQUIRED)

# Hack for unsupported blas vendors
if (NOT DEFINED ENV{SKIP_BLAS})
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()
