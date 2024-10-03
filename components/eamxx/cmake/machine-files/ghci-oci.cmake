include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(CMAKE_Fortran_FLAGS "-Wno-maybe-uninitialized -Wno-unused-dummy-argument -fallow-argument-mismatch"  CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "-fvisibility-inlines-hidden -fmessage-length=0 -Wno-use-after-free -Wno-unused-variable -Wno-maybe-uninitialized" CACHE STRING "" FORCE)

# TODO: figure out a better way to handle this, e.g.,
# TODO: --map-by ppr:1:node:pe=1 doesn't work with mpich,
# TODO: but -map-by core:1:numa:hwthread=1 may work well?
# TODO: this will need to be handled in EKAT at some point
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "-np")

# TODO: hack in place to get eamxx to recognize CPRNC
# TODO: See note in BuildCprnc.cmake...
set(ENV{CCSM_CPRNC} "/usr/local/packages/bin/cprnc")
