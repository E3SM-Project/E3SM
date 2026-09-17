# Mac arm64 (Apple Silicon) compiler overrides for the generic gnu.cmake macro.
# Included after gnu.cmake by Macros.cmake, so plain set() wins.
#
# Toolchain (all arm64-native):
#   MPI      : Homebrew Open MPI 5.0 wrappers (mpicc->Apple clang, mpif90->gfortran 16)
#   Serial C : Homebrew gcc-14   (MacPorts gcc is broken against the macOS 26 SDK)
#   Serial F : Homebrew gfortran-14
#
# NOTE: the conda gcc driver needs e3sm-build/bin on PATH for its vendored
# clang linker; the Homebrew toolchain has no such requirement.
set(SCC "/opt/homebrew/bin/gcc-14")
set(SCXX "/opt/homebrew/bin/g++-14")
set(SFC "/opt/homebrew/bin/gfortran-14")
set(MPIFC "/opt/homebrew/bin/mpif90")
set(MPICC "/opt/homebrew/bin/mpicc")
set(MPICXX "/opt/homebrew/bin/mpicxx")
set(E3SM_LINK_WITH_FORTRAN "TRUE")

# gnu.cmake maps Darwin's arm64 uname to -mcmodel=large, but GCC on aarch64
# does not implement code model 'large'. aarch64 gcc only supports 'small'
# (its default), matching the aarch64 branch of gnu.cmake.
string(REPLACE "-mcmodel=large" "-mcmodel=small" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
string(REPLACE "-mcmodel=large" "-mcmodel=small" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

# GCC >= 14 hard-errors on legacy constructs still present in old bundled
# code (mpi-serial, mct). Keep them as warnings.
string(APPEND CMAKE_C_FLAGS " -Wno-error=implicit-int -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types -Wno-error=int-conversion")

# The conda netcdf-c hides the legacy '_FillValue' attribute-name macro behind
# this option; scorpio (a pinned submodule) still uses it.
string(APPEND CPPDEFS " -DNETCDF_ENABLE_LEGACY_MACROS")