string(APPEND CONFIG_ARGS " --host=cray")
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")

if (COMP_NAME STREQUAL gptl)
	string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()

# required to resolve bshr_infnan_mod.F90 compile issue
string(APPEND CPPDEFS " -DCPRCRAY")

string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_HIP=OFF -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ARCH_AMD_GFX942=OFF -DKokkos_ARCH_AMD_GFX942_APU=OFF")

# resolve SPIO compile issue
string(APPEND SPIO_CMAKE_OPTS " -DPIO_ENABLE_TOOLS:BOOL=OFF")

# cray fortran outputs uppercase mod name, forces lowercase
string(APPEND CMAKE_Fortran_FLAGS " -ef")

string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")

if (COMP_NAME STREQUAL cpl)
	# auto-detection not working, tries to use AMD variant, we need Cray variant
	string(APPEND CMAKE_EXE_LINKER_FLAGS " -L/opt/cray/pe/libsci/24.11.0/CRAY/18.0/x86_64/lib -lsci_cray")
endif()

set(E3SM_LINK_WITH_FORTRAN "TRUE")

set(PIO_FILESYSTEM_HINTS "lustre")

# Keep the runtime search path attached to the binaries instead of relying on
# LD_LIBRARY_PATH in the user's shell or batch environment.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

foreach (_prefix
	"$ENV{NETCDF_C_PATH}"
	"$ENV{NETCDF_FORTRAN_PATH}"
	"$ENV{PNETCDF_PATH}"
	"$ENV{HDF5_ROOT}"
	"$ENV{TEMPESTREMAP_ROOT}"
	"$ENV{MOAB_ROOT}")
	list(APPEND CMAKE_BUILD_RPATH "${_prefix}/lib" "${_prefix}/lib64")
	list(APPEND CMAKE_INSTALL_RPATH "${_prefix}/lib" "${_prefix}/lib64")
endforeach()

# The MOAB coupler pulls in a MOAB build that is UBSan-instrumented in this
# environment, so the final executable must link the UBSan runtime.
if (DEFINED COMP_INTERFACE AND COMP_INTERFACE STREQUAL "moab")
	string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--no-as-needed -lubsan -Wl,--as-needed")
endif()

list(REMOVE_DUPLICATES CMAKE_BUILD_RPATH)
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,--enable-new-dtags")
string(APPEND CMAKE_SHARED_LINKER_FLAGS " -Wl,--enable-new-dtags")
