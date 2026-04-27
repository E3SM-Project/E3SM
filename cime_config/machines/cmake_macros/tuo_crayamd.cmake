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
