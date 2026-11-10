string(APPEND CONFIG_ARGS " --host=cray")
set(USE_CUDA "TRUE")
set(HOMME_QUAD_PREC FALSE CACHE BOOL "") # nvidia does not support QUAD (real*16)
string(APPEND CPPDEFS " -DGPU")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_SLASHPROC -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
string(APPEND CMAKE_Fortran_FLAGS " -tp=zen3")
string(APPEND CMAKE_C_FLAGS       " -tp=zen3")
string(APPEND CMAKE_CXX_FLAGS     " -tp=zen3")

# Determine if case uses EAMxx / Kokkos
if (NOT DEFINED COMP_ATM AND EXISTS "${CASEROOT}/env_case.xml")
  file(READ "${CASEROOT}/env_case.xml" _env_case_content)
  if (_env_case_content MATCHES "<entry id=\"COMP_ATM\"[^>]*value=\"([^\"]*)\"")
    set(COMP_ATM "${CMAKE_MATCH_1}")
  endif()
endif()

if (COMP_ATM STREQUAL "scream" OR USE_KOKKOS OR CAM_TARGET MATCHES "kokkos")
  set(CASE_USES_KOKKOS TRUE)
else()
  set(CASE_USES_KOKKOS FALSE)
endif()

if (CASE_USES_KOKKOS)
  # EAMxx / Kokkos path:
  # Must use -noacc so NVHPC/crayftn does not generate __acc_compiled conflicts with Kokkos/CUDA.
  # Remove -Mstack_arrays set by generic nvidiagpu.cmake since it is not safe with -noacc builds.
  # -Mnostackarrays (the old PGI negation flag) was removed in NVHPC 25.x; use
  # string(REPLACE) to strip -Mstack_arrays from the flags instead.
  string(REPLACE "-Mstack_arrays" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  string(APPEND CMAKE_C_FLAGS       " -noacc")
  string(APPEND CMAKE_CXX_FLAGS     " -noacc")
  string(APPEND CMAKE_Fortran_FLAGS " -noacc")
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ARCH_AMPERE80=On -DKokkos_ENABLE_CUDA=On -DKokkos_ENABLE_CUDA_LAMBDA=On -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=Off -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF -DCMAKE_CXX_FLAGS=-noacc -DCMAKE_C_FLAGS=-noacc")
  # Use KOKKOS_OPTIONS (above) to build the EKAT sharedlib instead of the ekat
  # machine file.  pm-gpu.cmake sets this to TRUE, but the ekat machine file for
  # pm-gpu lacks -noacc, so NVHPC embeds __acc_compiled in EKAT objects causing
  # linker failures.  KOKKOS_OPTIONS already supplies -noacc and the required
  # Kokkos arch/CUDA flags, so the machine file is redundant here.
  set(KOKKOS_USE_EKAT_MACH_FILE "FALSE")
  # NVHPC embeds __acc_compiled in every object file unless compiled with -noacc.
  # EKAT objects (e.g. libekat_yamlparser) may reference this symbol, so we must
  # satisfy it at link time via libacchost.
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_cuda -noacc -lacchost")
else()
  # Non-EAMxx / OpenACC path (e.g. MPAS OpenACC):
  string(APPEND CPPDEFS " -DMPAS_OPENACC")
  string(APPEND CMAKE_Fortran_FLAGS " -acc -gpu=cc80 -Minfo=accel")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -acc -gpu=cc80 -Minfo=accel")
endif()

if (COMP_NAME STREQUAL cice)
  # nvfortran -Mbounds produces false-positive subscript-out-of-range errors
  # for module-level fixed-size arrays (e.g. TLON) accessed inside OpenMP
  # parallel regions in CICE.  Disable bounds checking for CICE only.
  string(APPEND CMAKE_Fortran_FLAGS " -Mnobounds")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -Mnobounds")
endif()
string(APPEND CMAKE_CUDA_FLAGS " -ccbin CC -O2 -arch sm_80 --use_fast_math")
set(CMAKE_CUDA_ARCHITECTURES "80")
string(APPEND CMAKE_Fortran_FLAGS " -Mnofma")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Mnofma")
string(REPLACE "-Ktrap=fp" "-Ktrap=divz,ovf" CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
set(SCC "cc")
set(SCXX "CC")
set(SFC "ftn")
