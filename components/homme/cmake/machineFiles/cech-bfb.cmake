# CMake initial cache file for cech (NVIDIA GB200 Grace-Blackwell) for BFB testing
#
# NOTE: BFB-specific cache variables must be set BEFORE including cech.cmake
# (the first value set in the cache wins).
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")
SET (OPT_FLAGS "-O0" CACHE STRING "")
SET (ADD_Fortran_FLAGS "-i4 -Mstack_arrays -Mextend -byteswapio -Mflushz -DHAVE_IEEE_ARITHMETIC -Mallocatable=03 -DNO_R16 -traceback -DCPRNVIDIA -DFORTRANUNDERSCORE -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")
SET (ADD_C_FLAGS       "-DFORTRANUNDERSCORE -DNO_R16 -DCPRNVIDIA -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-DTHRUST_IGNORE_CUB_VERSION_CHECK -DFORTRANUNDERSCORE -DNO_R16 -DCPRNVIDIA -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")

# Prevent CMake's default NVHPC Debug flags (-g -O0 -Mbounds / -g -O0) from
# being appended after ADD_*_FLAGS when CMAKE_BUILD_TYPE=Debug (HOMMEBFB_D).
# SetCompilerFlags.cmake already prepends -g, and -Mbounds inserts pgf90_subchk64
# calls on every subscript causing a ~46x Fortran slowdown and Slurm timeouts.
SET (CMAKE_Fortran_FLAGS_DEBUG "" CACHE STRING "")
SET (CMAKE_C_FLAGS_DEBUG       "" CACHE STRING "")
SET (CMAKE_CXX_FLAGS_DEBUG     "" CACHE STRING "")

SET (HOMMEXX_EXEC_SPACE CUDA CACHE STRING "")
SET (HOMMEXX_CUDA_MIN_BLOCKS "3" CACHE STRING "")

SET (BUILD_HOMME_WITHOUT_PIOLIBRARY FALSE CACHE BOOL "")


SET (Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET (Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
SET (Kokkos_ENABLE_CUDA ON CACHE BOOL "")
SET (Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
SET (Kokkos_ARCH_BLACKWELL100 ON CACHE BOOL "")
SET (CMAKE_CUDA_ARCHITECTURES "100" CACHE STRING "")
SET (Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC OFF CACHE BOOL "")
SET (Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET (CXXLIB_SUPPORTED_CACHE FALSE CACHE BOOL "")

SET (ENABLE_OPENMP OFF CACHE BOOL "")
SET (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET (ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

SET (CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "")

SET (USE_NUM_PROCS 4 CACHE STRING "") # only default

GET_FILENAME_COMPONENT(CECH_AFFINITY_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/../../../../cime_config/machines/scripts/cech_set_affinity.sh" ABSOLUTE)
SET (MPIEXEC_EXECUTABLE "srun" CACHE FILEPATH "")
SET (MPIEXEC_PREFLAGS "--mpi=pmix -K -c 36 --cpu-bind=cores -m block:block ${CECH_AFFINITY_SCRIPT}" CACHE STRING "")
SET (HOMME_TESTING_TIMELIMIT 1800 CACHE STRING "")
