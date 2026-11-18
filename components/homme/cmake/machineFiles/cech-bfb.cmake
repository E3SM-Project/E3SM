# CMake initial cache file for cech (NVIDIA GB200 Grace-Blackwell) for BFB testing
#
# NOTE: BFB-specific cache variables must be set BEFORE including cech.cmake
# (the first value set in the cache wins).
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")
SET (OPT_FLAGS "-O0" CACHE STRING "")
SET (ADD_Fortran_FLAGS "-i4 -Mstack_arrays -Mextend -byteswapio -Mflushz -DHAVE_IEEE_ARITHMETIC -Mallocatable=03 -DNO_R16 -traceback -DCPRNVIDIA -DFORTRANUNDERSCORE -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")
SET (ADD_C_FLAGS       "-DFORTRANUNDERSCORE -DNO_R16 -DCPRNVIDIA -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-DTHRUST_IGNORE_CUB_VERSION_CHECK -DFORTRANUNDERSCORE -DNO_R16 -DCPRNVIDIA -tp=grace -noacc -Mnofma -O0 -Mnovect -Kieee" CACHE STRING "")

INCLUDE ("${CMAKE_CURRENT_LIST_DIR}/cech.cmake")
