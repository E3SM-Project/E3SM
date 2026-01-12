set(KOKKOS_USE_EKAT_MACH_FILE "TRUE")
string(APPEND CMAKE_C_FLAGS " -I/global/common/software/nersc/pe/conda-envs/24.1.0/python-3.11/nersc-python/include/python3.11")
string(APPEND SLIBS " -L/global/common/software/nersc/pe/conda-envs/24.1.0/python-3.11/nersc-python/lib -lpython3.11")
set(Python_EXECUTABLE "/global/common/software/nersc/pe/conda-envs/24.1.0/python-3.11/nersc-python/bin/python")
