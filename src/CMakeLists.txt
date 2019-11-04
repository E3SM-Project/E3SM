include(${CASEROOT}/Macros.cmake)

#
# General setup
#

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

set(CMAKE_C_COMPILER ${MPICC})
set(CMAKE_CXX_COMPILER ${MPICXX})
set(CMAKE_Fortran_COMPILER ${MPIFC})
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Set up CPPDEFS
set(FILE_OFFSET "-DOFFSET64BIT")
if (CPPDEFS)
  separate_arguments(CPPDEFS UNIX_COMMAND "${CPPDEFS}")
endif()
list(APPEND CPPDEFS "-DMPAS_NO_LOG_REDIRECT" "-DUSE_PIO2" "-DMPAS_NO_ESMF_INIT" "-DMPAS_ESM_SHR_CONST" "-DMPAS_PERF_MOD_TIMERS" "${MODEL_FORMULATION}" "${FILE_OFFSET}" "${ZOLTAN_DEFINE}" "-D_MPI" "-DMPAS_NAMELIST_SUFFIX=${NAMELIST_SUFFIX}" "-DMPAS_EXE_NAME=${EXE_NAME}")
if (DEBUG)
  list(APPEND CPPDEFS "-DMPAS_DEBUG")
endif()
if (compile_threaded)
  list(APPEND CPPDEFS "-DMPAS_OPENMP")
endif()

set(INCLUDES "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share" "${INSTALL_SHAREDPATH}/pio" "${PNETCDF_PATH}/include" "${CMAKE_CURRENT_SOURCE_DIR}/external/ezxml" "${CMAKE_BINARY_DIR}/framework" "${CMAKE_BINARY_DIR}/operators")
if (NETCDF_PATH)
  list(APPEND INCLUDES ${NETCDF_PATH}/include)
else()
  if (NETCDF_C_PATH)
    list(APPEND INCLUDES ${NETCDF_C_PATH}/include)
  endif()
  if (NETCDF_FORTRAN_PATH)
    list(APPEND INCLUDES ${NETCDF_FORTRAN_PATH}/include)
  endif()
endif()

if (USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/kokkos_generated_settings.cmake)
  string (REPLACE ";" " " KOKKOS_CXXFLAGS_STR "${KOKKOS_CXXFLAGS}")
  set(CXXFLAGS "${CXXFLAGS} ${KOKKOS_CXXFLAGS_STR}")
endif()

set(CMAKE_Fortran_FLAGS "${FFLAGS}")
set(CMAKE_C_FLAGS "${CFLAGS}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS}")

add_subdirectory(tools)

# Gather sources that are needed for all cores

# externals
set(COMMON_RAW_SOURCES external/ezxml/ezxml.c)

# framework
list(APPEND COMMON_RAW_SOURCES
  framework/mpas_kind_types.F
  framework/mpas_framework.F
  framework/mpas_timer.F
  framework/mpas_timekeeping.F
  framework/mpas_constants.F
  framework/mpas_attlist.F
  framework/mpas_hash.F
  framework/mpas_sort.F
  framework/mpas_block_decomp.F
  framework/mpas_block_creator.F
  framework/mpas_dmpar.F
  framework/mpas_abort.F
  framework/mpas_decomp.F
  framework/mpas_threading.F
  framework/mpas_io.F
  framework/mpas_io_streams.F
  framework/mpas_bootstrapping.F
  framework/mpas_io_units.F
  framework/mpas_stream_manager.F
  framework/mpas_stream_list.F
  framework/mpas_forcing.F
  framework/mpas_c_interfacing.F
  framework/random_id.c
  framework/pool_hash.c
  framework/mpas_derived_types.F
  framework/mpas_domain_routines.F
  framework/mpas_field_routines.F
  framework/mpas_pool_routines.F
  framework/xml_stream_parser.c
  framework/regex_matching.c
  framework/mpas_field_accessor.F
  framework/mpas_log.F
)

# operators
list(APPEND COMMON_RAW_SOURCES
  operators/mpas_vector_operations.F
  operators/mpas_matrix_operations.F
  operators/mpas_tensor_operations.F
  operators/mpas_rbf_interpolation.F
  operators/mpas_vector_reconstruction.F
  operators/mpas_spline_interpolation.F
  operators/mpas_tracer_advection_helpers.F
  operators/mpas_tracer_advection_mono.F
  operators/mpas_tracer_advection_std.F
  operators/mpas_geometry_utils.F
)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_utils.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)

add_library(common OBJECT)
target_compile_definitions(common PRIVATE ${CPPDEFS})
target_include_directories(common PRIVATE ${INCLUDES})

genf90_targets("${COMMON_RAW_SOURCES}" "${INCLUDES}" "${CPPDEFS}" "" "")
target_sources(common PRIVATE ${SOURCES})

foreach(CORE IN LISTS CORES)
  build_core(${CORE})
endforeach()
