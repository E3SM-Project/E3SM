set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

# From cime/share/util/CMakeLists.txt
set(share_util_path ${CIMEROOT}/src/share/util)
set(genf90_files ${share_util_path}/shr_infnan_mod.F90.in ${share_util_path}/shr_assert_mod.F90.in)

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_BINARY_DIR}
  share_genf90_sources)

sourcelist_to_parent(share_genf90_sources)

list(APPEND share_sources "${share_genf90_sources}")

list(APPEND share_sources
  ${share_util_path}/shr_file_mod.F90
  ${share_util_path}/shr_kind_mod.F90
  ${share_util_path}/shr_const_mod.F90
  ${share_util_path}/shr_sys_mod.F90
  ${share_util_path}/shr_log_mod.F90
  ${share_util_path}/shr_strconvert_mod.F90
  ${share_util_path}/shr_nl_mod.F90
  ${share_util_path}/shr_string_mod.F90
  ${share_util_path}/shr_timer_mod.F90
  ${share_util_path}/shr_abort_mod.F90
  ${share_util_path}/shr_mpi_mod.F90
  ${share_util_path}/shr_taskmap_mod.F90
  ${share_util_path}/../timing/perf_mod.F90
  ${share_util_path}/../timing/perf_utils.F90)

# end share/util

add_definitions(-DPCOLS=1 -DPSUBCOLS=1 -DPLEV=1)

set(P3_SRCS
  p3_f90.cpp
  p3_ic_cases.cpp
  micro_p3_iso_c.f90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/micro_p3.F90
  ${SCREAM_BASE_DIR}/../cam/src/control/cam_logfile.F90
  ${SCREAM_BASE_DIR}/../cam/src/utils/spmd_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/utils/cam_abortutils.F90
  ${SCREAM_BASE_DIR}/../cam/src/utils/namelist_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/control/units.F90
  ${SCREAM_BASE_DIR}/../cam/src/control/mpishorthand.F
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/micro_p3_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/control/physconst.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/ppgrid.F90
)

list(APPEND P3_SRCS ${share_sources})

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND P3_SRCS p3_functions_upwind.cpp p3_functions_table3.cpp p3_functions_find.cpp)
endif()

set(P3_HEADERS
  p3_f90.hpp
  p3_ic_cases.hpp
  p3_constants.hpp
  p3_functions_upwind_impl.hpp
  p3_functions_table3_impl.hpp
  p3_functions.hpp
  p3_functions_find_impl.hpp
)

# link_directories(${SCREAM_TPL_LIBRARY_DIRS} ${SCREAM_LIBRARY_DIRS})
add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(p3 PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
# target_link_libraries(p3 scream_share ${SCREAM_TPL_LIBRARIES})

add_subdirectory(tests)
