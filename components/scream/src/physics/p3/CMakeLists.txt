set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(P3_SRCS
  p3_f90.cpp
  p3_functions_f90.cpp
  p3_ic_cases.cpp
  micro_p3_iso_c.f90
  micro_p3_iso_f.f90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/micro_p3.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/micro_p3_utils.F90
  atmosphere_microphysics.cpp
  scream_p3_interface.F90
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND P3_SRCS
    p3_functions_upwind.cpp
    p3_functions_cloud_sed.cpp
    p3_functions_ice_sed.cpp
    p3_functions_ice_collection.cpp
    p3_functions_rain_sed.cpp
    p3_functions_table3.cpp
    p3_functions_table_ice.cpp
    p3_functions_dsd2.cpp
    p3_functions_find.cpp
    p3_functions_math.cpp
    p3_functions_update_prognostics.cpp
    p3_functions_autoconversion.cpp
    p3_functions_conservation.cpp
    p3_functions_back_to_cell_average.cpp
    p3_functions_cloud_rain_acc.cpp
    p3_functions_calc_rime_density.cpp
    p3_functions_cldliq_imm_freezing.cpp
    p3_functions_rain_imm_freezing.cpp
    p3_functions_droplet_self_coll.cpp
    p3_functions_update_prognostics.cpp
    p3_functions_impose_max_total_Ni.cpp)
endif()

# link_directories(${SCREAM_TPL_LIBRARY_DIRS} ${SCREAM_LIBRARY_DIRS})
add_library(p3 ${P3_SRCS})
target_include_directories(p3 PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(p3 SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(p3 PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(p3 scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
