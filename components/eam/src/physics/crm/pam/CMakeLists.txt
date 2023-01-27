
set(PAM_DRIVER_SRC 
    pam_driver.F90  
    pam_driver.cpp 
    pam_feedback.h  
    pam_radiation.h  
    pam_statistics.h  
    pam_state.h  
    params.F90)

add_library(pam_driver 
           ${PAM_DRIVER_SRC})

if (DEFINED PAM_CUDA)
  if (${USE_CUDA})
    # PAM will be CUDA-linked with device symbols resolved at library creation
    set_target_properties(pam_driver
                          PROPERTIES
                          LINKER_LANGUAGE
                          CUDA
                          CUDA_SEPARABLE_COMPILATION OFF
                          CUDA_RESOLVE_DEVICE_SYMBOLS ON
                         )
  endif()
endif()

set(PAM_DYCORE awfl)
set(PAM_MICRO  p3)
set(PAM_SGS    shoc)
set(PAM_RAD    forced)

add_subdirectory(external/pam_core)
add_subdirectory(external/physics)
add_subdirectory(external/dynamics)
add_subdirectory(external/externals/eigen)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   set(PAM_Fortran_FLAGS -g -real-size 64)
else()
   set(PAM_Fortran_FLAGS -g -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -ffixed-line-length-none )
endif()

target_include_directories(pam_driver PUBLIC 
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
      ${CMAKE_CURRENT_SOURCE_DIR}/external/pam_core
      ${CMAKE_CURRENT_SOURCE_DIR}/external/pam_core/modules
      ${CMAKE_CURRENT_BINARY_DIR}/external/pam_core/modules
      ${CMAKE_CURRENT_SOURCE_DIR}/external/pam_core/pam_interface
      ${CMAKE_CURRENT_BINARY_DIR}/external/pam_core/pam_interface
      ${CMAKE_CURRENT_SOURCE_DIR}/external/physics/micro/p3
      ${CMAKE_CURRENT_BINARY_DIR}/external/physics/micro/p3
      ${CMAKE_CURRENT_SOURCE_DIR}/external/physics/sgs/shoc
      ${CMAKE_CURRENT_BINARY_DIR}/external/physics/sgs/shoc
      ${CMAKE_CURRENT_SOURCE_DIR}/external/physics/radiation/forced
      ${CMAKE_CURRENT_BINARY_DIR}/external/physics/radiation/forced
      ${CMAKE_CURRENT_SOURCE_DIR}/external/dynamics/awfl
      ${CMAKE_CURRENT_BINARY_DIR}/external/dynamics/awfl
      #${YAMLCPP_SOURCE_DIR}/include
      #${YAMLCPP_BINARY_DIR}/include
      )

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(pam_driver)

set_target_properties(pam_driver PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

add_compile_options($<$<COMPILE_LANGUAGE:CXX>: -g >)

target_compile_options(pam_driver PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${PAM_Fortran_FLAGS}>
)

target_link_libraries(pam_driver pam_core physics dynamics eigen)

# set path for P3 lookup table data
#set (SCREAM_DATA_DIR ${SCREAM_INPUT_ROOT}/atm/scream CACHE PATH "" FORCE)

# Copy p3 lookup tables to local data directory
#FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
#CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
#               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
#CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
#               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)


#set (P3_TABLES
#  scream/tables/p3_lookup_table_1.dat-v4.1.1
#  scream/tables/mu_r_table_vals.dat${PRECISION_SUFFIX}
#  scream/tables/revap_table_vals.dat${PRECISION_SUFFIX}
#  scream/tables/vm_table_vals.dat${PRECISION_SUFFIX}
#  scream/tables/vn_table_vals.dat${PRECISION_SUFFIX}
#)

#include (ScreamUtils)
#foreach (file IN ITEMS ${P3_TABLES})
#  GetInputFile(${file})
#endforeach()
