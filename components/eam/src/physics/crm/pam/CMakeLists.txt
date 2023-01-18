
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

set(PAM_DYCORE awfl)
set(PAM_MICRO  p3)
set(PAM_SGS    shoc)
set(PAM_RAD    forced)

add_subdirectory(external/pam_core)
add_subdirectory(external/physics)
add_subdirectory(external/dynamics)
add_subdirectory(external/externals/eigen)

set(YAKL_CXX_FLAGS -g )

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

target_compile_options(pam_driver PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${PAM_Fortran_FLAGS}>
)

target_link_libraries(pam_driver pam_core physics dynamics eigen)

