
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

set(SCREAM_LIB_ONLY TRUE)
set(SCREAM_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../..)

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

string(FIND "${CAM_CONFIG_OPTS}" "-pam_dycor awfl" USE_AWFL)
string(FIND "${CAM_CONFIG_OPTS}" "-pam_dycor spam" USE_SPAM)

if (NOT USE_AWFL EQUAL -1) 
  set(PAM_DYCORE  awfl)
endif()

if (NOT USE_SPAM EQUAL -1) 
  set(PAM_DYCORE  spam)
endif()

set(PAM_MICRO   p3)
set(PAM_SGS     shoc)
set(PAM_RAD     forced)
set(PAM_SCREAM_USE_CXX True)

# removing this once the cime cmake_macros aren't using obsolete
# Kokkos settings like KOKKOS_OPTIONS.
if (KOKKOS_OPTIONS)
   unset(KOKKOS_OPTIONS)
endif()

if (PAM_SCREAM_USE_CXX)
  message(STATUS "*** building WITH C++ scream interface ***")
  target_compile_definitions(pam_driver PUBLIC "-DP3_CXX -DSHOC_CXX")
endif()

if ("${PAM_DYCORE}" STREQUAL "spam")
  set(PAMC_MODEL  extrudedmodel)
  set(PAMC_HAMIL  man)  # man = anelastic / mce_rho = moist compressible
  set(PAMC_THERMO constkappavirpottemp)
  set(PAMC_IO     none)
endif()

if (NOT PAM_SCREAM_USE_CXX)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
     set(PAM_Fortran_FLAGS -g )
  else()
     set(PAM_Fortran_FLAGS -g -ffree-line-length-none -ffixed-line-length-none )
  endif()
endif()


if (PAM_SCREAM_USE_CXX)
  set(EKAT_ENABLE_YAML_PARSER OFF CACHE BOOL "" FORCE)
  add_subdirectory(external/physics/scream_cxx_p3_shoc)
  add_subdirectory(external/physics/scream_cxx_interfaces)
  # add these here to allow driver to intialize P3
  include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3       )
  include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3/impl  )
  include_directories(${SCREAM_HOME}/components/eamxx/src/physics/share    )
  include_directories(${SCREAM_HOME}/components/eamxx/src                  )
endif()

add_subdirectory(external/pam_core)
add_subdirectory(external/physics)
add_subdirectory(external/dynamics)

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
      ${CMAKE_CURRENT_SOURCE_DIR}/external/dynamics/${PAM_DYCORE}
      ${CMAKE_CURRENT_BINARY_DIR}/external/dynamics/${PAM_DYCORE}
      )

if (PAM_SCREAM_USE_CXX)
target_include_directories(pam_driver PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/external/physics/scream_cxx_p3_shoc
      ${CMAKE_CURRENT_BINARY_DIR}/external/physics/scream_cxx_p3_shoc
      ${CMAKE_CURRENT_SOURCE_DIR}/external/physics/scream_cxx_interfaces
      ${CMAKE_CURRENT_BINARY_DIR}/external/physics/scream_cxx_interfaces
      )
endif()

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(pam_driver)

if (NOT PAM_SCREAM_USE_CXX)
  set_target_properties(pam_driver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
  )
endif()

add_compile_options($<$<COMPILE_LANGUAGE:CXX>: -g >)

target_compile_options(pam_driver PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${PAM_Fortran_FLAGS}>
)

if (PAM_SCREAM_USE_CXX)
  target_link_libraries(pam_driver pam_core physics dynamics pam_scream_cxx_interfaces ekat p3 shoc physics_share scream_share)
else()
  target_link_libraries(pam_driver pam_core physics dynamics )
endif()


