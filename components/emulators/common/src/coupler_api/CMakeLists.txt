add_library(coupler_api
  coupler_types.F90
  coupler_types.hpp
)

set(_coupler_api_mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")

set_target_properties(coupler_api PROPERTIES
  Fortran_MODULE_DIRECTORY "${_coupler_api_mod_dir}"
)

target_include_directories(coupler_api
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${_coupler_api_mod_dir}>
)

target_compile_features(coupler_api
  PUBLIC
    cxx_std_20
)

add_library(E3SM::coupler_api ALIAS coupler_api)
