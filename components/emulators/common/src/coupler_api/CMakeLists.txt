cmake_minimum_required(VERSION 3.28)
project(EmualtorCoupler LANGUAGES CXX)
add_library(coupler_api
  coupler_types.F90
  field_registry.cpp
  field_registry_c_api.cpp
  coupler.cpp
  coupler_driver.cpp
  process.cpp
)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(coupler_api PRIVATE
        $<$<CONFIG:Debug>:
            -Og
            -g3
            -fno-omit-frame-pointer
            -fno-optimize-sibling-calls
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wconversion
            -Wsign-conversion
            -Wformat=2
            -Wnull-dereference
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
        >
    )
 target_compile_definitions(coupler_api PRIVATE
        $<$<CONFIG:Debug>:_GLIBCXX_ASSERTIONS>
    )

    # PUBLIC propagates sanitizer runtime linking to consuming executables,
    # including when my_library is static.
    target_link_options(coupler_api PUBLIC
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    )
endif()

set(_coupler_api_mod_dir "${CMAKE_CURRENT_BINARY_DIR}/modules")

set_target_properties(coupler_api PROPERTIES
  Fortran_MODULE_DIRECTORY "${_coupler_api_mod_dir}"
)

target_include_directories(coupler_api
  PUBLIC
    $<BUILD_INTERFACE:${_coupler_api_mod_dir}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(coupler_api PUBLIC ekat::YamlParser)

target_compile_features(coupler_api PUBLIC cxx_std_20)

add_library(E3SM::coupler_api ALIAS coupler_api)
if(EMULATOR_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
