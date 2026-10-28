# enable_testing()
# if(EMULATOR_ENABLE_TESTS)
#   add_subdirectory(tests)
# endif()

set(CPL7_INTERFACE_CXX_SOURCES
  cpp/cpl7_component.cpp
  cpp/harness.cpp)

set(CPL7_INTERFACE_F90_SOURCES)

add_library(emulator_cpl7
          ${CPL7_INTERFACE_CXX_SOURCES}
          ${CPL7_INTERFACE_F90_SOURCES})

target_compile_features(emulator_cpl7 PUBLIC cxx_std_20)

set_target_properties(emulator_cpl7 PROPERTIES 
                      Fortran_MODULE_DIRECTORY 
                      ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(emulator_cpl7 PUBLIC
                          ${CMAKE_CURRENT_SOURCE_DIR}
                          ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
                          ${CMAKE_CURRENT_SOURCE_DIR}/modules
)


target_link_libraries(emulator_cpl7 PUBLIC E3SM::coupler_api emulator_common)
