enable_testing()
if(EMULATOR_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
# NOTE: for now, we assume an atmosphere emulator. Let's find a better place for this.
set(COMPONENT atm)
set(EMULATOR_NAME ace2)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran/emu_comp_mct.F90.in
  ${CMAKE_CURRENT_BINARY_DIR}/${COMPONENT}_comp_mct.F90
  @ONLY
)

set(CPL7_INTERFACE_CXX_SOURCES
  cpp/cpl7_component.cpp
  cpp/harness.cpp)

set(CPL7_INTERFACE_F90_SOURCES)

add_library(emulator_interface
          ${CPL7_INTERFACE_CXX_SOURCES}
          ${CPL7_INTERFACE_F90_SOURCES})

target_compile_features(emulator_interface PUBLIC cxx_std_20)

set_target_properties(emulator_interface PROPERTIES 
                      Fortran_MODULE_DIRECTORY 
                      ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(emulator_interface PUBLIC
                          ${CMAKE_CURRENT_SOURCE_DIR}
                          ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
                          ${CMAKE_CURRENT_SOURCE_DIR}/modules
)

target_link_libraries(emulator_interface PUBLIC E3SM::coupler_api emulator_common)
