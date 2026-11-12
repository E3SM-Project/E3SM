if(NOT EMULATOR_ENABLE_PYTHON)
  message(FATAL_ERROR "The SamudrACE demo needs -DEMULATOR_ENABLE_PYTHON=ON")
endif()

add_executable(samudrace_demo samudrace_demo.cpp)
target_link_libraries(samudrace_demo PRIVATE
  emulator_interface emulator_common emulator_inference)
