# Detect the installed nanobind package and import it into CMake
find_package(Python COMPONENTS Interpreter Development REQUIRED)
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind REQUIRED)
find_package(mpi4py REQUIRED)

nanobind_add_module(pyeamxx_ext pyeamxx_ext.cpp)
target_link_libraries(pyeamxx_ext PUBLIC
  mpi4py
  scream_share
  scream_io
  diagnostics
  eamxx_physics
  scream_test_support
)
