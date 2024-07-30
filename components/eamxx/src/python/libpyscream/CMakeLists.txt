find_package(pybind11 REQUIRED)
find_package(mpi4py REQUIRED)

pybind11_add_module(pyscream_ext pyscream_ext.cpp)
target_link_libraries(pyscream_ext PUBLIC
  mpi4py
  scream_share
  scream_io
  diagnostics
  eamxx_physics
  scream_test_support
)
