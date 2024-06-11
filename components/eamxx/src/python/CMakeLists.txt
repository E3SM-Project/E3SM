find_package(pybind11 REQUIRED)
find_package(mpi4py REQUIRED)

pybind11_add_module(pyeamxx pyeamxx.cpp)
target_link_libraries(pyeamxx PUBLIC mpi4py scream_share scream_io diagnostics eamxx_physics)
