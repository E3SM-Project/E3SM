find_package(pybind11 REQUIRED)

pybind11_add_module(pyeamxx pyeamxx.cpp)
target_link_libraries(pyeamxx PUBLIC scream_share scream_io diagnostics eamxx_physics)
