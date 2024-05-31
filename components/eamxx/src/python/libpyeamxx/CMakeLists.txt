find_package(pybind11 REQUIRED)

pybind11_add_module(libpyeamxx libpyeamxx.cpp)
target_link_libraries(libpyeamxx PUBLIC scream_share scream_io p3)
