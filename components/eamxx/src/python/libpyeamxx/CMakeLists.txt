find_package(pybind11 REQUIRED)

pybind11_add_module(libpyeamxx_ext libpyeamxx.cpp)
target_link_libraries(libpyeamxx_ext PUBLIC scream_share scream_io p3 scream_rrtmgp shoc spa cld_fraction)
