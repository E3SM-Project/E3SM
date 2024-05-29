find_package(pybind11 REQUIRED)

pybind11_add_module(pyscream pyscream.cpp)
target_link_libraries(pyscream PUBLIC scream_share scream_test_support scream_io)
