set (pybind11_DIR ${SCREAM_BASE_DIR}/../../externals/pybind11)
find_package(pybind11 REQUIRED
  HINTS ${pybind11_DIR})

pybind11_add_module(pyscream pyscream.cpp)
target_link_libraries(pyscream PUBLIC scream_share scream_io)
