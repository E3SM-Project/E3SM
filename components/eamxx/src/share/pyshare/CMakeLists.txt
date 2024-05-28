set (pybind11_DIR /home/lbertag/.conda/envs/bartgol/lib/python3.12/site-packages)
find_package(pybind11 REQUIRED
  HINTS ${pybind11_DIR})

pybind11_add_module(pyscream pyscream.cpp)
target_link_libraries(pyscream PUBLIC scream_share scream_io)
