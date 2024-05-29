# By now it should have already been found, but just for completeness
find_package(pybind11 REQUIRED)

pybind11_add_module(pyp3 pyp3.cpp)
target_link_libraries(pyp3 PUBLIC p3)
target_include_directories(pyp3 PRIVATE ${SCREAM_BASE_DIR}/src/share/pyshare)
