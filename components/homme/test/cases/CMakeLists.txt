# Setups for hand-run simulations. They are not ctest tests: CMake only copies their scripts and
# namelists in the build tree (see the CMakeLists.txt files in each case).
# Use 'make cases-execs' to build the executables they need.
# NOTE: the other folders in here (ASP, baroclinic, ...) contain standalone scripts, and are not
#       part of the build.

add_subdirectory (dcmip)
