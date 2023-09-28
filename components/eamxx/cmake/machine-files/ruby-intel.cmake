include(${CMAKE_CURRENT_LIST_DIR}/ruby.cmake)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/mkl/mkl-2022.1.0/lib/intel64/ -qmkl" CACHE STRING "" FORCE)
