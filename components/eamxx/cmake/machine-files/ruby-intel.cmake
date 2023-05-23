include(${CMAKE_CURRENT_LIST_DIR}/ruby.cmake)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/ -qmkl" CACHE STRING "" FORCE)
