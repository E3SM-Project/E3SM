include(${CMAKE_CURRENT_LIST_DIR}/quartz.cmake)
set(CMAKE_CXX_FLAGS "-w" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/tce/packages/gcc/gcc-8.3.1/rh/lib/gcc/x86_64-redhat-linux/8/" CACHE STRING "" FORCE)