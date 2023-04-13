include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include(${EKAT_MACH_FILES_PATH}/srun.cmake)
