set (F90_SRC rrtmgp_interface.F90)
set (CXX_SRC
    rrtmgp_interface.cpp
    mo_load_coefficients.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/examples/all-sky/mo_load_cloud_coefficients.cpp
)

# Add library for interface code
add_library(rrtmgp_interface STATIC ${CXX_SRC})

# Set compile flags for RRTMGP sources (handles CUDA flags when needed)
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(rrtmgp_interface)

# Libraries to link
target_link_libraries(rrtmgp_interface spio rrtmgp yakl)

# Where to find includes
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/rte)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/rrtmgp)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/cloud_optics)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/extensions/fluxes_byband)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/examples)
target_include_directories(rrtmgp_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/cpp/examples/all-sky)

# Set fortran compiler flags
set_source_files_properties(${F90_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${FFLAGS}")
