if (DEFINED ENV{SCREAM_FAKE_ONLY})
  set(SHARE_SRC
    scream_config.cpp
    scream_session.cpp
    util/scream_test_session.cpp
    util/scream_time_stamp.cpp
    )
else()
  set(SHARE_SRC
    scream_config.cpp
    scream_session.cpp
    util/scream_test_session.cpp
    util/scream_time_stamp.cpp
    )
endif()

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_SRC_DIR} ${SCREAM_BINARY_DIR}/src)
target_link_libraries(scream_share PUBLIC ekat)
set_target_properties(scream_share PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_compile_options(scream_share PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)
# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_librarires(scream_share PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(scream_share PUBLIC "--remove-duplicate-link-files")
  endif()
endif()
