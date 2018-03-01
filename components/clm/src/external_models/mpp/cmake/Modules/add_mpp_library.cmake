function(add_mpp_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${MPP_LIBRARIES})
  endif()
endfunction(add_mpp_library)

