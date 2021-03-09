function(add_emi_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${EMI_LIBRARIES})
  endif()
endfunction(add_emi_library)

