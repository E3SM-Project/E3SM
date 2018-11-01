function(add_pflotran_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${PFLOTRAN_LIBRARIES})
  endif()
endfunction(add_pflotran_library)

