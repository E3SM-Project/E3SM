function(add_betr_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
    target_link_libraries(${lib} ${BETR_LIBRARIES})
  endif()
endfunction(add_betr_library)

