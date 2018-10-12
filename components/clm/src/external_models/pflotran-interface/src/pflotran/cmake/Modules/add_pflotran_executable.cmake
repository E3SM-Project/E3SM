function(add_pflotran_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${PFLOTRAN_LIBRARIES})
endfunction(add_pflotran_executable)

