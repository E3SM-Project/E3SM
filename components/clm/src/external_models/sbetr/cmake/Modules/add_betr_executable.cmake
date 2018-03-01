function(add_betr_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${BETR_LIBRARIES})
endfunction(add_betr_executable)

