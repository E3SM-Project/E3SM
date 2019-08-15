function(add_emi_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${EMI_LIBRARIES})
endfunction(add_emi_executable)

