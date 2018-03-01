function(add_mpp_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${MPP_LIBRARIES})
endfunction(add_mpp_executable)

