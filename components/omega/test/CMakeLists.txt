
# Add a test
add_executable(OMEGA_TEST testdummy.cpp)
target_link_libraries(OMEGA_TEST ${OMEGA_LIB_NAME})
add_test(NAME OMEGA_TEST COMMAND OMEGA_TEST)

