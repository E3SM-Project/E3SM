
# Add tests

add_executable(testLogging.exe infra/LoggingTest.cpp)
target_link_libraries(testLogging.exe ${OMEGA_LIB_NAME})
target_include_directories(
	testLogging.exe
	PRIVATE
    "${OMEGA_SOURCE_DIR}/src/base"
	"${OMEGA_SOURCE_DIR}/src/infra"
)
add_test(NAME LOGGING_TEST COMMAND testLogging.exe)
set_tests_properties(LOGGING_TEST
    PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
