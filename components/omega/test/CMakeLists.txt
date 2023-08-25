# Omega Unit Tests

# Add data type test

set(_TestDataTypesName testDataTypes.exe)

add_executable(${_TestDataTypesName} base/DataTypesTest.cpp)
target_link_libraries(${_TestDataTypesName} ${OMEGA_LIB_NAME})
target_include_directories(
    ${_TestDataTypesName}
    PRIVATE
    "${OMEGA_SOURCE_DIR}/src/base"
    ${OMEGA_INCLUDE_DIRECTORIES}
)

target_link_directories(
    ${_TestDataTypesName}
    PRIVATE 
    ${OMEGA_LINK_DIRECTORIES}
)

target_link_options(
    ${_TestDataTypesName}
    PRIVATE 
    ${OMEGA_LINK_OPTIONS}
)

add_test(NAME DATA_TYPES_TEST COMMAND ${_TestDataTypesName})

# Add machine env test
# NOTE: disabled temporary due to not-launching mininum 8 ranks
#add_executable(testMachEnv.exe base/MachEnvTest.cpp)
#target_link_libraries(testMachEnv.exe ${OMEGA_LIB_NAME})
#target_include_directories(
#    testMachEnv.exe
#    PRIVATE
#    "${OMEGA_SOURCE_DIR}/src/base"
#)
#add_test(NAME MACHINE_ENV_TEST COMMAND testMachEnv.exe)

set(_TestLoggingName testLogging.exe)

add_executable(${_TestLoggingName} infra/LoggingTest.cpp)
target_link_libraries(${_TestLoggingName} ${OMEGA_LIB_NAME})
target_include_directories(
    ${_TestLoggingName}
	PRIVATE
    "${OMEGA_SOURCE_DIR}/src/base"
	"${OMEGA_SOURCE_DIR}/src/infra"
)

add_test(NAME LOGGING_TEST COMMAND testLogging.exe)

# set test properties
set_tests_properties(
	DATA_TYPES_TEST
#	MACHINE_ENV_TEST
	LOGGING_TEST 
    PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
