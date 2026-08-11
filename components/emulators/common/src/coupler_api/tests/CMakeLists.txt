enable_testing()

# set(COUPLER_API_TEST_SOURCES
#   test_coupler_api.F90
# )
#
# add_executable(test_coupler_api ${COUPLER_API_TEST_SOURCES})
#
# target_link_libraries(test_coupler_api PRIVATE
#   E3SM::coupler_api
# )
# add_test(NAME coupler_api_basic
#          COMMAND test_coupler_api)

add_executable(test_ekat_yaml
    test_yaml_fields.cpp
)

target_link_libraries(test_ekat_yaml
    PRIVATE
        Catch2
        E3SM::coupler_api
)
target_compile_definitions(test_ekat_yaml
    PRIVATE
        TEST_FIELDS_YAML="${CMAKE_CURRENT_SOURCE_DIR}/test_fields.yaml"
)

add_test(
    NAME test_ekat_yaml
    COMMAND test_ekat_yaml
)
