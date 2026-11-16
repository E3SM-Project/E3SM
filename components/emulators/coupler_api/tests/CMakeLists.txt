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


# CouplerDriver Test
add_executable(test_coupler_driver test_coupler_driver.cpp)
target_link_libraries(test_coupler_driver PRIVATE Catch2 E3SM::coupler_api)
target_compile_definitions(test_coupler_driver PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/input-data")
add_test(NAME test_coupler_driver COMMAND test_coupler_driver)

# MOAB ADAPTER TEST
add_executable(test_moab_adapter
  moab_test_main.cpp
  test_moab_adapter.cpp
)

target_link_libraries(test_moab_adapter PRIVATE
  E3SM::coupler_api
  MOAB
  MPI::MPI_C
  Catch2
)

foreach(nranks IN ITEMS 1 2)
  add_test(
    NAME moab_adapter_${nranks}
    COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${nranks}
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_moab_adapter>
            ${MPIEXEC_POSTFLAGS}
            "[moab]"
  )

  set_tests_properties(moab_adapter_${nranks} PROPERTIES
    PROCESSORS ${nranks}
    TIMEOUT 60
  )
endforeach()

add_executable(test_moab_transport
  moab_test_main.cpp
  test_moab_transport.cpp
)

target_link_libraries(test_moab_transport PRIVATE
  E3SM::coupler_api
  MOAB
  MPI::MPI_C
  Catch2
)

add_test(
  NAME moab_transport_6
  COMMAND ${MPIEXEC_EXECUTABLE}
          ${MPIEXEC_NUMPROC_FLAG} 6
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_moab_transport>
          ${MPIEXEC_POSTFLAGS}
          "[transport]"
)

set_tests_properties(moab_transport_6 PROPERTIES
  PROCESSORS 6
  TIMEOUT 120
)
