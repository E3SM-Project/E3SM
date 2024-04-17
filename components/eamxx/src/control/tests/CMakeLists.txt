# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_ONLY_GENERATE_BASELINES})
  include (ScreamUtils)

  # Unit test the ad
  CreateUnitTest(ad_ut "ad_tests.cpp"
    LIBS scream_control
    LABELS driver)

  # Copy yaml input file to run directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ad_tests.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/ad_tests.yaml COPYONLY)

endif()
