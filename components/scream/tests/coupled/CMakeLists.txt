# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  add_subdirectory(p3_shoc)
endif()
