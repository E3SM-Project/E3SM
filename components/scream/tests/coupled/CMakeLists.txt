# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  if (SCREAM_DOUBLE_PRECISION)
    add_subdirectory(shoc_cld_p3_rrtmgp)
  endif()
endif()
