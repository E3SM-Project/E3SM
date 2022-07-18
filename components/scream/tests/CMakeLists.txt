# Some tests for checking testing works
set(EAMxx_tests_IC_FILE_72lev "screami_unit_tests_ne4np4L72_20220718.nc")
set(EAMxx_tests_IC_FILE_128lev "screami_unit_tests_ne2np4L128_20220718.nc")

add_subdirectory(generic)

if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  # Testing individual atm processes
  add_subdirectory(uncoupled)

  # Testing multiple atm processes coupled together
  add_subdirectory(coupled)
endif()
