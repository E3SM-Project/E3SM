if (DEFINED ENV{SCREAM_FAKE_ONLY})
  # Some fake tests for testing infrastructure
  add_subdirectory(generic)
else()
  # Testing individual atm processes
  add_subdirectory(uncoupled)

  # Testing multiple atm processes coupled together
  add_subdirectory(coupled)
endif()
