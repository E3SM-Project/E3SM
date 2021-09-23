if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  # Testing individual atm processes
  add_subdirectory(uncoupled)

  # Testing multiple atm processes coupled together
  add_subdirectory(coupled)
endif()

# Some fake tests for testing infrastructure
add_subdirectory(generic)
