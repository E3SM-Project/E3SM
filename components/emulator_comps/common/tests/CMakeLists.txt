#==============================================================================
# Common Library Tests
#==============================================================================

message(STATUS "Configuring common library tests...")

# Test support library (used by all tests)
add_subdirectory(test_support)

# Unit tests
add_subdirectory(unit)

