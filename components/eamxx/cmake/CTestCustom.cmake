# Exception rules for ctest automatic error detection
# These are strings in the build output that ctest mistakenly
# interprets as errors
list(APPEND CTEST_CUSTOM_ERROR_EXCEPTION
  ".*error_handler.*"
  ".*spdlog/fmt/bundled/format.h.*"
)

# nvcc emits a pseudo source file name in debug line info; it does not exist on disk
list(APPEND CTEST_CUSTOM_COVERAGE_EXCLUDE
  ".*nvcc_internal.*"
)
