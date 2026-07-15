# Exception rules for ctest automatic error detection
# These are strings in the build output that ctest mistakenly
# interprets as errors
list(APPEND CTEST_CUSTOM_ERROR_EXCEPTION
  ".*error_handler.*"
  ".*spdlog/fmt/bundled/format.h.*"
)
