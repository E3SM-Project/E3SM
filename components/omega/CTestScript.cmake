cmake_minimum_required(VERSION 3.20)

ctest_start(Experimental)

ctest_build(
  RETURN_VALUE BuildRetval
  CAPTURE_CMAKE_ERROR BuildResult
)

ctest_test(
  RETURN_VALUE TestRetval
  CAPTURE_CMAKE_ERROR TestResult
)

ctest_submit(
  RETURN_VALUE SubmitRetval
  CAPTURE_CMAKE_ERROR SubmitResult
)
