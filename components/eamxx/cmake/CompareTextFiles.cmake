# Compare two text files

# Get the source and target file paths from the command line arguments
# NOTE: when using cmake -P <script> args, CMAKE_ARGVn will contain
#  n=0: cmake
#  n=1: -P
#  n=2: script name
#  n=3,..: args
set(SRC_FILE "${CMAKE_ARGV3}")
set(TGT_FILE "${CMAKE_ARGV4}")

# Print the files being compared
message(STATUS "Comparing files: ${SRC_FILE} and ${TGT_FILE}")

# Check if the source file exists
if(NOT EXISTS ${SRC_FILE})
    message(FATAL_ERROR "Error: Source file not found: ${SRC_FILE}")
endif()

# Check if the target file exists
if(NOT EXISTS ${TGT_FILE})
    message(FATAL_ERROR "Error: Target file not found: ${TGT_FILE}")
endif()

# Compare the files (note: works only on UNIX)
execute_process(
    COMMAND diff ${SRC_FILE} ${TGT_FILE}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE diff_output
    ERROR_VARIABLE diff_error
)

# Check the result of the comparison
if(result)
    message(STATUS "Files differ:")
    message(STATUS "${diff_output}")
    message(FATAL_ERROR "Comparison failed.")
else()
    message(STATUS "Files are identical.")
endif()
