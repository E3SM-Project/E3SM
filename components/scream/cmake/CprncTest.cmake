# Dump vars that are expected to be equal and those expected to be different
execute_process (
  COMMAND @CPRNC_BINARY@ @SRC_FILE@ @TGT_FILE@
  RESULT_VARIABLE cprnc_result
  OUTPUT_VARIABLE cprnc_output
  ERROR_VARIABLE  cprnc_stderr)

if (NOT cprnc_result EQUAL 0)
  string (CONCAT msg
          "Command\n"
          "  'cprnc @SRC_FILE@ @TGT_FILE@'\n"
          "returned non-zero exit code.\n")
  message ("${msg}")
  message (FATAL_ERROR "Aborting.")
endif()

# Search output for "IDENTICAL". -1 means it does not exist. Use REVERSE on
# the off chance that makes it faster, since "IDENTICAL", if it exists, is
# near the end of the output.
string (FIND "${cprnc_output}" "IDENTICAL" identical_pos REVERSE)

if (identical_pos EQUAL -1)
  string (CONCAT msg
          "Command\n"
          "  'cprnc @SRC_FILE@ @TGT_FILE@'\n"
          "reported differences between the files. Here's the output from cprnc:\n")
  message ("${msg}")
  message ("${cprnc_output}")
  message (FATAL_ERROR "Aborting.")
endif ()
