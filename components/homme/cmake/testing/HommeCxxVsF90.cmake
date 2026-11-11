# Compare the outputs of a F90 test and the corresponding C++ test. Usage:
#   cmake -DCPRNC=<cprnc> -DF90_DIR=<dir> -DF90_TEST_NAME=<name> -DCXX_DIR=<dir> -DCXX_TEST_NAME=<name>
#         -DFILES=<f1,f2,...> -P HommeCxxVsF90.cmake
# 1. The netcdf output files <dir>/movies/<f> must be bit-for-bit identical
# 2. For theta tests only, the diagnostics written in the log (<test>.log) must match
#    (for preqx, cxx runs do not push the results back to F90, so the cxx diagnostics are meaningless).

cmake_minimum_required (VERSION 3.16)

string (REPLACE "," ";" FILES "${FILES}")

# Run cprnc on each of the output files against the target
foreach (NC_FILE_NAME IN LISTS FILES)
  execute_process (COMMAND ${CPRNC} ${F90_DIR}/movies/${NC_FILE_NAME} ${CXX_DIR}/movies/${NC_FILE_NAME}
                   RESULT_VARIABLE ERROR_CODE
                   OUTPUT_VARIABLE cprnc_stdout
                   ERROR_VARIABLE cprnc_stderr)

  # Search stdout for "IDENTICAL". -1 means it does not exist. Use REVERSE on
  # the off chance that makes it faster, since "IDENTICAL", if it exists, is
  # near the end of the output.
  string (FIND "${cprnc_stdout}" "IDENTICAL" identical_pos REVERSE)

  if (identical_pos EQUAL -1)
    message ("${cprnc_stdout}")
    message (FATAL_ERROR "Test did not succeed. Netcdf outputs '${NC_FILE_NAME}' differ.")
  endif()
endforeach()

message ("CXX and F90 netcdf outputs match.\n")

# stdout part
string (FIND "${CXX_TEST_NAME}" "preqx" preqxtest)
if (preqxtest EQUAL -1)
  set (GREP_ARGS "  u     =")
  foreach (s
    "  v     ="
    "  w     ="
    "  T     ="
    "  mu    ="
    "  vTh_dp="
    "  dz(m) ="
    "  dp    ="
    "    fu  ="
    "    fv  ="
    "    ft  ="
    "    fq1 ="
    "min dz/w"
    " qv(  1)="
    "    TBOT="
    "      ps="
    "      M ="
    "  dry M ="
    "dQ1/dt"
    "Q  1,Q"
    "(Q-Q0)"
    "pmax_")
    set (GREP_ARGS "${GREP_ARGS}\\|${s}")
  endforeach()

  message ("grep_args are ${GREP_ARGS} -----------------------------------")

  execute_process (COMMAND grep ${GREP_ARGS} ${F90_DIR}/${F90_TEST_NAME}.log
                   RESULT_VARIABLE ERROR_CODE
                   OUTPUT_VARIABLE grep_f90_stdout
                   ERROR_VARIABLE  grep_f90_stderr)
  execute_process (COMMAND grep ${GREP_ARGS} ${CXX_DIR}/${CXX_TEST_NAME}.log
                   RESULT_VARIABLE ERROR_CODE
                   OUTPUT_VARIABLE grep_cxx_stdout
                   ERROR_VARIABLE  grep_cxx_stderr)

  if (NOT "${grep_f90_stdout}" STREQUAL "${grep_cxx_stdout}")
    message ("Test did not succeed. CXX and F90 diagnostics outputs differ:\nf90:\n${grep_f90_stdout}\ncxx:\n${grep_cxx_stdout}")
    message (FATAL_ERROR "Test did not succeed. Stdout outputs differ.")
  else ()
    message ("CXX and F90 diagnostics outputs match.\n")
  endif()
endif()
