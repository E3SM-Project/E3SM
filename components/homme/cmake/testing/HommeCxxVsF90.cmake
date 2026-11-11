# Compare the outputs of a F90 test and the corresponding C++ test. Usage:
#   cmake -DCPRNC=<cprnc> -DF90_DIR=<dir> -DCXX_DIR=<dir> -DFILES=<f1,f2,...> -P HommeCxxVsF90.cmake
# The netcdf output files <F90_DIR>/movies/<f> and <CXX_DIR>/movies/<f> must be bit-for-bit identical.

cmake_minimum_required (VERSION 3.16)

string (REPLACE "," ";" FILES "${FILES}")

set (nfails 0)
foreach (f IN LISTS FILES)
  set (f90_file ${F90_DIR}/movies/${f})
  set (cxx_file ${CXX_DIR}/movies/${f})
  foreach (file IN ITEMS ${f90_file} ${cxx_file})
    if (NOT EXISTS ${file})
      message (FATAL_ERROR "FAIL: file ${file} not found")
    endif()
  endforeach()
  execute_process (COMMAND ${CPRNC} ${f90_file} ${cxx_file}
                   OUTPUT_VARIABLE cprnc_out
                   ERROR_VARIABLE  cprnc_out)
  if (cprnc_out MATCHES "diff_test: the two files seem to be IDENTICAL")
    message ("PASS: ${f}: the F90 and C++ outputs are IDENTICAL")
  else()
    message ("${cprnc_out}")
    message ("FAIL: ${f}: the F90 and C++ outputs differ")
    math (EXPR nfails "${nfails} + 1")
  endif()
endforeach()

if (nfails GREATER 0)
  message (FATAL_ERROR "${nfails} file(s) differ between the F90 and C++ runs")
endif()
