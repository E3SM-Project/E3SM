# This utility checks output of the form
#   planar_conv>                          l2                    linf
#   planar_conv> Q 1  3.7612654153016604E-03  9.8620728873592828E-03
#   planar_conv> Q 2  3.9526251584022527E-03  8.2639849437505641E-03
# against caller-provided upper bounds on the l2 norms. Output of this sort
# appears on stdout at the end of some tracer transport tests. 'planar_conv' can
# be something else; it is caller-configurable. Usage:
#   cmake -DANCHOR=planar_conv -DFILENAME=out.txt -DUPPER_BOUNDS=3.9e-3,4.1e-3 -P HommeCheckErrorNorms.cmake
# See homme_add_error_norms_test in HommeTesting.cmake.

string (REPLACE "," ";" UPPER_BOUNDS "${UPPER_BOUNDS}")
execute_process (COMMAND grep "${ANCHOR}> Q" ${FILENAME} OUTPUT_VARIABLE grep_stdout)
separate_arguments (tokens NATIVE_COMMAND "${grep_stdout}")
list (LENGTH UPPER_BOUNDS nchk)
math (EXPR nchk "${nchk} - 1")
foreach (idx RANGE 0 ${nchk})
  math (EXPR tidx "3 + 5 * ${idx}")
  list (GET tokens ${tidx} e)
  list (GET UPPER_BOUNDS ${idx} ub)
  if (e GREATER ub)
    message (FATAL_ERROR "${e} should be <= ${ub} but is not")
  endif()
endforeach()
