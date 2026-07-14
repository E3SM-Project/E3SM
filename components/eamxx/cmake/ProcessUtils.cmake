# For process scripts (invoked with cmake -P), it is critical to minimize
# external dependencies because the usual include paths are not available.

###############################################################################
function(GetScriptPositionalArguments out_var)
###############################################################################
  # Useful for getting the "real" arguments to a cmake script (not including
  # cmake flags)
  set(real_args "")
  set(found_script FALSE)

  # Calculate total arguments passed to the command line
  math(EXPR max_index "${CMAKE_ARGC} - 1")

  foreach(i RANGE ${max_index})
    set(current_arg "${CMAKE_ARGV${i}}")

    # Everything after the .cmake file is a real positional argument
    if(found_script)
      list(APPEND real_args "${current_arg}")
    endif()

    # Mark that we found the script file itself
    if(current_arg MATCHES "\\.cmake$")
      set(found_script TRUE)
    endif()
  endforeach()

  # Propagate the filtered list back to the caller's scope
  set(${out_var} "${real_args}" PARENT_SCOPE)
endfunction()
