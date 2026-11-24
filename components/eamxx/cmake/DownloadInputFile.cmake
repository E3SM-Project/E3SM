# Script invoked at build time (cmake -P) by the targets created in GetInputFile.
# Required vars: INPUT_ROOT, REL_PATH, CHECK_INPUT_SCRIPT
set(full_path ${INPUT_ROOT}/atm/${REL_PATH})
if (NOT EXISTS ${full_path})
  execute_process(COMMAND ${CHECK_INPUT_SCRIPT} ${INPUT_ROOT} atm/${REL_PATH}
    RESULT_VARIABLE status
    OUTPUT_VARIABLE output
    ERROR_VARIABLE  output)
  if (NOT status EQUAL 0)
    message(FATAL_ERROR "Could not download ${REL_PATH}, output was: ${output}")
  endif()
endif()
