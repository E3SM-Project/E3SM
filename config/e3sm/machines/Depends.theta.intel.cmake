# Same as generic intel except remove fp-model fast

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    e3sm_remove_flags("${ITEM}" "-fp-model fast")
  endforeach()
endif()
