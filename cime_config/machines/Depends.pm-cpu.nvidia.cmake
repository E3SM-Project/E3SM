list(APPEND REDUCE_OPT_LIST
  homme/src/share/derivative_mod_base.F90
)

# Can use this flag to avoid internal compiler error for this file (with nvidia/21.11)
if (NOT DEBUG)
  foreach(ITEM IN LISTS REDUCE_OPT_LIST)
    e3sm_add_flags("${ITEM}" " -Mnovect")
  endforeach()
endif()
