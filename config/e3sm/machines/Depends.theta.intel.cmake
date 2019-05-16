# Same as generic intel except remove fp-model fast

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    get_property(ITEM_FLAGS SOURCE ${ITEM} PROPERTY COMPILE_FLAGS)
    string(REPLACE "-fp-model fast" "" ITEM_FLAGS "${ITEM_FLAGS}")
    set_property(SOURCE ${ITEM} PROPERTY COMPILE_FLAGS "${ITEM_FLAGS}")
  endforeach()
endif()
