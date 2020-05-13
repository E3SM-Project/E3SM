foreach(ITEM IN LISTS FILES_NEED_OPENMP_FLAGS)
  set_property(SOURCE ${ITEM} APPEND_STRING PROPERTY COMPILE_FLAGS " -g -O2 -WF,-D_OPENMP -qsmp=omp -qoffload ")
endforeach()
