foreach(ITEM IN LISTS FILES_NEED_OPENMP_FLAGS)
  e3sm_add_flags("${ITEM}" " -g -O2 -WF,-D_OPENMP -qsmp=omp -qoffload ")
endforeach()
