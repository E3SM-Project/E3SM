
set(CPPDEFS "${CPPDEFS} -DMPAS_OPENMP_OFFLOAD")
foreach(ITEM IN LISTS MPAS_ADD_ACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-qsmp -qoffload")
endforeach()

