set(REDOPT
  ../driver-mct/main/seq_io_mod.F90
  elm/src/biogeophys/BandDiagonalMod.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS REDOPT)
    e3sm_add_flags("${ITEM}" "-O1 -g")
    e3sm_remove_flags("${ITEM}" "-O2")
  endforeach()
endif()



