set(LND_F90
  BandDiagonalMod.F90)

foreach(ITEM IN LISTS LND_F90)
  e3sm_add_flags("elm/src/biogeophys/${ITEM}" "-O0 -g")
endforeach()
