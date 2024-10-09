# For this file, fixes non-BFB behavior of stealth feature on pm-cpu with -O2
set(NOOPT
  eam/src/physics/cam/zm_conv.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS NOOPT)
    e3sm_deoptimize_file("${ITEM}")
  endforeach()
endif()

# On pm-cpu (and muller-cpu), with gcc-native/12.3, we see hang with DEBUG runs of certain tests.
# https://github.com/E3SM-Project/E3SM/issues/6516
# Currently, we have pm-cpu using gcc/12.2.0 which does not have this issue, but using muller-cpu to test 12.3
# Turning off -O0 for these 2 files (by adding -O) at least avoids hang and will produce FPE in HOMME code
if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 12.3)
  if (DEBUG)

    set(ADJUST
      eam/src/dynamics/se/inidat.F90
      eam/src/dynamics/se/dyn_comp.F90
      )

    foreach(ITEM IN LISTS ADJUST)
      e3sm_add_flags("${ITEM}" "-O")
      #e3sm_add_flags("${ITEM}" "-DNDEBUG -O")
    endforeach()

  endif()
endif()
