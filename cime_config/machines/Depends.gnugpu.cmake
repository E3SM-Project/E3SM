list(APPEND ALLOW_INVALID_BOZ_LIST
  eam/src/control/cam_history.F90
  elm/src/biogeochem/MEGANFactorsMod.F90
)

if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  foreach(ITEM IN LISTS ALLOW_INVALID_BOZ_LIST)
    e3sm_add_flags("${ITEM}" "-fallow-invalid-boz") # avoids build error for integer, parameter :: gen_hash_key_offset = z'000053db'
  endforeach()
endif()

