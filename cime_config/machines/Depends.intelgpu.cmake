# https://github.com/E3SM-Project/E3SM/issues/8140
if (DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90" "-init=zero")
  endif()
endif()

# https://github.com/E3SM-Project/E3SM/issues/8036
# Turn off vectorization in one source file to avoid issue with Intel compiler
if (NOT DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/optics/quickbeam_optics.F90" "-no-vec")
  endif()
endif()

# Turn off fast-math for ifx 2025.3+
if (NOT DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.3")
    # ifx: clear out -freciprocal-math/unsafe-fp-math from the default fast model; source first clears them
    string(PREPEND CMAKE_Fortran_FLAGS " -fp-model=source")
  endif()
endif()

