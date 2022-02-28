# This macro adds interface flags to the input target, by probing the
# CACHE variables SCREAM_FLAGS and SCREAM_<LANG>_FLAGS. The former
# are added to the C, CXX, and Fortran interface compile options
# of the target, while the latter are only added to the interface
# compile options of the corresponding language. If LANG=LD, the
# macro sets interface link options on the target.

set (SCREAM_FLAGS "" CACHE STRING "Mach specific global flags")
set (SCREAM_F_FLAGS "" CACHE STRING "Mach specific Fortran flags")
set (SCREAM_C_FLAGS "" CACHE STRING "Mach specific C flags")
set (SCREAM_CXX_FLAGS "" CACHE STRING "Mach specific CXX flags")
set (SCREAM_LD_FLAGS "" CACHE STRING "Mach specific Linker flags")

macro (AddMachSpecificFlags targetName)
  if (SCREAM_FLAGS)
    target_compile_options (${targetName} PUBLIC ${SCREAM_FLAGS})
  endif()
  if (SCREAM_F_FLAGS)
    target_compile_options (${targetName} PUBLIC
      $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_F_FLAGS}>)
  endif()
  if (SCREAM_C_FLAGS)
    target_compile_options (${targetName} PUBLIC
      $<$<COMPILE_LANGUAGE:C>:${SCREAM_C_FLAGS}>)
  endif()
  if (SCREAM_CXX_FLAGS)
    target_compile_options (${targetName} PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:${SCREAM_CXX_FLAGS}>)
  endif()
  if (SCREAM_LD_FLAGS)
    target_link_options (${targetName} PUBLIC ${SCREAM_LD_FLAGS})
  endif()
endmacro()
