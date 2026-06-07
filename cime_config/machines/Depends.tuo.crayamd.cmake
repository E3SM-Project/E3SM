set(REDOPT_O0
	../driver-mct/main/prep_ocn_mod.F90
	../driver-mct/main/prep_atm_mod.F90
	../driver-mct/main/prep_lnd_mod.F90
	../driver-mct/main/prep_ice_mod.F90
	../driver-mct/main/prep_glc_mod.F90
	../driver-mct/main/prep_rof_mod.F90)

set(REDOPT_O1
	elm/src/biogeochem/FATESFireFactoryMod.F90
	elm/src/external_models/sbetr/src/driver/shared/BeTRSimulation.F90)

# Speed-oriented debug-symbols configuration.
# - DEBUG=TRUE in CIME implies very conservative compilation (`-O0 -g`) which can
#   significantly increase compile/link walltime for large C++/HIP builds.
# - This file is sourced by the CMake configure step, so we can safely adjust the
#   effective flags here without editing the machine macro files.
#
# Strategy:
# - In DEBUG builds, keep debug symbols but use `-Og` and reduced debug info for C/CXX.
# - In non-DEBUG builds, keep release optimizations but add minimal debug symbols.
function(_tuo_adjust_debug_symbols_flags)
	if (DEBUG)
		foreach(LANG IN ITEMS C CXX)
			set(VAR "CMAKE_${LANG}_FLAGS_DEBUG")
			if (DEFINED ${VAR})
				set(FLAGS "${${VAR}}")
				string(REGEX REPLACE "(^|[ \t])-O0([ \t]|$)" " -Og " FLAGS " ${FLAGS} ")
				string(REGEX REPLACE "(^|[ \t])-g([0-9]|line-tables-only)?([ \t]|$)" " " FLAGS "${FLAGS}")
				string(APPEND FLAGS " -g1")
				string(STRIP "${FLAGS}" FLAGS)
				set(${VAR} "${FLAGS}" PARENT_SCOPE)
			endif()
		endforeach()
	else()
		foreach(LANG IN ITEMS C CXX)
			set(VAR "CMAKE_${LANG}_FLAGS_RELEASE")
			if (DEFINED ${VAR})
				set(FLAGS "${${VAR}}")
				string(REGEX REPLACE "(^|[ \t])-g([0-9]|line-tables-only)?([ \t]|$)" " " FLAGS " ${FLAGS} ")
				string(APPEND FLAGS " -g1")
				string(STRIP "${FLAGS}" FLAGS)
				set(${VAR} "${FLAGS}" PARENT_SCOPE)
			endif()
		endforeach()

		# Cray Fortran: keep standard debug symbols in release builds (no -g1 equivalent).
		if (DEFINED CMAKE_Fortran_FLAGS_RELEASE)
			set(FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
			string(REGEX REPLACE "(^|[ \t])-g([ \t]|$)" " " FLAGS " ${FLAGS} ")
			string(APPEND FLAGS " -g")
			string(STRIP "${FLAGS}" FLAGS)
			set(CMAKE_Fortran_FLAGS_RELEASE "${FLAGS}" PARENT_SCOPE)
		endif()
	endif()
endfunction()

# This file is included after `${CASEROOT}/Macros.cmake` (via `common_setup.cmake`),
# so the flag variables are already populated and we can adjust them immediately.
_tuo_adjust_debug_symbols_flags()

if (NOT DEBUG)
	foreach(ITEM IN LISTS REDOPT_O0)
		e3sm_add_flags("${ITEM}" "-O0")
	endforeach()

	foreach(ITEM IN LISTS REDOPT_O1)
		e3sm_add_flags("${ITEM}" "-O1")
	endforeach()
endif()

