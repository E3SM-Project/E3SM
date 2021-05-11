# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

set(genf90_files dynVarMod.F90.in dynVarTimeInterpMod.F90.in dynVarTimeUninterpMod.F90.in)

process_genf90_source_list("${genf90_files}" ${CMAKE_CURRENT_SOURCE_DIR} clm_genf90_sources)

sourcelist_to_parent(clm_genf90_sources)

list(APPEND clm_sources "${clm_genf90_sources}")

list(APPEND clm_sources dynPriorWeightsMod.F90 dynTimeInfoMod.F90 dynLandunitAreaMod.F90 dynInitColumnsMod.F90)

sourcelist_to_parent(clm_sources)
