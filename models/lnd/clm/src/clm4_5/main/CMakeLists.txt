# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources clmtype.F90 clm_varctl.F90 clm_varcon.F90 clm_varpar.F90 decompMod.F90 surfrdUtilsMod.F90 subgridWeightsMod.F90)

sourcelist_to_parent(clm_sources)
