# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources GridcellType.F90 LandunitType.F90 ColumnType.F90 PatchType.F90 accumulMod.F90 clm_varctl.F90 clm_varcon.F90 clm_varpar.F90 clm_varsur.F90 lnd2glcMod.F90 column_varcon.F90 decompMod.F90 initSubgridMod.F90 landunit_varcon.F90 surfrdUtilsMod.F90 subgridWeightsMod.F90 abortutils.F90)

sourcelist_to_parent(clm_sources)
