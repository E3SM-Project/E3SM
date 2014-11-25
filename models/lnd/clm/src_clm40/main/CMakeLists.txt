# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources domainMod.F90)

sourcelist_to_parent(clm_sources)
