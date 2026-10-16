set(REDOPT_O0
	../driver-mct/main/prep_ocn_mod.F90
	../driver-mct/main/prep_atm_mod.F90
	../driver-mct/main/prep_lnd_mod.F90
	../driver-mct/main/prep_ice_mod.F90
	../driver-mct/main/prep_glc_mod.F90
	../driver-mct/main/prep_rof_mod.F90
	../driver-moab/main/prep_atm_mod.F90
	../driver-moab/main/prep_ocn_mod.F90
	mpas-seaice/src/model_forward/mpas_seaice_core_interface_structs.F)

set(REDOPT_O1
	elm/src/biogeochem/FATESFireFactoryMod.F90
	elm/src/external_models/sbetr/src/driver/shared/BeTRSimulation.F90
	${CMAKE_BINARY_DIR}/framework/mpas_io_streams.f90
	${CMAKE_BINARY_DIR}/framework/mpas_forcing.f90
	${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_forcing.f90
	${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_initialize.f90
	${CMAKE_BINARY_DIR}/core_seaice/analysis_members/mpas_seaice_regional_statistics.f90
	${CMAKE_BINARY_DIR}/core_seaice/analysis_members/mpas_seaice_time_series_stats.f90
	${CMAKE_BINARY_DIR}/core_seaice/analysis_members/mpas_seaice_conservation_check.f90)

set(REDOPT_BINARY_O0
	${CMAKE_BINARY_DIR}/framework/mpas_stream_manager.f90
	${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_icepack.f90
	${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver.f90
	${CMAKE_BINARY_DIR}/core_seaice/model_forward/mpas_seaice_core_interface_structs.f90)

if (NOT DEBUG)
	foreach(ITEM IN LISTS REDOPT_O0)
		e3sm_add_flags("${ITEM}" "-O0")
	endforeach()

	foreach(ITEM IN LISTS REDOPT_O1)
		e3sm_add_flags("${ITEM}" "-O1")
	endforeach()

	foreach(ITEM IN LISTS REDOPT_BINARY_O0)
		e3sm_add_flags("${ITEM}" "-O0")
	endforeach()
endif()
