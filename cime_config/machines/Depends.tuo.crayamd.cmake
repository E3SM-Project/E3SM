set(REDOPT_O0
	../driver-mct/main/prep_ocn_mod.F90
	../driver-mct/main/prep_atm_mod.F90
	../driver-mct/main/prep_lnd_mod.F90
	../driver-mct/main/prep_ice_mod.F90
	../driver-mct/main/prep_glc_mod.F90
	../driver-mct/main/prep_rof_mod.F90
	../driver-moab/main/prep_atm_mod.F90
	../driver-moab/main/prep_ocn_mod.F90
	components/eam/src/physics/rrtmgp/external/rte/kernels/mo_rte_solver_kernels.F90
	components/homme/src/share/compose/compose_slmm_departure_point.cpp
	components/homme/src/share/compose/compose_slmm_islmpi_adp.cpp
	components/homme/src/share/compose/compose_slmm_islmpi_calc_trajectory.cpp
	mpas-seaice/src/model_forward/mpas_seaice_core_interface_structs.F)

set(REDOPT_O1
	elm/src/biogeochem/FATESFireFactoryMod.F90
	../driver-mct/main/seq_hist_mod.F90
	../driver-mct/main/seq_io_mod.F90
	../driver-mct/main/seq_rest_mod.F90
	../driver-moab/main/seq_frac_mct.F90
	../driver-moab/main/seq_hist_mod.F90
	../driver-moab/main/seq_io_mod.F90
	../driver-moab/main/seq_rest_mod.F90
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
