cmake_minimum_required (VERSION 2.6.2)

project (pflotran Fortran)

# Version number
set (PFLOTRAN_VERSION_MAJOR 1)
set (PFLOTRAN_VERSION_MINOR 0)

#configure_file (
#  "${PROJECT_SOURCE_DIRECTORY}/pflotran_config.h.in"
#  "${PROJECT_BINARY_DIRECTORY}/pflotran_config.h"
#  )

include ($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/conf/PETScConfig.cmake)

#include ($ENV{PETSC_DIR}/conf/variables)
#include ($ENV{PETSC_DIR}/conf/rules)

include_directories(".")
include_directories($ENV{PETSC_DIR}/include)
#include_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)
include_directories(${PETSC_PACKAGE_INCLUDES})

list (APPEND SRC
	input.F90
	string.F90
	water_eos.F90
	gaseos_mod.F90
	co2eos.F90
	co2_sw_rtsafe.F90 
	co2_span_wagner.F90
	spline.F90
	co2_span_wagner_spline.F90
	co2_sw.F90
	utility.F90
	option.F90
	geometry.F90
	units.F90
	logging.F90
	hdf5_aux.F90
	global_aux.F90
	material_aux.F90
	matrix_block_aux.F90
	reactive_transport_aux.F90
	reaction_aux.F90
	thc_aux.F90
	richards_aux.F90
	general_aux.F90
	mphase_aux.F90
	miscible_aux.F90
	immis_aux.F90
	mfd_aux.F90
	flash2_aux.F90
	auxilliary.F90
	field.F90
	debug.F90
	waypoint.F90
	region.F90
	fluid.F90
	time_series.F90
	dataset_aux.F90
	dataset.F90
	saturation_function.F90
	material.F90
	pckr_mod.F90
	strata.F90
	connection.F90
	observation.F90
	solver.F90
	cf90bridge.F90
	unstructured_cell.F90
	unstructured_grid.F90
	structured_grid.F90
	condition.F90
	coupler.F90
	grid.F90
	matrix_buffer.F90
	hydrostatic.F90
	saturation.F90
	patch.F90
	mfd.F90
	level.F90
	discretization.F90
	database.F90
	database_hpt.F90
	velocity.F90
	realization.F90
	hdf5.F90
	global.F90
	mphase.F90
	immis.F90
	miscible.F90
	flash2.F90
	transport.F90
	reaction.F90
	reactive_transport.F90
	thc.F90
	richards.F90
	general.F90
	condition_control.F90
	general_grid.F90
	checkpoint.F90
	convergence.F90
	timestepper.F90
	simulation.F90
	stochastic_aux.F90
	stochastic.F90
	surface_field.F90
	surface_flow.F90
	surface_material.F90
	surface_realization.F90
	output.F90
	init.F90
	pflotran.F90
)

add_executable(pflotran ${SRC})

target_link_libraries(pflotran ${PETSC_PACKAGE_LIBS})
