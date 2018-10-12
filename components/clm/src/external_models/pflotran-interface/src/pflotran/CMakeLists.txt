cmake_minimum_required (VERSION 2.8.12)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
option(BUILD_SHARED_LIBS "Builds shared libraries [ON]." ON)

# For now, we disable shared libs on Macs.
if (APPLE)
  set(BUILD_SHARED_LIBS OFF)
endif()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

set(PETSC_ARCH $ENV{PETSC_ARCH})
set(PETSC_DIR $ENV{PETSC_DIR})
include($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/initial_cache_file.cmake)

# Make sure compilers are set. This must be done before enabling languages.
if (NOT CMAKE_C_COMPILER)
  if (NOT $ENV{CC} STREQUAL "")
    set(CMAKE_C_COMPILER $ENV{CC})
  else()
    set(CMAKE_C_COMPILER cc)
  endif()
endif()
if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS $ENV{CFLAGS})
endif()
if (NOT CMAKE_CXX_COMPILER)
  if (NOT $ENV{CXX} STREQUAL "")
    set(CMAKE_CXX_COMPILER $ENV{CXX})
  else()
    set(CMAKE_CXX_COMPILER c++)
  endif()
endif()
if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()
enable_language(C)
enable_language(CXX)
enable_language(Fortran)

# We declare the project here.
project (PFLOTRAN)

message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- CXX compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

if (BUILD_SHARED_LIBS)
  message("-- PFLOTRAN will be built as a shared library.")
else()
  message("-- PFLOTRAN will be built as a static library.")
endif()

set (PFLOTRAN_VERSION_MAJOR 1)
set (PFLOTRAN_VERSION_MINOR 0)
set (PFLOTRAN_VERSION_PATCH 0)
set (PFLOTRAN_VERSION "${PFLOTRAN_MAJOR_VERSION}.${PFLOTRAN_MINOR_VERSION}.${PFLOTRAN_PATCH_VERSION}")

# General C compiler flags.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")

  if (BUILD_SHARED_LIBS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
  endif()

  if (LINUX EQUAL 1)
    # Counter some of GCC's more recent stinginess on Linux.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
  endif()

elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration -fno-builtin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# Fortran compiler flags.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -ffree-line-length-0 -Wno-unused-variable -Wno-unused-parameter -DCPRGNU ")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRINTEL")
endif()


# Figure out the system type.
set(HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
endif ()

# Here we make sure CMake-installed binaries use the correct runpath, and
# that the path is not stripped during installation.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(PETSc)
if (NOT PETSC_FOUND)
  message(FATAL_ERROR "PETSc was not found.")
endif()
include_directories(${PETSC_INCLUDES})



# Include PETSc in the rpath.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PETSC_DIR}/${PETSC_ARCH}/lib")

# Other third-party libraries.
#add_subdirectory(3rd-party)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
#include_directories("${PROJECT_BINARY_DIR}/include")
#link_directories("${PROJECT_BINARY_DIR}/lib")

set(UTIL_SOURCES
	appleyard.F90
	generic.F90
	input_aux.F90
	logging.F90
	matrix_block_aux.F90
	option.F90
	option_flow.F90
	option_transport.F90
	output_aux.F90
	pflotran_constants.F90
	string.F90
	units.F90
	utility.F90
	toil_ims_derivs.F90
	variables.F90
)

set(EOS_SOURCES
	eos.F90
	eos_database.F90
	eos_gas.F90
	eos_oil.F90
	gas_eos_mod.F90
)

set(MODE_AUX_SOURCES
	auxiliary.F90
	auxvars_base.F90
	auxvars_flow.F90
	auxvars_flow_energy.F90
	auxvars_toil_ims.F90
	auxvars_towg.F90
	flash2_aux.F90
	general_aux.F90
	pm_base_aux.F90
	pm_toil_ims_aux.F90
	pm_towg_aux.F90
	immis_aux.F90
	miscible_aux.F90
	mphase_aux.F90
	pm_auxiliary.F90
	pm_base.F90
	pm_base_pointer.F90
	pm_bragflo.F90
	pm_general.F90
	pm_toil_ims.F90
	pm_towg.F90
	pm_flash2.F90
	pm_immis.F90
	pm_mphase.F90
	pm_miscible.F90
	pm_richards.F90
	pm_rt.F90
	pm_subsurface_flow.F90
	pm_surface_flow.F90
	pm_surface_th.F90
	pm_surface.F90
	pm_th.F90
	pm_ufd_decay.F90
	pm_waste_form.F90
	pm_wipp_flow.F90
	pm_wipp_srcsink.F90
	pm_ufd_biosphere.F90
	richards_aux.F90
	th_aux.F90
	wipp_flow_aux.F90
	inlinesurface_aux.F90
)

set(MODE_SOURCES
	bragflo.F90
	bragflo_common.F90
	flash2.F90
	general.F90
	general_common.F90
	toil_ims.F90
	towg.F90
	global.F90
	immis.F90
	miscible.F90
	mphase.F90
	reactive_transport.F90
	richards.F90
	richards_common.F90
	th.F90
	wipp_flow.F90
	wipp_flow_common.F90
	inlinesurface.F90
)

set(SHARD_MODE_AUX_SOURCES
	co2eos.F90
	co2_sw_rtsafe.F90
	eos_water.F90
	geometry.F90
	global_aux.F90
	material_aux.F90
	reactive_transport_aux.F90 
)

set(CHEM_SOURCES
	transport_constraint.F90
	reaction_database.F90
	reaction_database_aux.F90
	reaction_database_hpt.F90
	reaction_gas.F90
	reaction_gas_aux.F90
	reaction_immobile.F90
	reaction_immobile_aux.F90
	reaction_mineral.F90
	reaction_mineral_aux.F90
	reaction.F90
	reaction_aux.F90
	reaction_sandbox.F90
	reaction_sandbox_base.F90
)

set(GRID_SOURCES
	connection.F90
	communicator_base.F90
	communicator_structured.F90
	communicator_unstructured.F90
	discretization.F90
	dm_kludge.F90
	grid.F90
	grid_structured.F90
	grid_unstructured_cell.F90
	grid_unstructured_explicit.F90
	grid_unstructured_aux.F90
	grid_unstructured.F90
	grid_unstructured_polyhedra.F90
)

set(RELATIONS_SOURCES
	saturation.F90
	saturation_function.F90
	characteristic_curves.F90
	characteristic_curves_base.F90
	characteristic_curves_common.F90
	characteristic_curves_owg.F90
	characteristic_curves_wipp.F90
	mphase_pckr_mod.F90
)

set(PROPERTIES_SOURCES
	fluid.F90
	material.F90 
)

set(DATASET_SOURCES
	dataset_ascii.F90
	dataset_base.F90
	dataset_common_hdf5.F90
	dataset_global_hdf5.F90
	dataset_map_hdf5.F90
	dataset.F90
	dataset_gridded_hdf5.F90
	time_storage.F90 
)

set(SRCSINK_SOURCES
	srcsink_sandbox.F90
	srcsink_sandbox_base.F90
	srcsink_sandbox_mass_rate.F90
	srcsink_sandbox_downreg.F90
	srcsink_sandbox_wipp_gas.F90
	srcsink_sandbox_wipp_well.F90
)

set(IO_SOURCES
	hdf5.F90
	hdf5_aux.F90
	output.F90
	output_common.F90
	output_ekg.F90
	output_hdf5.F90
	output_observation.F90
	output_surface.F90
	output_tecplot.F90
	output_vtk.F90
)

set(MISC_SOURCES
	block_solve.F90
	block_tridiag.F90
	checkpoint.F90
	condition.F90
	condition_control.F90
	convergence.F90
	coupler.F90
	data_mediator.F90
	data_mediator_dataset.F90
	data_mediator_vec.F90
	data_mediator_base.F90
	debug.F90
	e4d_mat_inv.F90
	e4d_run.F90
	e4d_setup.F90
	e4d_vars.F90
	factory_hydrogeophysics.F90
	factory_pflotran.F90
	factory_subsurface.F90
	field.F90
	hydrostatic.F90
	hydrostatic_common.F90
	hydrostatic_multi_phase.F90
	init_common.F90
	init_subsurface.F90
	init_subsurface_flow.F90
	init_subsurface_transport.F90
	integral_flux.F90
	lookup_table.F90
	matrix_buffer.F90
	multisimulation.F90
	observation.F90
	patch.F90
	pmc_auxiliary.F90
	pmc_base.F90
	pmc_hydrogeophysics.F90
	pmc_subsurface.F90
	pmc_third_party.F90
	realization_base.F90
	realization_subsurface.F90
	region.F90
	regression.F90
	secondary_continuum_aux.F90
	secondary_continuum.F90
	simulation_aux.F90
	simulation_base.F90
	simulation_hydrogeophysics.F90
	simulation_subsurface.F90
	solver.F90
	spline.F90
	strata.F90
	timestepper_base.F90
	timestepper_BE.F90
	transport.F90
	uniform_velocity.F90
	waypoint.F90
	wipp.F90
	wrapper_hydrogeophysics.F90
)

set(WELL_SOURCES
	well_spec_base.F90
	well_base.F90
	well_flow.F90
	well_flow_energy.F90
	well_water_injector.F90
	well_oil_producer.F90
	well_toil_ims.F90
	well.F90
)

set(SURFACE_SOURCES
	checkpoint_surface.F90
	factory_surface.F90
	factory_surfsubsurface.F90
	init_surface.F90
	pmc_surface.F90
	simulation_surface.F90
	simulation_surfsubsurface.F90
	surface_auxiliary.F90
	surface_field.F90
	surface_flow.F90
	surface_global_aux.F90
	surface_global.F90
	surface_material.F90
	realization_surface.F90
	surface_th_aux.F90
	surface_th.F90
	timestepper_surface.F90
)

set(GEOMECH_SOURCES
	factory_geomechanics.F90
	gauss.F90
	geomechanics_auxiliary.F90
	geomechanics_realization.F90
	geomechanics_discretization.F90
	geomech_grid.F90
	geomech_grid_aux.F90
	geomechanics_condition.F90
	geomechanics_coupler.F90
	geomechanics_debug.F90
	geomechanics_field.F90
	geomechanics_force.F90
	geomechanics_global.F90
	geomechanics_global_aux.F90
	geomechanics_material.F90
	geomechanics_patch.F90
	geomechanics_region.F90
	geomechanics_regression.F90
	geomechanics_strata.F90
        geomechanics_subsurface_properties.F90
	output_geomechanics.F90
	pm_geomechanics_force.F90
	pmc_geomechanics.F90
	simulation_geomechanics.F90
	shape_function.F90
	timestepper_steady.F90
)

set(SOLVER_SOURCES
	solver.F90
	solver_cpr.F90
	preconditioner_cpr.F90
)

set (PFLOTRAN_SOURCES
    ${UTIL_SOURCES}
    ${EOS_SOURCES}
    ${MODE_AUX_SOURCES}
    ${MODE_SOURCES}
    ${SHARD_MODE_AUX_SOURCES}
    ${CHEM_SOURCES}
    ${GRID_SOURCES}
    ${RELATIONS_SOURCES}
    ${PROPERTIES_SOURCES}
    ${DATASET_SOURCES}
    ${SRCSINK_SOURCES}
    ${IO_SOURCES}
    ${MISC_SOURCES}
    ${WELL_SOURCES}
    ${SURFACE_SOURCES}
    ${GEOMECH_SOURCES}
    ${SOLVER_SOURCES}
    pflotran.F90
)


include(add_pflotran_library)
add_pflotran_library(pflotran ${PFLOTRAN_SOURCES})

set(PFLOTRAN_LIBRARIES pflotran;${PFLOTRAN_LIBRARIES})

link_libraries(pflotran)
link_libraries(${petsc_lib})

include(add_pflotran_executable)
add_pflotran_executable(pflotran.exe pflotran.F90)

# Unit testing.
enable_testing()

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")

#add_test(
#NAME unit_tests
#COMMAND make PFLOTRAN=${CMAKE_BINARY_DIR}/pflotran.exe test
#WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/unittests
#)

add_test(
NAME regression_tests
COMMAND make PFLOTRAN=${CMAKE_BINARY_DIR}/pflotran.exe test
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/../../regression_tests
)
