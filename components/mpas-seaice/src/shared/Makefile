.SUFFIXES: .F .o

OBJS = 	mpas_seaice_time_integration.o \
	mpas_seaice_advection_incremental_remap_tracers.o \
	mpas_seaice_advection_incremental_remap.o \
	mpas_seaice_advection_upwind.o \
	mpas_seaice_advection.o \
	mpas_seaice_velocity_solver.o \
	mpas_seaice_velocity_solver_weak.o \
	mpas_seaice_velocity_solver_variational.o \
	mpas_seaice_velocity_solver_wachspress.o \
	mpas_seaice_velocity_solver_pwl.o \
	mpas_seaice_velocity_solver_variational_shared.o \
	mpas_seaice_velocity_solver_constitutive_relation.o \
	mpas_seaice_forcing.o \
	mpas_seaice_initialize.o \
	mpas_seaice_testing.o \
	mpas_seaice_mesh.o \
	mpas_seaice_diagnostics.o \
	mpas_seaice_numerics.o \
	mpas_seaice_constants.o \
	mpas_seaice_column.o \
	mpas_seaice_diagnostics.o \
	mpas_seaice_error.o \
	mpas_seaice_mesh_pool.o \
	mpas_seaice_prescribed.o \
	mpas_seaice_special_boundaries.o

all: $(OBJS)

mpas_seaice_constants.o: 

mpas_seaice_error.o:

mpas_seaice_mesh_pool.o:

mpas_seaice_column.o: mpas_seaice_error.o

mpas_seaice_diagnostics.o: mpas_seaice_constants.o

mpas_seaice_mesh.o: mpas_seaice_constants.o

mpas_seaice_numerics.o: mpas_seaice_constants.o

mpas_seaice_testing.o: mpas_seaice_constants.o

mpas_seaice_velocity_solver_constitutive_relation.o: mpas_seaice_constants.o mpas_seaice_testing.o

mpas_seaice_forcing.o: mpas_seaice_constants.o mpas_seaice_mesh.o

mpas_seaice_velocity_solver_weak.o: mpas_seaice_constants.o mpas_seaice_testing.o mpas_seaice_velocity_solver_constitutive_relation.o mpas_seaice_mesh_pool.o

mpas_seaice_velocity_solver_variational_shared.o: mpas_seaice_constants.o

mpas_seaice_velocity_solver_wachspress.o: mpas_seaice_constants.o mpas_seaice_numerics.o mpas_seaice_mesh.o mpas_seaice_testing.o mpas_seaice_velocity_solver_variational_shared.o

mpas_seaice_velocity_solver_pwl.o: mpas_seaice_constants.o mpas_seaice_numerics.o mpas_seaice_mesh.o mpas_seaice_testing.o mpas_seaice_velocity_solver_variational_shared.o

mpas_seaice_velocity_solver_variational.o: mpas_seaice_constants.o mpas_seaice_velocity_solver_constitutive_relation.o mpas_seaice_velocity_solver_wachspress.o mpas_seaice_velocity_solver_pwl.o mpas_seaice_mesh_pool.o

mpas_seaice_velocity_solver.o: mpas_seaice_constants.o mpas_seaice_mesh.o mpas_seaice_testing.o mpas_seaice_velocity_solver_weak.o mpas_seaice_velocity_solver_constitutive_relation.o mpas_seaice_velocity_solver_variational.o mpas_seaice_diagnostics.o mpas_seaice_mesh_pool.o mpas_seaice_special_boundaries.o

mpas_seaice_advection_upwind.o: mpas_seaice_constants.o mpas_seaice_mesh.o

mpas_seaice_advection_incremental_remap_tracers.o: mpas_seaice_constants.o mpas_seaice_diagnostics.o

mpas_seaice_advection_incremental_remap.o: mpas_seaice_constants.o mpas_seaice_mesh.o mpas_seaice_advection_incremental_remap_tracers.o mpas_seaice_diagnostics.o mpas_seaice_error.o

mpas_seaice_advection.o: mpas_seaice_advection_upwind.o mpas_seaice_advection_incremental_remap.o

mpas_seaice_prescribed.o: mpas_seaice_constants.o mpas_seaice_column.o

mpas_seaice_time_integration.o: mpas_seaice_constants.o mpas_seaice_velocity_solver.o mpas_seaice_forcing.o mpas_seaice_advection.o mpas_seaice_diagnostics.o mpas_seaice_column.o mpas_seaice_prescribed.o mpas_seaice_special_boundaries.o

mpas_seaice_initialize.o: mpas_seaice_constants.o mpas_seaice_mesh.o mpas_seaice_velocity_solver.o mpas_seaice_testing.o mpas_seaice_forcing.o mpas_seaice_advection.o mpas_seaice_column.o mpas_seaice_forcing.o mpas_seaice_mesh_pool.o mpas_seaice_special_boundaries.o

mpas_seaice_core.o: mpas_seaice_constants.o mpas_seaice_time_integration.o mpas_seaice_velocity_solver.o mpas_seaice_forcing.o mpas_seaice_initialize.o mpas_seaice_column.o mpas_seaice_mesh_pool.o

mpas_seaice_core_interface.o: mpas_seaice_core.o

clean:
	$(RM) *.o *.i *.mod *.f90

FW = ../../../mpas-framework/src

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90

	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I$(FW)/framework -I$(FW)/operators -I$(FW)/external/esmf_time_f90 -I../column
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I$(FW)/framework -I$(FW)/operators -I$(FW)/external/esmf_time_f90 -I../column
endif
