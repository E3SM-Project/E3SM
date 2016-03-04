.SUFFIXES: .F .o

OBJS = 	mpas_cice_time_integration.o \
	mpas_cice_advection_incremental_remap_tracers.o \
	mpas_cice_advection_incremental_remap.o \
	mpas_cice_advection_upwind.o \
	mpas_cice_advection.o \
	mpas_cice_velocity_solver_unit_tests.o \
	mpas_cice_velocity_solver.o \
	mpas_cice_velocity_solver_weak.o \
	mpas_cice_velocity_solver_variational.o \
	mpas_cice_velocity_solver_wachspress.o \
	mpas_cice_velocity_solver_pwl.o \
	mpas_cice_velocity_solver_variational_shared.o \
	mpas_cice_velocity_solver_constitutive_relation.o \
	mpas_cice_forcing.o \
	mpas_cice_initialize.o \
	mpas_cice_testing.o \
	mpas_cice_debug.o \
	mpas_cice_unit_test.o \
	mpas_cice_mesh.o \
	mpas_cice_diagnostics.o \
	mpas_cice_numerics.o \
	mpas_cice_constants.o \
	mpas_cice_column.o \
	mpas_cice_history.o \
	mpas_cice_diagnostics.o

all: $(OBJS)

mpas_cice_constants.o: 

mpas_cice_column.o:

mpas_cice_history.o:

mpas_cice_diagnostics.o: mpas_cice_constants.o

mpas_cice_mesh.o: mpas_cice_constants.o

mpas_cice_diagnostics.o: mpas_cice_constants.o 

mpas_cice_numerics.o: mpas_cice_constants.o

mpas_cice_debug.o: mpas_cice_constants.o mpas_cice_mesh.o

mpas_cice_testing.o: mpas_cice_constants.o mpas_cice_debug.o

mpas_cice_velocity_solver_constitutive_relation.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o

mpas_cice_forcing.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_debug.o

mpas_cice_velocity_solver_weak.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_velocity_solver_constitutive_relation.o

mpas_cice_velocity_solver_variational_shared.o: mpas_cice_constants.o mpas_cice_debug.o

mpas_cice_velocity_solver_wachspress.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_velocity_solver_variational_shared.o

mpas_cice_velocity_solver_pwl.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_velocity_solver_variational_shared.o

mpas_cice_velocity_solver_variational.o: mpas_cice_constants.o mpas_cice_velocity_solver_constitutive_relation.o mpas_cice_debug.o mpas_cice_velocity_solver_wachspress.o mpas_cice_velocity_solver_pwl.o

mpas_cice_velocity_solver.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_weak.o mpas_cice_velocity_solver_constitutive_relation.o mpas_cice_debug.o mpas_cice_velocity_solver_variational.o

mpas_cice_velocity_solver_unit_tests.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_velocity_solver_weak.o mpas_cice_velocity_solver_variational.o

mpas_cice_advection_upwind.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_debug.o

mpas_cice_advection_incremental_remap_tracers.o: mpas_cice_constants.o 

mpas_cice_advection_incremental_remap.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_advection_incremental_remap_tracers.o

mpas_cice_advection.o: mpas_cice_advection_upwind.o mpas_cice_advection_incremental_remap.o

mpas_cice_unit_test.o: mpas_cice_constants.o mpas_cice_velocity_solver_unit_tests.o

mpas_cice_time_integration.o: mpas_cice_constants.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_unit_test.o mpas_cice_advection.o mpas_cice_diagnostics.o mpas_cice_column.o mpas_cice_history.o

mpas_cice_initialize.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_velocity_solver.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_forcing.o mpas_cice_advection.o mpas_cice_column.o mpas_cice_forcing.o

mpas_cice_core.o: mpas_cice_constants.o mpas_cice_time_integration.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_initialize.o mpas_cice_column.o

mpas_cice_core_interface.o: mpas_cice_core.o

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90

	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column
endif
