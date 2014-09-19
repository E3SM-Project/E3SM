.SUFFIXES: .F .o

OBJS = 	mpas_cice_mpas_core.o \
	mpas_cice_time_integration.o \
	mpas_cice_advection.o \
	mpas_cice_velocity_solver_unit_tests.o \
	mpas_cice_velocity_solver.o \
	mpas_cice_velocity_solver_weak.o \
	mpas_cice_velocity_solver_variational.o \
	mpas_cice_velocity_solver_wachspress.o \
	mpas_cice_velocity_solver_pwl.o \
	mpas_cice_velocity_solver_shared.o \
	mpas_cice_forcing.o \
	mpas_cice_therm_zerolayer.o \
	mpas_cice_initialize.o \
	mpas_cice_testing.o \
	mpas_cice_debug.o \
	mpas_cice_unit_test.o \
	mpas_cice_mesh.o \
	mpas_cice_diagnostics.o \
	mpas_cice_numerics.o \
	mpas_cice_constants.o \
	MPAS_forcing.o \
	mpas_cice_forcing_test.o \

all: core_cice

core_cice: $(OBJS)
	ar -ru libdycore.a $(OBJS)

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

mpas_cice_constants.o:

MPAS_forcing.o:

mpas_cice_forcing_test.o: MPAS_forcing.o

mpas_cice_mesh.o: mpas_cice_constants.o mpas_cice_forcing_test.o

mpas_cice_diagnostics.o: mpas_cice_constants.o 

mpas_cice_numerics.o: mpas_cice_constants.o

mpas_cice_debug.o: mpas_cice_constants.o

mpas_cice_testing.o: mpas_cice_constants.o mpas_cice_debug.o

mpas_cice_velocity_solver_shared.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o

mpas_cice_therm_zerolayer.o: mpas_cice_constants.o mpas_cice_numerics.o

mpas_cice_forcing.o: mpas_cice_constants.o mpas_cice_mesh.o

mpas_cice_velocity_solver_weak.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_velocity_solver_shared.o

mpas_cice_velocity_solver_variational.o: mpas_cice_constants.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o

mpas_cice_velocity_solver_wachspress.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_variational.o

mpas_cice_velocity_solver_pwl.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_variational.o

mpas_cice_velocity_solver.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_pwl.o mpas_cice_velocity_solver_wachspress.o mpas_cice_velocity_solver_weak.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o

mpas_cice_velocity_solver_unit_tests.o: mpas_cice_constants.o mpas_cice_mesh.o

mpas_cice_advection.o: 

mpas_cice_unit_test.o: mpas_cice_constants.o mpas_cice_velocity_solver_unit_tests.o

mpas_cice_time_integration.o: mpas_cice_constants.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_unit_test.o mpas_cice_advection.o mpas_cice_diagnostics.o

mpas_cice_initialize.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_velocity_solver.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_advection.o

mpas_cice_mpas_core.o: mpas_cice_constants.o mpas_cice_time_integration.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_initialize.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
