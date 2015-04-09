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
	mpas_cice_velocity_solver_variational_shared.o \
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

all: core_cice

core_cice: $(OBJS)
	ar -ru libdycore.a $(OBJS)

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.cice )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.cice stream_list.cice. listed )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

mpas_cice_constants.o:

mpas_cice_mesh.o: mpas_cice_constants.o

mpas_cice_diagnostics.o: mpas_cice_constants.o 

mpas_cice_numerics.o: mpas_cice_constants.o

mpas_cice_debug.o: mpas_cice_constants.o mpas_cice_mesh.o

mpas_cice_testing.o: mpas_cice_constants.o mpas_cice_debug.o

mpas_cice_velocity_solver_shared.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o

mpas_cice_therm_zerolayer.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_debug.o

mpas_cice_forcing.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_debug.o

mpas_cice_velocity_solver_weak.o: mpas_cice_constants.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_velocity_solver_shared.o

mpas_cice_velocity_solver_variational_shared.o: mpas_cice_constants.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o

mpas_cice_velocity_solver_wachspress.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_variational_shared.o

mpas_cice_velocity_solver_pwl.o: mpas_cice_constants.o mpas_cice_numerics.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_variational_shared.o

mpas_cice_velocity_solver_variational.o: mpas_cice_constants.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_wachspress.o mpas_cice_velocity_solver_pwl.o

mpas_cice_velocity_solver.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_testing.o mpas_cice_velocity_solver_weak.o mpas_cice_velocity_solver_shared.o mpas_cice_debug.o mpas_cice_velocity_solver_variational.o

mpas_cice_velocity_solver_unit_tests.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_velocity_solver_weak.o mpas_cice_velocity_solver_variational.o

mpas_cice_advection.o: mpas_cice_mesh.o mpas_cice_debug.o

mpas_cice_unit_test.o: mpas_cice_constants.o mpas_cice_velocity_solver_unit_tests.o

mpas_cice_time_integration.o: mpas_cice_constants.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_unit_test.o mpas_cice_advection.o mpas_cice_diagnostics.o

mpas_cice_initialize.o: mpas_cice_constants.o mpas_cice_mesh.o mpas_cice_velocity_solver.o mpas_cice_testing.o mpas_cice_debug.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_advection.o

mpas_cice_mpas_core.o: mpas_cice_constants.o mpas_cice_time_integration.o mpas_cice_velocity_solver.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o mpas_cice_initialize.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
endif
