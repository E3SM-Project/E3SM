.SUFFIXES: .F .o

OBJS = 	mpas_cice_mpas_core.o \
	mpas_cice_time_integration.o \
	mpas_cice_dynamics.o \
	mpas_cice_dynamics_hex.o \
	mpas_cice_dynamics_tri.o \
	mpas_cice_dyn_hex_weak.o \
	mpas_cice_dyn_hex_variational.o \
	mpas_cice_dyn_hex_pwl.o \
	mpas_cice_dyn_tri_weak.o \
	mpas_cice_dynamics_shared.o \
	mpas_cice_forcing.o \
	mpas_cice_therm_zerolayer.o \
	mpas_cice_testing.o

all: core_cice

core_cice: $(OBJS)
	ar -ru libdycore.a $(OBJS)

mpas_cice_testing.o:

mpas_cice_dynamics_shared.o: mpas_cice_testing.o

mpas_cice_therm_zerolayer.o: mpas_cice_dynamics_shared.o

mpas_cice_forcing.o: mpas_cice_dynamics_shared.o

mpas_cice_dyn_hex_weak.o: mpas_cice_testing.o mpas_cice_dynamics_shared.o

mpas_cice_dyn_hex_variational.o: mpas_cice_testing.o mpas_cice_dynamics_shared.o

mpas_cice_dyn_hex_pwl.o: mpas_cice_testing.o mpas_cice_dynamics_shared.o

mpas_cice_dynamics_hex.o: mpas_cice_testing.o mpas_cice_dyn_hex_pwl.o mpas_cice_dyn_hex_variational.o mpas_cice_dyn_hex_weak.o mpas_cice_dynamics_shared.o

mpas_cice_dyn_tri_weak.o: mpas_cice_testing.o mpas_cice_dynamics_shared.o

mpas_cice_dynamics_tri.o: mpas_cice_testing.o mpas_cice_dyn_tri_weak.o mpas_cice_dynamics_shared.o

mpas_cice_dynamics.o: mpas_cice_dynamics_tri.o mpas_cice_dynamics_hex.o mpas_cice_dynamics_shared.o

mpas_cice_time_integration.o: mpas_cice_dynamics.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o

mpas_cice_mpas_core.o: mpas_cice_time_integration.o mpas_cice_dynamics.o mpas_cice_forcing.o mpas_cice_therm_zerolayer.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
