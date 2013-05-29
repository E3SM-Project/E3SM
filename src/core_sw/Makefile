.SUFFIXES: .F .o

OBJS = 	mpas_sw_mpas_core.o \
        mpas_sw_test_cases.o \
	mpas_sw_advection.o \
	mpas_sw_time_integration.o \
	mpas_sw_global_diagnostics.o

all: core_sw

core_sw: $(OBJS)
	ar -ru libdycore.a $(OBJS)

mpas_sw_test_cases.o:

mpas_sw_advection.o:

mpas_sw_time_integration.o:

mpas_sw_global_diagnostics.o:

mpas_sw_mpas_core.o: mpas_sw_global_diagnostics.o mpas_sw_test_cases.o mpas_sw_time_integration.o mpas_sw_advection.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
