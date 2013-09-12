.SUFFIXES: .F .o

OBJS = 	mpas_li_mpas_core.o \
	mpas_li_time_integration.o \
	mpas_li_time_integration_fe.o \
	mpas_li_tendency.o \
	mpas_li_setup.o \
	mpas_li_velocity.o \
	mpas_li_sia.o 

all: core_li

core_sw: $(OBJS)
	ar -ru libdycore.a $(OBJS)

mpas_li_mpas_core.o: mpas_li_time_integration.o mpas_li_setup.o mpas_li_velocity.o mpas_li_diagnostic_vars.o

mpas_li_setup.o:

mpas_li_time_integration.o: mpas_li_time_integration_fe.o

mpas_li_time_integration_fe.o: mpas_li_velocity.o mpas_tendency.o mpas_li_diagnostic_vars.o

mpas_tendency.o: 

mpas_li_diagnostic_vars.o:

mpas_li_velocity.o: mpas_li_sia.o

mpas_li_sia.o:

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
