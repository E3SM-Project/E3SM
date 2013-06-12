.SUFFIXES: .F .o

OBJS = mpas_atm_time_integration.o \
       mpas_atm_advection.o

all: $(OBJS)

mpas_atm_time_integration.o: 

mpas_atm_advection.o: 

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
