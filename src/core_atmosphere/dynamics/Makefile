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
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
endif
