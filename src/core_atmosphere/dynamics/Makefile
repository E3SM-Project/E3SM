.SUFFIXES: .F .o

OBJS = mpas_atm_time_integration.o \
       mpas_atm_boundaries.o

all: $(OBJS)

mpas_atm_time_integration.o: mpas_atm_boundaries.o mpas_atm_iau.o

mpas_atm_boundaries.o:


clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
endif
