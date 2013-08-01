.SUFFIXES: .F .o

OBJS = mpas_subdriver.o \
       mpas.o

all: $(OBJS)

mpas_subdriver.o: 

mpas.o: mpas_subdriver.o

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
endif
