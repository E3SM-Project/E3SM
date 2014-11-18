.SUFFIXES: .F .o

OBJS = mpas_subdriver.o \
       mpas.o

all:
	($(MAKE) clean)
	($(MAKE) driver)

driver: $(OBJS)

mpas_subdriver.o: 

mpas.o: mpas_subdriver.o

clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../core_$(CORE) -I../external/esmf_time_f90
endif
