.SUFFIXES: .F .o

OBJS = mpas_oac_driver.o \
       mpas_oac_global_stats.o

all: shared_oac

shared_oac: $(OBJS)
	ar -ru libshared_oac.a $(OBJS)

mpas_oac_driver.o: mpas_oac_global_stats.o

mpas_oac_global_stats.o:

clean:
	$(RM) *.o *.mod *.f90 lib*.a 

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90 -I../../core_ocean/shared_ocn
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90 -I../../core_ocean/shared_ocn
endif
