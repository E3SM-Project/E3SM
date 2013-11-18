.SUFFIXES: .F .o

OBJS = mpas_oac_driver.o \
       mpas_oac_am_global_stats.o  

all: core_ocean_analysis

core_ocean_analysis: $(OBJS)
	ar -ru liboac.a $(OBJS)

mpas_oac_am_global_stats.o:

mpas_oac_driver.o: mpas_oac_am_global_stats.o

clean:
	$(RM) *.o *.mod *.f90 liboac.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../core_ocean/shared_ocean
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../core_ocean/shared_ocean
endif
