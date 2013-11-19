.SUFFIXES: .F .o

COREOBJS = mpas_oac_mpas_core.o

DRIVEROBJS = mpas_oac_driver.o \
       mpas_oac_am_global_stats.o  

all: sharedocean oacdriver oaccore

run_time_mode: oacdriver

sharedocean:
	( cd ../core_ocean/shared_ocean; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf ../core_ocean/shared_ocean/libsharedocn.a .

oacdriver: $(DRIVEROBJS)
	ar -ru liboac.a $(DRIVEROBJS)

oaccore: $(COREOBJS)
	ar -ru libdycore.a  $(COREOBJS) $(DRIVEROBJS) ../core_ocean/shared_ocean/*.o

mpas_oac_am_global_stats.o:

mpas_oac_driver.o: mpas_oac_am_global_stats.o

mpas_ocn_mpas_core.o: mpas_oac_driver.o

clean:
	$(RM) *.o *.mod *.f90 liboac.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -I../core_ocean/shared_ocean
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -I../core_ocean/shared_ocean
endif
