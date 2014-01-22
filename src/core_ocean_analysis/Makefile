.SUFFIXES: .F .o

OBJS = mpas_oac_mpas_core.o

all: oac_shared ocean_shared core_oac

ocean_shared:
	( cd ../core_ocean; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" ocean_shared ) 
	ln -sf ../core_ocean/ocn_init.mod .
	ln -sf ../core_ocean/ocn_diagnostics_routines.mod .

oac_shared: ocean_shared
	( cd shared_oac; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf shared_oac/*.mod .

core_oac: $(OBJS)
	ar -ru libdycore.a  $(OBJS) ../core_ocean/mpas_ocn_init.o ../core_ocean/mpas_ocn_diagnostics_routines.o shared_oac/*.o

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

mpas_ocn_mpas_core.o: ocean_shared oac_shared

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	(cd shared_oac; make clean)

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 
endif
