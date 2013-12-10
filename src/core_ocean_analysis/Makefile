.SUFFIXES: .F .o

OBJS = mpas_oac_mpas_core.o

all: shared_ocn shared_oac core_oac

shared_ocn:
	( cd ../core_ocean/shared_ocn; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf ../core_ocean/shared_ocn/libshared_ocn.a .
	ln -sf ../core_ocean/shared_ocn/*.mod .
	( cd shared_oac; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf shared_oac/libshared_oac.a libshared_oac.a
	ln -sf shared_oac/*.mod .

core_oac: $(OBJS)
	ar -ru libdycore.a  $(OBJS) ../core_ocean/shared_ocn/*.o

mpas_ocn_mpas_core.o: 

clean:
	$(RM) *.o *.mod *.f90 lib*.a 
	(cd shared_oac; make clean)
	(cd ../core_ocean/shared_ocn; make clean)

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -Ishared_oac -I../core_ocean/shared_ocn 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -Ishared_oac -I../core_ocean/shared_ocn 
endif
