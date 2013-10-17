.SUFFIXES: .F .o

PHYSICS=-DDO_PHYSICS
#PHYSICS=

OBJS = mpas_atm_mpas_core.o \
       mpas_atm_interp_diagnostics.o

all: physcore dycore atmcore

physcore:
	( cd physics; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DBL .)

dycore:
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: physcore dycore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o libphys/*.o

mpas_atm_mpas_core.o: dycore mpas_atm_interp_diagnostics.o

clean:
	( cd physics; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DBL )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
endif
