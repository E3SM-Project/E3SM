.SUFFIXES: .F .o

PHYSICS=-DDO_PHYSICS
#PHYSICS=

OBJS = mpas_atm_mpas_core.o \
       mpas_atm_time_integration.o \
       mpas_atm_advection.o \
       mpas_atm_interp_diagnostics.o

all: physcore core_hyd

physcore:
	( cd ../core_atmos_physics; make all )
	( mkdir phys; cd phys; ar -x ../../core_atmos_physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmos_physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmos_physics/physics_wrf/files/*DBL .)

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS) phys/*.o

mpas_atm_time_integration.o: 

mpas_atm_advection.o: 

mpas_atm_mpas_core.o: mpas_atm_advection.o mpas_atm_time_integration.o mpas_atm_interp_diagnostics.o

clean:
	( cd ../core_atmos_physics; make clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DBL )
	$(RM) -r phys
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../core_atmos_physics -I../core_atmos_physics/physics_wrf -I../external/esmf_time_f90
