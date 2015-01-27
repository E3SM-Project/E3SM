.SUFFIXES: .F .o

PHYSICS=-DDO_PHYSICS
#PHYSICS=

OBJS = mpas_atm_mpas_core.o \
       mpas_atm_interp_diagnostics.o

all: physcore dycore atmcore

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.atmosphere )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.atmosphere stream_list.atmosphere. )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

physcore:
	( cd physics; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DATA* .)

dycore:
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: physcore dycore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o libphys/*.o

mpas_atm_mpas_core.o: dycore mpas_atm_interp_diagnostics.o

clean:
	( cd physics; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DATA* )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
endif
