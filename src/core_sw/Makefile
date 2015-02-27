.SUFFIXES: .F .o

OBJS = 	mpas_sw_mpas_core.o \
        mpas_sw_test_cases.o \
        mpas_sw_advection.o \
        mpas_sw_time_integration.o \
        mpas_sw_global_diagnostics.o \
        mpas_sw_constants.o

all: core_sw

core_sw: $(OBJS)
	ar -ru libdycore.a $(OBJS)

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	(cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.sw )
	(cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.sw stream_list.sw. listed )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )



mpas_sw_constants.o:

mpas_sw_test_cases.o: mpas_sw_constants.o

mpas_sw_advection.o: mpas_sw_constants.o

mpas_sw_time_integration.o: mpas_sw_constants.o

mpas_sw_global_diagnostics.o: mpas_sw_constants.o

mpas_sw_mpas_core.o: mpas_sw_global_diagnostics.o mpas_sw_test_cases.o mpas_sw_time_integration.o mpas_sw_advection.o mpas_sw_constants.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
endif
