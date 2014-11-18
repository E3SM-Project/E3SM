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

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
endif
