.SUFFIXES: .F .c .o

OBJS = mpas_seaice_analysis_driver.o

MEMBERS = mpas_seaice_high_frequency_output.o mpas_seaice_temperatures.o mpas_seaice_thicknesses.o mpas_seaice_regional_statistics.o mpas_seaice_ridging_diagnostics.o mpas_seaice_conservation_check.o mpas_seaice_geographical_vectors.o mpas_seaice_ice_present.o mpas_seaice_time_series_stats.o mpas_seaice_load_balance.o mpas_seaice_maximum_ice_presence.o mpas_seaice_miscellaneous.o mpas_seaice_area_variables.o mpas_seaice_pond_diagnostics.o mpas_seaice_deactivate_unneeded_fields.o mpas_seaice_pointwise_stats.o mpas_seaice_unit_conversion.o mpas_seaice_ice_shelves.o

all: $(OBJS)

mpas_seaice_area_variables.o: mpas_seaice_deactivate_unneeded_fields.o

mpas_seaice_pond_diagnostics.o:

mpas_seaice_analysis_driver.o: $(MEMBERS)

clean:
	$(RM) *.o *.i *.mod *.f90

FW = ../../../mpas-framework/src

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 -I$(FW)/external/esmf_time_f90 -I$(FW)/framework -I$(FW)/operators -I../column -I../shared -I../model_forward $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F -I$(FW)/external/esmf_time_f90 -I$(FW)/framework -I$(FW)/operators -I../column -I../shared -I../model_forward $(CPPINCLUDES) $(FCINCLUDES)
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -I$(FW)/external/esmf_time_f90 -I$(FW)/framework -I$(FW)/operators -I../column -I../shared -I../model_forward -c $<
