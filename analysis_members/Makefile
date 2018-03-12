.SUFFIXES: .F .c .o

OBJS = mpas_cice_analysis_driver.o

MEMBERS = mpas_cice_high_frequency_output.o mpas_cice_temperatures.o mpas_cice_regional_statistics.o mpas_cice_ridging_diagnostics.o mpas_cice_conservation_check.o mpas_cice_geographical_vectors.o mpas_cice_ice_present.o mpas_cice_time_series_stats.o mpas_cice_load_balance.o mpas_cice_maximum_ice_presence.o mpas_cice_miscellaneous.o mpas_cice_area_variables.o mpas_cice_pond_diagnostics.o mpas_cice_deactivate_unneeded_fields.o mpas_cice_pointwise_stats.o mpas_cice_unit_conversion.o mpas_cice_ice_shelves.o

all: $(OBJS)

mpas_cice_area_variables.o: mpas_cice_deactivate_unneeded_fields.o

mpas_cice_pond_diagnostics.o:

mpas_cice_analysis_driver.o: $(MEMBERS)

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../column -I../shared -I../model_forward $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../column -I../shared -I../model_forward $(CPPINCLUDES) $(FCINCLUDES) 
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../column -I../shared -I../model_forward -c $<
