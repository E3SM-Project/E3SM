.SUFFIXES: .F .c .o

OBJS = mpas_ocn_analysis_driver.o

MEMBERS = mpas_ocn_global_stats.o \
          mpas_ocn_okubo_weiss.o \
          mpas_ocn_layer_volume_weighted_averages.o \
          mpas_ocn_surface_area_weighted_averages.o \
          mpas_ocn_water_mass_census.o \
          mpas_ocn_meridional_heat_transport.o \
          mpas_ocn_test_compute_interval.o \
          mpas_ocn_high_frequency_output.o \
          mpas_ocn_zonal_mean.o

all: $(OBJS)

mpas_ocn_analysis_driver.o: $(MEMBERS)

mpas_ocn_okubo_weiss.o: mpas_ocn_okubo_weiss_eigenvalues.o

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) 
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -c $<
