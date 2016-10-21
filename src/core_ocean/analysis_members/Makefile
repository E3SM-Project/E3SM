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
          mpas_ocn_zonal_mean.o \
          mpas_ocn_lagrangian_particle_tracking_interpolations.o \
          mpas_ocn_particle_list.o \
          mpas_ocn_lagrangian_particle_tracking_reset.o \
          mpas_ocn_lagrangian_particle_tracking.o \
          mpas_ocn_eliassen_palm.o \
          mpas_ocn_time_filters.o \
          mpas_ocn_mixed_layer_depths.o \
          mpas_ocn_pointwise_stats.o \
          mpas_ocn_debug_diagnostics.o \
          mpas_ocn_time_series_stats.o \
          mpas_ocn_regional_stats.o \
          mpas_ocn_rpn_calculator.o \
          mpas_ocn_transect_transport.o\
          mpas_ocn_eddy_product_variables.o

all: $(OBJS)

mpas_ocn_analysis_driver.o: $(MEMBERS)

mpas_ocn_okubo_weiss.o: mpas_ocn_okubo_weiss_eigenvalues.o

mpas_ocn_particle_list.o:

mpas_ocn_lagrangian_particle_tracking_reset.o:

mpas_ocn_lagrangian_particle_tracking.o: mpas_ocn_particle_list.o mpas_ocn_lagrangian_particle_tracking_interpolations.o mpas_ocn_lagrangian_particle_tracking_reset.o

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
