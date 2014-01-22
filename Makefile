.SUFFIXES: .F .o

CVMIX_REPO_ADDRESS=http://cvmix.googlecode.com/svn/trunk/src/shared

OBJS = mpas_ocn_mpas_core.o \
       mpas_ocn_init.o \
       mpas_ocn_thick_hadv.o \
       mpas_ocn_thick_vadv.o \
       mpas_ocn_thick_surface_flux.o \
       mpas_ocn_gm.o \
       mpas_ocn_vel_coriolis.o \
       mpas_ocn_vel_vadv.o \
       mpas_ocn_vel_hmix.o \
       mpas_ocn_vel_hmix_del2.o \
       mpas_ocn_vel_hmix_leith.o \
       mpas_ocn_vel_hmix_del4.o \
       mpas_ocn_vel_forcing.o \
       mpas_ocn_vel_forcing_windstress.o \
       mpas_ocn_vel_forcing_rayleigh.o \
       mpas_ocn_vel_pressure_grad.o \
       mpas_ocn_vmix.o \
       mpas_ocn_vmix_coefs_const.o \
       mpas_ocn_vmix_coefs_rich.o \
       mpas_ocn_vmix_coefs_tanh.o \
       mpas_ocn_vmix_cvmix.o \
       mpas_ocn_tendency.o \
       mpas_ocn_diagnostics.o \
       mpas_ocn_diagnostics_routines.o \
       mpas_ocn_thick_ale.o \
       mpas_ocn_tracer_hmix.o \
       mpas_ocn_tracer_hmix_del2.o \
       mpas_ocn_tracer_hmix_del4.o \
       mpas_ocn_tracer_advection.o \
	   mpas_ocn_tracer_short_wave_absorption.o \
	   mpas_ocn_tracer_short_wave_absorption_jerlov.o \
       mpas_ocn_high_freq_thickness_hmix_del2.o \
       mpas_ocn_tracer_surface_flux.o \
       mpas_ocn_time_integration.o \
       mpas_ocn_time_integration_rk4.o \
       mpas_ocn_time_integration_split.o \
       mpas_ocn_equation_of_state.o \
       mpas_ocn_equation_of_state_jm.o \
       mpas_ocn_equation_of_state_linear.o \
       mpas_ocn_global_diagnostics.o \
       mpas_ocn_test.o \
       mpas_ocn_constants.o \
       mpas_ocn_forcing.o \
       mpas_ocn_forcing_bulk.o \
       mpas_ocn_forcing_restoring.o \
       mpas_ocn_time_average.o \
       mpas_ocn_time_average_coupled.o \
       mpas_ocn_sea_ice.o

all: libcvmix oac_shared core_ocean

libcvmix:
	if [ ! -d cvmix ]; then \
		(svn checkout $(CVMIX_REPO_ADDRESS) cvmix) \
	fi
	if [ -d cvmix ]; then \
		(cd cvmix; svn update; make all FC="$(FC)" FFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	fi
	ln -sf cvmix/*.mod .

oac_shared: mpas_ocn_diagnostics_routines.o
	( cd ../core_ocean_analysis/shared_oac; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)" all ) 
	ln -sf ../core_ocean_analysis/shared_oac/*.mod .

core_ocean: libcvmix oac_shared $(OBJS) 
	ar -ru libdycore.a $(OBJS) cvmix/*.o ../core_ocean_analysis/shared_oac/*.o 

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

mpas_ocn_time_integration.o: mpas_ocn_time_integration_rk4.o mpas_ocn_time_integration_split.o

mpas_ocn_time_integration_rk4.o: mpas_ocn_tendency.o mpas_ocn_diagnostics.o mpas_ocn_time_average_coupled.o mpas_ocn_sea_ice.o

mpas_ocn_time_integration_split.o: mpas_ocn_tendency.o mpas_ocn_diagnostics.o mpas_ocn_time_average_coupled.o mpas_ocn_sea_ice.o

mpas_ocn_tendency.o: mpas_ocn_time_average.o mpas_ocn_high_freq_thickness_hmix_del2.o mpas_ocn_tracer_surface_flux.o mpas_ocn_thick_surface_flux.o mpas_ocn_tracer_short_wave_absorption.o

mpas_ocn_diagnostics.o: mpas_ocn_time_average.o mpas_ocn_thick_ale.o mpas_ocn_diagnostics_routines.o

mpas_ocn_thick_ale.o: 

mpas_ocn_global_diagnostics.o: 

mpas_ocn_time_average.o:

mpas_ocn_time_average_coupled.o: mpas_ocn_constants.o

mpas_ocn_thick_hadv.o:

mpas_ocn_thick_vadv.o:

mpas_ocn_thick_surface_flux.o: mpas_ocn_forcing.o

mpas_ocn_gm.o: 

mpas_ocn_vel_pressure_grad.o:

mpas_ocn_vel_vadv.o:

mpas_ocn_vel_hmix.o: mpas_ocn_vel_hmix_del2.o mpas_ocn_vel_hmix_leith.o mpas_ocn_vel_hmix_del4.o

mpas_ocn_vel_hmix_del2.o:

mpas_ocn_vel_hmix_leith.o:

mpas_ocn_vel_hmix_del4.o:

mpas_ocn_vel_forcing.o: mpas_ocn_vel_forcing_windstress.o mpas_ocn_vel_forcing_rayleigh.o mpas_ocn_forcing.o

mpas_ocn_vel_forcing_windstress.o:

mpas_ocn_vel_forcing_rayleigh.o:

mpas_ocn_vel_coriolis.o:

mpas_ocn_tracer_hmix.o: mpas_ocn_tracer_hmix_del2.o mpas_ocn_tracer_hmix_del4.o

mpas_ocn_tracer_hmix_del2.o:

mpas_ocn_tracer_hmix_del4.o:

mpas_ocn_tracer_advection.o:

mpas_ocn_high_freq_thickness_hmix_del2.o:

mpas_ocn_tracer_surface_flux.o: mpas_ocn_forcing.o

mpas_ocn_tracer_short_wave_absorption.o: mpas_ocn_tracer_short_wave_absorption_jerlov.o

mpas_ocn_tracer_short_wave_absorption_jerlov.o:

mpas_ocn_vmix.o: mpas_ocn_vmix_coefs_const.o mpas_ocn_vmix_coefs_rich.o mpas_ocn_vmix_coefs_tanh.o mpas_ocn_vmix_cvmix.o

mpas_ocn_vmix_coefs_const.o:

mpas_ocn_vmix_coefs_rich.o: mpas_ocn_equation_of_state.o

mpas_ocn_vmix_coefs_tanh.o:

mpas_ocn_vmix_cvmix.o: libcvmix

mpas_ocn_equation_of_state.o: mpas_ocn_equation_of_state_jm.o mpas_ocn_equation_of_state_linear.o

mpas_ocn_equation_of_state_jm.o:

mpas_ocn_equation_of_state_linear.o:

mpas_ocn_test.o: 

mpas_ocn_constants.o:

mpas_ocn_forcing.o: mpas_ocn_constants.o mpas_ocn_forcing_bulk.o mpas_ocn_forcing_restoring.o

mpas_ocn_forcing_bulk.o:

mpas_ocn_forcing_restoring.o:

mpas_ocn_sea_ice.o:

mpas_ocn_mpas_core.o: mpas_ocn_thick_hadv.o \
                      mpas_ocn_thick_vadv.o \
                      mpas_ocn_init.o \
                      mpas_ocn_thick_surface_flux.o \
                      mpas_ocn_gm.o \
                      mpas_ocn_vel_coriolis.o \
                      mpas_ocn_vel_vadv.o \
                      mpas_ocn_vel_hmix.o \
                      mpas_ocn_vel_hmix_del2.o \
                      mpas_ocn_vel_hmix_leith.o \
                      mpas_ocn_vel_hmix_del4.o \
                      mpas_ocn_vel_forcing.o \
                      mpas_ocn_vel_forcing_windstress.o \
                      mpas_ocn_vel_pressure_grad.o \
                      mpas_ocn_tracer_hmix.o \
                      mpas_ocn_tracer_hmix_del2.o \
                      mpas_ocn_tracer_hmix_del4.o \
                      mpas_ocn_high_freq_thickness_hmix_del2.o \
                      mpas_ocn_vmix.o \
                      mpas_ocn_vmix_coefs_const.o \
                      mpas_ocn_vmix_coefs_rich.o \
                      mpas_ocn_vmix_coefs_tanh.o \
                      mpas_ocn_vmix_cvmix.o \
                      mpas_ocn_tracer_advection.o \
                      mpas_ocn_tracer_surface_flux.o \
					  mpas_ocn_tracer_short_wave_absorption.o \
					  mpas_ocn_tracer_short_wave_absorption_jerlov.o \
                      mpas_ocn_tendency.o \
                      mpas_ocn_diagnostics.o \
                      mpas_ocn_thick_ale.o \
                      mpas_ocn_time_integration.o \
                      mpas_ocn_time_integration_rk4.o \
                      mpas_ocn_time_integration_split.o \
                      mpas_ocn_equation_of_state.o \
                      mpas_ocn_equation_of_state_jm.o \
                      mpas_ocn_equation_of_state_linear.o \
                      mpas_ocn_global_diagnostics.o \
                      mpas_ocn_test.o \
                      mpas_ocn_constants.o \
                      mpas_ocn_forcing.o \
                      mpas_ocn_forcing_bulk.o \
                      mpas_ocn_forcing_restoring.o \
                      mpas_ocn_time_average.o \
                      mpas_ocn_time_average_coupled.o \
					  mpas_ocn_sea_ice.o

clean:
	if [ -d cvmix ]; then \
		(cd cvmix; make clean) \
	fi
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	(cd ../core_ocean_analysis; make clean)
	(cd ../core_ocean_analysis/shared_oac; make clean)

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -I./cvmix/ -I../core_ocean_analysis/shared_oac 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90 -I./cvmix/ -I../core_ocean_analysis/shared_oac 
endif
