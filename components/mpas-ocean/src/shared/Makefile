.SUFFIXES: .F .o

OBJS = mpas_ocn_init_routines.o \
	   mpas_ocn_gm.o \
	   mpas_ocn_eddy_parameterization_helpers.o \
	   mpas_ocn_submesoscale_eddies.o \
	   mpas_ocn_diagnostics.o \
	   mpas_ocn_diagnostics_variables.o \
	   mpas_ocn_thick_ale.o \
	   mpas_ocn_equation_of_state.o \
	   mpas_ocn_equation_of_state_jm.o \
	   mpas_ocn_equation_of_state_linear.o \
	   mpas_ocn_equation_of_state_wright.o \
	   mpas_ocn_mesh.o \
	   mpas_ocn_thick_hadv.o \
	   mpas_ocn_thick_vadv.o \
	   mpas_ocn_thick_surface_flux.o \
	   mpas_ocn_vel_hadv_coriolis.o \
	   mpas_ocn_vel_vadv.o \
	   mpas_ocn_vel_hmix.o \
	   mpas_ocn_vel_hmix_del2.o \
	   mpas_ocn_vel_hmix_leith.o \
	   mpas_ocn_vel_hmix_del4.o \
	   mpas_ocn_vel_forcing.o \
	   mpas_ocn_vel_forcing_surface_stress.o \
	   mpas_ocn_vel_forcing_explicit_bottom_drag.o \
	   mpas_ocn_vel_pressure_grad.o \
	   mpas_ocn_vertical_remap.o \
	   mpas_ocn_vertical_regrid.o \
	   mpas_ocn_vmix.o \
	   mpas_ocn_vmix_coefs_redi.o \
	   mpas_ocn_vmix_cvmix.o \
	   mpas_ocn_vmix_gotm.o \
	   mpas_ocn_tendency.o \
	   mpas_ocn_tracer_hmix.o \
	   mpas_ocn_tracer_hmix_del2.o \
	   mpas_ocn_tracer_hmix_del4.o \
	   mpas_ocn_tracer_hmix_redi.o \
	   mpas_ocn_tracer_advection.o \
	   mpas_ocn_tracer_advection_mono.o \
	   mpas_ocn_tracer_advection_std.o \
	   mpas_ocn_tracer_advection_vert.o \
	   mpas_ocn_tracer_advection_shared.o \
	   mpas_ocn_tracer_nonlocalflux.o \
	   mpas_ocn_tracer_short_wave_absorption.o \
	   mpas_ocn_tracer_short_wave_absorption_jerlov.o \
	   mpas_ocn_tracer_short_wave_absorption_variable.o \
	   mpas_ocn_tracer_surface_restoring.o \
	   mpas_ocn_tracer_interior_restoring.o \
	   mpas_ocn_tracer_exponential_decay.o \
	   mpas_ocn_tracer_ideal_age.o \
	   mpas_ocn_tracer_TTD.o \
	   mpas_ocn_tracer_ecosys.o \
	   mpas_ocn_tracer_DMS.o \
	   mpas_ocn_tracer_MacroMolecules.o \
	   mpas_ocn_tracer_CFC.o \
	   mpas_ocn_high_freq_thickness_hmix_del2.o \
	   mpas_ocn_tracer_surface_flux_to_tend.o \
	   mpas_ocn_test.o \
	   mpas_ocn_constants.o \
	   mpas_ocn_config.o \
	   mpas_ocn_forcing.o \
	   mpas_ocn_surface_bulk_forcing.o \
	   mpas_ocn_surface_land_ice_fluxes.o \
	   mpas_ocn_effective_density_in_land_ice.o \
	   mpas_ocn_frazil_forcing.o \
	   mpas_ocn_tidal_forcing.o \
	   mpas_ocn_forcing_restoring.o \
	   mpas_ocn_time_average_coupled.o \
	   mpas_ocn_framework_forcing.o \
	   mpas_ocn_time_varying_forcing.o \
	   mpas_ocn_subgrid.o \
	   mpas_ocn_wetting_drying.o \
	   mpas_ocn_vel_tidal_potential.o \
	   mpas_ocn_vel_forcing_topographic_wave_drag.o \
	   mpas_ocn_transport_tests.o \
	   mpas_ocn_vel_self_attraction_loading.o \
	   mpas_ocn_vertical_advection.o \
           mpas_ocn_stokes_drift.o \
           mpas_ocn_manufactured_solution.o

all: $(OBJS)

mpas_ocn_init_routines.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics.o mpas_ocn_diagnostics_variables.o mpas_ocn_gm.o mpas_ocn_submesoscale_eddies.o mpas_ocn_forcing.o mpas_ocn_surface_land_ice_fluxes.o

mpas_ocn_tendency.o: mpas_ocn_high_freq_thickness_hmix_del2.o mpas_ocn_tracer_surface_restoring.o mpas_ocn_thick_surface_flux.o mpas_ocn_tracer_short_wave_absorption.o mpas_ocn_tracer_advection.o mpas_ocn_tracer_hmix.o mpas_ocn_tracer_nonlocalflux.o mpas_ocn_surface_bulk_forcing.o mpas_ocn_surface_land_ice_fluxes.o mpas_ocn_tracer_surface_flux_to_tend.o mpas_ocn_tracer_interior_restoring.o mpas_ocn_tracer_exponential_decay.o mpas_ocn_tracer_ideal_age.o mpas_ocn_tracer_TTD.o mpas_ocn_vmix.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_frazil_forcing.o mpas_ocn_tidal_forcing.o mpas_ocn_tracer_ecosys.o mpas_ocn_tracer_DMS.o mpas_ocn_tracer_MacroMolecules.o mpas_ocn_tracer_CFC.o mpas_ocn_diagnostics.o mpas_ocn_wetting_drying.o mpas_ocn_vel_self_attraction_loading.o  mpas_ocn_vel_tidal_potential.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o mpas_ocn_thick_hadv.o mpas_ocn_thick_vadv.o mpas_ocn_vel_hadv_coriolis.o mpas_ocn_vel_pressure_grad.o mpas_ocn_vel_vadv.o mpas_ocn_vel_hmix.o mpas_ocn_vel_forcing.o mpas_ocn_manufactured_solution.o

mpas_ocn_diagnostics.o: mpas_ocn_thick_ale.o mpas_ocn_equation_of_state.o mpas_ocn_gm.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o mpas_ocn_surface_land_ice_fluxes.o mpas_ocn_vertical_advection.o mpas_ocn_submesoscale_eddies.o mpas_ocn_subgrid.o

mpas_ocn_diagnostics_variables.o: mpas_ocn_config.o

mpas_ocn_mesh.o: mpas_ocn_config.o

mpas_ocn_thick_ale.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_time_average_coupled.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_tracer_ecosys.o mpas_ocn_tracer_DMS.o mpas_ocn_tracer_MacroMolecules.o mpas_ocn_diagnostics_variables.o

mpas_ocn_thick_hadv.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_thick_vadv.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_thick_surface_flux.o: mpas_ocn_forcing.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_gm.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_submesoscale_eddies.o

mpas_ocn_submesoscale_eddies.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o

mpas_ocn_eddy_parameterization_helpers.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o

mpas_ocn_vel_pressure_grad.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_vadv.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_hmix.o: mpas_ocn_vel_hmix_del2.o mpas_ocn_vel_hmix_leith.o mpas_ocn_vel_hmix_del4.o mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_vel_hmix_del2.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_hmix_leith.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_hmix_del4.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_forcing.o: mpas_ocn_vel_forcing_surface_stress.o mpas_ocn_vel_forcing_explicit_bottom_drag.o mpas_ocn_forcing.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_vel_forcing_topographic_wave_drag.o

mpas_ocn_vel_forcing_surface_stress.o: mpas_ocn_forcing.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_vel_forcing_topographic_wave_drag.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o mpas_ocn_forcing.o

mpas_ocn_vel_forcing_explicit_bottom_drag.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_forcing.o mpas_ocn_mesh.o

mpas_ocn_vel_hadv_coriolis.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_tracer_hmix.o: mpas_ocn_tracer_hmix_del2.o mpas_ocn_tracer_hmix_del4.o mpas_ocn_tracer_hmix_redi.o mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_hmix_del2.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_hmix_del4.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_advection.o: mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_tracer_advection_mono.o mpas_ocn_tracer_advection_std.o mpas_ocn_tracer_advection_vert.o mpas_ocn_tracer_advection_shared.o

mpas_ocn_tracer_advection_mono.o: mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o mpas_ocn_tracer_advection_vert.o mpas_ocn_tracer_advection_shared.o

mpas_ocn_tracer_advection_std.o: mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_tracer_advection_vert.o mpas_ocn_tracer_advection_shared.o

mpas_ocn_tracer_advection_vert.o: mpas_ocn_mesh.o mpas_ocn_config.o

mpas_ocn_tracer_advection_shared.o: mpas_ocn_constants.o mpas_ocn_mesh.o mpas_ocn_config.o

mpas_ocn_tracer_hmix_redi.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_high_freq_thickness_hmix_del2.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_tracer_nonlocalflux.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_short_wave_absorption.o: mpas_ocn_tracer_short_wave_absorption_jerlov.o mpas_ocn_tracer_short_wave_absorption_variable.o mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_short_wave_absorption_variable.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_framework_forcing.o

mpas_ocn_tracer_short_wave_absorption_jerlov.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_vmix.o: mpas_ocn_vmix_cvmix.o mpas_ocn_vmix_coefs_redi.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_vmix_gotm.o

mpas_ocn_vmix_cvmix.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o mpas_ocn_stokes_drift.o

mpas_ocn_vmix_coefs_redi.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o

mpas_ocn_vmix_gotm.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o

mpas_ocn_equation_of_state.o: mpas_ocn_equation_of_state_jm.o mpas_ocn_equation_of_state_linear.o mpas_ocn_equation_of_state_wright.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o

mpas_ocn_equation_of_state_jm.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_equation_of_state_linear.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_equation_of_state_wright.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_test.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o

mpas_ocn_constants.o: mpas_ocn_config.o

mpas_ocn_forcing.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_forcing_restoring.o mpas_ocn_diagnostics_variables.o

mpas_ocn_surface_bulk_forcing.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_equation_of_state.o

mpas_ocn_surface_land_ice_fluxes.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_equation_of_state.o mpas_ocn_diagnostics_variables.o

mpas_ocn_frazil_forcing.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o mpas_ocn_equation_of_state.o

mpas_ocn_tidal_forcing.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_equation_of_state.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o mpas_ocn_subgrid.o

mpas_ocn_transport_tests.o: mpas_ocn_config.o

mpas_ocn_effective_density_in_land_ice.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_forcing_restoring.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o

mpas_ocn_tracer_surface_restoring.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_framework_forcing.o mpas_ocn_diagnostics_variables.o

mpas_ocn_tracer_interior_restoring.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_exponential_decay.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_ideal_age.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_TTD.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_ecosys.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_framework_forcing.o mpas_ocn_diagnostics_variables.o

mpas_ocn_tracer_DMS.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_MacroMolecules.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_tracer_CFC.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_framework_forcing.o

mpas_ocn_tracer_surface_flux_to_tend.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_forcing.o

mpas_ocn_time_average_coupled.o: mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_framework_forcing.o:

mpas_ocn_time_varying_forcing.o:  mpas_ocn_framework_forcing.o mpas_ocn_diagnostics_variables.o mpas_ocn_constants.o mpas_ocn_config.o

mpas_ocn_wetting_drying.o: mpas_ocn_diagnostics.o mpas_ocn_gm.o mpas_ocn_diagnostics_variables.o mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_gm.o mpas_ocn_mesh.o

mpas_ocn_tidal_potential_forcing.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o mpas_ocn_diagnostics_variables.o

mpas_ocn_vel_self_attraction_loading.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o

mpas_ocn_vel_pressure_grad.o: mpas_ocn_constants.o

mpas_ocn_vel_tidal_potential.o: mpas_ocn_diagnostics_variables.o

mpas_ocn_vertical_advection.o: mpas_ocn_config.o

mpas_ocn_vertical_regrid.o: mpas_ocn_config.o mpas_ocn_constants.o mpas_ocn_mesh.o

mpas_ocn_vertical_remap.o: mpas_ocn_config.o mpas_ocn_constants.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o mpas_ocn_vertical_regrid.o

mpas_ocn_stokes_drift.o: mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_diagnostics_variables.o mpas_ocn_mesh.o

mpas_ocn_manufactured_solution.o:  mpas_ocn_constants.o mpas_ocn_config.o mpas_ocn_mesh.o 

mpas_ocn_subgrid.o: mpas_ocn_diagnostics_variables.o mpas_ocn_config.o mpas_ocn_constants.o

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -L../global/homes/k/knbarton/repos/repos/MPAS-Model/SAL/src/core_ocean/SHTNS/lib
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -L../global/homes/k/knbarton/repos/repos/MPAS-Model/SAL/src/core_ocean/SHTNS/lib
endif
