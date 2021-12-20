.SUFFIXES: .F .o

OBJS = \
        felix_dycore_interface.o \
        glad_constants.o \
        glad_initialise.o \
        glad_input_averages.o \
        glad_io.o \
        glad_main.o \
        glad_mbal_coupling.o \
        glad_mbal_io.o \
        glad_output_fluxes.o \
        glad_output_states.o \
        glad_restart_gcm.o \
        glad_timestep.o \
        glad_type.o \
        glam_grid_operators.o \
        glam_strs2.o \
        glam_velo.o \
        glc_broadcast.o \
        glc_communicate.o \
        glc_comp_mct.o \
        glc_constants.o \
        glc_coupling_flags.o \
        glc_cpl_indices.o \
        glc_ensemble.o \
        glc_ErrorMod.o \
        glc_exit_mod.o \
        glc_fields.o \
        glc_files.o \
        glc_FinalMod.o \
        glc_history.o \
        glc_import_export.o \
        glc_io.o \
        glc_InitMod.o \
        glc_indexing_info.o \
        glc_kinds_mod.o \
        glc_override_frac.o \
        glc_route_ice_runoff.o \
        glc_RunMod.o \
        glc_time_management.o \
        glide.o \
        glide_bwater.o \
        glide_calc_sigma.o \
        glide_diagnostics.o \
        glide_grid_operators.o \
        glide_ground.o \
        glide_io.o \
        glide_lithot.o \
        glide_lithot1d.o \
        glide_lithot3d.o \
        glide_mask.o \
        glide_nc_custom.o \
        glide_nonlin.o \
        glide_profile.o \
        glide_setup.o \
        glide_stop.o \
        glide_stress.o \
        glide_temp.o \
        glide_thck.o \
        glide_types.o \
        glide_velo.o \
        glimmer_config.o \
        glimmer_coordinates.o \
        glimmer_filenames.o \
        glimmer_global.o \
        glimmer_log.o \
        glimmer_map_CFproj.o \
        glimmer_map_init.o \
        glimmer_map_types.o \
        glimmer_ncdf.o \
        glimmer_ncio.o \
        glimmer_ncparams.o \
        glimmer_paramets.o \
        glimmer_physcon.o \
        glimmer_scales.o \
        glimmer_sparse.o \
        glimmer_sparse_slap.o \
        glimmer_sparse_type.o \
        glimmer_utils.o \
        glimmer_vers.o \
        glissade.o \
        glissade_basal_traction.o \
        glissade_calving.o \
        glissade_enthalpy.o \
        glissade_grid_operators.o \
        glissade_masks.o \
        glissade_remap.o \
        glissade_temp.o \
        glissade_test.o \
        glissade_therm.o \
        glissade_transport.o \
        glissade_velo.o \
        glissade_velo_higher.o \
        glissade_velo_higher_pcg.o \
        glissade_velo_higher_slap.o \
        glissade_velo_sia.o \
        history_tape_base.o \
        history_tape_coupler.o \
        history_tape_standard.o \
        isostasy.o \
        isostasy_elastic.o \
        isostasy_kelvin.o \
        mpi_mod.o \
        parallel_mpi.o \
        profile.o

all: $(OBJS)

felix_dycore_interface.o: \
        glide_mask.o \
        glide_types.o \
        glimmer_global.o \
        glimmer_log.o \
        glimmer_paramets.o \
        glimmer_physcon.o \
        glissade_grid_operators.o

glad_constants.o: \
        glimmer_global.o \
        glimmer_physcon.o

glad_initialise.o: \
        glad_constants.o \
        glad_io.o \
        glad_mbal_coupling.o \
        glad_mbal_io.o \
        glad_restart_gcm.o \
        glad_type.o \
        glide.o \
        glide_diagnostics.o \
        glide_nc_custom.o \
        glimmer_config.o \
        glimmer_coordinates.o \
        glimmer_global.o \
        glimmer_log.o \
        glimmer_ncio.o \
        glimmer_paramets.o \
        glissade.o \
        parallel_mpi.o

glad_input_averages.o: \
        glimmer_paramets.o

glad_io.o: \
        glad_type.o \
        glide_types.o \
        glimmer_map_types.o \
        glimmer_ncdf.o \
        glimmer_ncio.o \
        glimmer_paramets.o \
        glimmer_scales.o

glad_main.o: \
        glad_initialise.o \
        glad_timestep.o \
        glimmer_log.o \
        glad_output_states.o \
        glimmer_paramets.o \
        glimmer_utils.o \
        parallel_mpi.o

glad_mbal_coupling.o: \
        glimmer_config.o \
        glimmer_coordinates.o

glad_mbal_io.o:

glad_output_fluxes.o: \
        glide_types.o

glad_output_states.o:

glad_restart_gcm.o:

glad_timestep.o:

glad_type.o: \
        glad_input_averages.o \
        glad_mbal_coupling.o \
        glad_output_fluxes.o \
        glide_types.o \
        glimmer_ncparams.o

glam_grid_operators.o: \
        glide_grid_operators.o

glam_strs2.o: \
        glissade_basal_traction.o

glam_velo.o: \
        glide_stress.o

glc_broadcast.o: \
        glc_communicate.o \
        glc_kinds_mod.o

glc_communicate.o: \
        glc_kinds_mod.o \
        parallel_mpi.o

glc_comp_mct.o: \
        glad_main.o \
        glc_communicate.o \
        glc_constants.o \
        glc_coupling_flags.o \
        glc_cpl_indices.o \
        glc_ErrorMod.o \
        glc_fields.o \
        glc_FinalMod.o \
        glc_import_export.o \
        glc_indexing_info.o \
        glc_InitMod.o \
        glc_io.o \
        glc_RunMod.o

glc_constants.o: \
        glc_kinds_mod.o \
        glimmer_paramets.o

glc_coupling_flags.o: \
        glc_constants.o \
        glc_exit_mod.o \
        glc_kinds_mod.o \
        glc_route_ice_runoff.o

glc_cpl_indices.o: \
        glc_constants.o

glc_ensemble.o:

glc_ErrorMod.o:

glc_exit_mod.o: \
        glc_communicate.o \
        glc_constants.o \
        glc_kinds_mod.o 

glc_fields.o: \
        glad_main.o \
        glc_kinds_mod.o

glc_files.o: \
        glc_ensemble.o

glc_FinalMod.o:

glc_history.o: \
        history_tape_base.o \
        history_tape_coupler.o \
        history_tape_standard.o

glc_import_export.o: \
        glc_communicate.o \
        glc_constants.o \
        glc_cpl_indices.o \
        glc_fields.o \
        glc_indexing_info.o \
        glc_override_frac.o \
        glc_route_ice_runoff.o

glc_indexing_info.o: \
        glad_main.o

glc_InitMod.o: \
        glc_history.o \
        glc_io.o

glc_io.o:

glc_kinds_mod.o:

glc_override_frac.o: \
        glc_broadcast.o \
        glc_communicate.o \
        glc_constants.o \
        glc_exit_mod.o \
        glc_files.o \
        glc_kinds_mod.o \
        glc_time_management.o

glc_route_ice_runoff.o: \
        glc_constants.o \
        glc_exit_mod.o \
        glc_kinds_mod.o 

glc_RunMod.o:

glc_time_management.o:

glide.o: \
        glide_bwater.o \
        glide_diagnostics.o \
        glide_grid_operators.o \
        glide_ground.o \
        glide_io.o \
        glide_lithot.o \
        glide_mask.o \
        glide_profile.o \
        glide_setup.o \
        glide_stop.o \
        glide_temp.o \
        glide_thck.o \
        glide_velo.o \
        isostasy.o

glide_bwater.o: \
        glam_grid_operators.o \
        glide_grid_operators.o \
        glissade_grid_operators.o

glide_calc_sigma.o:

glide_diagnostics.o:

glide_grid_operators.o:

glide_ground.o:

glide_io.o:

glide_lithot.o: \
        glide_calc_sigma.o \
        glide_lithot1d.o \
        glide_lithot3d.o \
        glide_mask.o \
        glide_setup.o

glide_lithot1d.o: \
        glimmer_utils.o

glide_lithot3d.o: \
        glide_stop.o

glide_mask.o: \
        glide_types.o \
        glimmer_global.o \
        glimmer_physcon.o \
        parallel_mpi.o

glide_nc_custom.o:

glide_nonlin.o:

glide_profile.o:

glide_setup.o:

glide_stop.o:

glide_stress.o:

glide_temp.o:

glide_thck.o: \
        glide_nonlin.o \
        glide_velo.o \
        glimmer_sparse.o

glide_types.o: \
        glimmer_coordinates.o \
        glimmer_map_types.o \
        glimmer_ncdf.o \
        glimmer_paramets.o \
        glimmer_sparse_type.o \
        profile.o

glide_velo.o:

glimmer_config.o: \
        glimmer_log.o \
        parallel_mpi.o

glimmer_coordinates.o: \
        glimmer_global.o \
        glimmer_log.o

glimmer_filenames.o:

glimmer_global.o:

glimmer_log.o: \
        parallel_mpi.o

glimmer_map_CFproj.o: \
        glimmer_map_init.o

glimmer_map_init.o:

glimmer_map_types.o: \
        glimmer_global.o \
        glimmer_log.o \
        glimmer_physcon.o

glimmer_ncdf.o:

glimmer_ncio.o: \
        glimmer_map_CFproj.o

glimmer_ncparams.o: \
        glimmer_config.o \
        glimmer_filenames.o

glimmer_paramets.o: \
        glimmer_global.o \
        glimmer_physcon.o

glimmer_physcon.o: \
        glimmer_global.o

glimmer_scales.o:

glimmer_sparse.o: \
        glimmer_sparse_slap.o

glimmer_sparse_slap.o:

glimmer_sparse_type.o:

glimmer_utils.o:

glimmer_vers.o:

glissade.o: \
        felix_dycore_interface.o \
        glam_strs2.o \
        glam_velo.o \
        glissade_calving.o \
        glissade_enthalpy.o \
        glissade_masks.o \
        glissade_temp.o \
        glissade_test.o \
        glissade_therm.o \
        glissade_transport.o \
        glissade_velo.o \
        glissade_velo_higher.o

glissade_basal_traction.o:

glissade_calving.o: \
        glissade_masks.o

glissade_enthalpy.o:

glissade_grid_operators.o:

glissade_masks.o:

glissade_remap.o:

glissade_temp.o:

glissade_test.o: \
        glissade_transport.o

glissade_therm.o:

glissade_transport.o: \
        glissade_remap.o

glissade_velo.o: \
        glissade_velo_higher.o \
        glissade_velo_sia.o

glissade_velo_higher.o: \
        glissade_velo_higher_pcg.o \
        glissade_velo_higher_slap.o

glissade_velo_higher_pcg.o:

glissade_velo_higher_slap.o:

glissade_velo_sia.o:

history_tape_base.o: \
        glc_io.o

history_tape_coupler.o:

history_tape_standard.o:

isostasy.o: \
        isostasy_elastic.o

isostasy_elastic.o: \
        isostasy_kelvin.o

isostasy_kelvin.o:

mpi_mod.o:

parallel_mpi.o: \
        glimmer_global.o \
        mpi_mod.o

profile.o:

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
endif
