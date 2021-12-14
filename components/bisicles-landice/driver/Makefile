.SUFFIXES: .F .o

OBJS = \
       glad_constants.o \
       glad_initialise.o \
       glad_io.o \
       glad_main.o \
       glad_mbal_coupling.o \
       glad_restart_gcm.o \
       glad_timestep.o \
       glad_type.o \
       glc_communicate.o \
       glc_comp_mct.o \
       glc_constants.o \
       glc_coupling_flags.o \
       glc_cpl_indices.o \
       glc_exit_mod.o \
       glc_fields.o \
       glc_import_export.o \
       glc_indexing_info.o \
       glc_kinds_mod.o \
       glc_override_frac.o \
       glc_route_ice_runoff.o \
       glide.o \
       glide_diagnostics.o \
       glide_nc_custom.o \
       glide_types.o \
       glimmer_config.o \
       glimmer_coordinates.o \
       glimmer_global.o \
       glimmer_log.o \
       glimmer_map_types.o \
       glimmer_ncdf.o \
       glimmer_ncio.o \
       glimmer_paramets.o \
       glimmer_physcon.o \
       glimmer_scales.o \
       glimmer_util.o \
       glissade.o \
       mpi_mod.o \
       parallel_mpi.o

all: $(OBJS)

glad_constants.o: \
       glimmer_global.o \
       glimmer_physcon.o

glad_initialise.o: \
       glad_constants.o \
       glad_mbal_coupling.o \
       glad_io.o \
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
       glimmer_paramets.o \
       glimmer_util.o \
       parallel_mpi.o

glad_mbal_coupling.o: \
       glimmer_config.o \
       glimmer_coordinates.o

glad_restart_gcm.o:

glad_timestep.o:

glad_type.o:

glc_communicate.o: \
       glc_kinds_mod.o \
       parallel_mpi.o

glc_comp_mct.o: \
       glad_main.o \
       glc_communicate.o \
       glc_constants.o \
       glc_coupling_flags.o \
       glc_cpl_indices.o \
       glc_fields.o \
       glc_import_export.o \
       glc_indexing_info.o

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

glc_exit_mod.o: \
       glc_communicate.o \
       glc_constants.o \
       glc_kinds_mod.o 

glc_fields.o: \
       glad_main.o \
       glc_kinds_mod.o

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

glc_kinds_mod.o:

glc_override_frac.o: \
       glc_communicate.o \
       glc_constants.o \
       glc_exit_mod.o \
       glc_kinds_mod.o

glc_route_ice_runoff.o: \
       glc_constants.o \
       glc_exit_mod.o \
       glc_kinds_mod.o 

glide.o:

glide_diagnostics.o:

glide_nc_custom.o:

glide_types.o:

glimmer_config.o: \
       glimmer_log.o \
       parallel_mpi.o

glimmer_coordinates.o:

glimmer_global.o:

glimmer_log.o: \
       parallel_mpi.o

glimmer_map_types.o:

glimmer_ncdf.o:

glimmer_ncio.o:

glimmer_paramets.o: \
       glimmer_global.o \
       glimmer_physcon.o

glimmer_physcon.o: \
       glimmer_global.o

glimmer_scales.o:

glimmer_util.o:

glissade.o:

mpi_mod.o:

parallel_mpi.o: \
       glimmer_global.o \
       mpi_mod.o

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
