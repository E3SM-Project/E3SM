.SUFFIXES: .F .o

OBJS = \
        glad_main.o \
        glad_type.o \
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
        glc_import_export.o \
        glc_InitMod.o \
        glc_indexing_info.o \
        glc_kinds_mod.o \
        glc_override_frac.o \
        glc_route_ice_runoff.o \
        glc_RunMod.o \
        glc_time_management.o \
        glimmer_global.o \
        glimmer_log.o \
        glimmer_paramets.o \
        glimmer_physcon.o \
        mpi_mod.o \
        parallel_mpi.o \
        profile.o

all: $(OBJS)

glad_main.o: \
        glad_type.o \
        glimmer_paramets.o \
        parallel_mpi.o

glad_type.o:

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

glc_InitMod.o:

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

glimmer_global.o:

glimmer_log.o: \
        parallel_mpi.o

glimmer_paramets.o: \
        glimmer_global.o \
        glimmer_physcon.o

glimmer_physcon.o: \
        glimmer_global.o

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
