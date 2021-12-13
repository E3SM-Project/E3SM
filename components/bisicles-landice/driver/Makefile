.SUFFIXES: .F .o

OBJS = \
       mpi_mod.o \
       glimmer_global.o \
       parallel_mpi.o \
       glc_kinds_mod.o \
       glc_communicate.o \
       glimmer_physcon.o \
       glimmer_paramets.o \
       glc_constants.o \
       glc_exit_mod.o \
       glc_route_ice_runoff.o \
       glc_coupling_flags.o \
       glc_fields.o \
       glc_indexing_info.o \
       glc_override_frac.o \
       glc_comp_mct.o \
       glc_communicator.o \
       glc_cpl_indices.o \
       glc_mct_vars.o \
       glc_import_export.o

all: $(OBJS)

glc_comp_mct.o: \
       mpi_mod.o \
       glimmer_global.o \
       parallel_mpi.o \
       glc_kinds_mod.o \
       glc_communicate.o \
       glimmer_physcon.o \
       glimmer_paramets.o \
       glc_constants.o \
       glc_exit_mod.o \
       glc_route_ice_runoff.o \
       glc_coupling_flags.o \
       glc_fields.o \
       glc_indexing_info.o \
       glc_override_frac.o \
       glc_communicator.o \
       glc_cpl_indices.o \
       glc_mct_vars.o \
       glc_import_export.o

mpi_mod.o:

glimmer_global.o:

parallel_mpi.o:
 
glc_kinds_mod.o:

glc_communicate.o:

glimmer_physcon.o:

glimmer_paramets.o:

glc_constants.o:

glc_exit_mod.o:

glc_route_ice_runoff.o:

glc_coupling_flags.o:

glc_fields.o:

glc_indexing_info.o:

glc_override_frac.o:

glc_communicator.o:

glc_cpl_indices.o:

glc_mct_vars.o:

glc_import_export.o:

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

