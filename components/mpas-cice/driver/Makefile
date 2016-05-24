.SUFFIXES: .F .o

OBJS = ice_comp_mct.o \
       mpas_cice_cpl_indices.o \
       mpas_cice_mct_vars.o

ICE_SHARED_INCLUDES=-I../core_cice/column -I../core_cice/forward_model -I../core_cice/analysis_members -I../framework -I../operators

all: $(OBJS)

ice_comp_mct.o: mpas_cice_cpl_indices.o mpas_cice_mct_vars.o

mpas_cice_cpl_indices.o:

mpas_cice_mct_vars.o:

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
