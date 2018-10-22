.SUFFIXES: .F .o

OBJS = ocn_comp_mct.o \
       mpaso_cpl_indices.o \
	   mpaso_mct_vars.o

OCEAN_SHARED_INCLUDES=-I../core_ocean/mode_forward -I../core_ocean/shared -I../core_ocean/analysis_members -I../core_ocean/cvmix -I../framework -I../operators

all: $(OBJS)

ocn_comp_mct.o: mpaso_cpl_indices.o mpaso_mct_vars.o

mpaso_cpl_indices.o:

mpaso_mct_vars.o:

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
