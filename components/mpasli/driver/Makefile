.SUFFIXES: .F .o

OBJS = glc_comp_mct.o \
       glc_communicator.o \
       glc_cpl_indices.o \
       glc_mct_vars.o \
       glc_import_export.o

all: $(OBJS)

glc_comp_mct.o: glc_communicator.o glc_cpl_indices.o glc_mct_vars.o glc_import_export.o

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

