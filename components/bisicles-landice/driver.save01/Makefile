.SUFFIXES: .F .o

OBJS = glc_comp_mct_bisicles.o \
       glc_cpl_indices_bisicles.o \
       glc_mct_vars_bisicles.o \

all: $(OBJS)

glc_comp_mct_bisicles.o: glc_cpl_indices_bisicles.o glc_mct_vars_bisicles.o 

glc_cpl_indices_bisicles.o:

glc_mct_vars_bisicles.o:

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

