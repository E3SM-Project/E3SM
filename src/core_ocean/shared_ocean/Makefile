.SUFFIXES: .F .o

OBJS = mpas_ocn_diagnostics_routines.o

all: sharedocn

sharedocn: $(OBJS)
	ar -ru libsharedocn.a $(OBJS)



clean:
	$(RM) *.o *.mod *.f90 libsharedocn.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
endif
