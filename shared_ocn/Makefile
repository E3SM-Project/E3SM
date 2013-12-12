.SUFFIXES: .F .o

OBJS = mpas_ocn_diagnostics_routines.o \
       mpas_ocn_init.o 

all: shared_ocn

shared_ocn: $(OBJS)
	ar -ru libshared_ocn.a $(OBJS)

clean:
	$(RM) *.o *.mod *.f90 lib*.a 

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
endif
