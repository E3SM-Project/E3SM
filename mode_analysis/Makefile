.SUFFIXES: .F .o

OBJS = mpas_ocn_am_mpas_core.o

all: $(OBJS)

mpas_ocn_am_mpas_core.o:

clean:
	$(RM) *.o *.i *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES)
endif
