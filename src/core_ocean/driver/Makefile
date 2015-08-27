.SUFFIXES: .F .o

OBJS = mpas_ocn_core.o \
       mpas_ocn_core_interface.o

all: core_ocean

core_ocean: $(OBJS) 

mpas_ocn_core.o:

mpas_ocn_core_interface.o: mpas_ocn_core.o

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES)
endif
