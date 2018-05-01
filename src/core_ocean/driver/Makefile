.SUFFIXES: .F .o

OBJS = mpas_ocn_core.o \
       mpas_ocn_core_interface.o

all: core_ocean

core_ocean: $(OBJS) 

mpas_ocn_core.o:

mpas_ocn_core_interface.o: mpas_ocn_core.o

clean:
	$(RM) *.o *.mod *.f90

ifneq (,$(findstring CPRIBM,$(CPPFLAGS)))
FFLAGS_noSMP := $(filter-out -qsmp%,$(FFLAGS))
else
FFLAGS_noSMP := $(FFLAGS)
endif

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS_noSMP) -c $*.f90 $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS_noSMP) -c $*.F $(CPPINCLUDES) $(FCINCLUDES)
endif
