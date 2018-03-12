.SUFFIXES: .F .o

OBJS = 	mpas_cice_core.o \
	mpas_cice_core_interface.o

all: $(OBJS)

mpas_cice_core.o: 

mpas_cice_core_interface.o: mpas_cice_core.o

clean:
	$(RM) *.o *.i *.mod *.f90

ifneq (,$(findstring CPRIBM,$(CPPFLAGS)))
FFLAGS_noSMP := $(filter-out -qsmp%,$(FFLAGS))
else
FFLAGS_noSMP := $(FFLAGS)
endif

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90

	$(FC) $(FFLAGS_noSMP) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column -I../shared
else
	$(FC) $(CPPFLAGS) $(FFLAGS_noSMP) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column -I../shared
endif
