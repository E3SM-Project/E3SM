.SUFFIXES: .F .o

#
# Add new diagnostic modules to DIAGNOSTIC_MODULES
#
DIAGNOSTIC_MODULES = \
	mpas_atm_diagnostic_template.o


OBJS = mpas_atm_diagnostics_manager.o

all: $(DIAGNOSTIC_MODULS) $(OBJS)

mpas_atm_diagnostics_manager.o: $(DIAGNOSTIC_MODULES)


clean:
	$(RM) *.o *.mod *.f90
	@# Some Intel compilers generate *.i files; clean them up, too
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../dynamics -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../dynamics -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
endif
