.SUFFIXES:
.SUFFIXES: .F90 .o 

ifneq ('$(OBJECT_MODE)','32')
  export OBJECT_MODE=64
endif

SRC=grid_init_mod.F90 moist_init_mod.F90 prof_init_mod.F90 absor_mod.F90 \
    noise_mod.F90 surfflux_mod.F90 \
    gcrk_mod.F90 diag_mod.F90 velprd_mod.F90 rad_mod.F90 mpdat_mod.F90 \
    thermo_mod.F90 rain_mod.F90 \
    integxz_mod.F90 zero_mean_mod.F90 crcp_mod.F90 main.F90

OBJ=$(SRC:.F90=.o)


default:
	@echo "Usage: $(MAKE) [ sun | ibm | ibmdebug | blrts | intel | pgi]"

sun:
	(cd timers && $(MAKE) clean default CFLAGS='-DUNDERSCORE -O')
	@$(MAKE) FC=f77 LD=ar FFLAGS="-xO5 -xtypemap=real:64,double:64,integer:32 -r8const" crcp norm

ibmdebug: 
	(cd timers && $(MAKE) clean ibm)
	@$(MAKE)  FC=xlf90_r LD=ar LDFLAGS='-lmass -lmassv' FFLAGS="-qarch=auto  -g -qfullpath -qrealsize=8 -qsuffix=cpp=F90 -q$(OBJECT_MODE) -qsigtrap -qflttrap=en:inv:zero -qinitauto=7FF7FFFF -WF,-D_AIX -WF,-DTESTMODE" crcp norm

ibm: 
	(cd timers && $(MAKE) clean ibm)
	@$(MAKE) FC=xlf90_r LD=ar LDFLAGS='-lmass -lmassv' FFLAGS="-O3 -qdpc -g -qfullpath -qrealsize=8 -qsuffix=cpp=F90 -q$(OBJECT_MODE) -qdpc -WF,-D_AIX -WF,-DTESTMODE" crcp norm 


blrts: 
	(cd timers && $(MAKE) clean blrts)
	@$(MAKE) FC=blrts_xlf LD=ar LDFLAGS='-lmass -lmassv' FFLAGS="-O3 -qarch=440 -qautodbl=dbl4 -qdpc -qsuffix=cpp=F90 -WF,-D_BGL -WF,-DTESTMODE" crcp norm

intel:
	(cd timers && $(MAKE) clean default)
	@$(MAKE) FC=ifort LD=ar FFLAGS="-O3 -real_size 64" crcp norm

pgi:
	(cd timers && $(MAKE) clean default CFLAGS='-DUNDERSCORE -O')
	@$(MAKE) FC=pgf90 LD=ar FFLAGS="-Mfree -fast -Kieee -r8" crcp norm

$(OBJ): Makefile

.F90.o:
	$(FC) -c $(FFLAGS) $<

crcp: $(OBJ)
	$(FC) -o $@ $(FFLAGS) $(OBJ) timers/libtimers.a $(LDFLAGS)

norm: norm.F90
	$(FC) -o $@ $(FFLAGS) $<

clean:
	rm -f crcp norm $(OBJ) *.mod
	(cd timers && $(MAKE) clean)

ifeq ($(TARGET),blrts)
mpdat_mod.o : mpdat_mod.F90
	$(FC) -c $(FFLAGS)  -WF,-DUSE_FSEL -qlist -qsource  $<

gcrk_mod.o : gcrk_mod.F90
	$(FC) -c $(FFLAGS) -qhot=novector -qlist -qsource  $<

laplc_mod.o : laplc_mod.F90
	$(FC) -c $(FFLAGS) -qhot=novector -qlist -qsource  $<

thermo_mod.o : thermo_mod.F90
	$(FC) -c $(FFLAGS) -qhot=novector -qlist -qsource  $<
endif
