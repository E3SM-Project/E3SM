.SUFFIXES: .F90 .o

OBJS =  ice_colpkg.o \
	ice_kinds_mod.o \
	ice_warnings.o \
	ice_colpkg_shared.o \
	ice_constants_colpkg.o \
	ice_therm_shared.o \
	ice_orbital.o \
	ice_mushy_physics.o \
	ice_therm_mushy.o \
	ice_atmo.o \
	ice_age.o \
	ice_firstyear.o \
	ice_flux_colpkg.o \
	ice_meltpond_cesm.o \
	ice_meltpond_lvl.o \
	ice_meltpond_topo.o \
	ice_therm_vertical.o \
	ice_therm_bl99.o \
	ice_therm_0layer.o \
	ice_itd.o \
	ice_colpkg_tracers.o \
	ice_therm_itd.o \
	ice_shortwave.o \
	ice_mechred.o \
        ice_aerosol.o \
        ice_brine.o \
        ice_algae.o \
        ice_zbgc.o \
        ice_zbgc_shared.o

all: $(OBJS)

ifneq "$(ESM)" ""
ice_constants_colpkg.F90:
	cp constants/cesm/ice_constants_colpkg.F90 .
else
ice_constants_colpkg.F90:
	cp constants/cice/ice_constants_colpkg.F90 .
endif

ice_colpkg.o: ice_kinds_mod.o ice_constants_colpkg.o ice_warnings.o ice_colpkg_shared.o ice_therm_shared.o ice_orbital.o ice_atmo.o ice_age.o ice_firstyear.o ice_flux_colpkg.o ice_meltpond_cesm.o ice_meltpond_lvl.o ice_meltpond_topo.o ice_therm_vertical.o ice_itd.o ice_therm_itd.o ice_shortwave.o ice_mechred.o ice_colpkg_tracers.o ice_atmo.o ice_mushy_physics.o ice_zbgc.o ice_zbgc_shared.o  ice_aerosol.o ice_algae.o ice_brine.o ice_zsalinity.o

ice_kinds_mod.o:

ice_warnings.o: ice_kinds_mod.o

ice_constants_colpkg.o: ice_constants_colpkg.F90 ice_kinds_mod.o

ice_colpkg_shared.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_mushy_physics.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_therm_shared.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_orbital.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_therm_mushy.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_therm_shared.o ice_mushy_physics.o ice_colpkg_tracers.o

ice_atmo.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_tracers.o

ice_age.o: ice_kinds_mod.o

ice_firstyear.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_flux_colpkg.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_meltpond_cesm.o: ice_kinds_mod.o ice_constants_colpkg.o

ice_meltpond_lvl.o: ice_kinds_mod.o ice_constants_colpkg.o ice_therm_shared.o

ice_meltpond_topo.o: ice_kinds_mod.o ice_constants_colpkg.o ice_therm_shared.o

ice_therm_vertical.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_therm_shared.o ice_therm_bl99.o ice_therm_0layer.o ice_therm_mushy.o ice_mushy_physics.o

ice_therm_bl99.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_therm_shared.o

ice_therm_0layer.o: ice_kinds_mod.o ice_constants_colpkg.o ice_therm_bl99.o

ice_itd.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_tracers.o ice_therm_shared.o ice_zbgc_shared.o

ice_colpkg_tracers.o: ice_kinds_mod.o ice_constants_colpkg.o 

ice_therm_itd.o: ice_kinds_mod.o ice_constants_colpkg.o ice_itd.o ice_colpkg_tracers.o ice_therm_shared.o ice_therm_shared.o  ice_mushy_physics.o ice_zbgc.o ice_zbgc_shared.o

ice_shortwave.o: ice_kinds_mod.o ice_constants_colpkg.o ice_orbital.o

ice_mechred.o: ice_kinds_mod.o ice_constants_colpkg.o ice_itd.o ice_colpkg_tracers.o

ice_aerosol.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_colpkg_tracers.o ice_zbgc_shared.o 

ice_algae.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_colpkg_tracers.o ice_zbgc_shared.o ice_aerosol.o

ice_brine.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_colpkg_tracers.o ice_zbgc_shared.o ice_therm_mushy.o ice_mushy_physics.o ice_therm_shared.o

ice_zbgc.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o ice_colpkg_tracers.o ice_zbgc_shared.o ice_therm_shared.o ice_itd.o 

ice_zbgc_shared.o: ice_kinds_mod.o ice_constants_colpkg.o ice_colpkg_shared.o 


.F90.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90

	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES)
endif

clean:
	$(RM) *.o *.mod *.i90
	rm -f ice_constants_colpkg.F90
