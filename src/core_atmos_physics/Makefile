.SUFFIXES: .F .o

ifeq ($(CORE),hyd_atmos)
COREDEF = -Dhydrostatic_core
endif
ifeq ($(CORE),atmosphere)
COREDEF = -Dnon_hydrostatic_core 
endif
HYDROSTATIC = -Ddo_hydrostatic_pressure
#HYDROSTATIC =

dummy:
	echo "****** make non-hydrostatic core ******"

OBJS_init = \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

OBJS = \
	mpas_atmphys_driver_cloudiness.o      \
	mpas_atmphys_driver_convection_deep.o \
	mpas_atmphys_driver_gwdo.o            \
	mpas_atmphys_driver_lsm.o             \
	mpas_atmphys_driver_microphysics.o    \
	mpas_atmphys_driver_radiation_lw.o    \
	mpas_atmphys_driver_radiation_sw.o    \
	mpas_atmphys_driver_sfclayer.o        \
	mpas_atmphys_driver_pbl.o             \
	mpas_atmphys_driver.o                 \
	mpas_atmphys_camrad_init.o            \
	mpas_atmphys_control.o                \
	mpas_atmphys_date_time.o              \
	mpas_atmphys_init.o                   \
	mpas_atmphys_landuse.o                \
	mpas_atmphys_lsm_noahinit.o           \
	mpas_atmphys_manager.o                \
	mpas_atmphys_rrtmg_lwinit.o           \
	mpas_atmphys_rrtmg_swinit.o           \
	mpas_atmphys_todynamics.o             \
	mpas_atmphys_update_surface.o         \
	mpas_atmphys_update.o                 \
	mpas_atmphys_vars.o

ifeq ($(CORE),hyd_atmos)
	OBJS_dyn = mpas_atmphys_interface_hyd.o
else
	OBJS_dyn = mpas_atmphys_interface_nhyd.o
endif

all: lookup_tables core_physics_init core_physics_wrf core_dyn core_physics

lookup_tables:
	./checkout_data_files.sh

core_physics_wrf:
	(cd physics_wrf; make all COREDEF="$(COREDEF)")

core_dyn: $(OBJS_dyn)
	ar -ru libphys.a $(OBJS_dyn)

core_physics_init: $(OBJS_init)
	ar -ru libphys.a $(OBJS_init)

core_physics: $(OBJS)
	ar -ru libphys.a $(OBJS)

# DEPENDENCIES:
mpas_atmphys_driver_cloudines.o: \
	mpas_atmphys_driver_cloudiness.o    \
	mpas_atmphys_vars.o

mpas_atmphys_driver_convection_deep.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_utilities.o            \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_cu_kfeta.o     \
	./physics_wrf/module_cu_tiedtke.o

mpas_atmphys_driver_gwdo.o: \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_bl_gwdo.o

mpas_atmphys_driver_lsm.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_landuse.o              \
	mpas_atmphys_lsm_noahinit.o         \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_sf_noahdrv.o

mpas_atmphys_driver_pbl.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_bl_ysu.o

mpas_atmphys_driver_radiation_lw.o: \
	mpas_atmphys_driver_radiation_sw.o  \
	mpas_atmphys_camrad_init.o          \
	mpas_atmphys_constants.o            \
	mpas_atmphys_manager.o              \
	mpas_atmphys_rrtmg_lwinit.o         \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_ra_cam.o       \
	./physics_wrf/module_ra_rrtmg_lw.o

mpas_atmphys_driver_radiation_sw.o: \
	mpas_atmphys_camrad_init.o          \
	mpas_atmphys_constants.o            \
	mpas_atmphys_manager.o              \
	mpas_atmphys_rrtmg_swinit.o         \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_ra_cam.o       \
	./physics_wrf/module_ra_rrtmg_sw.o

mpas_atmphys_driver_sfclayer.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_vars.o                 \
	./physics_wrf/module_sf_sfclay.o

mpas_atmphys_camrad_init.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_utilities.o            \
	./physics_wrf/module_ra_cam_support.o

mpas_atmphys_control.o: \
	mpas_atmphys_utilities.o

mpas_atmphys_init.o: \
	mpas_atmphys_driver_convection_deep.o       \
	mpas_atmphys_driver_lsm.o                   \
	mpas_atmphys_driver_microphysics.o          \
	mpas_atmphys_driver_radiation_lw.o          \
	mpas_atmphys_driver_radiation_sw.o          \
	mpas_atmphys_driver_sfclayer.o              \
	mpas_atmphys_landuse.o

mpas_atmphys_landuse.o: \
	mpas_atmphys_utilities.o            \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahinit.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_utilities.o            \
	./physics_wrf/module_sf_noahlsm.o

mpas_atmphys_manager.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_vars.o                 \
	mpas_atmphys_update_surface.o       \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_lwinit.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_utilities.o            \
	./physics_wrf/module_ra_rrtmg_lw.o

mpas_atmphys_rrtmg_swinit.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_utilities.o            \
	./physics_wrf/module_ra_rrtmg_sw.o

mpas_atmphys_todynamics.o: \
	mpas_atmphys_constants.o

mpas_atmphys_update_surface.o: \
	mpas_atmphys_date_time.o            \
	mpas_atmphys_constants.o            \
	mpas_atmphys_landuse.o              \
	mpas_atmphys_vars.o

mpas_atmphys_update.o: \
	mpas_atmphys_driver_convection_deep.o

# CORE DEPENDENT DEPENDENCIES:
ifeq ($(CORE),hyd_atmos)

mpas_atmphys_interface_hyd.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_vars.o

module_driver_microphysics.o: \
	./physics_wrf/module_mp_kessler.o     \
	./physics_wrf/module_mp_wsm6.o        \
	mpas_atmphys_constants.o            \
	mpas_atmphys_interface_hyd.o        \
	mpas_atmphys_vars.o

mpas_atmphys_driver.o: \
	mpas_atmphys_driver_convection_deep.o       \
	mpas_atmphys_driver_lsm.o                   \
	mpas_atmphys_driver_pbl.o                   \
	mpas_atmphys_driver_sfclayer.o              \
	mpas_atmphys_constants.o            \
	mpas_atmphys_interface_hyd.o        \
	mpas_atmphys_vars.o

else

mpas_atmphys_interface_nhyd.o: \
	mpas_atmphys_constants.o            \
	mpas_atmphys_vars.o

mpas_atmphys_driver_microphysics.o: \
	./physics_wrf/module_mp_kessler.o   \
	./physics_wrf/module_mp_wsm6.o      \
	mpas_atmphys_constants.o            \
	mpas_atmphys_interface_nhyd.o       \
	mpas_atmphys_vars.o

mpas_atmphys_driver.o: \
	mpas_atmphys_driver_convection_deep.o \
	mpas_atmphys_driver_gwdo.o            \
	mpas_atmphys_driver_pbl.o             \
	mpas_atmphys_driver_radiation_lw.o    \
	mpas_atmphys_driver_radiation_sw.o    \
	mpas_atmphys_driver_sfclayer.o        \
	mpas_atmphys_constants.o              \
	mpas_atmphys_interface_nhyd.o         \
	mpas_atmphys_update.o                 \
	mpas_atmphys_vars.o

mpas_atmphys_update.o: \
	mpas_atmphys_driver_convection_deep.o \
	mpas_atmphys_vars.o

endif

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_wrf; make clean )

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES)  -DIWORDSIZE=4 -DRWORDSIZE=8 $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics_wrf -I../external/esmf_time_f90
