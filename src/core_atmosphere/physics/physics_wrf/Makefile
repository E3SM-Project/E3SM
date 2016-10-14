.SUFFIXES: .F .o

all: dummy physics_wrf

dummy:
	echo "****** compile physics_wrf ******"

OBJS = \
	libmassv.o                     \
	module_bl_gwdo.o               \
	module_bl_mynn.o               \
	module_bl_mynn.wrf-3.8.o       \
	module_bl_ysu.o                \
	module_cam_error_function.o    \
	module_cam_shr_kind_mod.o      \
	module_cam_support.o           \
	module_cu_gf.mpas.o            \
	module_cu_gf_momentum.mpas.o   \
	module_cu_tiedtke.wrf-3.8.1.o  \
	module_cu_ntiedtke.wrf-3.8.1.o \
	module_cu_kfeta.o              \
	module_cu_tiedtke.o            \
	module_mp_kessler.o            \
	module_mp_radar.o              \
	module_mp_thompson.o           \
	module_mp_thompson_cldfra3.o   \
	module_mp_wsm6.wrf-3.8.1.o     \
	module_ra_cam.o                \
	module_ra_cam_support.o        \
	module_ra_rrtmg_lw.o           \
	module_ra_rrtmg_sw.o           \
	module_sf_bem.o                \
	module_sf_bep.o                \
	module_sf_bep_bem.o            \
	module_sf_sfcdiags.o           \
	module_sf_mynn.o               \
	module_sf_mynn.wrf-3.8.o       \
	module_sf_noahdrv.o            \
	module_sf_noahlsm.o            \
	module_sf_oml.o                \
	module_sf_sfclay.o             \
	module_sf_urban.o

physics_wrf: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_bl_mynn.o: \
	module_cam_error_function.o \
	../mpas_atmphys_constants.o

module_bl_mynn.wrf3.8.o: \
	module_cam_error_function.o \
	../mpas_atmphys_constants.o \
	../mpas_atmphys_utilities.o

module_cam_support.o: \
	module_cam_shr_kind_mod.o  \
	../mpas_atmphys_utilities.o

module_cu_tiedtke.o: \
	../mpas_atmphys_constants.o

module_cu_tiedtke.wrf-3.8.1.o: \
	../mpas_atmphys_constants.o

module_cu_ntiedtke.wrf-3.8.1.o: \
	../mpas_atmphys_constants.o

module_mp_radar.o: \
	../mpas_atmphys_functions.o \
	../mpas_atmphys_utilities.o

module_mp_thompson.o: \
	module_mp_radar.o \
	../mpas_atmphys_functions.o \
	../mpas_atmphys_utilities.o

module_ra_cam.o: \
	module_cam_support.o \
	module_ra_cam_support.o \
	../mpas_atmphys_utilities.o

module_ra_cam_support.o: \
	../mpas_atmphys_utilities.o

module_ra_rrtmg_lw.o: \
	../mpas_atmphys_constants.o

module_ra_rrtmg_sw.o: \
	../mpas_atmphys_constants.o

module_sf_bep.o: \
	module_sf_urban.o

module_sf_bep_bem.o: \
	module_sf_bem.o \
	module_sf_urban.o

module_sf_mynn.o: \
	module_bl_mynn.o \
	module_sf_sfclay.o \
	../mpas_atmphys_constants.o

module_sf_mynn_wrf3.8.o: \
	module_bl_mynn_wrf-3.8.o \
	module_sf_sfclay.o \
	../mpas_atmphys_constants.o

module_sf_noahdrv.o: \
	module_sf_bem.o \
	module_sf_bep.o \
	module_sf_bep_bem.o \
	module_sf_noahlsm.o \
	module_sf_urban.o

module_sf_noahlsm.o: \
	../mpas_atmphys_constants.o \
	../mpas_atmphys_utilities.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I.. -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I.. -I../../../external/esmf_time_f90
endif
