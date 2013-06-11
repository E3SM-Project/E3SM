.SUFFIXES: .F .o

all: dummy physics_wrf

dummy:
	echo "****** compile physics_wrf ******"

OBJS = \
	libmassv.o                 \
	module_bl_gwdo.o           \
	module_bl_ysu.o            \
	module_cam_shr_kind_mod.o  \
	module_cam_support.o       \
	module_cu_kfeta.o          \
	module_cu_tiedtke.o        \
	module_mp_kessler.o        \
	module_mp_radar.o          \
	module_mp_wsm6.o           \
	module_ra_cam.o            \
	module_ra_cam_support.o    \
	module_ra_rrtmg_lw.o       \
	module_ra_rrtmg_sw.o       \
	module_sf_bem.o            \
	module_sf_bep.o            \
	module_sf_bep_bem.o        \
	module_sf_noahdrv.o        \
	module_sf_noahlsm.o        \
	module_sf_sfclay.o         \
	module_sf_urban.o

physics_wrf: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_cam_support.o: \
	module_cam_shr_kind_mod.o  \
	../mpas_atmphys_utilities.o

module_mp_radar.o: \
	../mpas_atmphys_utilities.o

module_mp_wsm6.o: \
	libmassv.o \
	module_mp_radar.o

module_ra_cam.o: \
	module_cam_support.o \
	module_ra_cam_support.o \
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

module_sf_noahdrv.o: \
 	module_sf_bem.o \
	module_sf_bep.o \
	module_sf_bep_bem.o \
	module_sf_noahlsm.o \
	module_sf_urban.o

module_sf_noahlsm.o: \
	../mpas_atmphys_constants.o

clean:
	$(RM) *.f90 *.o *.mod

.F.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..
