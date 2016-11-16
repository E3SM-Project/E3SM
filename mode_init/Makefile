.SUFFIXES: .F .o

OBJS = mpas_ocn_init_mode.o

UTILS = mpas_ocn_init_spherical_utils.o \
        mpas_ocn_init_vertical_grids.o \
        mpas_ocn_init_cell_markers.o \
        mpas_ocn_init_interpolation.o \
	mpas_ocn_init_ssh_and_landIcePressure.o

TEST_CASES = mpas_ocn_init_baroclinic_channel.o \
             mpas_ocn_init_lock_exchange.o \
             mpas_ocn_init_internal_waves.o \
             mpas_ocn_init_overflow.o \
             mpas_ocn_init_cvmix_WSwSBF.o \
             mpas_ocn_init_iso.o \
             mpas_ocn_init_soma.o \
             mpas_ocn_init_ziso.o \
             mpas_ocn_init_sub_ice_shelf_2D.o \
             mpas_ocn_init_periodic_planar.o \
             mpas_ocn_init_ecosys_column.o \
             mpas_ocn_init_sea_mount.o \
             mpas_ocn_init_global_ocean.o \
             mpas_ocn_init_isomip.o \
             mpas_ocn_init_isomip_plus.o
             #mpas_ocn_init_TEMPLATE.o

all: init_mode

init_mode: $(UTILS) $(TEST_CASES) $(OBJS)

mpas_ocn_init_mode.o: $(UTILS) $(TEST_CASES)

mpas_ocn_init_cell_markers.o:

mpas_ocn_init_interpolation.o:

mpas_ocn_init_ssh_and_landIcePressure.o: mpas_ocn_init_interpolation.o mpas_ocn_init_vertical_grids.o

mpas_ocn_init_spherical_utils.o:

mpas_ocn_init_vertical_grids.o:

mpas_ocn_init_seaSurfaceHeightAndPressure.o:

mpas_ocn_init_baroclinic_channel.o: $(UTILS)

mpas_ocn_init_iso.o: $(UTILS)

mpas_ocn_init_soma.o: $(UTILS)

mpas_ocn_init_lock_exchange.o: $(UTILS)

mpas_ocn_init_internal_waves.o: $(UTILS)

mpas_ocn_init_overflow.o: $(UTILS)

mpas_ocn_init_global_ocean.o: $(UTILS)

mpas_ocn_init_sub_ice_shelf_2D.o: $(UTILS)

mpas_ocn_init_cvmix_WSwSBF.o: $(UTILS)

mpas_ocn_init_periodic_planar.o: $(UTILS)

mpas_ocn_init_ecosys_column.o: $(UTILS)

mpas_ocn_init_sea_mount.o: $(UTILS)

mpas_ocn_init_isomip.o: $(UTILS)

mpas_ocn_init_isomip_plus.o: $(UTILS)

mpas_ocn_init_ziso.o: $(UTILS)

#mpas_ocn_init_TEMPLATE.o: $(UTILS)

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
