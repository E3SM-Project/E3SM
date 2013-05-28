.SUFFIXES: .F .o

OBJS = mpas_init_atm_mpas_core.o \
       mpas_init_atm_cases.o \
       mpas_atm_advection.o \
       mpas_init_atm_read_met.o \
       mpas_init_atm_llxy.o \
       mpas_init_atm_bitarray.o \
       mpas_init_atm_queue.o \
       mpas_init_atm_hinterp.o \
       mpas_init_atm_static.o \
       mpas_init_atm_surface.o \
       read_geogrid.o \
       mpas_atmphys_date_time.o \
       mpas_atmphys_initialize_real.o \
       mpas_atmphys_utilities.o

all: core_hyd

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS)

mpas_init_atm_cases.o: \
	read_geogrid.o \
	mpas_atm_advection.o \
	mpas_init_atm_read_met.o \
	mpas_init_atm_llxy.o \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_static.o \
	mpas_init_atm_surface.o \
	mpas_atmphys_initialize_real.o

mpas_init_atm_hinterp.o: mpas_init_atm_queue.o mpas_init_atm_bitarray.o

mpas_advection.o: 

mpas_init_atm_read_met.o:

read_geogrid.o:

mpas_init_atm_llxy.o:

mpas_init_atm_mpas_core.o: mpas_advection.o mpas_init_atm_cases.o

mpas_init_atm_static.o: \
	mpas_atm_advection.o \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_llxy.o \
	mpas_atmphys_utilities.o

mpas_init_atm_surface.o: \
	mpas_init_atm_hinterp.o  \
	mpas_init_atm_llxy.o     \
	mpas_init_atm_read_met.o

mpas_atmphys_initialize_real.o:  \
	mpas_init_atm_surface.o  \
	mpas_atmphys_date_time.o \
	mpas_atmphys_utilities.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators  -I../external/esmf_time_f90

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
