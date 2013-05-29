.SUFFIXES: .F .o

OBJS = mpas_kind_types.o \
       mpas_framework.o \
       mpas_timer.o \
       mpas_timekeeping.o \
       mpas_configure.o \
       mpas_constants.o \
       mpas_dmpar_types.o \
       mpas_attlist.o \
       mpas_grid_types.o \
       mpas_hash.o \
       mpas_sort.o \
       mpas_block_decomp.o \
	   mpas_block_creator.o \
       mpas_dmpar.o \
       mpas_io.o \
       mpas_io_streams.o \
       mpas_io_input.o \
       mpas_io_output.o \
       streams.o

all: framework

framework: $(OBJS)
	ar -ru libframework.a $(OBJS)

mpas_framework.o: mpas_dmpar.o mpas_io_input.o mpas_io_output.o mpas_io.o mpas_grid_types.o mpas_configure.o mpas_timer.o

mpas_configure.o: mpas_dmpar.o

mpas_constants.o: mpas_kind_types.o

mpas_dmpar_types.o : mpas_kind_types.o

mpas_attlist.o: mpas_kind_types.o

mpas_grid_types.o: mpas_kind_types.o mpas_dmpar_types.o mpas_attlist.o

mpas_dmpar.o: mpas_sort.o streams.o mpas_kind_types.o mpas_grid_types.o mpas_hash.o

mpas_sort.o: mpas_kind_types.o

mpas_timekeeping.o: mpas_kind_types.o

mpas_timer.o: mpas_kind_types.o

mpas_block_decomp.o: mpas_grid_types.o mpas_hash.o mpas_configure.o

mpas_block_creator.o: mpas_dmpar.o mpas_hash.o mpas_sort.o mpas_configure.o

mpas_io.o: mpas_dmpar_types.o

mpas_io_streams.o: mpas_attlist.o mpas_grid_types.o mpas_timekeeping.o mpas_io.o

mpas_io_input.o: mpas_grid_types.o mpas_dmpar.o mpas_block_decomp.o mpas_block_creator.o mpas_sort.o mpas_configure.o mpas_timekeeping.o mpas_io_streams.o

mpas_io_output.o: mpas_grid_types.o mpas_dmpar.o mpas_sort.o mpas_configure.o mpas_io_streams.o

clean:
	$(RM) *.o *.mod *.f90 libframework.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../external/esmf_time_f90

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
