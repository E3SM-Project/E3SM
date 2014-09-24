.SUFFIXES: .F .o

DEPS := $(shell find ../core_$(CORE)/ -type f -name "*.xml" ! -name "*processed.xml")

OBJS = mpas_kind_types.o \
       mpas_framework.o \
       mpas_timer.o \
       mpas_timekeeping.o \
       mpas_configure.o \
       mpas_packages.o \
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
       mpas_io_units.o \
       mpas_stream_manager.o \
       mpas_stream_list.o \
       mpas_c_interfacing.o \
       random_id.o \
       streams.o \
       pool_hash.o \
       xml_stream_parser.o \
       ../registry/ezxml/ezxml.o

all: framework $(DEPS)

framework: $(OBJS)
	ar -ru libframework.a $(OBJS)

mpas_configure.o: mpas_dmpar.o mpas_io_units.o $(DEPS)

mpas_packages.o: $(DEPS)

mpas_framework.o: mpas_dmpar.o \
                  mpas_io_input.o \
                  mpas_io_output.o \
                  mpas_io.o \
                  mpas_grid_types.o \
                  mpas_configure.o \
                  mpas_timer.o \
                  mpas_sort.o \
                  mpas_io_units.o \
                  mpas_packages.o \
                  mpas_stream_manager.o \
                  mpas_c_interfacing.o

mpas_constants.o: mpas_kind_types.o mpas_io_units.o

mpas_dmpar_types.o : mpas_kind_types.o mpas_io_units.o

mpas_attlist.o: mpas_kind_types.o mpas_io_units.o

mpas_grid_types.o: mpas_kind_types.o mpas_dmpar_types.o mpas_attlist.o mpas_io_units.o mpas_packages.o mpas_io_units.o pool_hash.o mpas_timekeeping.o pool_subroutines.inc duplicate_field_array.inc duplicate_field_scalar.inc $(DEPS)

mpas_dmpar.o: mpas_sort.o streams.o mpas_kind_types.o mpas_grid_types.o mpas_hash.o mpas_io_units.o

mpas_sort.o: mpas_kind_types.o mpas_io_units.o

mpas_timekeeping.o: mpas_kind_types.o mpas_io_units.o

mpas_timer.o: mpas_kind_types.o mpas_io_units.o

mpas_block_decomp.o: mpas_grid_types.o mpas_hash.o mpas_configure.o mpas_io_units.o

mpas_block_creator.o: mpas_dmpar.o mpas_hash.o mpas_sort.o mpas_configure.o mpas_io_units.o $(DEPS)

mpas_io.o: mpas_dmpar_types.o mpas_io_units.o

mpas_io_streams.o: mpas_attlist.o mpas_grid_types.o mpas_timekeeping.o mpas_io.o mpas_io_units.o $(DEPS)

mpas_io_input.o: mpas_grid_types.o mpas_dmpar.o mpas_block_decomp.o mpas_block_creator.o mpas_sort.o mpas_configure.o mpas_timekeeping.o mpas_io_streams.o mpas_io_units.o mpas_stream_manager.o random_id.o $(DEPS)

mpas_io_output.o: mpas_grid_types.o mpas_dmpar.o mpas_sort.o mpas_configure.o mpas_io_streams.o mpas_io_units.o $(DEPS)

mpas_io_units.o:

mpas_stream_list.o: mpas_grid_types.o mpas_kind_types.o mpas_io_units.o mpas_io_streams.o mpas_timekeeping.o

mpas_stream_manager.o: mpas_io_streams.o mpas_timekeeping.o mpas_grid_types.o mpas_io_units.o mpas_kind_types.o mpas_c_interfacing.o mpas_stream_list.o mpas_dmpar.o

xml_stream_parser.o: xml_stream_parser.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -I../registry -c xml_stream_parser.c

clean:
	$(RM) *.o *.mod *.f90 libframework.a
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) -I../inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../inc -I../external/esmf_time_f90
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
