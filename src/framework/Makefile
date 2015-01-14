.SUFFIXES: .F .c .o

DEPS := $(shell find ../core_$(CORE)/ -type f -name "*.xml" ! -name "*processed.xml")

OBJS = mpas_kind_types.o \
       mpas_framework.o \
       mpas_timer.o \
       mpas_timekeeping.o \
       mpas_configure.o \
       mpas_constants.o \
       mpas_attlist.o \
       mpas_hash.o \
       mpas_sort.o \
       mpas_block_decomp.o \
       mpas_block_creator.o \
       mpas_dmpar.o \
       mpas_decomp.o \
       mpas_io.o \
       mpas_io_streams.o \
       mpas_bootstrapping.o \
       mpas_io_units.o \
       mpas_stream_manager.o \
       mpas_stream_list.o \
       mpas_forcing.o \
       mpas_c_interfacing.o \
       random_id.o \
       streams.o \
       timer.o \
       pool_hash.o \
       mpas_derived_types.o \
       mpas_domain_routines.o \
       mpas_field_routines.o \
       mpas_pool_routines.o \
       xml_stream_parser.o \

all: framework $(DEPS)

framework: $(OBJS)
	ar -ru libframework.a $(OBJS) ../external/ezxml/ezxml.o

mpas_configure.o: mpas_dmpar.o mpas_io_units.o mpas_pool_routines.o $(DEPS)

mpas_framework.o: mpas_dmpar.o \
                  mpas_io.o \
                  mpas_derived_types.o \
                  mpas_domain_routines.o \
                  mpas_field_routines.o \
                  mpas_pool_routines.o \
                  mpas_configure.o \
                  mpas_timer.o \
                  mpas_sort.o \
                  mpas_io_units.o \
                  mpas_stream_manager.o \
                  mpas_c_interfacing.o

mpas_constants.o: mpas_kind_types.o

mpas_attlist.o: mpas_kind_types.o mpas_io_units.o mpas_derived_types.o

mpas_derived_types.o: mpas_kind_types.o mpas_constants.o

mpas_domain_routines.o: mpas_derived_types.o mpas_pool_routines.o

mpas_field_routines.o: mpas_derived_types.o

mpas_pool_routines.o: mpas_derived_types.o mpas_field_routines.o

mpas_decomp.o: mpas_derived_types.o mpas_stream_manager.o

mpas_hash.o : mpas_derived_types.o

mpas_dmpar.o: mpas_sort.o streams.o mpas_kind_types.o mpas_derived_types.o mpas_hash.o mpas_io_units.o

mpas_sort.o: mpas_kind_types.o mpas_io_units.o

mpas_timekeeping.o: mpas_kind_types.o mpas_io_units.o mpas_derived_types.o

mpas_timer.o: mpas_kind_types.o mpas_io_units.o mpas_dmpar.o

mpas_block_decomp.o: mpas_derived_types.o mpas_hash.o mpas_configure.o mpas_io_units.o

mpas_block_creator.o: mpas_dmpar.o mpas_hash.o mpas_sort.o mpas_configure.o mpas_io_units.o mpas_block_decomp.o mpas_stream_manager.o mpas_decomp.o $(DEPS)

mpas_io.o: mpas_dmpar.o mpas_io_units.o mpas_attlist.o

mpas_io_streams.o: mpas_attlist.o mpas_derived_types.o mpas_timekeeping.o mpas_io.o mpas_io_units.o mpas_pool_routines.o $(DEPS)

mpas_bootstrapping.o: mpas_derived_types.o mpas_dmpar.o mpas_block_decomp.o mpas_block_creator.o mpas_sort.o mpas_configure.o mpas_timekeeping.o mpas_io_streams.o mpas_io_units.o mpas_stream_manager.o random_id.o $(DEPS)

mpas_io_units.o: mpas_kind_types.o

mpas_stream_list.o: mpas_derived_types.o mpas_kind_types.o mpas_io_units.o mpas_io_streams.o mpas_timekeeping.o

mpas_stream_manager.o: mpas_io_streams.o mpas_timekeeping.o mpas_derived_types.o mpas_io_units.o mpas_kind_types.o mpas_c_interfacing.o mpas_stream_list.o mpas_dmpar.o mpas_io.o

mpas_forcing.o: mpas_derived_types.o mpas_timekeeping.o mpas_io_streams.o mpas_stream_manager.o

xml_stream_parser.o: xml_stream_parser.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -I../external/ezxml -c xml_stream_parser.c

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
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -I../external/ezxml -c $<
