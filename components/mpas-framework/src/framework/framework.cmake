# framework
list(APPEND CPPDEFS "-Dcoupled")
list(APPEND COMMON_RAW_SOURCES
  framework/mpas_kind_types.F
  framework/mpas_framework.F
  framework/mpas_timer.F
  framework/mpas_timekeeping.F
  framework/mpas_constants.F
  framework/mpas_attlist.F
  framework/mpas_hash.F
  framework/mpas_sort.F
  framework/mpas_block_decomp.F
  framework/mpas_block_creator.F
  framework/mpas_dmpar.F
  framework/mpas_abort.F
  framework/mpas_decomp.F
  framework/mpas_threading.F
  framework/mpas_io.F
  framework/mpas_io_streams.F
  framework/mpas_bootstrapping.F
  framework/mpas_io_units.F
  framework/mpas_stream_manager.F
  framework/mpas_stream_list.F
  framework/mpas_forcing.F
  framework/mpas_c_interfacing.F
  framework/random_id.c
  framework/pool_hash.c
  framework/mpas_derived_types.F
  framework/mpas_domain_routines.F
  framework/mpas_field_routines.F
  framework/mpas_pool_routines.F
  framework/xml_stream_parser.c
  framework/regex_matching.c
  framework/mpas_field_accessor.F
  framework/mpas_log.F
  framework/mpas_global_sum_mod.F
)
