module seaice_core_interface_structs

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units
   use seaice_core
   use seaice_analysis_driver
   use mpas_log, only: mpas_log_write
   
   implicit none
   public
   
   contains

#include "../inc/structs_and_variables.inc"

#include "../inc/namelist_call.inc"

#include "../inc/namelist_defines.inc"

end module seaice_core_interface_structs
