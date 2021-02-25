module shr_nuopc_scalars_mod

  !----------------------------------------------------------------------------
  ! scalars
  !----------------------------------------------------------------------------

  implicit none
  public

  character(len=*) , parameter :: flds_scalar_name = "cpl_scalars"

  integer, parameter :: flds_scalar_index_nx                 = 1
  integer, parameter :: flds_scalar_index_ny                 = 2
  integer, parameter :: flds_scalar_index_precip_fact        = 3
  integer, parameter :: flds_scalar_index_nextsw_cday        = 4
  integer, parameter :: flds_scalar_index_dead_comps         = 5
  integer, parameter :: flds_scalar_index_rofice_present     = 6  ! does rof have iceberg coupling on
  integer, parameter :: flds_scalar_index_flood_present      = 7  ! does rof have flooding on
  integer, parameter :: flds_scalar_index_iceberg_prognostic = 10 ! does ice model support icebergs
  integer, parameter :: flds_scalar_index_glclnd_present     = 11 ! does glc have land coupling fields on
  integer, parameter :: flds_scalar_index_glcocn_present     = 12 ! does glc have ocean runoff on
  integer, parameter :: flds_scalar_index_glcice_present     = 13 ! does glc have iceberg coupling on
  integer, parameter :: flds_scalar_index_glc_valid_input    = 14 ! does glc have is valid accumulated data being sent to it?
                                                                  ! (only valid if glc_prognostic is .true.)
  integer, parameter :: flds_scalar_index_glc_coupled        = 15 ! does glc send fluxes to other components
                                                                  ! (only relevant if glc_present is .true.)
  integer, parameter :: flds_scalar_num                      = 15

end module shr_nuopc_scalars_mod
