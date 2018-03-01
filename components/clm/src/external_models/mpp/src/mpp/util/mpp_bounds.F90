module mpp_bounds

  !
  ! !PUBLIC TYPES:
  implicit none
  !
  save

  integer, public :: bounds_proc_begg      = 0
  integer, public :: bounds_proc_endg      = 0
  integer, public :: bounds_proc_begc      = 0
  integer, public :: bounds_proc_endc      = 0
  integer, public :: bounds_proc_begg_all  = 0
  integer, public :: bounds_proc_endg_all  = 0
  integer, public :: bounds_proc_begc_all  = 0
  integer, public :: bounds_proc_endc_all  = 0
  integer, public :: nclumps               = 0

  ! !PUBLIC MEMBER FUNCTIONS:
  public mpp_bounds_init_proc_bounds
  public mpp_bounds_init_clump
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine mpp_bounds_init_proc_bounds(begg_val, endg_val, begg_all_val, endg_all_val, &
      begc_val, endc_val, begc_all_val, endc_all_val)

    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent (in) :: begg_val, endg_val
    integer, intent (in) :: begg_all_val, endg_all_val
    integer, intent (in) :: begc_val, endc_val
    integer, intent (in) :: begc_all_val, endc_all_val

    bounds_proc_begg = begg_val
    bounds_proc_endg = endg_val
    bounds_proc_begc = begc_val
    bounds_proc_endc = endc_val

    bounds_proc_begg_all = begg_all_val
    bounds_proc_endg_all = endg_all_val
    bounds_proc_begc_all = begc_all_val
    bounds_proc_endc_all = endc_all_val

  end subroutine mpp_bounds_init_proc_bounds

  !------------------------------------------------------------------------------
  subroutine mpp_bounds_init_clump(nclumps_val)

    !
    ! !DESCRIPTION:
    ! Initialize module variables 
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent (in) :: nclumps_val

    nclumps = nclumps_val

  end subroutine mpp_bounds_init_clump

end module mpp_bounds
