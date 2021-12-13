module glc_override_frac

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides functionality to allow overriding the ice fractions (i.e., the
  ! ice_covered field) and topographic heights that are sent to the coupler

  ! !USES:

  use glc_kinds_mod
  
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_glc_frac_overrides ! initialize stuff in this module, including reading the namelist
  public :: frac_overrides_enabled  ! return true if overrides are enabled, false otherwise
  public :: do_frac_overrides       ! do all overrides

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: read_namelist        ! read namelist setting options in this module
  private :: apply_increase_frac
  private :: apply_decrease_frac
  private :: apply_rearrange_freq
  private :: time_since_baseline  ! return time (days) since the baseline for adjustments
  
  
  ! !PRIVATE DATA MEMBERS:
  logical(log_kind) :: enable_frac_overrides  ! whether the overrides in this module are enabled for this run
  integer(int_kind) :: override_delay         ! time delay before beginning any overrides (days)
  real(r8)          :: decrease_frac          ! fractional decrease per day (should be positive)
  real(r8)          :: increase_frac          ! fractional increase per day
  integer(int_kind) :: rearrange_freq         ! frequency (days) at which we rearrange elevation classes

  ! Assumed maximum topographic height. It's okay for heights to go above this value, but
  ! they may not be handled exactly as desired by the overrides here.
  real(r8), parameter :: max_height = 3500._r8
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_glc_frac_overrides
    !
    ! !DESCRIPTION:
    ! Initialize stuff in this module, including reading the namelist
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'init_glc_frac_overrides'
    !-----------------------------------------------------------------------
    
    call read_namelist

  end subroutine init_glc_frac_overrides
  

  !-----------------------------------------------------------------------
  subroutine read_namelist
    !
    ! !DESCRIPTION:
    ! Read namelist setting options in this module
    !
    ! !USES:
    use glc_files      , only : nml_filename
    use glc_constants  , only : stdout, nml_in, blank_fmt, ndelim_fmt
    use glc_communicate, only : my_task, master_task
    use glc_broadcast  , only : broadcast_scalar
    use glc_exit_mod   , only : exit_glc, sigAbort
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer(int_kind) :: nml_error   ! namelist i/o error flag

    namelist /glc_override_nml/ enable_frac_overrides, override_delay, decrease_frac, &
         increase_frac, rearrange_freq
    
    character(len=*), parameter :: subname = 'read_namelist'
    !-----------------------------------------------------------------------
    
    ! Initialize namelist inputs
    enable_frac_overrides = .false.
    override_delay = 0
    decrease_frac = 0._r8
    increase_frac = 0._r8
    rearrange_freq = 0

    ! Read namelist
    if (my_task == master_task) then
       open(nml_in, file=nml_filename, status='old', iostat=nml_error)
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       endif
       do while (nml_error > 0)
          read(nml_in, nml=glc_override_nml,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    end if

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call exit_glc(sigAbort,'ERROR reading glc_override_nml')
    end if

    ! Write namelist settings
    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' GLC Override Frac:'
       write(stdout,blank_fmt)
       write(stdout,*) ' glc_override_nml namelist settings:'
       write(stdout,blank_fmt)
       write(stdout, glc_override_nml)
    end if

    ! Send namelist settings to all procs
    call broadcast_scalar(enable_frac_overrides, master_task)
    call broadcast_scalar(override_delay, master_task)
    call broadcast_scalar(decrease_frac, master_task)
    call broadcast_scalar(increase_frac, master_task)
    call broadcast_scalar(rearrange_freq, master_task)

  end subroutine read_namelist

  !-----------------------------------------------------------------------
  function frac_overrides_enabled() result(enabled)
    !
    ! !DESCRIPTION:
    ! Return true if glc fraction overrides are enabled in this run
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    logical :: enabled  ! function result
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'frac_overrides_enabled'
    !-----------------------------------------------------------------------
    
    enabled = enable_frac_overrides

  end function frac_overrides_enabled

  !-----------------------------------------------------------------------
  subroutine do_frac_overrides(ice_covered, topo, ice_sheet_grid_mask)
    !
    ! !DESCRIPTION:
    ! Do all overrides of glc fraction
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: ice_covered(:,:)
    real(r8), intent(inout) :: topo(:,:)
    real(r8), intent(in)    :: ice_sheet_grid_mask(:,:)
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'do_frac_overrides'
    !-----------------------------------------------------------------------

    call apply_increase_frac(ice_covered, topo, ice_sheet_grid_mask)
    call apply_decrease_frac(ice_covered, topo, ice_sheet_grid_mask)
    call apply_rearrange_freq(topo, ice_sheet_grid_mask)

  end subroutine do_frac_overrides

  !-----------------------------------------------------------------------
  subroutine apply_increase_frac(ice_covered, topo, ice_sheet_grid_mask)
    !
    ! !DESCRIPTION:
    ! Apply increase_frac to ice_covered
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: ice_covered(:,:)
    real(r8), intent(in)    :: topo(:,:)
    real(r8), intent(in)    :: ice_sheet_grid_mask(:,:)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: increase_topo_threshold

    character(len=*), parameter :: subname = 'apply_increase_frac'
    !-----------------------------------------------------------------------
    
    if (time_since_baseline() > 0) then
       ! When time_since_baseline * increase_frac is 0, we'll set elevations >= max_height
       !   to ice_covered = 1
       ! When time_since_baseline * increase_frac is 1, we'll set all elevations >= 0 to
       !   ice_covered = 1
       ! In between those times, we'll use an intermediate threshold
       increase_topo_threshold = (1._r8 - time_since_baseline() * increase_frac) * max_height
       increase_topo_threshold = max(increase_topo_threshold, 0._r8)
       increase_topo_threshold = min(increase_topo_threshold, max_height)

       where (ice_sheet_grid_mask > 0._r8 .and. topo >= increase_topo_threshold)
          ice_covered = 1._r8
       end where
    end if

  end subroutine apply_increase_frac

  !-----------------------------------------------------------------------
  subroutine apply_decrease_frac(ice_covered, topo, ice_sheet_grid_mask)
    !
    ! !DESCRIPTION:
    ! Apply decrease_frac to ice_covered
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: ice_covered(:,:)
    real(r8), intent(in)    :: topo(:,:)
    real(r8), intent(in)    :: ice_sheet_grid_mask(:,:)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: decrease_topo_threshold
    
    character(len=*), parameter :: subname = 'apply_decrease_frac'
    !-----------------------------------------------------------------------

    if (time_since_baseline() > 0) then
       ! When time_since_baseline * decrease_frac is 0, we'll set elevations < 0 to
       !   ice_covered = 0
       ! When time_since_baseline * decrease_frac is 1, we'll set all elevations < max_height
       !   to ice_covered = 0
       ! In between those times, we'll use an intermediate threshold
       decrease_topo_threshold = (time_since_baseline() * decrease_frac) * max_height
       decrease_topo_threshold = max(decrease_topo_threshold, 0._r8)
       decrease_topo_threshold = min(decrease_topo_threshold, max_height)

       where (ice_sheet_grid_mask > 0._r8 .and. topo < decrease_topo_threshold)
          ice_covered = 0._r8
       end where
    end if
    
  end subroutine apply_decrease_frac

  !-----------------------------------------------------------------------
  subroutine apply_rearrange_freq(topo, ice_sheet_grid_mask)
    !
    ! !DESCRIPTION:
    ! Apply rearrange_freq to topographic heights.
    !
    ! Example: if rearrange_freq = 3, then for the first 3 days, there will be no
    ! rearrangement, for the next 3 days (days 4-6) the topographic heights will be
    ! rearranged, for the next 3 days (days 7-9) the topographic heights will be back to
    ! normal, etc.
    !
    ! If rearrange_frac is <= 0, no rearrangement is done.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(inout) :: topo(:,:)
    real(r8), intent(in)    :: ice_sheet_grid_mask(:,:)
    !
    ! !LOCAL VARIABLES:
    integer(int_kind) :: num_intervals ! number of intervals of rearrange_freq

    character(len=*), parameter :: subname = 'apply_rearrange_freq'
    !-----------------------------------------------------------------------

    if (time_since_baseline() > 0 .and. rearrange_freq > 0) then
       ! num_intervals will be 0 in the first interval, 1 in the second, etc.
       num_intervals = time_since_baseline() / rearrange_freq

       ! Rearrange topographic heights half the time
       if (modulo(num_intervals, 2) == 1) then
          where(ice_sheet_grid_mask > 0._r8)
             topo = max(0._r8, max_height - topo)
          end where
       end if
    end if
    
  end subroutine apply_rearrange_freq


  !-----------------------------------------------------------------------
  integer(int_kind) function time_since_baseline()
    !
    ! !DESCRIPTION:
    ! Return time (days) since the baseline for adjustments (based on override_delay)
    !
    ! !USES:
    use glc_time_management, only : elapsed_days_init_date
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'time_since_baseline'
    !-----------------------------------------------------------------------
    
    time_since_baseline = elapsed_days_init_date - override_delay

  end function time_since_baseline




end module glc_override_frac
