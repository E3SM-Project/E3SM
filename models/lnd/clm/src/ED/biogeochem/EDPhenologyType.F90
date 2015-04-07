module EDPhenologyType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module holds routines dealing with phenology in ED.  The primary use
  ! is to hold extract and accumulate routines
  
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_cal_mod       , only : calParams
  use shr_const_mod     , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use accumulMod        , only : update_accum_field, extract_accum_field, accumResetVal
  use clm_varctl        , only : iulog
  use clm_time_manager  , only : get_nstep, get_step_size
  !
  ! !USES:
  implicit none
  private
  !
  type, public :: ed_phenology_type
     !
     ! change these to allocatable
     ! add a rbuf variable that is a part of this type
     !
     real(r8), pointer :: ED_GDD_patch               (:)   ! ED Phenology growing degree days.
     ! This (phen_cd_status_patch?) could and should be site-level. RF
     integer , pointer :: phen_cd_status_patch       (:)   ! ED Phenology cold deciduous status
     character(10)     :: accString = 'ED_GDD0'
     real(r8)          :: checkRefVal = 26._r8

   contains

     ! Public procedures
     procedure, public  :: accumulateAndExtract
     procedure, public  :: init
     procedure, public  :: initAccVars
     procedure, public  :: initAccBuffer
     procedure, public  :: clean

     ! Private procedures
     procedure, private :: initAllocate
     procedure, private :: initHistory

  end type ed_phenology_type
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine accumulateAndExtract( this, bounds,     &
       t_ref2m_patch,    &
       gridcell, latdeg, &
       day, month, secs )
    !
    ! start formal argument list --
    ! group formal (dummy) arguments by use/similarity
    !
    class(ed_phenology_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds                       ! beginning and ending pft index
    ! data arguments
    real(r8)                , intent(in)     :: t_ref2m_patch(bounds%begp: ) ! patch 2 m height surface air temperature (K)
    ! arguments for the grid
    integer                 , intent(in)     :: gridcell(bounds%begp: )  ! gridcell
    real(r8)                , intent(in)     :: latdeg(bounds%begg: )    ! latitude (degrees)
    ! time related arguments
    integer                 , intent(in)     :: day           ! day
    integer                 , intent(in)     :: month         ! month
    integer                 , intent(in)     :: secs          ! secs
    !
    ! -- end formal argument list
    !

    !
    ! local variables
    !
    ! update_accum_field expects a pointer, can't make this an allocatable
    real(r8), pointer :: rbufslp(:)           ! temporary single level - pft level
    integer           :: g, p                 ! local index for gridcell and pft
    integer           :: ier                  ! error code
    integer           :: m                    ! local month variable

    allocate(rbufslp(bounds%begp:bounds%endp), stat=ier)
    if (ier/=0) then
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract GDD0 for ED
    do p = bounds%begp,bounds%endp

       g = gridcell(p)

       if (latdeg(g) >= 0._r8) then
          m = calParams%january
       else
          m = calParams%june
       endif

       ! FIX(RF,032414) - is this accumulation a bug in the normal phenology code,
       ! as it means to count from november but ctually counts from january?
       if ( month==m .and. day==calParams%firstDayOfMonth .and. secs==get_step_size() ) then
          rbufslp(p) = accumResetVal ! reset ED_GDD
       else
          rbufslp(p) = max(0._r8, min(this%checkRefVal, t_ref2m_patch(p)-SHR_CONST_TKFRZ)) &
               * get_step_size()/SHR_CONST_CDAY
       end if

       if( this%phen_cd_status_patch(p) == 2 ) then ! we have over-counted past the maximum possible  range
          rbufslp(p) = accumResetVal !don't understand how this doens't make it negative, but it doesn't. RF
       endif

       if( latdeg(g) >= 0._r8 .and. month >= calParams%july ) then !do not accumulate in latter half of year.
          rbufslp(p) = accumResetVal
       endif

       if( latdeg(g) < 0._r8 .and. month < calParams%june ) then !do not accumulate in earlier half of year.
          rbufslp(p) = accumResetVal
       endif

    end do

    call update_accum_field  ( trim(this%accString), rbufslp, get_nstep() )
    call extract_accum_field ( trim(this%accstring), this%ED_GDD_patch, get_nstep() )

    deallocate(rbufslp)

  end subroutine accumulateAndExtract

  !---------------------------------------------------------------------
  subroutine clean( this )
    !
    ! !DESCRIPTION:
    ! clean up memory
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(ed_phenology_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    deallocate(this%ED_GDD_patch)
    deallocate(this%phen_cd_status_patch)

  end subroutine clean

  subroutine init(this, bounds)

    class(ed_phenology_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds

    call this%initAllocate ( bounds )
    call this%initHistory ()

  end subroutine init

  !------------------------------------------------------------------------
  subroutine initAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(ed_phenology_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    !------------------------------------------------------------------------

    allocate(this%ED_GDD_patch         (bounds%begp:bounds%endp)) ; this%ED_GDD_patch         (:) = 0.0_r8
    allocate(this%phen_cd_status_patch (bounds%begp:bounds%endp)) ; this%phen_cd_status_patch (:) = 0

  end subroutine initAllocate

  !------------------------------------------------------------------------
  subroutine initHistory(this)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(Ed_phenology_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    call hist_addfld1d (fname=trim(this%accString), units='deg C',  &
         avgflag='A', long_name='ED phenology growing degree days', &
         ptr_patch=this%ED_GDD_patch, set_lake=0._r8, set_urb=0._r8)

  end subroutine initHistory

  !-----------------------------------------------------------------------
  subroutine initAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation.
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field].
    ! To accumulate a field, it must first be defined in subroutine [initAccVars]
    ! and then accumulated by calls to [updateAccVars].
    ! Four types of accumulations are possible:
    !   o average over time interval
    !   o running mean over time interval
    !   o running accumulation over time interval
    !  Time average fields are only valid at the end of the averaging interval.
    ! Running means are valid once the length of the simulation exceeds the
    ! averaging interval. Accumulated fields are continuously accumulated.
    ! The trigger value "-99999." resets the accumulation to zero.
    !
    ! !USES
    use accumulMod       , only : init_accum_field
    !
    !  !ARGUMENTS:
    class(ed_phenology_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds

    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    call init_accum_field (name=this%accString, units='K', &
         desc='growing degree-days base 0C from planting', accum_type='runaccum', accum_period=huge(1), &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine initAccBuffer

  !-----------------------------------------------------------------------
  subroutine initAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES
    !
    ! !ARGUMENTS:
    class(ed_phenology_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------

    allocate(rbufslp(bounds%begp:bounds%endp), stat=ier)
    if (ier/=0) then
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    call extract_accum_field (this%accString, rbufslp, get_nstep())
    this%ED_GDD_patch(bounds%begp:bounds%endp) = rbufslp(bounds%begp:bounds%endp)

    deallocate(rbufslp)

  end subroutine initAccVars

end module EDPhenologyType
