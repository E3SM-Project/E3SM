module FATESFireDataMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for FATES to obtain fire inputs from data
  !
  ! !USES:
  use shr_kind_mod  , only: r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod   , only: errmsg => shr_log_errMsg
  use abortutils    , only: endrun
  use elm_varctl    , only: iulog
  use decompMod     , only: bounds_type
  use FatesFireBase , only: fates_fire_base_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fates_fire_data_type
  !
  type, extends(fates_fire_base_type) :: fates_fire_data_type
      ! !PRIVATE MEMBER DATA:
      real(r8), private, pointer :: lnfm24(:) ! Daily avg lightning by grid cell (#/km2/hr)

      contains
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: need_lightning_and_popdens
      procedure, public :: GetLight24     ! Return 24-hour averaged lightning data
      procedure, public :: GetGDP         ! Return the global gdp data
      procedure, public :: InitAccBuffer  ! Initialize accumulation processes
      procedure, public :: InitAccVars  ! Initialize accumulation variables
      procedure, public :: UpdateAccVars  ! Update/extract accumulations vars

  end type fates_fire_data_type

  character(len=*), parameter, private :: sourcefile = __FILE__

contains

  !------------------------------------------------------------------------
  function need_lightning_and_popdens(this)
    ! !ARGUMENTS:
    class(fates_fire_data_type), intent(in) :: this
    logical :: need_lightning_and_popdens  ! function result
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'need_lightning_and_popdens'
    !-----------------------------------------------------------------------

    need_lightning_and_popdens = .true.
  end function need_lightning_and_popdens

  !-----------------------------------------------------------------------
  function GetLight24( this ) result(lnfm24)
    !
    ! !DESCRIPTION: Get the 24-hour averaged lightning data
    ! !USES
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    real(r8), pointer :: lnfm24(:)
    !---------------------------------------------------------------------
    lnfm24 => this%lnfm24
    !---------------------------------------------------------------------
  end function
  
  !-----------------------------------------------------------------------
  function GetGDP( this ) result(gdp)
    !
    ! !DESCRIPTION: Get the global gross domestic product data
    ! !USES
    
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    real(r8), pointer :: gdp(:)
    
    !---------------------------------------------------------------------
    gdp => this%gdp_lf_col
    !---------------------------------------------------------------------
  end function GetGDP
  
  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds )
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES
    use elm_varcon, only : spval
    use accumulMod, only : init_accum_field
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: ier
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    allocate(this%lnfm24(begg:endg), stat=ier)
    if (ier/=0) then
       call endrun(msg="allocation error for lnfm24"//&
            errmsg(sourcefile, __LINE__))
    endif
    this%lnfm24(:) = spval
    call init_accum_field (name='lnfm24', units='strikes/km2/hr', &
         desc='24hr average of lightning strikes',  accum_type='runmean', &
         accum_period=-1, subgrid_type='gridcell', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart
    ! file is read in and the accumulation buffer is obtained)
    !
    ! !USES
    use accumulMod       , only : extract_accum_field
    use elm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslg(:)  ! temporary
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslg(begg:endg), stat=ier)
    if (ier/=0) then
       call endrun(msg="allocation error for rbufslg"//&
            errmsg(sourcefile, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('lnfm24', rbufslg, nstep)
    this%lnfm24(begg:endg) = rbufslg(begg:endg)

    deallocate(rbufslg)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use elm_time_manager, only: get_nstep
    use accumulMod      , only: update_accum_field, extract_accum_field
    use abortutils      , only: endrun
    !
    ! !ARGUMENTS:
    class(fates_fire_data_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: dtime                 ! timestep size [seconds]
    integer :: nstep                 ! timestep number
    integer :: ier                   ! error status
    integer :: begg, endg
    real(r8), pointer :: rbufslg(:)  ! temporary single level - gridcell level
    !---------------------------------------------------------------------

    begg = bounds%begg; endg = bounds%endg

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level gridcell field

    allocate(rbufslg(begg:endg), stat=ier)
    if (ier/=0) then
       write(iulog,*)'UpdateAccVars allocation error for rbuf1dg'
       call endrun(msg=errmsg(sourcefile, __LINE__))
    endif

    ! Accumulate and extract lnfm24
    rbufslg(begg:endg) = this%forc_lnfm(begg:endg)
    call update_accum_field  ('lnfm24', rbufslg, nstep)
    call extract_accum_field ('lnfm24', this%lnfm24, nstep)

    deallocate(rbufslg)

  end subroutine UpdateAccVars

end module FATESFireDataMod
