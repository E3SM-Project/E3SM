module AnnualFluxDribbler

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! Defines a class for handling fluxes that are generated once per year (e.g., due to
  ! transient landcover changes that happen at the year boundary), but are meant to be
  ! dribbled in evenly throughout the year.
  !
  ! This assumes that the once-per-year fluxes are generated on the first timestep of the
  ! year. Any flux given on the first timestep of the year is dribbled evenly for every
  ! timestep of the coming year. Any flux given on other timesteps is applied entirely in
  ! the current timestep. (Note that, if there is a combination of an annual flux and an
  ! every-time-step flux, with both combined in the same delta term, then, on the first
  ! timestep of the year, the every-time-step flux generated on that timestep will be
  ! dribbled over the year rather than applied in that timestep.)
  !
  ! NOTE(wjs, 2016-08-30) If we change the glc coupling time to be more frequent, then
  ! we'll need to make this more dynamic: e.g., for coupling every 73 days (5 times per
  ! year), we'd need to dribble fluxes over the next 73 days.
  !
  ! Typical usage:
  !
  !   - call mydribbler%set_curr_delta every time step
  !
  !     This must be called every timestep, even if the delta is currently zero, in order
  !     to zero out any existing stored delta. This can (and generally should) even be
  !     called when it isn't the first timestep of the year. For deltas that are non-zero
  !     at times other than the first timestep of the year, they will simply be passed on
  !     to the output flux in get_curr_flux, making for easier handling by the client.
  !
  !   - call mydribbler%get_curr_flux every time step, AFTER set_curr_delta
  !
  !     This will get the current flux for this timestep, which is the sum of (1) the
  !     dribbled flux from the last start-of-year timestep, and (2) the current timestep's
  !     flux, based on the delta passed in to set_curr_delta in this timestep, if this is
  !     not the start-of-year timestep.
  !
  !     Alternatively, you can call mydribbler%get_dribbled_delta, if you need the result as
  !     a delta over the time step rather than as a per-second flux.
  !
  ! And, for the sake of checking conservation:
  !
  !   - To get gridcell water (or whatever) content at the start of the time step:
  !
  !     call mydribbler%get_amount_left_to_dribble_beg
  !
  !   - To get gridcell water (or whatever) content at the end of the time step:
  !
  !     call mydribbler%get_amount_left_to_dribble_end
  !
  !   These both return the pseudo-state representing how much of the original delta
  !   still needs to be dribbled. The 'beg' version includes the amount left to dribble
  !   in the current time step; the 'end' version does not.
  !  
  !
  ! !USES:
  use elm_varctl       , only : iulog
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_beg, get_end
  use decompMod        , only : BOUNDS_SUBGRID_GRIDCELL, BOUNDS_SUBGRID_PATCH
  use elm_varcon       , only : secspday, nameg, namep
  use clm_time_manager , only : get_days_per_year, get_step_size_real, is_beg_curr_year
  use clm_time_manager , only : get_curr_yearfrac, get_prev_yearfrac, get_prev_date
  use clm_time_manager , only : is_first_step
  !
  implicit none
  private

  ! Compiler support for allocatable characters isn't fully robust (particularly for
  ! pgi), so using a max lengths for now
  !
  ! (If we used allocatable characters, these max lengths could be removed
  integer, parameter :: name_maxlen = 128
  integer, parameter :: units_maxlen = 64
  integer, parameter :: subgrid_maxlen = 64

  ! !PUBLIC TYPES:

  type, public :: annual_flux_dribbler_type
     private
     ! Metadata
     character(len=name_maxlen) :: name
     character(len=units_maxlen) :: units

     ! Whether this dribbler allows non-zero deltas on time steps other than the first
     ! time step of the year
     logical :: allows_non_annual_delta

     ! Which subgrid level this dribbler is operating at, stored in various ways
     character(len=subgrid_maxlen) :: dim1name
     character(len=subgrid_maxlen) :: name_subgrid
     integer :: bounds_subgrid_level

     ! Annual amount to dribble in over the year
     real(r8), pointer :: amount_to_dribble(:)

     ! Amount from the current timestep to pass through to the flux, if this isn't the
     ! first timestep of the year
     real(r8), pointer :: amount_from_this_timestep(:)
   contains
     ! Public infrastructure methods
     procedure, public :: Restart
     procedure, public :: Clean

     ! Public science methods
     procedure, public :: set_curr_delta  ! Set the delta state for this time step
     procedure, public :: get_curr_flux   ! Get the current flux for this time step
     procedure, public :: get_dribbled_delta  ! Similar to get_curr_flux, but gets result as a delta rather than a per-second flux
     procedure, public :: get_amount_left_to_dribble_beg  ! Get the pseudo-state representing the amount that still needs to be dribbled in this and future time steps
     procedure, public :: get_amount_left_to_dribble_end  ! Get the pseudo-state representing the amount that still needs to be dribbled in just future time steps

     ! Private methods
     procedure, private :: allocate_and_initialize_data
     procedure, private :: set_metadata
     procedure, private :: get_amount_left_to_dribble
  end type annual_flux_dribbler_type

  public :: annual_flux_dribbler_gridcell  ! Creates an annual_flux_dribbler_type object at the gridcell-level
  public :: annual_flux_dribbler_patch     ! Creates an annual_flux_dribbler_type object at the patch-level

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Factory methods
  !
  ! For now, there are only factory methods for gridcell-level and patch-level. But
  ! adding the ability to work at other levels is as easy as adding another factory
  ! method like this (along with some variables in the 'only' clauses of the 'use'
  ! statements).
  ! ========================================================================

  !-----------------------------------------------------------------------
  function annual_flux_dribbler_gridcell(bounds, name, units, allows_non_annual_delta) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates an annual_flux_dribbler_type object at the gridcell-level
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(annual_flux_dribbler_type) :: this   ! function result
    type(bounds_type), intent(in)   :: bounds
    character(len=*) , intent(in)   :: name   ! name of this object, used for i/o
    character(len=*) , intent(in)   :: units  ! units metadata - should be state units, not flux (i.e., NOT per-second)

    ! If allows_non_annual_delta is .false., then an error check is performed for each
    ! call to set_curr_delta, ensuring that the delta is 0 at all times other than the
    ! first time step of the year. This is just provided as a convenient sanity check -
    ! to ensure that the code is behaving as expected. (However, non-zero deltas are
    ! always allowed on the first step of the run.)
    !
    ! If allows_non_annual_delta is not provided, it is assumed to be .true.
    logical, intent(in), optional :: allows_non_annual_delta
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'annual_flux_dribbler_gridcell'
    !-----------------------------------------------------------------------

    this%dim1name = 'gridcell'
    this%name_subgrid = nameg
    this%bounds_subgrid_level = BOUNDS_SUBGRID_GRIDCELL

    call this%allocate_and_initialize_data(bounds)
    call this%set_metadata(name, units, allows_non_annual_delta)

  end function annual_flux_dribbler_gridcell

  !-----------------------------------------------------------------------
  function annual_flux_dribbler_patch(bounds, name, units, allows_non_annual_delta) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Creates an annual_flux_dribbler_type object at the patch-level
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(annual_flux_dribbler_type) :: this   ! function result
    type(bounds_type), intent(in)   :: bounds
    character(len=*) , intent(in)   :: name   ! name of this object, used for i/o
    character(len=*) , intent(in)   :: units  ! units metadata - should be state units, not flux (i.e., NOT per-second)

    ! If allows_non_annual_delta is .false., then an error check is performed for each
    ! call to set_curr_delta, ensuring that the delta is 0 at all times other than the
    ! first time step of the year. This is just provided as a convenient sanity check -
    ! to ensure that the code is behaving as expected.
    !
    ! If allows_non_annual_delta is not provided, it is assumed to be .true.
    logical, intent(in), optional :: allows_non_annual_delta
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'annual_flux_dribbler_patch'
    !-----------------------------------------------------------------------

    this%dim1name = 'pft'
    this%name_subgrid = namep
    this%bounds_subgrid_level = BOUNDS_SUBGRID_PATCH

    call this%allocate_and_initialize_data(bounds)
    call this%set_metadata(name, units, allows_non_annual_delta)

  end function annual_flux_dribbler_patch

  ! ========================================================================
  ! Public methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine set_curr_delta(this, bounds, delta)
    !
    ! !DESCRIPTION:
    ! Sets the delta state for this time step. Note that the delta is specified just as
    ! the change in state - NOT as a flux (per-second) quantity.
    !
    ! This must be called every timestep, even if the deltas are currently 0, in order to
    ! zero out any existing stored delta. This can (and generally should) even be called
    ! when it isn't the first timestep of the year. For deltas that are non-zero at times
    ! other than the first timestep of the year, they will simply be passed on to the
    ! output flux in get_curr_flux, making for easier handling by the client. (i.e., this
    ! class handles the addition of the dribbled flux and the current flux for you.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: delta( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    integer :: beg_index, end_index
    integer :: i
    integer :: yr, mon, day, tod

    character(len=*), parameter :: subname = 'set_curr_delta'
    !-----------------------------------------------------------------------

    beg_index = lbound(delta, 1)
    end_index = get_end(bounds, this%bounds_subgrid_level)
    SHR_ASSERT_ALL((ubound(delta) == (/end_index/)), errMsg(sourcefile, __LINE__))

    if (is_beg_curr_year()) then
       do i = beg_index, end_index
          this%amount_to_dribble(i) = delta(i)

          ! On the first timestep of the year, we don't have any pass-through flux. Need
          ! to zero out any previously-set amount_from_this_timestep.
          this%amount_from_this_timestep(i) = 0._r8
       end do
    else
       do i = beg_index, end_index
          this%amount_from_this_timestep(i) = delta(i)
       end do
       if (.not. this%allows_non_annual_delta .and. .not. is_first_step()) then
          do i = beg_index, end_index
             if (this%amount_from_this_timestep(i) /= 0._r8) then
                write(iulog,*) subname//' ERROR: found unexpected non-zero delta mid-year'
                write(iulog,*) 'Dribbler name: ', trim(this%name)
                write(iulog,*) 'i, delta = ', i, this%amount_from_this_timestep(i)
                call get_prev_date(yr, mon, day, tod)
                write(iulog,*) 'Start of time step date (yr, mon, day, tod) = ', &
                     yr, mon, day, tod
                write(iulog,*) 'This indicates that some non-zero flux was generated at a time step'
                write(iulog,*) 'other than the first time step of the year, which this dribbler was told not to expect.'
                write(iulog,*) 'If this non-zero mid-year delta is expected, then you can suppress this error'
                write(iulog,*) 'by setting allows_non_annual_delta to .true. when constructing this dribbler.'
                call endrun(decomp_index=i, elmlevel=this%name_subgrid, &
                     msg=subname//': found unexpected non-zero delta mid-year: ' // &
                     errMsg(sourcefile, __LINE__))
             end if
          end do
       end if
    end if

  end subroutine set_curr_delta

  !-----------------------------------------------------------------------
  subroutine get_curr_flux(this, bounds, flux)
    !
    ! !DESCRIPTION:
    ! Gets the current flux for this timestep, and stores it in the flux argument.
    !
    ! This should be called AFTER set_curr_delta is called for the given timestep.
    !
    ! This will get the current flux for this timestep, which is the sum of (1) the
    ! dribbled flux from the last start-of-year timestep, and (2) the current timestep's
    ! flux, based on the delta passed in to set_curr_delta in this timestep, if this is
    ! not the start-of-year timestep.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(out) :: flux( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    integer :: beg_index, end_index
    integer :: i
    real(r8) :: secs_per_year
    real(r8) :: dtime
    real(r8) :: flux_from_dribbling
    real(r8) :: flux_from_this_timestep

    character(len=*), parameter :: subname = 'get_curr_flux'
    !-----------------------------------------------------------------------

    beg_index = lbound(flux, 1)
    end_index = get_end(bounds, this%bounds_subgrid_level)
    SHR_ASSERT_ALL((ubound(flux) == (/end_index/)), errMsg(sourcefile, __LINE__))

    secs_per_year = get_days_per_year() * secspday
    dtime = get_step_size_real()

    do i = beg_index, end_index
       flux_from_dribbling = this%amount_to_dribble(i) / secs_per_year
       flux_from_this_timestep = this%amount_from_this_timestep(i) / dtime
       flux(i) = flux_from_dribbling + flux_from_this_timestep
    end do

  end subroutine get_curr_flux

  !-----------------------------------------------------------------------
  subroutine get_dribbled_delta(this, bounds, delta)
    !
    ! !DESCRIPTION:
    ! Gets the current delta for this timestep, and stores it in the delta argument.
    !
    ! This is similar to get_curr_flux, but returns the total, dribbled delta over this
    ! timestep, rather than a per-second flux. See documentation in get_curr_flux for
    ! more usage details.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(out) :: delta( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    integer :: beg_index, end_index
    integer :: i
    real(r8) :: dtime
    real(r8), allocatable :: flux(:)

    character(len=*), parameter :: subname = 'get_dribbled_delta'
    !-----------------------------------------------------------------------

    beg_index = lbound(delta, 1)
    end_index = get_end(bounds, this%bounds_subgrid_level)
    SHR_ASSERT_ALL((ubound(delta) == (/end_index/)), errMsg(sourcefile, __LINE__))

    allocate(flux(beg_index:end_index))

    call this%get_curr_flux(bounds, flux(beg_index:end_index))

    dtime = get_step_size_real()
    do i = beg_index, end_index
       delta(i) = flux(i) * dtime
    end do

  end subroutine get_dribbled_delta


  !-----------------------------------------------------------------------
  subroutine get_amount_left_to_dribble_beg(this, bounds, amount_left_to_dribble)
    !
    ! !DESCRIPTION:
    ! Get the pseudo-state representing the amount that still needs to be dribbled in
    ! this and future time steps. This represents the pseudo-state before this time
    ! step's dribbling flux has been removed. (This behavior is regardless of whether
    ! get_curr_flux has been called already this time step.)
    !
    ! As a special case, this returns 0 in the first time step of the year, because we
    ! haven't created this year's dribbling pool as of the beginning of this time step.
    !
    ! i.e., if we imagined that the total amount to dribble was added to a state
    ! variable, and then this state variable was updated each time step as the flux
    ! dribbles out, then this subroutine gives the amount left in that state. (However,
    ! the actual implementation doesn't explicitly track this state, which is why we
    ! refer to it as a pseudo-state.)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(out) :: amount_left_to_dribble( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: yearfrac

    character(len=*), parameter :: subname = 'get_amount_left_to_dribble_beg'
    !-----------------------------------------------------------------------

    yearfrac = get_prev_yearfrac()
    call this%get_amount_left_to_dribble(bounds, yearfrac, amount_left_to_dribble)

  end subroutine get_amount_left_to_dribble_beg


  !-----------------------------------------------------------------------
  subroutine get_amount_left_to_dribble_end(this, bounds, amount_left_to_dribble)
    !
    ! !DESCRIPTION:
    ! Gets the pseudo-state representing the amount that still needs to be dribbled in
    ! future time steps. This represents the pseudo-state after this time step's dribbling
    ! flux has been removed. i.e., this includes the amount that will be dribbled starting
    ! with the *next* time step, through the end of this year. So this will return 0 on
    ! the last time step of the year. (This behavior is regardless of whether
    ! get_curr_flux has been called already this time step.)
    !
    ! See documentation of get_amount_left_to_dribble_beg for more details.
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(out) :: amount_left_to_dribble( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: yearfrac

    character(len=*), parameter :: subname = 'get_amount_left_to_dribble_end'
    !-----------------------------------------------------------------------

    yearfrac = get_curr_yearfrac()
    call this%get_amount_left_to_dribble(bounds, yearfrac, amount_left_to_dribble)

  end subroutine get_amount_left_to_dribble_end


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: restname  ! name of field on restart file
    logical :: readvar

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------

    restname = trim(this%name) // '_amt_to_dribble'
    call restartvar(ncid=ncid, flag=flag, varname=restname, xtype=ncd_double, &
         dim1name = this%dim1name, &
         long_name = 'total amount to dribble over the year for ' // trim(this%name), &
         units = trim(this%units), &
         interpinic_flag = 'interp', &
         readvar = readvar, &
         data = this%amount_to_dribble)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory associated with this object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------

    deallocate(this%amount_to_dribble)
    deallocate(this%amount_from_this_timestep)

  end subroutine Clean

  ! ========================================================================
  ! Private methods
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine allocate_and_initialize_data(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate arrays in this object and set them to initial values
    !
    ! Assumes this%bounds_subgrid_level is already set
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer :: beg_index, end_index

    character(len=*), parameter :: subname = 'allocate_and_initialize_data'
    !-----------------------------------------------------------------------

    beg_index = get_beg(bounds, this%bounds_subgrid_level)
    end_index = get_end(bounds, this%bounds_subgrid_level)

    allocate(this%amount_to_dribble(beg_index:end_index))
    this%amount_to_dribble(beg_index:end_index) = 0._r8

    allocate(this%amount_from_this_timestep(beg_index:end_index))
    this%amount_from_this_timestep(beg_index:end_index) = 0._r8

  end subroutine allocate_and_initialize_data

  !-----------------------------------------------------------------------
  subroutine set_metadata(this, name, units, allows_non_annual_delta)
    !
    ! !DESCRIPTION:
    ! Set metadata in this object
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(inout) :: this
    character(len=*) , intent(in)   :: name   ! name of this object, used for i/o
    character(len=*) , intent(in)   :: units  ! units metadata - should be state units, not flux (i.e., NOT per-second)

    ! If allows_non_annual_delta is .false., then an error check is performed for each
    ! call to set_curr_delta, ensuring that the delta is 0 at all times other than the
    ! first time step of the year. This is just provided as a convenient sanity check -
    ! to ensure that the code is behaving as expected.
    !
    ! If allows_non_annual_delta is not provided, it is assumed to be .true.
    logical, intent(in), optional :: allows_non_annual_delta
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_metadata'
    !-----------------------------------------------------------------------

    if (len_trim(name) > name_maxlen) then
       write(iulog,*) subname // ': name too long'
       write(iulog,*) trim(name) // ' exceeds max length: ', name_maxlen
       call endrun(msg=subname // ': name too long: ' // &
            errMsg(sourcefile, __LINE__))
    end if
    this%name = trim(name)

    if (len_trim(units) > units_maxlen) then
       write(iulog,*) subname // ': units too long'
       write(iulog,*) trim(units) // ' exceeds max length: ', units_maxlen
       call endrun(msg=subname // ': units too long: ' // &
            errMsg(sourcefile, __LINE__))
    end if
    this%units = trim(units)

    if (present(allows_non_annual_delta)) then
       this%allows_non_annual_delta = allows_non_annual_delta
    else
       this%allows_non_annual_delta = .true.
    end if

  end subroutine set_metadata

  !-----------------------------------------------------------------------
  subroutine get_amount_left_to_dribble(this, bounds, yearfrac, amount_left_to_dribble)
    !
    ! !DESCRIPTION:
    ! Helper method shared by get_amount_left_to_dribble_beg and
    ! get_amount_left_to_dribble_end. Returns amount left to dribble as of a given
    ! yearfrac.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(annual_flux_dribbler_type), intent(in) :: this
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in)  :: yearfrac
    real(r8), intent(out) :: amount_left_to_dribble( get_beg(bounds, this%bounds_subgrid_level) : )
    !
    ! !LOCAL VARIABLES:
    integer :: beg_index, end_index
    integer :: i

    character(len=*), parameter :: subname = 'get_amount_left_to_dribble'
    !-----------------------------------------------------------------------

    beg_index = lbound(amount_left_to_dribble, 1)
    end_index = get_end(bounds, this%bounds_subgrid_level)
    SHR_ASSERT_ALL((ubound(amount_left_to_dribble) == (/end_index/)), errMsg(sourcefile, __LINE__))

    do i = beg_index, end_index
       if (yearfrac < 1.e-15_r8) then
          ! last time step of year; we'd like this to be given a yearfrac of 1 rather than
          ! 0 in this case; since it's given as 0, we need to handle it specially
          amount_left_to_dribble(i) = 0._r8
       else
          amount_left_to_dribble(i) = this%amount_to_dribble(i) * (1._r8 - yearfrac)
       end if
    end do

  end subroutine get_amount_left_to_dribble


end module AnnualFluxDribbler
