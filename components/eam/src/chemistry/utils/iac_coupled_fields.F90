module iac_coupled_fields
  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Converts the airborne co2 coupled from iac (gcam) into a field in
  ! the pbuf, interpolated onto the right vertical and temporal grid.
  ! The coupler should take care of getting it on the horizontal grid.
  ! Follows the structure of aircraft_emit.F90, as it should replace
  ! the co2 flux from the aircraft input file.  I'm attempting to make
  ! this compatible with aircraft_emit.F90, so you could use aircraft
  ! for other species and iac for co2, so I want to keep things as
  ! similar as possible.
  !
  ! Authors: Tim Shippert, May 2023 and Balwinder Singh
  !-----------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use time_manager,     only: get_nstep
  use shr_kind_mod,     only: r8 => shr_kind_r8, cxx =>SHR_KIND_CXX
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use time_manager,     only: set_time_float_from_date, get_curr_date

  implicit none
  private
  save

  public :: iac_coupled_fields_register
  public :: iac_coupled_fields_init
  public :: iac_coupled_fields_adv
  public :: iac_coupled_timeinterp

  ! Public data types
  public iac_vertical_emiss_t             ! Data structure for airhi and airlo

  ! Public data
  public iac_vertical_emiss
  logical, public, protected :: iac_present=.false.

  !------------------------------------------------
  ! This is the chunked and columnated data from the iac coupling
  !------------------------------------------------
  type iac_vertical_emiss_t
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     real(r8), allocatable :: fco2_low_height(:)   ! co2 flux from iac (low alt)
     real(r8), allocatable :: fco2_high_height(:)  ! co2 flux from iac (high alt)
  end type iac_vertical_emiss_t

  character(len=16), parameter, public :: iac_co2_name = 'iac_co2'
  integer           :: iac_co2_pbuf_ndx = -1

  ! This is our internal data structure to hold the chunked airlo and airhi
  type(iac_vertical_emiss_t), pointer:: iac_vertical_emiss(:)

contains

  subroutine iac_coupled_fields_register()

    !------------------------------------------------------------------
    ! **** Add the iac coupled fields the physics buffer ****
    ! called by: phys_register (in physpkg.F90)
    !------------------------------------------------------------------
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use ppgrid,         only: pver, pcols

    integer            :: idx

    ! FIXMEB: Hardwire it true for now - I don't know if this is an infodata or namelist
    ! variable or some other way of setting this for non-iac coupled runs.
    iac_present=.false.
    if (.not. iac_present) return

    ! Register the co2 field with the pbuf
    call pbuf_add_field(iac_co2_name,'physpkg',dtype_r8,(/pcols,pver/),idx)
  end subroutine iac_coupled_fields_register

  subroutine iac_coupled_fields_init(state)
    !-------------------------------------------------------------------
    ! **** Initialize the iac_co2 data handling ****
    ! called by: phys_init (in physpkg.F90)
    !-------------------------------------------------------------------
    use physics_buffer,   only: pbuf_get_index
    use ppgrid,           only: begchunk, endchunk, pcols
    use phys_grid,        only: get_ncols_p
    use physics_types,    only: physics_state

    implicit none

    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)

    !Local vars
    integer :: ichunk        ! chunk index
    integer :: ierror   ! Error code
    character(len=cxx)  :: err_str

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return

    !get pbuf index to store the field in pbuf
    iac_co2_pbuf_ndx = pbuf_get_index(iac_co2_name,ierror)

    if(ierror < 0 ) then
       write(err_str,*)'ERROR:failed to get pbuf index for species:',trim(iac_co2_name),' errorcode is:',ierror,',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif

    ! Allocate and initialize the coupled data arrays
    allocate (iac_vertical_emiss(begchunk:endchunk), stat=ierror)

    if ( ierror /= 0 )then
      write(err_str,*) 'ERROR: Failed to allocate iac_vertical_emiss, allocation error: ', ierror,',',errmsg(__FILE__, __LINE__)
      call endrun(err_str)
    end if

    do ichunk = begchunk,endchunk
       iac_vertical_emiss(ichunk)%lchnk = ichunk
       iac_vertical_emiss(ichunk)%ncol  = get_ncols_p(ichunk)

       allocate (iac_vertical_emiss(ichunk)%fco2_low_height(pcols), stat=ierror)
       if ( ierror /= 0 ) then
         write(err_str,*) 'ERROR: Failed to allocate fco2_low_height, allocation error: ', ierror,',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
       endif

       allocate (iac_vertical_emiss(ichunk)%fco2_high_height(pcols), stat=ierror)
       if ( ierror /= 0 ) then
         write(err_str,*) 'ERROR: Failed to allocate fco2_high_height, allocation error: ', ierror,',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
       endif
       iac_vertical_emiss(ichunk)%fco2_low_height (:) = 0._r8
       iac_vertical_emiss(ichunk)%fco2_high_height (:) = 0._r8
    end do
  end subroutine iac_coupled_fields_init


  subroutine iac_coupled_fields_adv( state, pbuf2d)
    !-------------------------------------------------------------------
    ! Update pbuf2d with the iac coupled fields (just co2, for now)
    ! called by: phys_timestep_init at every time step (in physpkg.F90)
    !-------------------------------------------------------------------

    use perf_mod,     only: t_startf, t_stopf
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use ppgrid,       only: pcols, pver
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    character(len=cxx)  :: err_str
    integer, parameter :: HIGH_LAYER = 11000 !Height at which fco2_high should be released [m]
    integer, parameter :: iac_low_height_option = 1 !FIXMEB: This should be a namelist option

    integer  :: ichunk, ncol, klev, icol, high_lev
    real(r8), pointer :: tmpptr(:,:) ! pointer for storing CO2 values [kg/m2/s]
    real(r8) :: denominator, dist_amount

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return
    call t_startf('iac_coupled_fields_adv')


    do ichunk = begchunk,endchunk
       ncol = state(ichunk)%ncol !number of columns in this chunk
       pbuf_chnk => pbuf_get_chunk(pbuf2d, ichunk)
       call pbuf_get_field(pbuf_chnk, iac_co2_pbuf_ndx, tmpptr )! get tmpptr point to the iac field

       ! In the following do-loop, we will populate tmpptr with the CO2 values read in from
       ! iac_vertical_emiss. iac_vertical_emiss has CO2 at two heights, the high height and the low
       ! height. The value at high height are assigned to the layer at HIGH_LAYER. There are several
       ! options for distributing the low height values. The select-case block implements all those
       ! options

       ! Loop over columns in this chunk
       do icol = 1, ncol
          ! init all vertical layers to 0
          tmpptr(icol,:) = 0._r8

          ! Find the layer containing HIGH_LAYER value
          do klev = 1, pver
            if (state(ichunk)%zi(icol, klev+1) > HIGH_LAYER) cycle
            high_lev = klev !found the level at HIGH_LAYER
            exit
          enddo
          !assign the fco2_high_height to the "high_lev" level
          tmpptr(icol,high_lev) = iac_vertical_emiss(ichunk)%fco2_high_height(icol)

          !For the fco2_low_height, we have several options in the following select case statement
          select case (iac_low_height_option)
             case (1)
              !Distribut the fco2_low_height equally among all the levels below "high_lev"

              !Divide by zero check
              denominator = (pver - high_lev+1)
              if (denominator < 1.e-11 .and. denominator > -1.e-11) then !Checking if denominator is zero
                write(err_str,*) 'ERROR: possible divide by zero, denominator is: ',denominator ,',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
              end if
              !Amount of CO2 to distribute
              dist_amount = iac_vertical_emiss(ichunk)%fco2_low_height(icol)/(pver - high_lev+1)
              do klev = high_lev+1, pver
                tmpptr(icol,klev) = dist_amount
              enddo

             !Default case
             case default
                write(err_str,*) 'ERROR: Invalid fco2_low_height distribution option : ',iac_low_height_option,',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
          end select
       end do
    enddo

    call t_stopf('iac_coupled_fields_adv')
  end subroutine iac_coupled_fields_adv

  subroutine iac_coupled_timeinterp(lower_bound, upper_bound, time_interp_frac) !output
    use cam_logfile,          only: iulog
    !-------------------------------------------------------------------
    ! The IAC variables are passed at monthly boundaries, we linearly
    ! interpolate these variables in time. This subroutine computes the
    ! time interpolation fraction for the IAC component variables.
    !-------------------------------------------------------------------

    !--Arguments
    integer,  intent (out) :: lower_bound, upper_bound ! Bounds for interpolation in "mid_mon_num_days" array below
    real(r8), intent (out) :: time_interp_frac         ! fraction between boundaries to interp

    !--Local variables
    character(len=cxx)  :: err_str

    ! Parameters
    integer, parameter   :: tot_mon_in_a_year = 12  ! Total # of months in a year

    ! Number of days from the start of the year at mid month of each month
    real(r8), parameter  :: mid_mon_num_days(tot_mon_in_a_year) = [15.0_r8,   &
         45.0_r8, 74.0_r8, 105.0_r8, 135.0_r8, 166.0_r8, 196.0_r8, 227.0_r8, &
         258.0_r8, 288.0_r8, 319.0_r8, 349.0_r8]

    integer  :: year, month ! Year and month of the current simulation time
    integer  :: day, secs   ! Day and time of day (in seconds past 0Z) of the current date
    real(r8) :: num_days_in_current_year !Number of days in the current year
    real(r8) :: curr_model_time ! Current model time in days (can be fractional)
    real(r8) :: cyclic_curr_model_time ! ! Current model time in days starting from the current year (can be fractional)
    real(r8) :: day_at_lower_bnd, day_at_upper_bnd ! Number of days into the current year at bounds
    real(r8) :: model_time                         ! Model time in days used for interpolation

    !If IAC component is not active, return
    if (.not. iac_present) return

    !Get the year, month, day and secs of this time step
    !(all arguments of the call below are outputs)
    call get_curr_date(year, month, day, secs ) !outputs


    !Get current model time in the number of "days" counted from the start of the simulation
    !or RUN_START_DATE. Output from the following call is "curr_mdl_time", which is in days,
    !it can be fraction based on how far we are in a given day
    call set_time_float_from_date( curr_model_time, & !output
         year, month, day, secs ) !input

    !Get the number of days in the current year
    num_days_in_current_year = get_num_days_in_year(year)

    !Convert curr_model_time to cyclic_curr_model_time. cyclic_curr_model_time
    !is the number of days into the simulation starting
    !from the current year. This is done so that we can reuse "mid_mon_num_days"
    !for every year to find the bounds for the interpolation
    !(it might be an overkill to compute it every time step but we do it to protect
    !against calendars with leap years)
    cyclic_curr_model_time = mod(curr_model_time, num_days_in_current_year)

    !Find the lower bound for the interpolation
    !Lower bound of an interval in mid_mon_num_days array based on cyclic_curr_model_time
    call findplb(mid_mon_num_days, tot_mon_in_a_year, cyclic_curr_model_time, & !input
         lower_bound) !output

    !Day of the year at the lower bound from the "mid_mon_num_days" array
    day_at_lower_bnd = mid_mon_num_days(lower_bound)

    !The upper bound is next to the lower bound
    upper_bound = lower_bound + 1

    model_time = cyclic_curr_model_time ! store the model time in days (it can be fractional)

    !We have a special case for dates between Dec 16th to Jan 15th. In this case, the lower
    !bound should be 12 (i.e., around december 15th) and the upper bound should be 1 (i.e., around Jan 15th).
    !We identify this special case by "lower_bound == tot_mon_in_a_year" if condition (i.e., the lower bound is
    !the last month - [December, i.e., 12] of the year)
    !We set upper bound to 1 (January) and recompute "model_time" for this special case.
    if (lower_bound == tot_mon_in_a_year) then

       !We are in last half of December (16th to 31st) or first half Jan (1st to 15th)
       upper_bound = 1 ! set upper bound to the month of Jan

       !The upper bound is now the January of the "next" year, so the number of days should be
       !more than "num_days_in_current_year". Upper bound number of  days is now 15 + 365 = 380
       day_at_upper_bnd = mid_mon_num_days(upper_bound) + num_days_in_current_year

       ! If we are between Jan 1st and Jan 15th of the "next year", we need to recompute
       ! model_time so that it is in the "next" year. Therefore, we add "num_days_in_current_year"
       ! to the cyclic_curr_model_time
       !NOTE: For december 16th to 31st, the model_time computed above should be used
       if (cyclic_curr_model_time < day_at_lower_bnd) then
          !FIXME: we should use "cyclic_curr_model_time <= mid_mon_num_days(1)" or explicit
          !month and day instead of "cyclic_curr_model_time < day_at_lower_bnd" to be clear
          model_time = cyclic_curr_model_time + num_days_in_current_year
       endif
    else
       !For all other cases, we can simply find days at upper bounds from mid_mon_num_days
       day_at_upper_bnd = mid_mon_num_days(upper_bound)
    endif

    !fraction of time into the current model time
    time_interp_frac = (model_time-day_at_lower_bnd)/(day_at_upper_bnd-day_at_lower_bnd)

    !Sanity check 1: tfrac should be between 0 and 1
    if(time_interp_frac>1._r8 .or. time_interp_frac<0._r8) then
       write(err_str,*)'ERROR:time_interp_frac should be between 0.0 and 1.0, time_interp_frac is:',time_interp_frac,',',errmsg(__FILE__, __LINE__)
       write(iulog,*)'ERROR: time_interp_frac out of bounds, it should be between 0 and 1'
       write(iulog,*)'model_time:',model_time
       write(iulog,*)'day_at_lower_bnd:',day_at_lower_bnd
       write(iulog,*)'day_at_upper_bnd:',day_at_upper_bnd
       write(iulog,*)'upper_bound:',upper_bound
       write(iulog,*)'curr_model_time:',curr_model_time
       write(iulog,*)'lower_bound:',lower_bound
       write(iulog,*)'year:', year
       write(iulog,*)'month:',month
       write(iulog,*)'day:',day
       write(iulog,*)'secs:', secs
       call endrun(trim(err_str))
    endif

    !Sanity check 2: At monthly boundaries (i.e. start of the next mid month), the time_interp_frac should be zero
    ! (Adds 1 to mid_mon_num_days(month) in the following if condition to capture start of the next mid month interval)
    if (day == mid_mon_num_days(month)+1 .and. secs == 0) then
       if (time_interp_frac > tiny(time_interp_frac) ) then
          write(err_str,*)'ERROR:time_interp_frac should be 0.0 at mid month days, time_interp_frac:',&
               time_interp_frac,', day:',day,', month:',month,', mid_mon_num_days:',mid_mon_num_days(month) &
               ,errmsg(__FILE__, __LINE__)
          call endrun(trim(err_str))
       endif
    endif
  end subroutine iac_coupled_timeinterp

  function get_num_days_in_year(year_in) result(days_in_curr_yr)
    !--------------------------------------------------------------------------------
    !Compute the total number of days in current year
    !--------------------------------------------------------------------------------
    !*Logic*: We compute the total number of days at 01/01 (MM/DD) for the next year
    !and this year. The difference between the two is the total number of days
    !in the current year
    !--------------------------------------------------------------------------------
    !input
    integer, intent(in) :: year_in

    !output (return value)
    real(r8) :: days_in_curr_yr

    !local
    real(r8) :: num_days_next_year, num_days_this_year

    !number of days on Jan 1st the next year
    call set_time_float_from_date( num_days_next_year, & !out
         year_in+1, 1, 1, 0 ) !in (year, month, day, time)

    !number of days on Jan 1st this year
    call set_time_float_from_date( num_days_this_year,  & !out
         year_in, 1, 1, 0 ) !in (year, month, day, time)

    !number of days in the current year (return value)
    days_in_curr_yr = num_days_next_year-num_days_this_year

  end function get_num_days_in_year

end module iac_coupled_fields
