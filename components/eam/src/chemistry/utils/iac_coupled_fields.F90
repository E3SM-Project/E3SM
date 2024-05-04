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
  logical, protected :: iac_present=.false.

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

  !Number of days in a year can be a real number
  real(r8)          :: days_in_one_yr

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

    ! BFIXME: Hardwire it true for now - I don't know if this is an infodata or namelist
    ! variable or some other way of setting this for non-iac coupled runs.
    iac_present=.true.
    if (.not. iac_present) return

    ! Register the co2 field with the pbuf
    call pbuf_add_field(iac_co2_name,'physpkg',dtype_r8,(/pcols,pver/),idx)!FIXME:B what is this idx, do we need it, do we need to capture the error?

    if (masterproc)write(102,*)"Rgstr iac_co2 PBUF;iac_coupled_fields_register ts:", get_nstep()

  end subroutine iac_coupled_fields_register

  subroutine iac_coupled_fields_init(state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Initialize the iac_co2 data handling ****
    ! called by: phys_init (in physpkg.F90)
    !-------------------------------------------------------------------
    use physics_buffer,   only: pbuf_get_index
    use physics_buffer,   only: pbuf_add_field, dtype_r8, physics_buffer_desc
    use ppgrid,           only: begchunk, endchunk, pcols, pver
    use phys_grid,        only: get_ncols_p, phys_grid_initialized
    use physics_types,    only: physics_state

    implicit none

    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !Local vars
    integer :: ichunk        ! chunk index
    integer :: ierror   ! Error code
    character(len=cxx)  :: err_str

    real(r8) :: first_year, second_year

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return

    !Find the number of days in a year by subtracting 
    !number of days between year 1 and year 2.

    !number of days in the second year
    call set_time_float_from_date( second_year, & !out
         2, 1, 1, 0 ) !in (year, month, day, time)
   
    !number of days in the first year 
    call set_time_float_from_date( first_year,  & !out
         1, 1, 1, 0 ) !in (year, month, day, time)

    !number of days in a year
    days_in_one_yr = second_year-first_year


    !get pbuf index to store the field in pbuf
    iac_co2_pbuf_ndx = pbuf_get_index(iac_co2_name,ierror)

    if(ierror < 0 ) then
       write(err_str,*)'ERROR:failed to get pbuf index for species:',iac_co2_name,' errorcode is:',ierror,',',errmsg(__FILE__, __LINE__)
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

    if (masterproc)write(102,*)"Alloc iac_vertical_emiss-hh-ll; iac_coupled_fields_init-ts:", get_nstep()

  end subroutine iac_coupled_fields_init


  subroutine iac_coupled_fields_adv( state, pbuf2d)
    !-------------------------------------------------------------------
    ! Update pbuf2d with the iac coupled fields (just co2, for now)
    ! called by: phys_timestep_init at every time step (in physpkg.F90)
    !-------------------------------------------------------------------

    use perf_mod,     only: t_startf, t_stopf
    use tracer_data,  only: advance_trcdata
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use ppgrid,       only: pcols, pver
    use string_utils, only: to_lower, GLC
    use cam_history,  only: outfld
    use physconst,    only: mwdry       ! molecular weight dry air ~ kg/kmole
    use physconst,    only: boltz                ! J/K/molecule
    ! C.-C. Chen
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer  :: ind, ichunk, ncol, klev, icol, caseid, m, pbuf_ndx
    real(r8) :: to_mmr(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return
    call t_startf('iac_coupled_fields_adv')


    do ichunk = begchunk,endchunk
       ncol = state(ichunk)%ncol
       pbuf_chnk => pbuf_get_chunk(pbuf2d, ichunk)
       call pbuf_get_field(pbuf_chnk, iac_co2_pbuf_ndx, tmpptr )

       ! So tmpptr(col,h) is where we put our data.
       ! We have coupled "low" and "high" aircraft fields from iac.  The structure
       ! of the typical aircraft files used in aircraft_emit.F90 to prescribe co2
       ! generally have some positive values at ~11km and either values all 0.0
       ! below or some positive values for every level below.  To mimic this
       ! structure, we will put the "airhi" values at the 11km layer and evenly spread
       ! the "airlo" values across all layers from the surface up to 11km.  To get a
       ! finer vertical structure, we'll need to either use a standard file as a
       ! template to scale to total values, or calculate and couple more fields from
       ! iac.

       ! Loop over columns in this chunk
       do icol = 1, ncol
          ! init all vertical layers to 0
          tmpptr(icol,:) = 0._r8

          ! Find the layer containing 11 km
          do klev = 1, pver
             ! Check upper layer boundary
             if (state(ichunk)%zi(icol, klev+1) < 11000) cycle
             tmpptr(icol,klev) = iac_vertical_emiss(ichunk)%fco2_high_height(klev)

             ! Layers 1 to klev-1 now get an even distribution of the co2 in airlo.
             tmpptr(icol,1:klev-1) = iac_vertical_emiss(ichunk)%fco2_low_height(klev)/(klev-1)
             !if (ichunk==begchunk .and. klev==65 .and. icol ==1 .and. masterproc)write(102,*)"Called iac_coupled_fields_adv loop:",get_nstep(),":",iac_vertical_emiss(ichunk)%fco2_low_height(klev), iac_vertical_emiss(ichunk)%fco2_high_height(klev)
             write(104,*)"Called iac_coupled_fields_adv loop:",get_nstep(),":",iac_vertical_emiss(ichunk)%fco2_low_height(klev), iac_vertical_emiss(ichunk)%fco2_high_height(klev)

             ! Exit the vertical loop
             exit
          end do
       end do
    enddo

    call t_stopf('iac_coupled_fields_adv')
    if (masterproc)write(102,*)"Asgn PBUF iac_co2 v-intrp; iac_coupled_fields_adv ts:", get_nstep()
  end subroutine iac_coupled_fields_adv

  subroutine iac_coupled_timeinterp(mon_spec, day_spec, tod_spec, & !input
       lower_bound, upper_bound, time_interp_frac) !output

    !-------------------------------------------------------------------
    !  This subroutine takes the current month, day, and dayfrac and
    !  calculates the boundary month indices and the interpolation
    !  fraction between them.  This is called once at the top of
    !  atm_import, and we do the interpolation in that faction for each
    !  of co2sfc, co2airlo, and co2airhi, for each column.
    !-------------------------------------------------------------------

    !arguments
    integer,  intent (in)  :: mon_spec   ! simulation month
    integer,  intent (in)  :: day_spec   ! Simulation day
    integer,  intent (in)  :: tod_spec   ! Simulation time of day (secs)

    integer,  intent (out) :: lower_bound, upper_bound ! Bounds for interpolation in "mid_mon_num_days" array
    real(r8), intent (out) :: time_interp_frac         ! fraction between boundaries to interp

    ! Dumb
    integer, parameter   :: tot_mon_in_a_year = 12  ! Total # of months in a year

    ! Number of days from the start of the year at mid month of each month
    real(r8), parameter  :: mid_mon_num_days(tot_mon_in_a_year) = [15.0_r8,   &
         45.0_r8, 74.0_r8, 105.0_r8, 135.0_r8, 166.0_r8, 196.0_r8, 227.0_r8, &
         258.0_r8, 288.0_r8, 319.0_r8, 349.0_r8]

    integer  :: year, month ! Year and month of the current date
    integer  :: day, secs   ! Day and time of day (in seconds past 0Z) of the current date
    real(r8) :: curr_model_time                    ! Current model time in days (can be fractional)
    real(r8) :: day_at_lower_bnd, day_at_upper_bnd ! Number of days at bounds
    real(r8) :: model_time                         ! Model time in days used for interpolation

    if (.not. iac_present) return

    call get_curr_date(year, month, day, secs )


    !Get current model time in "days".
    !Output from the following call is "curr_mdl_time", which is in days, it can be fraction
    !based on how far we are in a given day
    call set_time_float_from_date( curr_model_time, & !output
         year , month, day, secs ) !input

    !Find the lower bound of an interval in mid_mon_num_days array based on curr_model_time
    call findplb(mid_mon_num_days, tot_mon_in_a_year, curr_model_time, & !input
         lower_bound) !output

    !Day of the year at the lower bound
    day_at_lower_bnd = mid_mon_num_days(lower_bound) ! m is good

    !The upper bound is next to the lower bound
    upper_bound = lower_bound + 1

    model_time = curr_model_time ! store the model time in days (it can be fractional)

    !if the lower bound is the last month (December), set lower bound to 1 (January)
    ! and recompute "model_time".
    if (lower_bound == tot_mon_in_a_year) then

       !We are in last half of December or first half Jan
       upper_bound = 1 ! set upper bound to the month of Jan

       !Let's assume we are in the Jan of the 2nd year
       !Upper bound number of  days is now 15 + 365 = 380
       day_at_upper_bnd = mid_mon_num_days(upper_bound) + days_in_one_yr

       !if we are in the first half of Jan, we need to add days (in one year), so that
       !it is Jan of 2nd year
       if (curr_model_time < day_at_lower_bnd) then
          model_time = curr_model_time + days_in_one_yr
       endif
    else
       !if we are in Dec 16th to 31st range
       day_at_upper_bnd = mid_mon_num_days(upper_bound)
    endif

    !fraction of time into the current model time
    time_interp_frac = (model_time-day_at_lower_bnd)/(day_at_upper_bnd-day_at_lower_bnd)


  end subroutine iac_coupled_timeinterp

end module iac_coupled_fields
