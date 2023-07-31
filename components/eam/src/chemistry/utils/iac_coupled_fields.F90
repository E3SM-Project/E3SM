module iac_coupled_fields
  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Converts the airborne CO2 coupled from iac (gcam) into a field in
  ! the pbuf, interpolated onto the right vertical and temporal grid.
  ! The coupler should take care of getting it on the horizontal grid.
  ! Follows the structure of aircraft_emit.F90, as it should replace
  ! the CO2 flux from the aircraft input file.  I'm attempting to make
  ! this compatible with aircraft_emit.F90, so you could use aircraft
  ! for other species and iac for CO2, so I want to keep things as
  ! similar as possible.
  !
  ! Authors: Tim Shippert, May 2023
  !-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8, cx =>SHR_KIND_CX, cl =>SHR_KIND_CL, &
       cs =>SHR_KIND_CS, cxx =>SHR_KIND_CXX
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use tracer_data,      only: trfld, trfile
  use cam_logfile,      only: iulog
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use input_data_utils, only: time_coordinate

  implicit none
  private
  save

  public :: iac_coupled_fields_register
  public :: iac_coupled_fields_init
  public :: iac_coupled_fields_adv
  public :: iac_coupled_timeinterp

  ! Public data types
  public iac_coupled_data_t             ! Data structure for airhi and airlo

  ! Public data
  public iac_present
  public iac_coupled_data
  logical :: iac_present

  !------------------------------------------------
  ! This is the chunked and columnated data from the iac coupling
  !------------------------------------------------
  type iac_coupled_data_t
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     real(r8), allocatable :: fco2_low_iac(:)   ! co2 flux from iac (low alt)
     real(r8), allocatable :: fco2_high_iac(:)  ! co2 flux from iac (high alt)
  end type iac_coupled_data_t

  logical  :: cam_outfld = .false.
  integer, parameter :: huge_int = huge(1)
  character(len=16) :: iac_CO2_name = 'iac_CO2'
  integer           :: iac_CO2_pbuf_ndx = -1
  real(r8) :: lev_bnds(2), time_bnds(1)

  ! This is our internal data structure to hold the chunked airlo and airhi
  type(iac_coupled_data_t), pointer:: iac_coupled_data(:)

contains

  subroutine iac_coupled_fields_register()

    !------------------------------------------------------------------
    ! **** Add the iac coupled fields the physics buffer ****
    ! called by: phys_register (in physpkg.F90)
    !------------------------------------------------------------------
    use physics_buffer, only: pbuf_add_field, dtype_r8
    !use tracer_data,    only: incr_filename
    !use constituents,   only: cnst_add
    use ppgrid,         only: pver, pcols

    integer            :: i,idx, mm, ind, n

    ! Hardwire it true for now - I don't know if this is an infodata or namelist
    ! variable or some other way of setting this for non-iac coupled runs.
    iac_present=.true.
    if (.not. iac_present) return

    ! Register the CO2 field with the pbuf
    call pbuf_add_field(iac_CO2_name,'physpkg',dtype_r8,(/pcols,pver/),idx)

    ! If we add as a tracer, make a call to cnst_add() here
  end subroutine iac_coupled_fields_register

  subroutine iac_coupled_fields_init(state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Initialize the iac_CO2 data handling ****
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
    integer :: c        ! chunk index
    integer :: ierror   ! Error code
    character(len=cxx)  :: err_str

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return

    !get pbuf index to store the field in pbuf
    iac_CO2_pbuf_ndx = pbuf_get_index(iac_CO2_name,ierror)

    if(ierror < 0 ) then
       write(err_str,*)'failed to get pbuf index for species:',iac_CO2_name,' errorcode is:',ierror,',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif

    ! Allocate and initialize the coupled data arrays
    allocate (iac_coupled_data(begchunk:endchunk), stat=ierror)
 
    if ( ierror /= 0 )then
       write(iulog,*) 'Allocation error: ', ierror
       call endrun('IAC_COUPLED_FIELDS_INIT error: allocation error')
    end if

    do c = begchunk,endchunk
       iac_coupled_data(c)%lchnk = c
       iac_coupled_data(c)%ncol  = get_ncols_p(c)

       allocate (iac_coupled_data(c)%fco2_low_iac(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('IAC_COUPLED_FIELDS_INIT error: allocation error fco2_low_iac')

       allocate (iac_coupled_data(c)%fco2_high_iac(pcols), stat=ierror)
       if ( ierror /= 0 ) call endrun('IAC_COUPLED_FIELDS_INIT error: allocation error fco2_high_iac')

       iac_coupled_data(c)%fco2_low_iac (:) = 0._r8
       iac_coupled_data(c)%fco2_high_iac (:) = 0._r8
    end do

  end subroutine iac_coupled_fields_init


  subroutine iac_coupled_fields_adv( state, pbuf2d)
    !-------------------------------------------------------------------
    ! Update pbuf2d with the iac coupled fields (just CO2, for now)
    ! called by: phys_timestep_init (in physpkg.F90)
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

    integer  :: ind, c, ncol, i, icol, caseid, m, pbuf_ndx
    real(r8) :: to_mmr(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)
    character(len = cs) :: units_spc

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return
    call t_startf('iac_coupled_fields_adv')

    ! Here is where the magic happens - we've time interpolated on import, and the
    ! surfexch fields are properly chunked and on the column grid.  So this is
    ! mostly about converrting airlo and airhi to a vertical distribution, and
    ! stuffing the values into the pbuf.

    ! 
    !$OMP PARALLEL DO PRIVATE (IC, NCOL, tmpptr, pbuf_chnk)
    do c = begchunk,endchunk
       ncol = state(c)%ncol
       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
       call pbuf_get_field(pbuf_chnk, iac_CO2_pbuf_ndx, tmpptr )

       ! So tmpptr(col,h) is where we put our data.
       ! We have coupled "low" and "high" aircraft fields from iac.  The structure
       ! of the typical aircraft files used in aircraft_emit.F90 to prescribe CO2
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
          do i=1,pver
             ! Check upper layer boundary
             if (state(c)%zi(icol, i+1) < 11000) cycle 

             ! i is now the layer with 11km, which gets the airhi value
             tmpptr(icol,i) = iac_coupled_data(c)%fco2_high_iac(i)

             ! Layers 1 to i-1 now get an even distribution of the co2 in airlo.  
             tmpptr(icol,1:i-1) = iac_coupled_data(c)%fco2_low_iac(i)/(i-1)

             ! Exit the vertical loop
             exit
          end do
       end do
    enddo

    call t_stopf('iac_coupled_fields_adv')
  end subroutine iac_coupled_fields_adv

  subroutine iac_coupled_timeinterp(mon_spec, day_spec, tod_spec, &
       b1, b2, tfrac)
    !-------------------------------------------------------------------
    !  This subroutine takes the current month, day, and dayfrac and
    !  calculates the boundary month indeces and the interpolation
    !  fraction between them.  This is called once at the top of
    !  atm_import, and we do the interpolation in that faction for each
    !  of co2sfc, co2airlo, and co2airhi, for each column.
    !-------------------------------------------------------------------

    !args
    integer, intent (in)            :: mon_spec   ! (simulation month)
    integer, intent (in)            :: day_spec   ! (simulation day)
    integer, intent (in)            :: tod_spec   ! (simulation time of day (secs)

    integer, intent (out)           :: b1  ! Lower boundary month index for interp
    integer, intent (out)           :: b2  ! Upper boundary month index for interp

    real(r8), intent (out)          :: tfrac ! fraction between boundaries to interp

    ! Dumb
    integer, dimension(12)          :: mon_length, mon_mid
    integer                         :: d0  ! number of days since b1 (zero-offset)

    real(r8)                        :: denom

    if (.not. iac_present) return

    ! Our coupled values represent the "middle" of the month, so we'll use midnight
    ! on the end of the mon_mid as the representative value
    mon_length= (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    mon_mid= (/15, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15/)

    ! Use cyclical seasons, which means we use the end of this December to
    ! interpolate the start of this January, and the other way around

    if (day_spec <= mon_mid(mon_spec)) then 
       ! b1 is prior month, b2 is this month
       if (mon_spec > 1) then 
          b1=mon_spec-1 
       else 
          b1 = 12
       endif
          
       b2=mon_spec

       ! d0 is the whole number of days since b1.  So, if we are in the first half
       ! of b2, we add number of days in last half of b1 to day_spec-1 (to
       ! zero_offset from b1) 
       d0 = day_spec - 1 + mon_length(b1)-mon_mid(b1) 
    else
       ! b1 is this month, b2 is next month
       b1=mon_spec
       if (mon_spec < 12) then 
          b2=mon_spec+1 
       else 
          b2=1
       endif

       ! For d0, we are in second half of b1, so subtract number of days in first
       ! half of b1 from day_spec-1, following similar logic as above
       d0 = day_spec - 1 - mon_mid(b1)
    endif

    ! Tfrac is the distance between b1 and b2 we want to go for our interp value 
    ! Thus, our denominator is the number days in the second half of b1 plus the
    ! number of days in the first half of b2
    denom = mon_length(b1)-mon_mid(b1) + mon_mid(b2)

    ! tfrac is then d0 plus dayfrac for today, divided by our denom
    tfrac = (d0 + tod_spec/86400.) / denom

    ! Voila
  end subroutine iac_coupled_timeinterp

end module iac_coupled_fields
