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
  use shr_kind_mod,     only: r8 => shr_kind_r8, cxx =>SHR_KIND_CXX
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use phys_control,     only: iac_present

  implicit none
  private
  save

  public :: iac_coupled_fields_register
  public :: iac_coupled_fields_init
  public :: iac_coupled_fields_adv

  ! Public data types
  public iac_vertical_emiss_t             ! Data structure for airhi and airlo

  ! Public data
  public iac_vertical_emiss

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


  subroutine iac_coupled_fields_adv( state, pbuf)
    !-------------------------------------------------------------------
    ! Update pbuf with the iac coupled fields (just co2, for now)
    ! called by: tphysac at every time step (in physpkg.F90)
    !-------------------------------------------------------------------

    use perf_mod,     only: t_startf, t_stopf
    use physics_types,only: physics_state
    use ppgrid,       only: pver
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    type(physics_state), intent(in)    :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    character(len=cxx)  :: err_str
    integer, parameter :: HIGH_LAYER = 11000 !Height at which fco2_high should be released [m]
    integer, parameter :: iac_low_height_option = 1 !FIXMEB: This should be a namelist option

    integer  :: ncol, klev, icol, high_lev
    real(r8), pointer :: tmpptr(:,:) ! pointer for storing CO2 values [kg/m2/s]
    real(r8) :: denominator, dist_amount

    !------------------------------------------------------------------
    ! Return if no iac coupling
    !------------------------------------------------------------------
    if (.not. iac_present) return
    call t_startf('iac_coupled_fields_adv')

    ncol = state%ncol
    call pbuf_get_field(pbuf, iac_co2_pbuf_ndx, tmpptr )! get tmpptr point to the iac field

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
        if (state%zi(icol, klev+1) > HIGH_LAYER) cycle
            high_lev = klev !found the level at HIGH_LAYER
            exit
          enddo
          !assign the fco2_high_height to the "high_lev" level
      tmpptr(icol,high_lev) = iac_vertical_emiss(state%lchnk)%fco2_high_height(icol)

          !For the fco2_low_height, we have several options in the following select case statement
          select case (iac_low_height_option)
             case (1)
              !Distribut the fco2_low_height equally among all the levels below "high_lev"

              !Divide by zero check
          denominator = (pver - high_lev)
          if (denominator < 1.e-11 .and. denominator > -1.e-11) then !Checking if denominator is near zero or a very tinay number
                write(err_str,*) 'ERROR: possible divide by zero, denominator is: ',denominator ,',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
              end if

              !Amount of CO2 to distribute
          dist_amount = iac_vertical_emiss(state%lchnk)%fco2_low_height(icol)/(pver - high_lev)
              do klev = high_lev+1, pver
                tmpptr(icol,klev) = dist_amount
              enddo

             !Default case
             case default
                write(err_str,*) 'ERROR: Invalid fco2_low_height distribution option : ',iac_low_height_option,',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
          end select
    end do ! End of loop over columns

    call t_stopf('iac_coupled_fields_adv')

  end subroutine iac_coupled_fields_adv

end module iac_coupled_fields
