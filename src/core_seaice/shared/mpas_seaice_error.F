!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_error

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_stream_manager
  use mpas_log

  implicit none

  private
  save

  ! column critical error codes
  integer, parameter, public :: &
       SEAICE_ERROR_COL_VERT_THERM = 10, &
       SEAICE_ERROR_COL_ITD_THERM  = 11, &
       SEAICE_ERROR_COL_RIDGING    = 12, &
       SEAICE_ERROR_COL_BGC        = 13

  ! incremental remap critical error codes
  integer, parameter, public :: &
       SEAICE_ERROR_IR_NEG_AREA    = 20, &
       SEAICE_ERROR_IR_NEG_MASS    = 21, &
       SEAICE_ERROR_IR_MONO        = 22, &
       SEAICE_ERROR_IR_CONS        = 23

  ! public routines
  public :: &
       seaice_critical_error_write_block, &
       seaice_check_critical_error

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_critical_error_write_block
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_critical_error_write_block(domain, block)

    type(domain_type), intent(in) :: &
         domain

    type(block_type), intent(in) :: &
         block

    ! write out block streams
    call mpas_stream_mgr_block_write(&
         domain % streamManager, &
         writeBlock=block, &
         streamID='abort_block', &
         forceWriteNow=.true.)

  end subroutine seaice_critical_error_write_block

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_check_critical_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_check_critical_error(domain, ierr)

    use mpas_dmpar, only: mpas_dmpar_max_int

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    integer, intent(in) :: &
         ierr

    logical, pointer :: &
         config_full_abort_write

    integer :: &
         ierrmax

    call MPAS_pool_get_config(domain % configs, "config_full_abort_write", config_full_abort_write)

    if (config_full_abort_write) then

       ! find if anyone failed
       call mpas_dmpar_max_int(domain % dminfo, ierr, ierrmax)
       if (ierrmax > 0) then

          ! make final full write of abort stream
          call MPAS_stream_mgr_write(domain % streamManager, streamID="abort", forceWriteNow=.true.)

       endif
    endif

    ! check if a critical error occured
    if (ierr > 0) then

       ! kill the model
       call mpas_log_write("Runtime error $i: "//trim(error_string(ierr)), &
            messageType=MPAS_LOG_CRIT, intArgs=(/ierr/))

    endif

  end subroutine seaice_check_critical_error

  !-----------------------------------------------------------------------

  function error_string(ierr) result(errorStr)

    character(len=strKIND) :: errorStr

    integer, intent(in) :: ierr

    select case(ierr)
    case(SEAICE_ERROR_COL_VERT_THERM)
       errorStr = "Column: Vertical Thermodynamics"
    case(SEAICE_ERROR_COL_ITD_THERM)
       errorStr = "Column: ITD thermodynamics"
    case(SEAICE_ERROR_COL_RIDGING)
       errorStr = "Column: Ridging"
    case(SEAICE_ERROR_COL_BGC)
       errorStr = "Column: BGC"
    case(SEAICE_ERROR_IR_NEG_AREA)
       errorStr = "IR: Negative area"
    case(SEAICE_ERROR_IR_NEG_MASS)
       errorStr = "IR: Negative mass"
    case(SEAICE_ERROR_IR_MONO)
       errorStr = "IR: Monotonicity violation"
    case(SEAICE_ERROR_IR_CONS)
       errorStr = "IR: Conservation violation"
    case default
       errorStr = "Unknown error code"
    end select

  end function error_string

  !-----------------------------------------------------------------------

end module seaice_error
