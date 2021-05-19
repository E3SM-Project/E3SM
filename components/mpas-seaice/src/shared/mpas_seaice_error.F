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

  subroutine seaice_critical_error_write_block(domain, block, abortFlag)

    type(domain_type), intent(in) :: &
         domain

    type(block_type), intent(in) :: &
         block

    logical, intent(in) :: &
         abortFlag

    if (abortFlag) then

       ! write out block streams
       call mpas_stream_mgr_block_write(&
            domain % streamManager, &
            writeBlock=block, &
            streamID='abort_block', &
            forceWriteNow=.true.)

    endif ! abortFlag

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

  subroutine seaice_check_critical_error(domain, abortFlag)

    use mpas_dmpar, only: mpas_dmpar_max_int

    type(domain_type), intent(in) :: &
         domain !< Input/Output:

    logical, intent(in) :: &
         abortFlag

    logical, pointer :: &
         config_full_abort_write

    integer :: &
         ierr, &
         ierrmax

    call MPAS_pool_get_config(domain % configs, "config_full_abort_write", config_full_abort_write)
    if (config_full_abort_write) then

       ! find if anyone failed
       ierr = 0
       if (abortFlag) ierr = 1
       call mpas_dmpar_max_int(domain % dminfo, ierr, ierrmax)
       if (ierrmax > 0) then

          ! make final full write of abort stream
          call MPAS_stream_mgr_write(domain % streamManager, streamID="abort", forceWriteNow=.true.)

       endif
    endif

    ! check if a critical error occured
    if (abortFlag) call mpas_log_write("Runtime error", messageType=MPAS_LOG_CRIT)

  end subroutine seaice_check_critical_error

  !-----------------------------------------------------------------------

end module seaice_error
