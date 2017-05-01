!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_error

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_stream_manager
  use mpas_log

  implicit none

  private
  save

  public :: &
       cice_critical_error, &
       cice_finalize_critical_error

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_critical_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_critical_error(domain, block, ierr)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    integer, intent(inout) :: &
         ierr

    ! set error code
    ierr = 1

    ! write out block streams
    call mpas_stream_mgr_block_write(domain % streamManager, writeBlock=block, streamID='abort_block', forceWriteNow=.true.)

  end subroutine cice_critical_error

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_finalize_critical_error
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date   May 1st 2017
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_finalize_critical_error(domain, ierr)

    use mpas_dmpar, only: mpas_dmpar_max_int

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    integer, intent(inout) :: &
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
       call mpas_log_write("Critical error encountered", messageType=MPAS_LOG_CRIT)

    endif

  end subroutine cice_finalize_critical_error

  !-----------------------------------------------------------------------

end module cice_error
