module scream_abortutils

!-------------------------------------------------
!Utilities to stop the model in case of
!catastrophic errors
!-------------------------------------------------
implicit none
private

!public subroutines
public :: endscreamrun

contains

  subroutine endscreamrun (msg)
    !-------------------------------------------------
    ! This subroutine will print the optional message
    ! received via optional arg "msg" and stop
    ! the simulation
    !-------------------------------------------------

#ifdef SPMD
    !Modules to use when SCREAM run with MPI under E3SM
    use micro_p3_utils, only:iulog_e3sm
    use cam_abortutils, only: endrun
#endif

    implicit none

    !intent-ins
    character(len=*), intent(in), optional :: msg

#ifdef SPMD
    !Call E3SM's utility to stop the model
    !for model runs with MPI
    if(present(msg)) then
       call endrun(msg)
    else
       call endrun()
    endif
#else
    !Stop the model when run in non-MPI mode
    write(*,*)'ERROR:'
    if(present(msg)) write(*,*)trim(adjustl(msg))

    stop
#endif

  end subroutine endscreamrun

end module scream_abortutils
