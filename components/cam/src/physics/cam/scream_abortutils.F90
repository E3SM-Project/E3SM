module scream_abortutils

!-------------------------------------------------
!Utility to stop the model in case of
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
#endif

    implicit none

    !intent-ins
    character(len=*), intent(in), optional :: msg
#ifdef SPMD
#include "mpif.h"
    integer:: ierr
#endif


#ifdef SPMD
    !for model runs with MPI
    if(present(msg)) then
       write(iulog_e3sm,*)msg
    else
       write(iulog_e3sm,*)'ERROR: Aborting...'
    endif
    call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
#else
    !Stop the model when run in non-MPI mode
    write(*,*)'ERROR: Aborting...'
    if(present(msg)) write(*,*)trim(adjustl(msg))

    call abort()
#endif

  end subroutine endscreamrun

end module scream_abortutils
