module mpp_abortutils

  private
  save

  public :: endrun

  interface endrun
     module procedure endrun_vanilla
  end interface

CONTAINS

  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg) 

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    use mpp_varctl  , only: iulog
    !
    !
    ! !ARGUMENTS:
    implicit none
#include <mpif.h>         ! mpi library include file
    character(len=*), intent(in), optional :: msg    ! string to be printed
    integer :: rc, ierr
    !-----------------------------------------------------------------------

    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    !call shr_sys_abort()
    !call print_error_to_logs("ERROR", local_string)

    call sys_backtrace()

    !call shr_mpi_initialized(flag)

    !call shr_mpi_abort()
    call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)


    ! A compiler's abort method may print a backtrace or do other nice
    ! things, but in fact we can rarely leverage this, because MPI_Abort
    ! usually sends SIGTERM to the process, and we don't catch that signal.
    call abort()

  end subroutine endrun_vanilla

!-------------------------------------------------------------------------------
  subroutine sys_backtrace()

    ! This routine uses compiler-specific facilities to print a backtrace to
    ! error_unit (standard error, usually unit 0).

   use, intrinsic :: iso_fortran_env, only: error_unit

#if defined(CPRIBM)

    ! This theoretically should be in xlfutility, but using it from that
    ! module doesn't seem to always work.
    interface
       subroutine xl_trbk()
       end subroutine xl_trbk
    end interface

    call xl__trbk()

#elif defined(CPRGNU) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8 ))

    ! gfortran 4.8 and later implement this intrinsic. We explicitly call it
    ! out as such to make sure that it really is available, just in case the
    ! CPP logic above screws up.
    intrinsic :: backtrace

    call backtrace()

#elif defined(CPRINTEL)

    ! tracebackqq uses optional arguments, so *must* have an explicit
    ! interface.
    use ifcore, only: tracebackqq
    
    ! An exit code of -1 is a special value that prevents this subroutine
    ! from aborting the run.
    call tracebackqq(user_exit_code=-1)

#else

    ! Currently we have no means to request a backtrace from the NAG runtime,
    ! even though it is capable of emitting backtraces itself, if you use the
    ! "-gline" option.
    
    ! Similarly, PGI has a -traceback option, but no user interface for
    ! requesting a backtrace to be printed.

#endif

    flush(error_unit)

  end subroutine sys_backtrace

end module mpp_abortutils
