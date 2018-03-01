module babortutils

  !-----------------------------------------------------------------------
  ! !MODULE: abortutils
  !
  ! !DESCRIPTION:
  ! Abort the model for abnormal termination
  !-----------------------------------------------------------------------

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  save

  public :: endrun
  interface endrun
     module procedure endrun_vanilla
     module procedure endrun_globalindex
  end interface

CONTAINS

  !-----------------------------------------------------------------------
  subroutine endrun_vanilla(msg)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Abort the model for abnormal termination
    !
    !
    ! !ARGUMENTS:
    use betr_ctrl      , only : iulog => biulog
    use bshr_sys_mod   , only: shr_sys_abort
    implicit none
    character(len=*), intent(in), optional :: msg    ! string to be printed
    !-----------------------------------------------------------------------

    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    call shr_sys_abort()

  end subroutine endrun_vanilla

  !-----------------------------------------------------------------------
  subroutine endrun_globalindex(decomp_index, clmlevel, msg)

    !-----------------------------------------------------------------------
    ! Description:
    ! Abort the model for abnormal termination
    !
    use betr_ctrl      , only : iulog  => biulog
    use bshr_sys_mod   , only: shr_sys_abort
    !
    ! Arguments:
    implicit none
    integer          , intent(in)           :: decomp_index
    character(len=*) , intent(in)           :: clmlevel
    character(len=*) , intent(in), optional :: msg    ! string to be printed
    !
    ! Local Variables:
    integer :: igrc, ilun, icol
    !-----------------------------------------------------------------------

    write(6,*)'calling getglobalwrite with decomp_index= ',decomp_index,' and clmlevel= ',trim(clmlevel)

    if (present (msg)) then
       write(iulog,*)'ENDRUN:', msg
    else
       write(iulog,*)'ENDRUN: called without a message string'
    end if

    call shr_sys_abort()

  end subroutine endrun_globalindex

end module babortutils
