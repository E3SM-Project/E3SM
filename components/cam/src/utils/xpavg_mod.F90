module xpavg_mod
  !----------------------------------------------------------------------- 
  !BOP
  ! !ROUTINE:  xpaxg --- Average a scalar latitude field
  !
  ! !INTERFACE:
  implicit none
  private
  public :: xpavg
contains
  subroutine xpavg(p, im)

    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none

    ! !INPUT PARAMETERS:
    integer, intent(in) ::  im

    ! !INPUT/OUTPUT PARAMETERS:
    real(r8), intent(inout):: p(:)

    !
    ! !DESCRIPTION:
    !   This routine determines the average of the scalar latitude field p
    !   and then sets all the values in p to the average.
    !
    ! !REVISION HISTORY: 
    !   ??.??.??    Lin?       Creation
    !   01.03.26    Sawyer     Added ProTeX documentation
    !   07.04.06    Edwards    moved to utils, placed in module form
    !
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter ::  D0_0                    =  0.0_r8

    integer :: i
    real(r8) :: sum1

    sum1 = D0_0
    do i=1,im
       sum1 = sum1 + p(i)
    enddo
    sum1 = sum1 / im

    do i=1,im
       p(i) = sum1
    enddo

    return
    !EOC
  end subroutine xpavg
  !-----------------------------------------------------------------------
end module xpavg_mod
