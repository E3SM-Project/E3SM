module quadraticMod

  use abortutils  ,   only: endrun
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_log_mod ,   only: errMsg => shr_log_errMsg
  use elm_varctl  ,   only: iulog
  use elm_varcon  ,   only : spval

  implicit none

  public :: quadratic

contains

  pure subroutine quadratic (a, b, c, r1, r2)
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     !==============================================================================!
     ! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
     ! Computing (Cambridge University Press, Cambridge), pp. 145.
     !
     ! !REVISION HISTORY:
     ! 4/5/10: Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
     !
     ! !USES:
     implicit none
     !$acc routine seq
     !
     ! !ARGUMENTS:
     real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
     real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
     !
     ! !LOCAL VARIABLES:
     real(r8) :: q                        ! Temporary term for quadratic solution
     !------------------------------------------------------------------------------

     if (b >= 0._r8) then
        q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
     else
        q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
     end if

     r1 = q / a
     if (q /= 0._r8) then
        r2 = c / q
     else
        r2 = spval
     end if

  end subroutine quadratic

end module quadraticMod
