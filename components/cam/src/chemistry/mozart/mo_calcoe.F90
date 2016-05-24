
      module mo_calcoe

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: calcoe

      contains

      subroutine calcoe( c, xz, tt, adjin, adjcoe )
!-----------------------------------------------------------------------------
!   parameters:
!   adjcoe - real(r8), coross section adjust coefficients (in and out)
!   c(5,28)-polynomal coef
!   tt     -nomarlized temperature
!-----------------------------------------------------------------------------*

      use ppgrid,   only : pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)    :: adjin
      real(r8), intent(in)    :: tt
      real(r8), intent(in)    :: c(5)
      real(r8), intent(in)    :: xz(pverp)
      real(r8), intent(inout) :: adjcoe(:)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer  :: k
      real(r8) :: x

      do k = 1,pverp
	 x = xz(k)
         adjcoe(k) = adjin * (1._r8 + .01_r8*(c(1) + x*(c(2) + x*(c(3) + x*(c(4) + x*c(5))))))
      end do
        
      end subroutine calcoe

      end module mo_calcoe
