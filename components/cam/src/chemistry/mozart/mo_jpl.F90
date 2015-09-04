
      module mo_jpl

      private
      public :: jpl

      contains

      subroutine jpl( rate, m, factor, ko, kinf, ncol )
!-----------------------------------------------------------------
!        ... Calculate JPL troe rate
!-----------------------------------------------------------------
      
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------
      integer, intent(in)   ::   ncol
      real(r8), intent(in)  ::   factor
      real(r8), intent(in)  ::   ko(ncol)
      real(r8), intent(in)  ::   kinf(ncol)
      real(r8), intent(in)  ::   m(ncol)
      real(r8), intent(out) ::   rate(ncol)

!-----------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------
      real(r8)  ::  xpo(ncol)

      xpo(:)  = ko(:) * m(:) / kinf(:)
      rate(:) = ko(:) / (1._r8 + xpo(:))
      xpo(:)  = log10( xpo(:) )
      xpo(:)  = 1._r8 / (1._r8 + xpo(:)*xpo(:))
      rate(:) = rate(:) * factor**xpo(:)

      end subroutine jpl

      end module mo_jpl
