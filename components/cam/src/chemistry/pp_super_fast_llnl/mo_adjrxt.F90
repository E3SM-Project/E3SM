      module mo_adjrxt
      private
      public :: adjrxt
      contains
      subroutine adjrxt( rate, inv, m, ncol )
      use ppgrid, only : pver
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in) :: inv(ncol,pver,nfs)
      real(r8), intent(in) :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
      real(r8) :: im(ncol)
      do k = 1,pver
         rate(:,k, 14) = rate(:,k, 14) * inv(:,k, 1)
         rate(:,k, 15) = rate(:,k, 15) * inv(:,k, 5)
         rate(:,k, 23) = rate(:,k, 23) * inv(:,k, 4)
         rate(:,k, 26) = rate(:,k, 26) * inv(:,k, 1)
         rate(:,k, 7) = rate(:,k, 7) * m(:,k)
         rate(:,k, 8) = rate(:,k, 8) * m(:,k)
         rate(:,k, 9) = rate(:,k, 9) * m(:,k)
         rate(:,k, 10) = rate(:,k, 10) * m(:,k)
         rate(:,k, 11) = rate(:,k, 11) * m(:,k)
         rate(:,k, 12) = rate(:,k, 12) * m(:,k)
         rate(:,k, 13) = rate(:,k, 13) * m(:,k)
         rate(:,k, 14) = rate(:,k, 14) * m(:,k)
         rate(:,k, 16) = rate(:,k, 16) * m(:,k)
         rate(:,k, 17) = rate(:,k, 17) * m(:,k)
         rate(:,k, 18) = rate(:,k, 18) * m(:,k)
         rate(:,k, 19) = rate(:,k, 19) * m(:,k)
         rate(:,k, 20) = rate(:,k, 20) * m(:,k)
         rate(:,k, 21) = rate(:,k, 21) * m(:,k)
         rate(:,k, 22) = rate(:,k, 22) * m(:,k)
         rate(:,k, 24) = rate(:,k, 24) * m(:,k)
         rate(:,k, 25) = rate(:,k, 25) * m(:,k)
         rate(:,k, 26) = rate(:,k, 26) * m(:,k)
         rate(:,k, 27) = rate(:,k, 27) * m(:,k)
         rate(:,k, 28) = rate(:,k, 28) * m(:,k)
         rate(:,k, 29) = rate(:,k, 29) * m(:,k)
         rate(:,k, 30) = rate(:,k, 30) * m(:,k)
      end do
      end subroutine adjrxt
      end module mo_adjrxt
