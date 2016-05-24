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
         rate(:,k, 3) = rate(:,k, 3) * inv(:,k, 6)
         rate(:,k, 4) = rate(:,k, 4) * inv(:,k, 6)
         rate(:,k, 5) = rate(:,k, 5) * inv(:,k, 6)
         rate(:,k, 6) = rate(:,k, 6) * inv(:,k, 6)
         rate(:,k, 7) = rate(:,k, 7) * inv(:,k, 7)
         rate(:,k, 8) = rate(:,k, 8) * inv(:,k, 6)
         im(:) = 1._r8 / m(:,k)
         rate(:,k, 2) = rate(:,k, 2) * inv(:,k, 8) * inv(:,k, 8) * im(:)
      end do
      end subroutine adjrxt
      end module mo_adjrxt
