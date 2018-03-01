




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
      real(r8) :: im(ncol,pver)


         rate(:,:, 3) = rate(:,:, 3) * inv(:,:, 6)
         rate(:,:, 4) = rate(:,:, 4) * inv(:,:, 6)
         rate(:,:, 5) = rate(:,:, 5) * inv(:,:, 6)
         rate(:,:, 6) = rate(:,:, 6) * inv(:,:, 6)
         rate(:,:, 7) = rate(:,:, 7) * inv(:,:, 7)
         im(:,:) = 1._r8 / m(:,:)
         rate(:,:, 2) = rate(:,:, 2) * inv(:,:, 8) * inv(:,:, 8) * im(:,:)

      end subroutine adjrxt

      end module mo_adjrxt
