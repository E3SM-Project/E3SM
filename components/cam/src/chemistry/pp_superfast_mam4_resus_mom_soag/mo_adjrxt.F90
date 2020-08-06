













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
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: inv(ncol,pver,nfs)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!--------------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------------
      real(r8) :: im(ncol,pver)


         rate(:,:, 14) = rate(:,:, 14) * inv(:,:, 1)
         rate(:,:, 15) = rate(:,:, 15) * inv(:,:, 6)
         rate(:,:, 23) = rate(:,:, 23) * inv(:,:, 4)
         rate(:,:, 29) = rate(:,:, 29) * inv(:,:, 5)
         rate(:,:,  7) = rate(:,:,  7) * m(:,:)
         rate(:,:,  8) = rate(:,:,  8) * m(:,:)
         rate(:,:,  9) = rate(:,:,  9) * m(:,:)
         rate(:,:, 10) = rate(:,:, 10) * m(:,:)
         rate(:,:, 11) = rate(:,:, 11) * m(:,:)
         rate(:,:, 12) = rate(:,:, 12) * m(:,:)
         rate(:,:, 13) = rate(:,:, 13) * m(:,:)
         rate(:,:, 14) = rate(:,:, 14) * m(:,:)
         rate(:,:, 16) = rate(:,:, 16) * m(:,:)
         rate(:,:, 17) = rate(:,:, 17) * m(:,:)
         rate(:,:, 18) = rate(:,:, 18) * m(:,:)
         rate(:,:, 19) = rate(:,:, 19) * m(:,:)
         rate(:,:, 20) = rate(:,:, 20) * m(:,:)
         rate(:,:, 21) = rate(:,:, 21) * m(:,:)
         rate(:,:, 22) = rate(:,:, 22) * m(:,:)
         rate(:,:, 24) = rate(:,:, 24) * m(:,:)
         rate(:,:, 25) = rate(:,:, 25) * m(:,:)
         rate(:,:, 26) = rate(:,:, 26) * m(:,:)
         rate(:,:, 27) = rate(:,:, 27) * m(:,:)
         rate(:,:, 28) = rate(:,:, 28) * m(:,:)

      end subroutine adjrxt

      end module mo_adjrxt
