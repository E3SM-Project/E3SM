













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


         rate(:,:, 20) = rate(:,:, 20) * inv(:,:, 4)
         rate(:,:, 21) = rate(:,:, 21) * inv(:,:, 5)
         rate(:,:, 22) = rate(:,:, 22) * inv(:,:, 6)
         rate(:,:, 24) = rate(:,:, 24) * inv(:,:, 1)
         rate(:,:, 25) = rate(:,:, 25) * inv(:,:, 5)
         rate(:,:, 26) = rate(:,:, 26) * inv(:,:, 6)
         rate(:,:, 29) = rate(:,:, 29) * inv(:,:, 1)
         rate(:,:, 38) = rate(:,:, 38) * inv(:,:, 1)
         rate(:,:, 40) = rate(:,:, 40) * inv(:,:, 1)
         rate(:,:, 50) = rate(:,:, 50) * inv(:,:, 1)
         rate(:,:, 51) = rate(:,:, 51) * inv(:,:, 1)
         rate(:,:, 52) = rate(:,:, 52) * inv(:,:, 1)
         rate(:,:, 53) = rate(:,:, 53) * inv(:,:, 1)
         rate(:,:, 23) = rate(:,:, 23) * m(:,:)
         rate(:,:, 24) = rate(:,:, 24) * m(:,:)
         rate(:,:, 27) = rate(:,:, 27) * m(:,:)
         rate(:,:, 28) = rate(:,:, 28) * m(:,:)
         rate(:,:, 29) = rate(:,:, 29) * m(:,:)
         rate(:,:, 30) = rate(:,:, 30) * m(:,:)
         rate(:,:, 31) = rate(:,:, 31) * m(:,:)
         rate(:,:, 32) = rate(:,:, 32) * m(:,:)
         rate(:,:, 33) = rate(:,:, 33) * m(:,:)
         rate(:,:, 34) = rate(:,:, 34) * m(:,:)
         rate(:,:, 35) = rate(:,:, 35) * m(:,:)
         rate(:,:, 36) = rate(:,:, 36) * m(:,:)
         rate(:,:, 37) = rate(:,:, 37) * m(:,:)
         rate(:,:, 38) = rate(:,:, 38) * m(:,:)
         rate(:,:, 39) = rate(:,:, 39) * m(:,:)
         rate(:,:, 40) = rate(:,:, 40) * m(:,:)
         rate(:,:, 41) = rate(:,:, 41) * m(:,:)
         rate(:,:, 42) = rate(:,:, 42) * m(:,:)
         rate(:,:, 43) = rate(:,:, 43) * m(:,:)
         rate(:,:, 44) = rate(:,:, 44) * m(:,:)
         rate(:,:, 45) = rate(:,:, 45) * m(:,:)
         rate(:,:, 46) = rate(:,:, 46) * m(:,:)
         rate(:,:, 47) = rate(:,:, 47) * m(:,:)
         rate(:,:, 48) = rate(:,:, 48) * m(:,:)
         rate(:,:, 49) = rate(:,:, 49) * m(:,:)
         rate(:,:, 50) = rate(:,:, 50) * m(:,:)
         rate(:,:, 51) = rate(:,:, 51) * m(:,:)
         rate(:,:, 52) = rate(:,:, 52) * m(:,:)
         rate(:,:, 53) = rate(:,:, 53) * m(:,:)
         rate(:,:, 57) = rate(:,:, 57) * m(:,:)
         rate(:,:, 58) = rate(:,:, 58) * m(:,:)
         rate(:,:, 59) = rate(:,:, 59) * m(:,:)
         rate(:,:, 60) = rate(:,:, 60) * m(:,:)
         rate(:,:, 61) = rate(:,:, 61) * m(:,:)
         rate(:,:, 62) = rate(:,:, 62) * m(:,:)
         rate(:,:, 63) = rate(:,:, 63) * m(:,:)
         rate(:,:, 64) = rate(:,:, 64) * m(:,:)
         rate(:,:, 65) = rate(:,:, 65) * m(:,:)
         rate(:,:, 66) = rate(:,:, 66) * m(:,:)
         rate(:,:, 67) = rate(:,:, 67) * m(:,:)
         rate(:,:, 68) = rate(:,:, 68) * m(:,:)
         rate(:,:, 69) = rate(:,:, 69) * m(:,:)
         rate(:,:, 70) = rate(:,:, 70) * m(:,:)
         rate(:,:, 71) = rate(:,:, 71) * m(:,:)
         rate(:,:, 72) = rate(:,:, 72) * m(:,:)
         rate(:,:, 73) = rate(:,:, 73) * m(:,:)
         rate(:,:, 74) = rate(:,:, 74) * m(:,:)
         rate(:,:, 75) = rate(:,:, 75) * m(:,:)
         rate(:,:, 76) = rate(:,:, 76) * m(:,:)
         rate(:,:, 77) = rate(:,:, 77) * m(:,:)
         rate(:,:, 78) = rate(:,:, 78) * m(:,:)
         rate(:,:, 79) = rate(:,:, 79) * m(:,:)
         rate(:,:, 80) = rate(:,:, 80) * m(:,:)
         rate(:,:, 81) = rate(:,:, 81) * m(:,:)
         rate(:,:, 82) = rate(:,:, 82) * m(:,:)
         rate(:,:, 83) = rate(:,:, 83) * m(:,:)
         rate(:,:, 85) = rate(:,:, 85) * m(:,:)
         rate(:,:, 86) = rate(:,:, 86) * m(:,:)
         rate(:,:, 87) = rate(:,:, 87) * m(:,:)
         rate(:,:, 88) = rate(:,:, 88) * m(:,:)

      end subroutine adjrxt

      end module mo_adjrxt
