











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


         rate(:,:, 23) = rate(:,:, 23) * inv(:,:, 4)
         rate(:,:, 24) = rate(:,:, 24) * inv(:,:, 5)
         rate(:,:, 27) = rate(:,:, 27) * inv(:,:, 1)
         rate(:,:, 28) = rate(:,:, 28) * inv(:,:, 5)
         rate(:,:, 32) = rate(:,:, 32) * inv(:,:, 1)
         rate(:,:, 41) = rate(:,:, 41) * inv(:,:, 1)
         rate(:,:, 43) = rate(:,:, 43) * inv(:,:, 1)
         rate(:,:, 53) = rate(:,:, 53) * inv(:,:, 1)
         rate(:,:, 54) = rate(:,:, 54) * inv(:,:, 1)
         rate(:,:, 55) = rate(:,:, 55) * inv(:,:, 1)
         rate(:,:, 56) = rate(:,:, 56) * inv(:,:, 1)
         rate(:,:, 57) = rate(:,:, 57) * inv(:,:, 1)
         rate(:,:, 58) = rate(:,:, 58) * inv(:,:, 1)
         rate(:,:, 59) = rate(:,:, 59) * inv(:,:, 1)
         rate(:,:, 92) = rate(:,:, 92) * inv(:,:, 9)
         rate(:,:, 93) = rate(:,:, 93) * inv(:,:, 9)
         rate(:,:, 94) = rate(:,:, 94) * inv(:,:, 9)
         rate(:,:, 95) = rate(:,:, 95) * inv(:,:, 9)
         rate(:,:, 96) = rate(:,:, 96) * inv(:,:, 9)
         rate(:,:, 97) = rate(:,:, 97) * inv(:,:, 7)
         rate(:,:, 98) = rate(:,:, 98) * inv(:,:, 7)
         rate(:,:, 99) = rate(:,:, 99) * inv(:,:, 8)
         rate(:,:,100) = rate(:,:,100) * inv(:,:, 8)
         rate(:,:,101) = rate(:,:,101) * inv(:,:, 9)
         rate(:,:,102) = rate(:,:,102) * inv(:,:, 9)
         rate(:,:,103) = rate(:,:,103) * inv(:,:, 9)
         rate(:,:,104) = rate(:,:,104) * inv(:,:, 9)
         rate(:,:, 25) = rate(:,:, 25) * m(:,:)
         rate(:,:, 26) = rate(:,:, 26) * m(:,:)
         rate(:,:, 27) = rate(:,:, 27) * m(:,:)
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
         rate(:,:, 54) = rate(:,:, 54) * m(:,:)
         rate(:,:, 55) = rate(:,:, 55) * m(:,:)
         rate(:,:, 56) = rate(:,:, 56) * m(:,:)
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
         rate(:,:, 84) = rate(:,:, 84) * m(:,:)
         rate(:,:, 85) = rate(:,:, 85) * m(:,:)
         rate(:,:, 86) = rate(:,:, 86) * m(:,:)
         rate(:,:, 88) = rate(:,:, 88) * m(:,:)
         rate(:,:, 89) = rate(:,:, 89) * m(:,:)
         rate(:,:, 90) = rate(:,:, 90) * m(:,:)
         rate(:,:, 91) = rate(:,:, 91) * m(:,:)

      end subroutine adjrxt

      end module mo_adjrxt
