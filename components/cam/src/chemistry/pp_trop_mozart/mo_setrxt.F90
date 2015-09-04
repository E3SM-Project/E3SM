
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,45) = 2.2e-10_r8
      rate(:,:,46) = 1.1e-10_r8
      rate(:,:,55) = 1.8e-12_r8
      rate(:,:,57) = 4.9e-11_r8
      rate(:,:,58) = 6.7e-11_r8
      rate(:,:,63) = 3.5e-12_r8
      rate(:,:,73) = 1.5e-10_r8
      rate(:,:,80) = 9.e-12_r8
      rate(:,:,84) = 4.5e-13_r8
      rate(:,:,93) = 1.e-14_r8
      rate(:,:,98) = 2.e-13_r8
      rate(:,:,99) = 6.8e-14_r8
      rate(:,:,107) = 1e-12_r8
      rate(:,:,108) = 4.e-14_r8
      rate(:,:,111) = 1.e-11_r8
      rate(:,:,112) = 1.1e-11_r8
      rate(:,:,113) = 7.e-13_r8
      rate(:,:,131) = 6.8e-13_r8
      rate(:,:,134) = 3.e-12_r8
      rate(:,:,135) = 5.4e-11_r8
      rate(:,:,142) = 3.5e-12_r8
      rate(:,:,149) = 2.4e-12_r8
      rate(:,:,153) = 1.4e-11_r8
      rate(:,:,156) = 2.4e-12_r8
      rate(:,:,164) = 2.4e-12_r8
      rate(:,:,167) = 1.4e-11_r8
      rate(:,:,170) = 5.e-12_r8
      rate(:,:,177) = 4.5e-11_r8
      rate(:,:,181) = 2.40e-12_r8
      rate(:,:,188) = 3.e-12_r8
      rate(:,:,189) = 1.e-11_r8
      rate(:,:,199) = 2.1e-6_r8
      rate(:,:,203) = 7.1e-6_r8
      rate(:,:,209) = 7.1e-6_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,42) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,43) = 2.1e-11_r8 * exp( 115._r8 * itemp(:,:) )
      rate(:,:,44) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,47) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,48) = 2.2e-11_r8 * exp( 120._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,49) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,78) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,100) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,120) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,125) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,130) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,140) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,145) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,151) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,168) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,192) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,198) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,50) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,51) = 1.e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,53) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,54) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,59) = 3.5e-12_r8 * exp_fac(:,:)
      rate(:,:,60) = 3e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:,61) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:,62) = 1.2e-13_r8 * exp( -2450._r8 * itemp(:,:) )
      rate(:,:,68) = 1.5e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,70) = 1.3e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,72) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,74) = 2.8e-12_r8 * exp_fac(:,:)
      rate(:,:,127) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,75) = 5.e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,76) = 1.9e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,77) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,79) = 6.0e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,83) = 7.3e-12_r8 * exp( -620._r8 * itemp(:,:) )
      rate(:,:,85) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      rate(:,:,86) = 2.4e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,87) = 2.6e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,88) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,97) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,119) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,123) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,128) = 8.6e-13_r8 * exp_fac(:,:)
      rate(:,:,139) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,144) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,150) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,157) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,165) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,182) = 8.e-13_r8 * exp_fac(:,:)
      rate(:,:,191) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,197) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,91) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,92) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,118) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,122) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,136) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,138) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,143) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,148) = 4.4e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,196) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,94) = 1.6e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,95) = 8.7e-12_r8 * exp( -1070._r8 * itemp(:,:) )
      rate(:,:,96) = 2.6e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,101) = 5.6e-12_r8 * exp_fac(:,:)
      rate(:,:,103) = 8.1e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,102) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,116) = 6.5e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,105) = 4.3e-13_r8 * exp_fac(:,:)
      rate(:,:,171) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,106) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,110) = 2.5e-12_r8 * exp_fac(:,:)
      rate(:,:,129) = 7.1e-13_r8 * exp_fac(:,:)
      rate(:,:,172) = 2.0e-12_r8 * exp_fac(:,:)
      rate(:,:,114) = 6.9e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,117) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,121) = 1.0e-11_r8 * exp( -665._r8 * itemp(:,:) )
      rate(:,:,124) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,132) = 8.4e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,133) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,178) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,137) = 2.3e-12_r8 * exp( -170._r8 * itemp(:,:) )
      rate(:,:,146) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,147) = 1.05e-14_r8 * exp( -2000._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,152) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,166) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,183) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,154) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,155) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,162) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,163) = 1.3e-13_r8 * exp_fac(:,:)
      rate(:,:,169) = 5.3e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 2.7e-12_r8 * exp_fac(:,:)
      rate(:,:,158) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,159) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,160) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,179) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,161) = 4.4e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,173) = 4.6e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,184) = 1.30e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,185) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,187) = 1.7e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,193) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,194) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,195) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,205) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,207) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,208) = 1.7e-12_r8 * exp( -710._r8 * itemp(:,:) )
      rate(:,:,212) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 6.9e-31_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,56), m, .6_r8, ko, kinf, n )

      ko(:,:) = 2.e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**.7_r8
      call jpl( rate(1,1,64), m, .6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,66), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.0e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,69), m, .6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,81), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.5e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2._r8)
      call jpl( rate(1,1,89), m, .6_r8, ko, kinf, n )

      ko(:,:) = 1.e-28_r8 * itemp(:,:)**.8_r8
      kinf(:,:) = 8.8e-12_r8
      call jpl( rate(1,1,90), m, .6_r8, ko, kinf, n )

      ko(:,:) = 8.5e-29_r8 * itemp(:,:)**6.5_r8
      kinf(:,:) = 1.1e-11_r8 * itemp(:,:)
      call jpl( rate(1,1,104), m, .6_r8, ko, kinf, n )

      ko(:,:) = 8.e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.e-11_r8
      call jpl( rate(1,1,115), m, .5_r8, ko, kinf, n )

      ko(:,:) = 8.e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.e-11_r8
      call jpl( rate(1,1,141), m, .5_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,211), m, 0.8_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
