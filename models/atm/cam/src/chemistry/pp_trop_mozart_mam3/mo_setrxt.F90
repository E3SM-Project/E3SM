
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

      rate(:,:,47) = 1.20e-10_r8
      rate(:,:,49) = 1.20e-10_r8
      rate(:,:,52) = 1.31e-10_r8
      rate(:,:,53) = 3.50e-11_r8
      rate(:,:,54) = 9.00e-12_r8
      rate(:,:,58) = 7.20e-11_r8
      rate(:,:,59) = 6.90e-12_r8
      rate(:,:,60) = 1.60e-12_r8
      rate(:,:,64) = 1.80e-12_r8
      rate(:,:,67) = 1.80e-12_r8
      rate(:,:,86) = 1.00e-11_r8
      rate(:,:,87) = 2.20e-11_r8
      rate(:,:,88) = 3.50e-12_r8
      rate(:,:,105) = 4.50e-13_r8
      rate(:,:,114) = 1.00e-14_r8
      rate(:,:,117) = 7.00e-13_r8
      rate(:,:,120) = 2.00e-13_r8
      rate(:,:,121) = 6.80e-14_r8
      rate(:,:,130) = 1.00e-12_r8
      rate(:,:,131) = 1.00e-11_r8
      rate(:,:,132) = 1.15e-11_r8
      rate(:,:,135) = 4.00e-14_r8
      rate(:,:,152) = 3.00e-12_r8
      rate(:,:,155) = 6.80e-13_r8
      rate(:,:,156) = 5.40e-11_r8
      rate(:,:,168) = 2.40e-12_r8
      rate(:,:,171) = 1.40e-11_r8
      rate(:,:,174) = 5.00e-12_r8
      rate(:,:,186) = 2.40e-12_r8
      rate(:,:,190) = 1.40e-11_r8
      rate(:,:,192) = 2.40e-12_r8
      rate(:,:,194) = 3.50e-12_r8
      rate(:,:,195) = 4.50e-11_r8
      rate(:,:,202) = 2.40e-12_r8
      rate(:,:,212) = 3.00e-12_r8
      rate(:,:,213) = 1.00e-11_r8
      rate(:,:,220) = 2.10e-6_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,44) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,45) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,46) = 1.65e-12_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,48) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,50) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,51) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,55) = 7.70e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:,57) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,61) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,112) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,139) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,144) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,157) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,161) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,185) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,198) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,209) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,217) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,62) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,63) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:,66) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,68) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,69) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,104) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,122) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,142) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,151) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,163) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,172) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,188) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,200) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,211) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,219) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,70) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,72) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,183) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,74) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,76) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,125) = 8.10e-12_r8 * exp_fac(:,:)
      rate(:,:,77) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:,78) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:,80) = 1.20e-13_r8 * exp( -2450._r8 * itemp(:,:) )
      rate(:,:,85) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,90) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,92) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,95) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,96) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,97) = 3.40e-11_r8 * exp( -1600._r8 * itemp(:,:) )
      rate(:,:,98) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,99) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,148) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,100) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,101) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,102) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,103) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,106) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,107) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,108) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,113) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,119) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,140) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,145) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,149) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,162) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,169) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,187) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,193) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,199) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,203) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,210) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,218) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,110) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,115) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,116) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,118) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      rate(:,:,123) = 4.63e-12_r8 * exp( 350._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,124) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,137) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,127) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,175) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,128) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,129) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,150) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,176) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,133) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,138) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,141) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,143) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,153) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,154) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,196) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,158) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,159) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,160) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,164) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,197) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,165) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,166) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,167) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,173) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,191) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,201) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,170) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,189) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,204) = 5.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,177) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,178) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,182) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,184) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,205) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,206) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,208) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,214) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,215) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,216) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,225) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,227) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,228) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,56), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,65), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,73), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,75), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,79), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,81), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,83), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,89), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,94), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,109), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,111), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,126), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,136), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,181), m, 0.5_r8, ko, kinf, n )

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
