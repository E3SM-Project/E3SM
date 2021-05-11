
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

      rate(:,:,90) = 1.20e-10_r8
      rate(:,:,91) = 2.02e-10_r8
      rate(:,:,92) = 1.204e-10_r8
      rate(:,:,93) = 1.50e-10_r8
      rate(:,:,94) = 9.75e-11_r8
      rate(:,:,95) = 1.50e-11_r8
      rate(:,:,96) = 7.20e-11_r8
      rate(:,:,97) = 1.794e-10_r8
      rate(:,:,98) = 1.628e-10_r8
      rate(:,:,99) = 2.84e-10_r8
      rate(:,:,100) = 1.674e-10_r8
      rate(:,:,101) = 9.60e-11_r8
      rate(:,:,102) = 4.10e-11_r8
      rate(:,:,103) = 1.012e-10_r8
      rate(:,:,104) = 1.20e-10_r8
      rate(:,:,105) = 4.49e-10_r8
      rate(:,:,106) = 2.57e-10_r8
      rate(:,:,107) = 1.31e-10_r8
      rate(:,:,108) = 3.50e-11_r8
      rate(:,:,109) = 9.00e-12_r8
      rate(:,:,110) = 1.20e-10_r8
      rate(:,:,111) = 1.50e-10_r8
      rate(:,:,112) = 1.20e-10_r8
      rate(:,:,116) = 7.20e-11_r8
      rate(:,:,117) = 6.90e-12_r8
      rate(:,:,118) = 1.60e-12_r8
      rate(:,:,122) = 1.80e-12_r8
      rate(:,:,125) = 1.80e-12_r8
      rate(:,:,149) = 1.00e-11_r8
      rate(:,:,150) = 2.20e-11_r8
      rate(:,:,151) = 3.50e-12_r8
      rate(:,:,176) = 1.70e-13_r8
      rate(:,:,223) = 4.50e-13_r8
      rate(:,:,233) = 1.00e-14_r8
      rate(:,:,236) = 7.00e-13_r8
      rate(:,:,239) = 2.00e-13_r8
      rate(:,:,240) = 6.80e-14_r8
      rate(:,:,249) = 1.00e-12_r8
      rate(:,:,250) = 1.00e-11_r8
      rate(:,:,251) = 1.15e-11_r8
      rate(:,:,254) = 4.00e-14_r8
      rate(:,:,271) = 3.00e-12_r8
      rate(:,:,274) = 6.80e-13_r8
      rate(:,:,275) = 5.40e-11_r8
      rate(:,:,287) = 2.40e-12_r8
      rate(:,:,290) = 1.40e-11_r8
      rate(:,:,293) = 5.00e-12_r8
      rate(:,:,305) = 2.40e-12_r8
      rate(:,:,309) = 1.40e-11_r8
      rate(:,:,311) = 2.40e-12_r8
      rate(:,:,313) = 3.50e-12_r8
      rate(:,:,314) = 4.50e-11_r8
      rate(:,:,321) = 2.40e-12_r8
      rate(:,:,331) = 3.00e-12_r8
      rate(:,:,332) = 1.00e-11_r8
      rate(:,:,339) = 2.10e-6_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,83) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,85) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      rate(:,:,86) = 3.30e-11_r8 * exp( 55._r8 * itemp(:,:) )
      rate(:,:,87) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,88) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,89) = 4.63e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,113) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,134) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,115) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,119) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,231) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,258) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,263) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,276) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,280) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,304) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,317) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,328) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,336) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,120) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,121) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,124) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,126) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,127) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,194) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,222) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,241) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,261) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,265) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,270) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,282) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,291) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,307) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,319) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,330) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,338) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,128) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,130) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,302) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,132) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      rate(:,:,133) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,135) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,136) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,137) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,139) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,158) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,163) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,244) = 8.10e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,140) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,195) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,141) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,143) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,169) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,148) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,153) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,155) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,156) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,157) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,159) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,160) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,161) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,162) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,164) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,185) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,193) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,165) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,167) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,192) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,166) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,170) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,171) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,174) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,175) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,177) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,178) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,179) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,206) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,180) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,181) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,182) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,183) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,184) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,208) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,188) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,191) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,196) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,197) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,198) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,199) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,200) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,201) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,204) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,215) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,202) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,203) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,205) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,207) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,209) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,210) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,213) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,214) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,216) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,217) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,267) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,219) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,220) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,221) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,224) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,225) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,226) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,232) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,238) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,259) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,264) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,268) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,281) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,288) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,306) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,312) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,318) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,322) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,329) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,337) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,227) = 7.20e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,229) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,234) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,235) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,237) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      rate(:,:,242) = 4.63e-12_r8 * exp( 350._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,243) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,256) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,246) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,294) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,247) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,248) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,269) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,295) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,252) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,257) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,260) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,262) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,272) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,273) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,315) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,277) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,278) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,279) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,283) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,316) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,284) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,285) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,286) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,292) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,310) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,320) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,289) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,308) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,323) = 5.e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,296) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,297) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,301) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,303) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,324) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,325) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,327) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,333) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,334) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,335) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,344) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,346) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,347) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,114), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,123), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,131), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,138), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,142), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,144), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,146), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,152), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,168), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,172), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,189), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,212), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,228), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,230), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,245), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,255), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,300), m, 0.5_r8, ko, kinf, n )

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
