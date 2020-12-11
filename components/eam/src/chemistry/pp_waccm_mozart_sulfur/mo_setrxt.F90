
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

      rate(:,:,83) = 8.00e-14_r8
      rate(:,:,84) = 3.90e-17_r8
      rate(:,:,87) = 4.20e-13_r8
      rate(:,:,88) = 8.50e-2_r8
      rate(:,:,89) = 1.30e-16_r8
      rate(:,:,91) = 1.00e-20_r8
      rate(:,:,92) = 2.58e-04_r8
      rate(:,:,99) = 1.20e-10_r8
      rate(:,:,100) = 1.70e-10_r8
      rate(:,:,101) = 1.20e-10_r8
      rate(:,:,102) = 1.50e-10_r8
      rate(:,:,103) = 7.20e-11_r8
      rate(:,:,104) = 2.84e-10_r8
      rate(:,:,105) = 1.80e-10_r8
      rate(:,:,106) = 9.60e-11_r8
      rate(:,:,107) = 4.10e-11_r8
      rate(:,:,108) = 1.125e-10_r8
      rate(:,:,109) = 3.00e-11_r8
      rate(:,:,110) = 7.50e-12_r8
      rate(:,:,111) = 1.10e-10_r8
      rate(:,:,112) = 1.50e-10_r8
      rate(:,:,113) = 1.50e-10_r8
      rate(:,:,114) = 5.00e-12_r8
      rate(:,:,115) = 7.00e-13_r8
      rate(:,:,130) = 1.00e-11_r8
      rate(:,:,131) = 2.20e-11_r8
      rate(:,:,132) = 3.50e-12_r8
      rate(:,:,140) = 5.80e-16_r8
      rate(:,:,147) = 7.20e-11_r8
      rate(:,:,148) = 6.90e-12_r8
      rate(:,:,149) = 1.60e-12_r8
      rate(:,:,153) = 1.80e-12_r8
      rate(:,:,156) = 1.80e-12_r8
      rate(:,:,181) = 1.70e-13_r8
      rate(:,:,210) = 6.60E-11_r8
      rate(:,:,211) = 2.30E-12_r8
      rate(:,:,212) = 1.20E-11_r8
      rate(:,:,216) = 1.40E-11_r8
      rate(:,:,217) = 2.80E-11_r8
      rate(:,:,218) = 5.70E-11_r8
      rate(:,:,219) = 1.90E-12_r8
      rate(:,:,243) = 9.e-10_r8
      rate(:,:,244) = 1.e-10_r8
      rate(:,:,245) = 4.4e-10_r8
      rate(:,:,246) = 4.e-10_r8
      rate(:,:,247) = 2.e-10_r8
      rate(:,:,248) = 1.e-12_r8
      rate(:,:,249) = 6.e-11_r8
      rate(:,:,250) = 5.e-16_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,81) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,85) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,86) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,90) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,93) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,94) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,95) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,96) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,97) = 6.70e-11_r8 * exp_fac(:,:)
      rate(:,:,98) = 4.70e-11_r8 * exp_fac(:,:)
      rate(:,:,116) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,117) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,118) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,171) = 2.70e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,120) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,152) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,191) = 1.70e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,121) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,200) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,122) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,124) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,174) = 3.00e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170._r8 * itemp(:,:) )
      rate(:,:,129) = 1.50e-11_r8 * exp_fac(:,:)
      rate(:,:,164) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,134) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,136) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,137) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,138) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,139) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,157) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,199) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,141) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,142) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,206) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,150) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,151) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,155) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,158) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,160) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,161) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,162) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,163) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,165) = 4.10e-11_r8 * exp( -450._r8 * itemp(:,:) )
      rate(:,:,166) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,167) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,168) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      rate(:,:,169) = 7.40e-12_r8 * exp( 270._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,170) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,190) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,198) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,172) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,197) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,175) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,176) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,179) = 2.60e-12_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,180) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,182) = 2.50e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,183) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,184) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,187) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,189) = 1.70e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -330._r8 * itemp(:,:) )
      rate(:,:,185) = 1.20e-12_r8 * exp_fac(:,:)
      rate(:,:,221) = 1.30E-12_r8 * exp_fac(:,:)
      rate(:,:,186) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,188) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,192) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,193) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,196) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,195) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,201) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,202) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,203) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,204) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,205) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,207) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,208) = 2.10E-11_r8 * exp( -2200.0_r8 * itemp(:,:) )
      rate(:,:,209) = 1.10E-13_r8 * exp( -1200.0_r8 * itemp(:,:) )
      rate(:,:,213) = 2.70E-11_r8 * exp( 335._r8 * itemp(:,:) )
      rate(:,:,214) = 1.25E-13_r8 * exp( -2190.0_r8 * itemp(:,:) )
      rate(:,:,215) = 3.40E-12_r8 * exp( -1100.0_r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,119), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,123), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,125), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,127), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,133), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,143), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( rate(1,1,145), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,154), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,173), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 2.0e-12_r8 * itemp(:,:)**2.4_r8
      call jpl( rate(1,1,177), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,194), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 3.3E-31_r8 * itemp(:,:)**4.3_r8
      kinf(:,:) = 1.60E-12_r8
      call jpl( rate(1,1,220), m, 0.6_r8, ko, kinf, n )

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

      rate(:,:kbot,83) = 8.00e-14_r8
      rate(:,:kbot,84) = 3.90e-17_r8
      rate(:,:kbot,88) = 8.50e-2_r8
      rate(:,:kbot,89) = 1.30e-16_r8
      rate(:,:kbot,91) = 1.00e-20_r8
      rate(:,:kbot,92) = 2.58e-04_r8
      rate(:,:kbot,114) = 5.00e-12_r8
      rate(:,:kbot,115) = 7.00e-13_r8
      rate(:,:kbot,148) = 6.90e-12_r8
      rate(:,:kbot,244) = 1.e-10_r8
      rate(:,:kbot,245) = 4.4e-10_r8
      rate(:,:kbot,246) = 4.e-10_r8
      rate(:,:kbot,247) = 2.e-10_r8
      rate(:,:kbot,248) = 1.e-12_r8
      rate(:,:kbot,249) = 6.e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,81) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,85) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,86) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,90) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,93) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,94) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,95) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,116) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,117) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,120) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,152) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,121) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,122) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,146) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,150) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:kbot,151) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,157) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,158) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)







      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,145) = wrk(:,:)






      end subroutine setrxt_hrates

      end module mo_setrxt
