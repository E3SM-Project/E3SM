
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

      rate(:,:,78) = 8.00e-14_r8
      rate(:,:,79) = 3.90e-17_r8
      rate(:,:,82) = 4.20e-13_r8
      rate(:,:,83) = 8.50e-2_r8
      rate(:,:,84) = 1.30e-16_r8
      rate(:,:,86) = 1.00e-20_r8
      rate(:,:,87) = 2.58e-04_r8
      rate(:,:,94) = 1.20e-10_r8
      rate(:,:,95) = 1.70e-10_r8
      rate(:,:,96) = 1.20e-10_r8
      rate(:,:,97) = 1.50e-10_r8
      rate(:,:,98) = 7.20e-11_r8
      rate(:,:,99) = 2.84e-10_r8
      rate(:,:,100) = 1.80e-10_r8
      rate(:,:,101) = 9.60e-11_r8
      rate(:,:,102) = 4.10e-11_r8
      rate(:,:,103) = 1.125e-10_r8
      rate(:,:,104) = 3.00e-11_r8
      rate(:,:,105) = 7.50e-12_r8
      rate(:,:,106) = 1.10e-10_r8
      rate(:,:,107) = 1.50e-10_r8
      rate(:,:,108) = 1.50e-10_r8
      rate(:,:,109) = 5.00e-12_r8
      rate(:,:,110) = 7.00e-13_r8
      rate(:,:,125) = 1.00e-11_r8
      rate(:,:,126) = 2.20e-11_r8
      rate(:,:,127) = 3.50e-12_r8
      rate(:,:,135) = 5.80e-16_r8
      rate(:,:,142) = 7.20e-11_r8
      rate(:,:,143) = 6.90e-12_r8
      rate(:,:,144) = 1.60e-12_r8
      rate(:,:,148) = 1.80e-12_r8
      rate(:,:,151) = 1.80e-12_r8
      rate(:,:,176) = 1.70e-13_r8
      rate(:,:,231) = 9.e-10_r8
      rate(:,:,232) = 1.e-10_r8
      rate(:,:,233) = 4.4e-10_r8
      rate(:,:,234) = 4.e-10_r8
      rate(:,:,235) = 2.e-10_r8
      rate(:,:,236) = 1.e-12_r8
      rate(:,:,237) = 6.e-11_r8
      rate(:,:,238) = 5.e-16_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,76) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,80) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,81) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,85) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,88) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,89) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,90) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,91) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,92) = 6.70e-11_r8 * exp_fac(:,:)
      rate(:,:,93) = 4.70e-11_r8 * exp_fac(:,:)
      rate(:,:,111) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,112) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,113) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,166) = 2.70e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,115) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,147) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 1.70e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,116) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,195) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,117) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,119) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,169) = 3.00e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 170._r8 * itemp(:,:) )
      rate(:,:,124) = 1.50e-11_r8 * exp_fac(:,:)
      rate(:,:,159) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,129) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,131) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,132) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,133) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,134) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,152) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,194) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,136) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,137) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,209) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,141) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,145) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:,146) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,150) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,153) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,155) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,156) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,157) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,158) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,160) = 4.10e-11_r8 * exp( -450._r8 * itemp(:,:) )
      rate(:,:,161) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,162) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,163) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      rate(:,:,164) = 7.40e-12_r8 * exp( 270._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,165) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,185) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,193) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,167) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,192) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,170) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,171) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,174) = 2.60e-12_r8 * exp( -350._r8 * itemp(:,:) )
      rate(:,:,175) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,177) = 2.50e-12_r8 * exp( -130._r8 * itemp(:,:) )
      rate(:,:,178) = 3.00e-12_r8 * exp( -500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,179) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,182) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,184) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,180) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,181) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,183) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,187) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,188) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,191) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,196) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,197) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,202) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,204) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,206) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,207) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,208) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,210) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,114), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,118), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,120), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,122), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,128), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,138), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( rate(1,1,140), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,149), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,168), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 2.0e-12_r8 * itemp(:,:)**2.4_r8
      call jpl( rate(1,1,172), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,189), m, 0.6_r8, ko, kinf, n )

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

      rate(:,:kbot,78) = 8.00e-14_r8
      rate(:,:kbot,79) = 3.90e-17_r8
      rate(:,:kbot,83) = 8.50e-2_r8
      rate(:,:kbot,84) = 1.30e-16_r8
      rate(:,:kbot,86) = 1.00e-20_r8
      rate(:,:kbot,87) = 2.58e-04_r8
      rate(:,:kbot,109) = 5.00e-12_r8
      rate(:,:kbot,110) = 7.00e-13_r8
      rate(:,:kbot,143) = 6.90e-12_r8
      rate(:,:kbot,232) = 1.e-10_r8
      rate(:,:kbot,233) = 4.4e-10_r8
      rate(:,:kbot,234) = 4.e-10_r8
      rate(:,:kbot,235) = 2.e-10_r8
      rate(:,:kbot,236) = 1.e-12_r8
      rate(:,:kbot,237) = 6.e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,76) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,80) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,81) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,85) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,88) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,89) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,90) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,111) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,112) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,115) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,147) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,116) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,117) = 5.20e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,141) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,145) = 2.20e-11_r8 * exp( 120._r8 * itemp(:,:) )
      rate(:,:kbot,146) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,152) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,153) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)







      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 4.7e-11_r8 * itemp(:,:)**0.2_r8
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,140) = wrk(:,:)





      end subroutine setrxt_hrates

      end module mo_setrxt
