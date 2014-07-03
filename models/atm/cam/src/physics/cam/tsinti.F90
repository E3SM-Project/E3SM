
subroutine tsinti (tmeltx, latvapx, rairx, stebolx, laticex)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize surface temperature calculation constants
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: L. Buja
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_control_mod, only: rair, stebol, snwedp, snwedp, latice, tmelt, latvap
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8) tmeltx         ! Melting temperature of snow and ice
   real(r8) latvapx        ! Latent heat of vaporization
   real(r8) rairx          ! Gas constant for dry air
   real(r8) stebolx        ! Stefan-Boltzmann constant
   real(r8) laticex        ! latent heat of fusion
!
!-----------------------------------------------------------------------
!
   latice = laticex    ! Latent heat of fusion at 0'C = 3.336e5 J/Kg
   tmelt  = tmeltx
   latvap = latvapx
   rair   = rairx
   stebol = stebolx
   snwedp = 10.0_r8       ! 10:1 Snow:water equivalent depth factor
!
   return
end subroutine tsinti

