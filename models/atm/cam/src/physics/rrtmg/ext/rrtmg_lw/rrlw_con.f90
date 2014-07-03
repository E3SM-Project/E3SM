      module rrlw_con

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw constants

! Initial version: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! fluxfac:  real   : radiance to flux conversion factor 
! heatfac:  real   : flux to heating rate conversion factor
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
! grav   :  real   : acceleration of gravity (m/s2)
! planck :  real   : planck constant
! boltz  :  real   : boltzman constant
! clight :  real   : speed of light
! avogad :  real   : avogadro's constant 
! alosmt :  real   : 
! gascon :  real   : gas constant
! radcn1 :  real   : 
! radcn2 :  real   : 
!------------------------------------------------------------------

      real(kind=r8) :: fluxfac, heatfac
      real(kind=r8) :: oneminus, pi, grav
      real(kind=r8) :: planck, boltz, clight
      real(kind=r8) :: avogad, alosmt, gascon
      real(kind=r8) :: radcn1, radcn2

      end module rrlw_con

