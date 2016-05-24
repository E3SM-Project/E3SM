      module rrlw_cld

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw cloud property coefficients

! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! abscld1:  real   : 
! absice0:  real   : 
! absice1:  real   : 
! absice2:  real   : 
! absice3:  real   : 
! absliq0:  real   : 
! absliq1:  real   : 
!------------------------------------------------------------------

      real(kind=r8) :: abscld1
      real(kind=r8) , dimension(2) :: absice0
      real(kind=r8) , dimension(2,5) :: absice1
      real(kind=r8) , dimension(43,16) :: absice2
      real(kind=r8) , dimension(46,16) :: absice3
      real(kind=r8) :: absliq0
      real(kind=r8) , dimension(58,16) :: absliq1

      end module rrlw_cld

