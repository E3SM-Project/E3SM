      module rrlw_kg12

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no12 = 16

      real(kind=r8) :: fracrefao(no12,9)
      real(kind=r8) :: kao(9,5,13,no12)
      real(kind=r8) :: selfrefo(10,no12)
      real(kind=r8) :: forrefo(4,no12)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer, parameter :: ng12 = 8

      real(kind=r8) :: fracrefa(ng12,9)
      real(kind=r8) :: ka(9,5,13,ng12) ,absa(585,ng12)
      real(kind=r8) :: selfref(10,ng12)
      real(kind=r8) :: forref(4,ng12)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg12
