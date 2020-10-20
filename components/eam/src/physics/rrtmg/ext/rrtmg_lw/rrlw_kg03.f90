      module rrlw_kg03

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
!fracrefbo: real
! kao     : real     
! kbo     : real     
! kao_mn2o: real     
! kbo_mn2o: real     
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer, parameter :: no3  = 16

      real(kind=r8) :: fracrefao(no3,10) ,fracrefbo(no3,5)
      real(kind=r8) :: kao(9,5,13,no3)
      real(kind=r8) :: kbo(5,5,13:59,no3)
      real(kind=r8) :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
      real(kind=r8) :: selfrefo(10,no3)
      real(kind=r8) :: forrefo(4,no3)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
!fracrefb : real
! ka      : real     
! kb      : real     
! ka_mn2o : real     
! kb_mn2o : real     
! selfref : real     
! forref  : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer, parameter :: ng3  = 16

      real(kind=r8) :: fracrefa(ng3,10) ,fracrefb(ng3,5)
      real(kind=r8) :: ka(9,5,13,ng3)  ,absa(585,ng3)
      real(kind=r8) :: kb(5,5,13:59,ng3),absb(1175,ng3)
      real(kind=r8) :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
      real(kind=r8) :: selfref(10,ng3)
      real(kind=r8) :: forref(4,ng3)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      end module rrlw_kg03


