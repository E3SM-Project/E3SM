      module rrlw_kg07

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! kao_mco2: real     
! kbo_mco2: real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no7  = 16

      real(kind=r8) , dimension(no7) :: fracrefbo
      real(kind=r8) :: fracrefao(no7,9)
      real(kind=r8) :: kao(9,5,13,no7)
      real(kind=r8) :: kbo(5,13:59,no7)
      real(kind=r8) :: kao_mco2(9,19,no7)
      real(kind=r8) :: kbo_mco2(19,no7)
      real(kind=r8) :: selfrefo(10,no7)
      real(kind=r8) :: forrefo(4,no7)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! ka_mco2 : real     
! kb_mco2 : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer, parameter :: ng7  = 12

      real(kind=r8) , dimension(ng7) :: fracrefb
      real(kind=r8) :: fracrefa(ng7,9)
      real(kind=r8) :: ka(9,5,13,ng7) ,absa(585,ng7)
      real(kind=r8) :: kb(5,13:59,ng7),absb(235,ng7)
      real(kind=r8) :: ka_mco2(9,19,ng7)
      real(kind=r8) :: kb_mco2(19,ng7)
      real(kind=r8) :: selfref(10,ng7)
      real(kind=r8) :: forref(4,ng7)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg07
