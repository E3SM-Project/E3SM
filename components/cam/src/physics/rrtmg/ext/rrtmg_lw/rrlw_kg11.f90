      module rrlw_kg11

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
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
! kao_mo2 : real     
! kbo_mo2 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no11 = 16

      real(kind=r8) , dimension(no11) :: fracrefao
      real(kind=r8) , dimension(no11) :: fracrefbo

      real(kind=r8) :: kao(5,13,no11)
      real(kind=r8) :: kbo(5,13:59,no11)
      real(kind=r8) :: kao_mo2(19,no11)
      real(kind=r8) :: kbo_mo2(19,no11)
      real(kind=r8) :: selfrefo(10,no11)
      real(kind=r8) :: forrefo(4,no11)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
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
! ka_mo2  : real     
! kb_mo2  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer, parameter :: ng11 = 8

      real(kind=r8) , dimension(ng11) :: fracrefa
      real(kind=r8) , dimension(ng11) :: fracrefb

      real(kind=r8) :: ka(5,13,ng11)   , absa(65,ng11)
      real(kind=r8) :: kb(5,13:59,ng11), absb(235,ng11)
      real(kind=r8) :: ka_mo2(19,ng11)
      real(kind=r8) :: kb_mo2(19,ng11)
      real(kind=r8) :: selfref(10,ng11)
      real(kind=r8) :: forref(4,ng11)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg11
