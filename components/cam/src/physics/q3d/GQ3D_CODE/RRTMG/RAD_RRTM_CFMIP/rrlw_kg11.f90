      module rrlw_kg11

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
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

      integer(kind=im), parameter :: no11 = 16

      real(kind=rb) , dimension(no11) :: fracrefao
      real(kind=rb) , dimension(no11) :: fracrefbo

      real(kind=rb) :: kao(5,13,no11)
      real(kind=rb) :: kbo(5,13:59,no11)
      real(kind=rb) :: kao_mo2(19,no11)
      real(kind=rb) :: kbo_mo2(19,no11)
      real(kind=rb) :: selfrefo(10,no11)
      real(kind=rb) :: forrefo(4,no11)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
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

      integer(kind=im), parameter :: ng11 = 8

      real(kind=rb) , dimension(ng11) :: fracrefa
      real(kind=rb) , dimension(ng11) :: fracrefb

      real(kind=rb) :: ka(5,13,ng11)   , absa(65,ng11)
      real(kind=rb) :: kb(5,13:59,ng11), absb(235,ng11)
      real(kind=rb) :: ka_mo2(19,ng11)
      real(kind=rb) :: kb_mo2(19,ng11)
      real(kind=rb) :: selfref(10,ng11)
      real(kind=rb) :: forref(4,ng11)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg11
