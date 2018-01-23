      module rrlw_kg09

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! kao_mn2o: real     
! kbo_mn2o: real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no9  = 16

      real(kind=rb) , dimension(no9) :: fracrefbo

      real(kind=rb) :: fracrefao(no9,9)
      real(kind=rb) :: kao(9,5,13,no9)
      real(kind=rb) :: kbo(5,13:59,no9)
      real(kind=rb) :: kao_mn2o(9,19,no9)
      real(kind=rb) :: kbo_mn2o(19,no9)
      real(kind=rb) :: selfrefo(10,no9)
      real(kind=rb) :: forrefo(4,no9)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! ka_mn2o : real     
! kb_mn2o : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng9  = 12

      real(kind=rb) , dimension(ng9) :: fracrefb
      real(kind=rb) :: fracrefa(ng9,9)
      real(kind=rb) :: ka(9,5,13,ng9) ,absa(585,ng9)
      real(kind=rb) :: kb(5,13:59,ng9) ,absb(235,ng9)
      real(kind=rb) :: ka_mn2o(9,19,ng9)
      real(kind=rb) :: kb_mn2o(19,ng9)
      real(kind=rb) :: selfref(10,ng9)
      real(kind=rb) :: forref(4,ng9)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg09
