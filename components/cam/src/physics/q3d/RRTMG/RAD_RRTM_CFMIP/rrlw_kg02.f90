      module rrlw_kg02

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no2  = 16

      real(kind=rb) :: fracrefao(no2)   , fracrefbo(no2)
      real(kind=rb) :: kao(5,13,no2)
      real(kind=rb) :: kbo(5,13:59,no2)
      real(kind=rb) :: selfrefo(10,no2) , forrefo(4,no2)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
!
! refparam: real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng2  = 12

      real(kind=rb) :: fracrefa(ng2)  , fracrefb(ng2)
      real(kind=rb) :: ka(5,13,ng2)   , absa(65,ng2)
      real(kind=rb) :: kb(5,13:59,ng2), absb(235,ng2)
      real(kind=rb) :: selfref(10,ng2), forref(4,ng2)

      real(kind=rb) :: refparam(13)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg02


