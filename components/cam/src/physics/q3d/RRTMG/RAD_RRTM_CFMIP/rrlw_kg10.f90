      module rrlw_kg10

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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

      integer(kind=im), parameter :: no10 = 16

      real(kind=rb) , dimension(no10) :: fracrefao
      real(kind=rb) , dimension(no10) :: fracrefbo

      real(kind=rb) :: kao(5,13,no10)
      real(kind=rb) :: kbo(5,13:59,no10)
      real(kind=rb) :: selfrefo(10,no10)
      real(kind=rb) :: forrefo(4,no10)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=im), parameter :: ng10 = 6

      real(kind=rb) , dimension(ng10) :: fracrefa
      real(kind=rb) , dimension(ng10) :: fracrefb

      real(kind=rb) :: ka(5,13,ng10)   , absa(65,ng10)
      real(kind=rb) :: kb(5,13:59,ng10), absb(235,ng10)
      real(kind=rb) :: selfref(10,ng10)
      real(kind=rb) :: forref(4,ng10)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg10
