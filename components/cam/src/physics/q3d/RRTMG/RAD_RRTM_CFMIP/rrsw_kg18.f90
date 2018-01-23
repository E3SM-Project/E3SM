      module rrsw_kg18

      use parkind ,only : im => kind_im, rb => kind_rb
      use parrrsw, only : ng18

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 18
! band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real
!sfluxrefo: real     
!-----------------------------------------------------------------

      integer(kind=im), parameter :: no18 = 16

      real(kind=rb) :: kao(9,5,13,no18)
      real(kind=rb) :: kbo(5,13:59,no18)
      real(kind=rb) :: selfrefo(10,no18), forrefo(3,no18)
      real(kind=rb) :: sfluxrefo(no18,9)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 18
! band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! kb      : real     
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
! sfluxref: real     
!-----------------------------------------------------------------

      real(kind=rb) :: ka(9,5,13,ng18), absa(585,ng18)
      real(kind=rb) :: kb(5,13:59,ng18), absb(235,ng18)
      real(kind=rb) :: selfref(10,ng18), forref(3,ng18)
      real(kind=rb) :: sfluxref(ng18,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg18

