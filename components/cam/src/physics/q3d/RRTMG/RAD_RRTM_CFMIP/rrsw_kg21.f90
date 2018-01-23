      module rrsw_kg21

      use parkind ,only : im => kind_im, rb => kind_rb
      use parrrsw, only : ng21

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 21
! band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
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

      integer(kind=im), parameter :: no21 = 16

      real(kind=rb) :: kao(9,5,13,no21)
      real(kind=rb) :: kbo(5,5,13:59,no21)
      real(kind=rb) :: selfrefo(10,no21), forrefo(4,no21)
      real(kind=rb) :: sfluxrefo(no21,9)

      real(kind=rb) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 21
! band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
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

      real(kind=rb) :: ka(9,5,13,ng21), absa(585,ng21)
      real(kind=rb) :: kb(5,5,13:59,ng21), absb(1175,ng21)
      real(kind=rb) :: selfref(10,ng21), forref(4,ng21)
      real(kind=rb) :: sfluxref(ng21,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

      end module rrsw_kg21

