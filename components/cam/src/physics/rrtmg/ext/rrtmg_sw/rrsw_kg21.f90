      module rrsw_kg21

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng21

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 21
! band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
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

      integer, parameter :: no21 = 16

      real(kind=r8) :: kao(9,5,13,no21)
      real(kind=r8) :: kbo(5,5,13:59,no21)
      real(kind=r8) :: selfrefo(10,no21), forrefo(4,no21)
      real(kind=r8) :: sfluxrefo(no21,9)

      integer :: layreffr
      real(kind=r8) :: rayl, strrat

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 21
! band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
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

      real(kind=r8) :: ka(9,5,13,ng21), absa(585,ng21)
      real(kind=r8) :: kb(5,5,13:59,ng21), absb(1175,ng21)
      real(kind=r8) :: selfref(10,ng21), forref(4,ng21)
      real(kind=r8) :: sfluxref(ng21,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

      end module rrsw_kg21

