      module rrsw_kg22

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng22

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 22
! band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
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

      integer, parameter :: no22 = 16

      real(kind=r8) :: kao(9,5,13,no22)
      real(kind=r8) :: kbo(5,13:59,no22)
      real(kind=r8) :: selfrefo(10,no22), forrefo(3,no22)
      real(kind=r8) :: sfluxrefo(no22,9)

      integer :: layreffr
      real(kind=r8) :: rayl, strrat

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 22
! band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
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

      real(kind=r8) :: ka(9,5,13,ng22), absa(585,ng22)
      real(kind=r8) :: kb(5,13:59,ng22), absb(235,ng22)
      real(kind=r8) :: selfref(10,ng22), forref(3,ng22)
      real(kind=r8) :: sfluxref(ng22,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg22

