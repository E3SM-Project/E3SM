      module rrsw_kg16

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng16

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      integer, parameter :: no16 = 16

      real(kind=r8) :: kao(9,5,13,no16)
      real(kind=r8) :: kbo(5,13:59,no16)
      real(kind=r8) :: selfrefo(10,no16), forrefo(3,no16)
      real(kind=r8) :: sfluxrefo(no16)

      integer :: layreffr
      real(kind=r8) :: rayl, strrat1

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      real(kind=r8) :: ka(9,5,13,ng16) , absa(585,ng16)
      real(kind=r8) :: kb(5,13:59,ng16), absb(235,ng16)
      real(kind=r8) :: selfref(10,ng16), forref(3,ng16)
      real(kind=r8) :: sfluxref(ng16)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg16

