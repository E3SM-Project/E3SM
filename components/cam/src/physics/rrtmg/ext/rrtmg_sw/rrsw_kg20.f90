      module rrsw_kg20

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng20

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
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
! absch4o : real     
!-----------------------------------------------------------------

      integer, parameter :: no20 = 16

      real(kind=r8) :: kao(5,13,no20)
      real(kind=r8) :: kbo(5,13:59,no20)
      real(kind=r8) :: selfrefo(10,no20), forrefo(4,no20)
      real(kind=r8) :: sfluxrefo(no20)
      real(kind=r8) :: absch4o(no20)

      integer :: layreffr
      real(kind=r8) :: rayl 

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 20
! band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
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
! absch4  : real     
!-----------------------------------------------------------------

      real(kind=r8) :: ka(5,13,ng20), absa(65,ng20)
      real(kind=r8) :: kb(5,13:59,ng20), absb(235,ng20)
      real(kind=r8) :: selfref(10,ng20), forref(4,ng20)
      real(kind=r8) :: sfluxref(ng20)
      real(kind=r8) :: absch4(ng20)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg20

