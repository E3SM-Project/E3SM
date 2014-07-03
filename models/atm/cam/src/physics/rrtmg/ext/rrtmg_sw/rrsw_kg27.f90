      module rrsw_kg27

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng27

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
! kbo     : real     
!sfluxrefo: real     
! raylo   : real     
!-----------------------------------------------------------------

      integer, parameter :: no27 = 16

      real(kind=r8) :: kao(5,13,no27)
      real(kind=r8) :: kbo(5,13:59,no27)
      real(kind=r8) :: sfluxrefo(no27)
      real(kind=r8) :: raylo(no27)

      integer :: layreffr
      real(kind=r8) :: scalekur

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
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
! sfluxref: real     
! rayl    : real     
!-----------------------------------------------------------------

      real(kind=r8) :: ka(5,13,ng27), absa(65,ng27)
      real(kind=r8) :: kb(5,13:59,ng27), absb(235,ng27)
      real(kind=r8) :: sfluxref(ng27)
      real(kind=r8) :: rayl(ng27)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg27

