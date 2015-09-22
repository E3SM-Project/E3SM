      module rrsw_kg28

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng28

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 28
! band 28: 38000-50000 cm-1 (low - o3, o2; high - o3, o2)
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
!-----------------------------------------------------------------

      integer, parameter :: no28 = 16

      real(kind=r8) :: kao(9,5,13,no28)
      real(kind=r8) :: kbo(5,5,13:59,no28)
      real(kind=r8) :: sfluxrefo(no28,5)

      integer :: layreffr
      real(kind=r8) :: rayl, strrat

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 28
! band 28: 38000-50000 cm-1 (low - o3, o2; high - o3, o2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! kb      : real     
! sfluxref: real     
!-----------------------------------------------------------------

      real(kind=r8) :: ka(9,5,13,ng28), absa(585,ng28)
      real(kind=r8) :: kb(5,5,13:59,ng28), absb(1175,ng28)
      real(kind=r8) :: sfluxref(ng28,5)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

      end module rrsw_kg28

