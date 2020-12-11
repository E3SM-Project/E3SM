      module rrsw_kg23

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb
      use parrrsw, only : ng23

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 23
! band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
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

      integer, parameter :: no23 = 16

      real(kind=r8) :: kao(5,13,no23)
      real(kind=r8) :: selfrefo(10,no23), forrefo(3,no23)
      real(kind=r8) :: sfluxrefo(no23)
      real(kind=r8) :: raylo(no23)

      integer :: layreffr
      real(kind=r8) :: givfac

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 23
! band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
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

      real(kind=r8) :: ka(5,13,ng23), absa(65,ng23)
      real(kind=r8) :: selfref(10,ng23), forref(3,ng23)
      real(kind=r8) :: sfluxref(ng23), rayl(ng23)

      equivalence (ka(1,1,1),absa(1,1))

      end module rrsw_kg23

