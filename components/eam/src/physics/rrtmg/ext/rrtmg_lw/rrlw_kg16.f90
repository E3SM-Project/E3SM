      module rrlw_kg16

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no16 = 16

      real(kind=r8) , dimension(no16) :: fracrefbo

      real(kind=r8) :: fracrefao(no16,9)
      real(kind=r8) :: kao(9,5,13,no16)
      real(kind=r8) :: kbo(5,13:59,no16)
      real(kind=r8) :: selfrefo(10,no16)
      real(kind=r8) :: forrefo(4,no16)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! kb      : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer, parameter :: ng16 = 2

      real(kind=r8) , dimension(ng16) :: fracrefb

      real(kind=r8) :: fracrefa(ng16,9)
      real(kind=r8) :: ka(9,5,13,ng16) ,absa(585,ng16)
      real(kind=r8) :: kb(5,13:59,ng16), absb(235,ng16)
      real(kind=r8) :: selfref(10,ng16)
      real(kind=r8) :: forref(4,ng16)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrlw_kg16

