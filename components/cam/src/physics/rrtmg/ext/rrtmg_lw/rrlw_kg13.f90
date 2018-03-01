      module rrlw_kg13

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mco2: real     
! kao_mco : real     
! kbo_mo3 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no13 = 16

      real(kind=r8) , dimension(no13) :: fracrefbo

      real(kind=r8) :: fracrefao(no13,9)
      real(kind=r8) :: kao(9,5,13,no13)
      real(kind=r8) :: kao_mco2(9,19,no13)
      real(kind=r8) :: kao_mco(9,19,no13)
      real(kind=r8) :: kbo_mo3(19,no13)
      real(kind=r8) :: selfrefo(10,no13)
      real(kind=r8) :: forrefo(4,no13)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mco2 : real     
! ka_mco  : real     
! kb_mo3  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer, parameter :: ng13 = 4

      real(kind=r8) , dimension(ng13) :: fracrefb

      real(kind=r8) :: fracrefa(ng13,9)
      real(kind=r8) :: ka(9,5,13,ng13) ,absa(585,ng13)
      real(kind=r8) :: ka_mco2(9,19,ng13)
      real(kind=r8) :: ka_mco(9,19,ng13)
      real(kind=r8) :: kb_mo3(19,ng13)
      real(kind=r8) :: selfref(10,ng13)
      real(kind=r8) :: forref(4,ng13)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg13
