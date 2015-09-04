      module rrlw_kg15

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, r8

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mn2 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer, parameter :: no15 = 16

      real(kind=r8) :: fracrefao(no15,9)
      real(kind=r8) :: kao(9,5,13,no15)
      real(kind=r8) :: kao_mn2(9,19,no15)
      real(kind=r8) :: selfrefo(10,no15)
      real(kind=r8) :: forrefo(4,no15)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mn2  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer, parameter :: ng15 = 2

      real(kind=r8) :: fracrefa(ng15,9)
      real(kind=r8) :: ka(9,5,13,ng15) ,absa(585,ng15)
      real(kind=r8) :: ka_mn2(9,19,ng15)
      real(kind=r8) :: selfref(10,ng15)
      real(kind=r8) :: forref(4,ng15)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg15
