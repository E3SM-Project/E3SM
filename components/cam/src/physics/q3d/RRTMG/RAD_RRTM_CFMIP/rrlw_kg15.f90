      module rrlw_kg15

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
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

      integer(kind=im), parameter :: no15 = 16

      real(kind=rb) :: fracrefao(no15,9)
      real(kind=rb) :: kao(9,5,13,no15)
      real(kind=rb) :: kao_mn2(9,19,no15)
      real(kind=rb) :: selfrefo(10,no15)
      real(kind=rb) :: forrefo(4,no15)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
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

      integer(kind=im), parameter :: ng15 = 2

      real(kind=rb) :: fracrefa(ng15,9)
      real(kind=rb) :: ka(9,5,13,ng15) ,absa(585,ng15)
      real(kind=rb) :: ka_mn2(9,19,ng15)
      real(kind=rb) :: selfref(10,ng15)
      real(kind=rb) :: forref(4,ng15)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg15
