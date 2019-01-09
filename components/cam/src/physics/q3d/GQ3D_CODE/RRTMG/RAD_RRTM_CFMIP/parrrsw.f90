
      module parrrsw

      use parkind ,only : im => kind_im, rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw main parameters
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! mxlay  :  integer: maximum number of layers
! mg     :  integer: number of original g-intervals per spectral band
! nbndsw :  integer: number of spectral bands
! naerec :  integer: number of aerosols (iaer=6, ecmwf aerosol option)
! ngptsw :  integer: total number of reduced g-intervals for rrtmg_lw
! ngNN   :  integer: number of reduced g-intervals per spectral band
! ngsNN  :  integer: cumulative number of g-intervals per band
!------------------------------------------------------------------

      integer(kind=im), parameter :: mxlay  = 203    !jplay, klev
      integer(kind=im), parameter :: mg     = 16     !jpg
      integer(kind=im), parameter :: nbndsw = 14     !jpsw, ksw
      integer(kind=im), parameter :: naerec  = 6     !jpaer
      integer(kind=im), parameter :: mxmol  = 38
      integer(kind=im), parameter :: nstr   = 2
      integer(kind=im), parameter :: nmol   = 7
! Use for 112 g-point model   
      integer(kind=im), parameter :: ngptsw = 112    !jpgpt
! Use for 224 g-point model   
!      integer(kind=im), parameter :: ngptsw = 224   !jpgpt

! may need to rename these - from v2.6
      integer(kind=im), parameter :: jpband   = 29
      integer(kind=im), parameter :: jpb1     = 16   !istart
      integer(kind=im), parameter :: jpb2     = 29   !iend

      integer(kind=im), parameter :: jmcmu    = 32
      integer(kind=im), parameter :: jmumu    = 32
      integer(kind=im), parameter :: jmphi    = 3
      integer(kind=im), parameter :: jmxang   = 4
      integer(kind=im), parameter :: jmxstr   = 16
! ^

! Use for 112 g-point model   
      integer(kind=im), parameter :: ng16 = 6
      integer(kind=im), parameter :: ng17 = 12
      integer(kind=im), parameter :: ng18 = 8
      integer(kind=im), parameter :: ng19 = 8
      integer(kind=im), parameter :: ng20 = 10
      integer(kind=im), parameter :: ng21 = 10
      integer(kind=im), parameter :: ng22 = 2
      integer(kind=im), parameter :: ng23 = 10
      integer(kind=im), parameter :: ng24 = 8
      integer(kind=im), parameter :: ng25 = 6
      integer(kind=im), parameter :: ng26 = 6
      integer(kind=im), parameter :: ng27 = 8
      integer(kind=im), parameter :: ng28 = 6
      integer(kind=im), parameter :: ng29 = 12

      integer(kind=im), parameter :: ngs16 = 6
      integer(kind=im), parameter :: ngs17 = 18
      integer(kind=im), parameter :: ngs18 = 26
      integer(kind=im), parameter :: ngs19 = 34
      integer(kind=im), parameter :: ngs20 = 44
      integer(kind=im), parameter :: ngs21 = 54
      integer(kind=im), parameter :: ngs22 = 56
      integer(kind=im), parameter :: ngs23 = 66
      integer(kind=im), parameter :: ngs24 = 74
      integer(kind=im), parameter :: ngs25 = 80
      integer(kind=im), parameter :: ngs26 = 86
      integer(kind=im), parameter :: ngs27 = 94
      integer(kind=im), parameter :: ngs28 = 100
      integer(kind=im), parameter :: ngs29 = 112

! Use for 224 g-point model   
!      integer(kind=im), parameter :: ng16 = 16
!      integer(kind=im), parameter :: ng17 = 16
!      integer(kind=im), parameter :: ng18 = 16
!      integer(kind=im), parameter :: ng19 = 16
!      integer(kind=im), parameter :: ng20 = 16
!      integer(kind=im), parameter :: ng21 = 16
!      integer(kind=im), parameter :: ng22 = 16
!      integer(kind=im), parameter :: ng23 = 16
!      integer(kind=im), parameter :: ng24 = 16
!      integer(kind=im), parameter :: ng25 = 16
!      integer(kind=im), parameter :: ng26 = 16
!      integer(kind=im), parameter :: ng27 = 16
!      integer(kind=im), parameter :: ng28 = 16
!      integer(kind=im), parameter :: ng29 = 16

!      integer(kind=im), parameter :: ngs16 = 16
!      integer(kind=im), parameter :: ngs17 = 32
!      integer(kind=im), parameter :: ngs18 = 48
!      integer(kind=im), parameter :: ngs19 = 64
!      integer(kind=im), parameter :: ngs20 = 80
!      integer(kind=im), parameter :: ngs21 = 96
!      integer(kind=im), parameter :: ngs22 = 112
!      integer(kind=im), parameter :: ngs23 = 128
!      integer(kind=im), parameter :: ngs24 = 144
!      integer(kind=im), parameter :: ngs25 = 160
!      integer(kind=im), parameter :: ngs26 = 176
!      integer(kind=im), parameter :: ngs27 = 192
!      integer(kind=im), parameter :: ngs28 = 208
!      integer(kind=im), parameter :: ngs29 = 224

! Source function solar constant
      real(kind=rb), parameter :: rrsw_scon = 1.36822e+03     ! W/m2
 
      end module parrrsw


