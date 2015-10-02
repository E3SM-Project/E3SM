
      module parrrsw

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind ,only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw main parameters
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
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

! Settings for single column mode.
! For GCM use, set nlon to number of longitudes, and
! mxlay to number of model layers
      integer, parameter :: mxlay  = 203    !jplay, klev
      integer, parameter :: mg     = 16     !jpg
      integer, parameter :: nbndsw = 14     !jpsw, ksw
      integer, parameter :: naerec  = 6     !jpaer
      integer, parameter :: mxmol  = 38
      integer, parameter :: nstr   = 2
      integer, parameter :: nmol   = 7
! Use for 112 g-point model   
      integer, parameter :: ngptsw = 112    !jpgpt
! Use for 224 g-point model   
!      integer, parameter :: ngptsw = 224   !jpgpt

! may need to rename these - from v2.6
      integer, parameter :: jpband   = 29
      integer, parameter :: jpb1     = 16   !istart
      integer, parameter :: jpb2     = 29   !iend

      integer, parameter :: jmcmu    = 32
      integer, parameter :: jmumu    = 32
      integer, parameter :: jmphi    = 3
      integer, parameter :: jmxang   = 4
      integer, parameter :: jmxstr   = 16
! ^

! Use for 112 g-point model   
      integer, parameter :: ng16 = 6
      integer, parameter :: ng17 = 12
      integer, parameter :: ng18 = 8
      integer, parameter :: ng19 = 8
      integer, parameter :: ng20 = 10
      integer, parameter :: ng21 = 10
      integer, parameter :: ng22 = 2
      integer, parameter :: ng23 = 10
      integer, parameter :: ng24 = 8
      integer, parameter :: ng25 = 6
      integer, parameter :: ng26 = 6
      integer, parameter :: ng27 = 8
      integer, parameter :: ng28 = 6
      integer, parameter :: ng29 = 12

      integer, parameter :: ngs16 = 6
      integer, parameter :: ngs17 = 18
      integer, parameter :: ngs18 = 26
      integer, parameter :: ngs19 = 34
      integer, parameter :: ngs20 = 44
      integer, parameter :: ngs21 = 54
      integer, parameter :: ngs22 = 56
      integer, parameter :: ngs23 = 66
      integer, parameter :: ngs24 = 74
      integer, parameter :: ngs25 = 80
      integer, parameter :: ngs26 = 86
      integer, parameter :: ngs27 = 94
      integer, parameter :: ngs28 = 100
      integer, parameter :: ngs29 = 112

! Use for 224 g-point model   
!      integer, parameter :: ng16 = 16
!      integer, parameter :: ng17 = 16
!      integer, parameter :: ng18 = 16
!      integer, parameter :: ng19 = 16
!      integer, parameter :: ng20 = 16
!      integer, parameter :: ng21 = 16
!      integer, parameter :: ng22 = 16
!      integer, parameter :: ng23 = 16
!      integer, parameter :: ng24 = 16
!      integer, parameter :: ng25 = 16
!      integer, parameter :: ng26 = 16
!      integer, parameter :: ng27 = 16
!      integer, parameter :: ng28 = 16
!      integer, parameter :: ng29 = 16

!      integer, parameter :: ngs16 = 16
!      integer, parameter :: ngs17 = 32
!      integer, parameter :: ngs18 = 48
!      integer, parameter :: ngs19 = 64
!      integer, parameter :: ngs20 = 80
!      integer, parameter :: ngs21 = 96
!      integer, parameter :: ngs22 = 112
!      integer, parameter :: ngs23 = 128
!      integer, parameter :: ngs24 = 144
!      integer, parameter :: ngs25 = 160
!      integer, parameter :: ngs26 = 176
!      integer, parameter :: ngs27 = 192
!      integer, parameter :: ngs28 = 208
!      integer, parameter :: ngs29 = 224

! Source function solar constant
      real(kind=r8), parameter :: rrsw_scon = 1.36822e+03     ! W/m2
 
      end module parrrsw


