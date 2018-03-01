
      module parrrtm

!      use parkind ,only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw main parameters
!
! Initial version:  JJMorcrette, ECMWF, Jul 1998
! Revised: MJIacono, AER, Jun 2006
! Revised: MJIacono, AER, Aug 2007
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! mxlay  :  integer: maximum number of layers
! mg     :  integer: number of original g-intervals per spectral band
! nbndlw :  integer: number of spectral bands
! maxxsec:  integer: maximum number of cross-section molecules
!                    (e.g. cfcs)
! maxinpx:  integer: 
! ngptlw :  integer: total number of reduced g-intervals for rrtmg_lw
! ngNN   :  integer: number of reduced g-intervals per spectral band
! ngsNN  :  integer: cumulative number of g-intervals per band
!------------------------------------------------------------------

      integer, parameter :: mxlay  = 203
      integer, parameter :: mg     = 16
      integer, parameter :: nbndlw = 16
      integer, parameter :: maxxsec= 4
      integer, parameter :: mxmol  = 38
      integer, parameter :: maxinpx= 38
      integer, parameter :: nmol   = 7
! Use for 140 g-point model 
      integer, parameter :: ngptlw = 140
! Use for 256 g-point model 
!      integer, parameter :: ngptlw = 256

! Use for 140 g-point model
      integer, parameter :: ng1  = 10
      integer, parameter :: ng2  = 12
      integer, parameter :: ng3  = 16
      integer, parameter :: ng4  = 14
      integer, parameter :: ng5  = 16
      integer, parameter :: ng6  = 8
      integer, parameter :: ng7  = 12
      integer, parameter :: ng8  = 8
      integer, parameter :: ng9  = 12
      integer, parameter :: ng10 = 6
      integer, parameter :: ng11 = 8
      integer, parameter :: ng12 = 8
      integer, parameter :: ng13 = 4
      integer, parameter :: ng14 = 2
      integer, parameter :: ng15 = 2
      integer, parameter :: ng16 = 2

      integer, parameter :: ngs1  = 10
      integer, parameter :: ngs2  = 22
      integer, parameter :: ngs3  = 38
      integer, parameter :: ngs4  = 52
      integer, parameter :: ngs5  = 68
      integer, parameter :: ngs6  = 76
      integer, parameter :: ngs7  = 88
      integer, parameter :: ngs8  = 96
      integer, parameter :: ngs9  = 108
      integer, parameter :: ngs10 = 114
      integer, parameter :: ngs11 = 122
      integer, parameter :: ngs12 = 130
      integer, parameter :: ngs13 = 134
      integer, parameter :: ngs14 = 136
      integer, parameter :: ngs15 = 138

! Use for 256 g-point model
!      integer, parameter :: ng1  = 16
!      integer, parameter :: ng2  = 16
!      integer, parameter :: ng3  = 16
!      integer, parameter :: ng4  = 16
!      integer, parameter :: ng5  = 16
!      integer, parameter :: ng6  = 16
!      integer, parameter :: ng7  = 16
!      integer, parameter :: ng8  = 16
!      integer, parameter :: ng9  = 16
!      integer, parameter :: ng10 = 16
!      integer, parameter :: ng11 = 16
!      integer, parameter :: ng12 = 16
!      integer, parameter :: ng13 = 16
!      integer, parameter :: ng14 = 16
!      integer, parameter :: ng15 = 16
!      integer, parameter :: ng16 = 16

!      integer, parameter :: ngs1  = 16
!      integer, parameter :: ngs2  = 32
!      integer, parameter :: ngs3  = 48
!      integer, parameter :: ngs4  = 64
!      integer, parameter :: ngs5  = 80
!      integer, parameter :: ngs6  = 96
!      integer, parameter :: ngs7  = 112
!      integer, parameter :: ngs8  = 128
!      integer, parameter :: ngs9  = 144
!      integer, parameter :: ngs10 = 160
!      integer, parameter :: ngs11 = 176
!      integer, parameter :: ngs12 = 192
!      integer, parameter :: ngs13 = 208
!      integer, parameter :: ngs14 = 224
!      integer, parameter :: ngs15 = 240
!      integer, parameter :: ngs16 = 256

      end module parrrtm
