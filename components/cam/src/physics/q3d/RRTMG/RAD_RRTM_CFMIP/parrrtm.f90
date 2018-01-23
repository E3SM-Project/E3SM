
      module parrrtm

      use parkind ,only : im => kind_im

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw main parameters
!
! Initial version:  JJMorcrette, ECMWF, Jul 1998
! Revised: MJIacono, AER, Jun 2006
! Revised: MJIacono, AER, Aug 2007
! Revised: MJIacono, AER, Aug 2008
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

      integer(kind=im), parameter :: mxlay  = 203
      integer(kind=im), parameter :: mg     = 16
      integer(kind=im), parameter :: nbndlw = 16
      integer(kind=im), parameter :: maxxsec= 4
      integer(kind=im), parameter :: mxmol  = 38
      integer(kind=im), parameter :: maxinpx= 38
      integer(kind=im), parameter :: nmol   = 7
! Use for 140 g-point model 
      integer(kind=im), parameter :: ngptlw = 140
! Use for 256 g-point model 
!      integer(kind=im), parameter :: ngptlw = 256

! Use for 140 g-point model
      integer(kind=im), parameter :: ng1  = 10
      integer(kind=im), parameter :: ng2  = 12
      integer(kind=im), parameter :: ng3  = 16
      integer(kind=im), parameter :: ng4  = 14
      integer(kind=im), parameter :: ng5  = 16
      integer(kind=im), parameter :: ng6  = 8
      integer(kind=im), parameter :: ng7  = 12
      integer(kind=im), parameter :: ng8  = 8
      integer(kind=im), parameter :: ng9  = 12
      integer(kind=im), parameter :: ng10 = 6
      integer(kind=im), parameter :: ng11 = 8
      integer(kind=im), parameter :: ng12 = 8
      integer(kind=im), parameter :: ng13 = 4
      integer(kind=im), parameter :: ng14 = 2
      integer(kind=im), parameter :: ng15 = 2
      integer(kind=im), parameter :: ng16 = 2

      integer(kind=im), parameter :: ngs1  = 10
      integer(kind=im), parameter :: ngs2  = 22
      integer(kind=im), parameter :: ngs3  = 38
      integer(kind=im), parameter :: ngs4  = 52
      integer(kind=im), parameter :: ngs5  = 68
      integer(kind=im), parameter :: ngs6  = 76
      integer(kind=im), parameter :: ngs7  = 88
      integer(kind=im), parameter :: ngs8  = 96
      integer(kind=im), parameter :: ngs9  = 108
      integer(kind=im), parameter :: ngs10 = 114
      integer(kind=im), parameter :: ngs11 = 122
      integer(kind=im), parameter :: ngs12 = 130
      integer(kind=im), parameter :: ngs13 = 134
      integer(kind=im), parameter :: ngs14 = 136
      integer(kind=im), parameter :: ngs15 = 138

! Use for 256 g-point model
!      integer(kind=im), parameter :: ng1  = 16
!      integer(kind=im), parameter :: ng2  = 16
!      integer(kind=im), parameter :: ng3  = 16
!      integer(kind=im), parameter :: ng4  = 16
!      integer(kind=im), parameter :: ng5  = 16
!      integer(kind=im), parameter :: ng6  = 16
!      integer(kind=im), parameter :: ng7  = 16
!      integer(kind=im), parameter :: ng8  = 16
!      integer(kind=im), parameter :: ng9  = 16
!      integer(kind=im), parameter :: ng10 = 16
!      integer(kind=im), parameter :: ng11 = 16
!      integer(kind=im), parameter :: ng12 = 16
!      integer(kind=im), parameter :: ng13 = 16
!      integer(kind=im), parameter :: ng14 = 16
!      integer(kind=im), parameter :: ng15 = 16
!      integer(kind=im), parameter :: ng16 = 16

!      integer(kind=im), parameter :: ngs1  = 16
!      integer(kind=im), parameter :: ngs2  = 32
!      integer(kind=im), parameter :: ngs3  = 48
!      integer(kind=im), parameter :: ngs4  = 64
!      integer(kind=im), parameter :: ngs5  = 80
!      integer(kind=im), parameter :: ngs6  = 96
!      integer(kind=im), parameter :: ngs7  = 112
!      integer(kind=im), parameter :: ngs8  = 128
!      integer(kind=im), parameter :: ngs9  = 144
!      integer(kind=im), parameter :: ngs10 = 160
!      integer(kind=im), parameter :: ngs11 = 176
!      integer(kind=im), parameter :: ngs12 = 192
!      integer(kind=im), parameter :: ngs13 = 208
!      integer(kind=im), parameter :: ngs14 = 224
!      integer(kind=im), parameter :: ngs15 = 240
!      integer(kind=im), parameter :: ngs16 = 256

      end module parrrtm
