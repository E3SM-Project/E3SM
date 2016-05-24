!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_init.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:13 $

      module rrtmg_sw_init

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use rrsw_wvn
      use rrtmg_sw_setcoef, only: swatmref

      implicit none

      contains

! **************************************************************************
      subroutine rrtmg_sw_ini
! **************************************************************************
!
!  Original version:   Michael J. Iacono; February, 2004
!  Revision for F90 formatting:  M. J. Iacono, July, 2006
!
!  This subroutine performs calculations necessary for the initialization
!  of the shortwave model.  Lookup tables are computed for use in the SW
!  radiative transfer, and input absorption coefficient data for each
!  spectral band are reduced from 224 g-point intervals to 112.
! **************************************************************************

      use parrrsw, only : mg, nbndsw, ngptsw
      use rrsw_tbl, only: ntbl, tblint, pade, bpade, tau_tbl, exp_tbl
      use rrsw_vsn, only: hvrini, hnamini

! ------- Local -------

      integer :: ibnd, igc, ig, ind, ipr
      integer :: igcsm, iprsm
      integer :: itr

      real(kind=r8) :: wtsum, wtsm(mg)
      real(kind=r8) :: tfn

! ------- Definitions -------
!     Arrays for 10000-point look-up tables:
!     TAU_TBL  Clear-sky optical depth 
!     EXP_TBL  Exponential lookup table for transmittance
!     PADE     Pade approximation constant (= 0.278)
!     BPADE    Inverse of the Pade approximation constant
!

      hvrini = '$Revision: 1.2 $'

! Initialize model data
      call swdatinit
      call swcmbdat              ! g-point interval reduction data
      call swaerpr               ! aerosol optical properties
      call swcldpr               ! cloud optical properties
      call swatmref              ! reference MLS profile
      call sw_kgb16              ! molecular absorption coefficients
      call sw_kgb17
      call sw_kgb18
      call sw_kgb19
      call sw_kgb20
      call sw_kgb21
      call sw_kgb22
      call sw_kgb23
      call sw_kgb24
      call sw_kgb25
      call sw_kgb26
      call sw_kgb27
      call sw_kgb28
      call sw_kgb29

! Define exponential lookup tables for transmittance. Tau is
! computed as a function of the tau transition function, and transmittance 
! is calculated as a function of tau.  All tables are computed at intervals 
! of 0.0001.  The inverse of the constant used in the Pade approximation to 
! the tau transition function is set to bpade.

      exp_tbl(0) = 1.0_r8
      exp_tbl(ntbl) = 0.0_r8
      bpade = 1.0_r8 / pade
      do itr = 1, ntbl-1
         tfn = float(itr) / float(ntbl)
         tau_tbl = bpade * tfn / (1._r8 - tfn)
         exp_tbl(itr) = exp(-tau_tbl)
      enddo

! Perform g-point reduction from 16 per band (224 total points) to
! a band dependent number (112 total points) for all absorption
! coefficient input data and Planck fraction input data.
! Compute relative weighting for new g-point combinations.

      igcsm = 0
      do ibnd = 1,nbndsw
         iprsm = 0
         if (ngc(ibnd).lt.mg) then
            do igc = 1,ngc(ibnd)
               igcsm = igcsm + 1
               wtsum = 0.
               do ipr = 1, ngn(igcsm)
                  iprsm = iprsm + 1
                  wtsum = wtsum + wt(iprsm)
               enddo
               wtsm(igc) = wtsum
            enddo
            do ig = 1, ng(ibnd+15)
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = wt(ig)/wtsm(ngm(ind))
            enddo
         else
            do ig = 1, ng(ibnd+15)
               igcsm = igcsm + 1
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = 1.0_r8
            enddo
         endif
      enddo

! Reduce g-points for absorption coefficient data in each LW spectral band.

      call cmbgb16s
      call cmbgb17
      call cmbgb18
      call cmbgb19
      call cmbgb20
      call cmbgb21
      call cmbgb22
      call cmbgb23
      call cmbgb24
      call cmbgb25
      call cmbgb26
      call cmbgb27
      call cmbgb28
      call cmbgb29

      end subroutine rrtmg_sw_ini

!***************************************************************************
      subroutine swdatinit
!***************************************************************************

! --------- Modules ----------

      use rrsw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2 
      use rrsw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave
      use shr_const_mod, only: shr_const_avogad
      use physconst,     only: cday, gravit, cpair
      use rrsw_vsn

      save 
 
! Shortwave spectral band limits (wavenumbers)
      wavenum1(:) = (/2600._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
                      8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,  820._r8/)
      wavenum2(:) = (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
                     12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)
      delwave(:) =  (/ 650._r8,  750._r8,  650._r8,  500._r8, 1000._r8, 1550._r8,  350._r8, &
                      4800._r8, 3150._r8, 6650._r8, 6350._r8, 9000._r8,12000._r8, 1780._r8/)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/9,9,9,9,1,9,9,1,9,1,0,1,9,1/)
      nspb(:) = (/1,5,1,1,1,5,1,0,1,0,0,1,5,1/)

! Use constants set in CAM for consistency
      grav = gravit
      avogad = shr_const_avogad * 1.e-3_r8

!     Heatfac is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
!      heatfac = 8.4391_r8

!     Modified values for consistency with CAM3:
!        =  (9.80616)(86400)(1e-5)/(1.00464)
!      heatfac = 8.43339130434_r8

!     Calculate heatfac directly from CAM constants:
      heatfac = grav * cday * 1.e-5_r8 / (cpair * 1.e-3_r8)

!    Constants from NIST 01/11/2002

!      grav = 9.8066_r8
      planck = 6.62606876e-27_r8
      boltz = 1.3806503e-16_r8
      clight = 2.99792458e+10_r8
!      avogad = 6.02214199e+23_r8
      alosmt = 2.6867775e+19_r8
      gascon = 8.31447200e+07_r8
      radcn1 = 1.191042722e-12_r8
      radcn2 = 1.4387752_r8

!
!     units are generally cgs
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!          radcn1 = 2.*planck*clight*clight*1.e-07
!          radcn2 = planck*clight/boltz

      end subroutine swdatinit

!***************************************************************************
      subroutine swcmbdat
!***************************************************************************

      use rrsw_wvn, only: ngc, ngs, ngn, ngb, ngm, wt

      save
 
! ------- Definitions -------
!     Arrays for the g-point reduction from 224 to 112 for the 16 LW bands:
!     This mapping from 224 to 112 points has been carefully selected to 
!     minimize the effect on the resulting fluxes and cooling rates, and
!     caution should be used if the mapping is modified.  The full 224
!     g-point set can be restored with ngpt=224, ngc=16*16, ngn=224*1., etc.
!     ngpt    The total number of new g-points
!     ngc     The number of new g-points in each band
!     ngs     The cumulative sum of new g-points for each band
!     ngm     The index of each new g-point relative to the original
!             16 g-points for each band.  
!     ngn     The number of original g-points that are combined to make
!             each new g-point in each band.
!     ngb     The band index for each new g-point.
!     wt      RRTM weights for 16 g-points.

! Use this set for 112 quadrature point (g-point) model
! ------- Data statements -------
      ngc(:) = (/ 6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8, 6,12 /)
      ngs(:) = (/ 6,18,26,34,44,54,56,66,74,80,86,94,100,112 /)
      ngm(:) = (/ 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 16
                  1,2,3,4,5,6,6,7,8,8,9,10,10,11,12,12, &      ! band 17
                  1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 18
                  1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 19
                  1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 20
                  1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 21
                  1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 22
                  1,1,2,2,3,4,5,6,7,8,9,9,10,10,10,10, &       ! band 23
                  1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 24
                  1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 25
                  1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 26
                  1,2,3,4,5,6,7,7,7,7,8,8,8,8,8,8, &           ! band 27
                  1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 28
                  1,2,3,4,5,5,6,6,7,7,8,8,9,10,11,12 /)        ! band 29
      ngn(:) = (/ 2,2,2,2,4,4, &                               ! band 16
                  1,1,1,1,1,2,1,2,1,2,1,2, &                   ! band 17
                  1,1,1,1,2,2,4,4, &                           ! band 18
                  1,1,1,1,2,2,4,4, &                           ! band 19
                  1,1,1,1,1,1,1,1,2,6, &                       ! band 20
                  1,1,1,1,1,1,1,1,2,6, &                       ! band 21
                  8,8, &                                       ! band 22
                  2,2,1,1,1,1,1,1,2,4, &                       ! band 23
                  2,2,2,2,2,2,2,2, &                           ! band 24
                  1,1,2,2,4,6, &                               ! band 25
                  1,1,2,2,4,6, &                               ! band 26
                  1,1,1,1,1,1,4,6, &                           ! band 27
                  1,1,2,2,4,6, &                               ! band 28
                  1,1,1,1,2,2,2,2,1,1,1,1 /)                   ! band 29
      ngb(:) = (/ 16,16,16,16,16,16, &                         ! band 16
                  17,17,17,17,17,17,17,17,17,17,17,17, &       ! band 17
                  18,18,18,18,18,18,18,18, &                   ! band 18
                  19,19,19,19,19,19,19,19, &                   ! band 19
                  20,20,20,20,20,20,20,20,20,20, &             ! band 20
                  21,21,21,21,21,21,21,21,21,21, &             ! band 21
                  22,22, &                                     ! band 22
                  23,23,23,23,23,23,23,23,23,23, &             ! band 23
                  24,24,24,24,24,24,24,24, &                   ! band 24
                  25,25,25,25,25,25, &                         ! band 25
                  26,26,26,26,26,26, &                         ! band 26
                  27,27,27,27,27,27,27,27, &                   ! band 27
                  28,28,28,28,28,28, &                         ! band 28
                  29,29,29,29,29,29,29,29,29,29,29,29 /)       ! band 29

! Use this set for full 224 quadrature point (g-point) model
! ------- Data statements -------
!      ngc(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /)
!      ngs(:) = (/ 16,32,48,64,80,96,112,128,144,160,176,192,208,224 /)
!      ngm(:) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 16
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 17
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 18
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 19
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 20
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 21
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 22
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 23
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 24
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 25
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 26
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 27
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 28
!                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 /)    ! band 29
!      ngn(:) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 16
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 17
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 18
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 19
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 20
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 21
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 22
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 23
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 24
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 25
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 26
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 27
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 28
!                  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)           ! band 29
!      ngb(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16, &   ! band 16
!                  17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17, &   ! band 17
!                  18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18, &   ! band 18
!                  19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19, &   ! band 19
!                  20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20, &   ! band 20
!                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21, &   ! band 21
!                  22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22, &   ! band 22
!                  23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23, &   ! band 23
!                  24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24, &   ! band 24
!                  25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25, &   ! band 25
!                  26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26, &   ! band 26
!                  27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27, &   ! band 27
!                  28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28, &   ! band 28
!                  29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29 /)   ! band 29


      wt(:) =  (/ 0.1527534276_r8, 0.1491729617_r8, 0.1420961469_r8, &
                  0.1316886544_r8, 0.1181945205_r8, 0.1019300893_r8, &
                  0.0832767040_r8, 0.0626720116_r8, 0.0424925000_r8, &
                  0.0046269894_r8, 0.0038279891_r8, 0.0030260086_r8, &
                  0.0022199750_r8, 0.0014140010_r8, 0.0005330000_r8, &
                  0.0000750000_r8 /)

      end subroutine swcmbdat

!***************************************************************************
      subroutine cmbgb16s
!***************************************************************************
!
!  Original version:       MJIacono; July 1998
!  Revision for RRTM_SW:   MJIacono; November 2002
!  Revision for RRTMG_SW:  MJIacono; December 2003
!  Revision for F90 reformatting:  MJIacono; July 2006
!
!  The subroutines CMBGB16->CMBGB29 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 14 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMG_SW_INIT.  Solar source 
!  function data in array SFLUXREF are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTMG_SW.
!
!  band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg16, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(1)
                  sumk = 0.
                  do ipr = 1, ngn(igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(1)
         sumf = 0.
         do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumf = sumf + sfluxrefo(iprsm)
         enddo
         sfluxref(igc) = sumf
      enddo

      end subroutine cmbgb16s

!***************************************************************************
      subroutine cmbgb17
!***************************************************************************
!
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg17, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(2)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(1)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+16)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(2)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(1)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+16)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(2)
            sumf = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb17

!***************************************************************************
      subroutine cmbgb18
!***************************************************************************
!
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg18, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(3)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(3)
               sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+32)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb18

!***************************************************************************
      subroutine cmbgb19
!***************************************************************************
!
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg19, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(4)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(4)
               sumk = 0.
               do ipr = 1, ngn(ngs(3)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+48)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb19

!***************************************************************************
      subroutine cmbgb20
!***************************************************************************
!
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg20, only : kao, kbo, selfrefo, forrefo, sfluxrefo, absch4o, &
                            ka, kb, selfref, forref, sfluxref, absch4

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(5)
               sumk = 0.
               do ipr = 1, ngn(ngs(4)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+64)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(5)
               sumk = 0.
               do ipr = 1, ngn(ngs(4)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+64)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(5)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + sfluxrefo(iprsm)
            sumf2 = sumf2 + absch4o(iprsm)*rwgt(iprsm+64)
         enddo
         sfluxref(igc) = sumf1
         absch4(igc) = sumf2
      enddo

      end subroutine cmbgb20

!***************************************************************************
      subroutine cmbgb21
!***************************************************************************
!
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg21, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(6)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(5)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+80)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(6)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(5)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+80)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(6)
            sumf = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb21

!***************************************************************************
      subroutine cmbgb22
!***************************************************************************
!
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg22, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            ka, kb, selfref, forref, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(7)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(6)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(7)
               sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(7)
            sumf = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb22

!***************************************************************************
      subroutine cmbgb23
!***************************************************************************
!
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg23, only : kao, selfrefo, forrefo, sfluxrefo, raylo, &
                            ka, selfref, forref, sfluxref, rayl

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(8)
               sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(8)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + sfluxrefo(iprsm)
            sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+112)
         enddo
         sfluxref(igc) = sumf1
         rayl(igc) = sumf2
      enddo

      end subroutine cmbgb23

!***************************************************************************
      subroutine cmbgb24
!***************************************************************************
!
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg24, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            abso3ao, abso3bo, raylao, raylbo, &
                            ka, kb, selfref, forref, sfluxref, &
                            abso3a, abso3b, rayla, raylb

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2, sumf3


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(9)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(8)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(9)
               sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,3
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(9)
         sumf1 = 0.
         sumf2 = 0.
         sumf3 = 0.
         do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + raylbo(iprsm)*rwgt(iprsm+128)
            sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+128)
            sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+128)
         enddo
         raylb(igc) = sumf1
         abso3a(igc) = sumf2
         abso3b(igc) = sumf3
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(9)
            sumf1 = 0.
            sumf2 = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumf1 = sumf1 + sfluxrefo(iprsm,jp)
               sumf2 = sumf2 + raylao(iprsm,jp)*rwgt(iprsm+128)
            enddo
            sfluxref(igc,jp) = sumf1
            rayla(igc,jp) = sumf2
         enddo
      enddo

      end subroutine cmbgb24

!***************************************************************************
      subroutine cmbgb25
!***************************************************************************
!
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg25, only : kao, sfluxrefo, &
                            abso3ao, abso3bo, raylo, &
                            ka, sfluxref, &
                            abso3a, abso3b, rayl

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2, sumf3, sumf4


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(10)
         sumf1 = 0.
         sumf2 = 0.
         sumf3 = 0.
         sumf4 = 0.
         do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + sfluxrefo(iprsm)
            sumf2 = sumf2 + abso3ao(iprsm)*rwgt(iprsm+144)
            sumf3 = sumf3 + abso3bo(iprsm)*rwgt(iprsm+144)
            sumf4 = sumf4 + raylo(iprsm)*rwgt(iprsm+144)
         enddo
         sfluxref(igc) = sumf1
         abso3a(igc) = sumf2
         abso3b(igc) = sumf3
         rayl(igc) = sumf4
      enddo

      end subroutine cmbgb25

!***************************************************************************
      subroutine cmbgb26
!***************************************************************************
!
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg26, only : sfluxrefo, raylo, &
                            sfluxref, rayl

! ------- Local -------
      integer :: igc, ipr, iprsm
      real(kind=r8) :: sumf1, sumf2


      iprsm = 0
      do igc = 1,ngc(11)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + raylo(iprsm)*rwgt(iprsm+160)
            sumf2 = sumf2 + sfluxrefo(iprsm)
         enddo
         rayl(igc) = sumf1
         sfluxref(igc) = sumf2
      enddo

      end subroutine cmbgb26

!***************************************************************************
      subroutine cmbgb27
!***************************************************************************
!
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg27, only : kao, kbo, sfluxrefo, raylo, &
                            ka, kb, sfluxref, rayl

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(12)
               sumk = 0.
               do ipr = 1, ngn(ngs(11)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+176)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(12)
               sumk = 0.
               do ipr = 1, ngn(ngs(11)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+176)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(12)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(11)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + sfluxrefo(iprsm)
            sumf2 = sumf2 + raylo(iprsm)*rwgt(iprsm+176)
         enddo
         sfluxref(igc) = sumf1
         rayl(igc) = sumf2
      enddo

      end subroutine cmbgb27

!***************************************************************************
      subroutine cmbgb28
!***************************************************************************
!
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg28, only : kao, kbo, sfluxrefo, &
                            ka, kb, sfluxref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(13)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(12)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(13)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(12)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+192)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(13)
            sumf = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumf = sumf + sfluxrefo(iprsm,jp)
            enddo
            sfluxref(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb28

!***************************************************************************
      subroutine cmbgb29
!***************************************************************************
!
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)
!-----------------------------------------------------------------------

      use rrsw_wvn, only : ngc, ngs, ngn, wt, rwgt
      use rrsw_kg29, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
                            absh2oo, absco2o, &
                            ka, kb, selfref, forref, sfluxref, &
                            absh2o, absco2

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm
      real(kind=r8) :: sumk, sumf1, sumf2, sumf3


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(14)
         sumf1 = 0.
         sumf2 = 0.
         sumf3 = 0.
         do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumf1 = sumf1 + sfluxrefo(iprsm)
            sumf2 = sumf2 + absco2o(iprsm)*rwgt(iprsm+208)
            sumf3 = sumf3 + absh2oo(iprsm)*rwgt(iprsm+208)
         enddo
         sfluxref(igc) = sumf1
         absco2(igc) = sumf2
         absh2o(igc) = sumf3
      enddo

      end subroutine cmbgb29

!***************************************************************************
      subroutine swaerpr
!***************************************************************************

! Purpose: Define spectral aerosol properties for six ECMWF aerosol types
! as used in the ECMWF IFS model (see module rrsw_aer.F90 for details)
!
! Original: Defined for rrtmg_sw 14 spectral bands, JJMorcrette, ECMWF Feb 2003
! Revision: Reformatted for consistency with rrtmg_lw, MJIacono, AER, Jul 2006

      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya

      save

      rsrtaua( 1, :) = (/ &
        0.10849_r8, 0.66699_r8, 0.65255_r8, 0.11600_r8, 0.06529_r8, 0.04468_r8/)
      rsrtaua( 2, :) = (/ &
        0.10849_r8, 0.66699_r8, 0.65255_r8, 0.11600_r8, 0.06529_r8, 0.04468_r8/)
      rsrtaua( 3, :) = (/ &
        0.20543_r8, 0.84642_r8, 0.84958_r8, 0.21673_r8, 0.28270_r8, 0.10915_r8/)
      rsrtaua( 4, :) = (/ &
        0.20543_r8, 0.84642_r8, 0.84958_r8, 0.21673_r8, 0.28270_r8, 0.10915_r8/)
      rsrtaua( 5, :) = (/ &
        0.20543_r8, 0.84642_r8, 0.84958_r8, 0.21673_r8, 0.28270_r8, 0.10915_r8/)
      rsrtaua( 6, :) = (/ &
        0.20543_r8, 0.84642_r8, 0.84958_r8, 0.21673_r8, 0.28270_r8, 0.10915_r8/)
      rsrtaua( 7, :) = (/ &
        0.20543_r8, 0.84642_r8, 0.84958_r8, 0.21673_r8, 0.28270_r8, 0.10915_r8/)
      rsrtaua( 8, :) = (/ &
        0.52838_r8, 0.93285_r8, 0.93449_r8, 0.53078_r8, 0.67148_r8, 0.46608_r8/)
      rsrtaua( 9, :) = (/ &
        0.52838_r8, 0.93285_r8, 0.93449_r8, 0.53078_r8, 0.67148_r8, 0.46608_r8/)
      rsrtaua(10, :) = (/ &
        1.69446_r8, 1.11855_r8, 1.09212_r8, 1.72145_r8, 1.03858_r8, 1.12044_r8/)
      rsrtaua(11, :) = (/ &
        1.69446_r8, 1.11855_r8, 1.09212_r8, 1.72145_r8, 1.03858_r8, 1.12044_r8/)
      rsrtaua(12, :) = (/ &
        1.69446_r8, 1.11855_r8, 1.09212_r8, 1.72145_r8, 1.03858_r8, 1.12044_r8/)
      rsrtaua(13, :) = (/ &
        1.69446_r8, 1.11855_r8, 1.09212_r8, 1.72145_r8, 1.03858_r8, 1.12044_r8/)
      rsrtaua(14, :) = (/ &
        0.10849_r8, 0.66699_r8, 0.65255_r8, 0.11600_r8, 0.06529_r8, 0.04468_r8/)
 
      rsrpiza( 1, :) = (/ &
        .5230504_r8, .7868518_r8, .8531531_r8, .4048149_r8, .8748231_r8, .2355667_r8/)
      rsrpiza( 2, :) = (/ &
        .5230504_r8, .7868518_r8, .8531531_r8, .4048149_r8, .8748231_r8, .2355667_r8/)
      rsrpiza( 3, :) = (/ &
        .8287144_r8, .9949396_r8, .9279543_r8, .6765051_r8, .9467578_r8, .9955938_r8/)
      rsrpiza( 4, :) = (/ &
        .8287144_r8, .9949396_r8, .9279543_r8, .6765051_r8, .9467578_r8, .9955938_r8/)
      rsrpiza( 5, :) = (/ &
        .8287144_r8, .9949396_r8, .9279543_r8, .6765051_r8, .9467578_r8, .9955938_r8/)
      rsrpiza( 6, :) = (/ &
        .8287144_r8, .9949396_r8, .9279543_r8, .6765051_r8, .9467578_r8, .9955938_r8/)
      rsrpiza( 7, :) = (/ &
        .8287144_r8, .9949396_r8, .9279543_r8, .6765051_r8, .9467578_r8, .9955938_r8/)
      rsrpiza( 8, :) = (/ &
        .8970131_r8, .9984940_r8, .9245594_r8, .7768385_r8, .9532763_r8, .9999999_r8/)
      rsrpiza( 9, :) = (/ &
        .8970131_r8, .9984940_r8, .9245594_r8, .7768385_r8, .9532763_r8, .9999999_r8/)
      rsrpiza(10, :) = (/ &
        .9148907_r8, .9956173_r8, .7504584_r8, .8131335_r8, .9401905_r8, .9999999_r8/)
      rsrpiza(11, :) = (/ &
        .9148907_r8, .9956173_r8, .7504584_r8, .8131335_r8, .9401905_r8, .9999999_r8/)
      rsrpiza(12, :) = (/ &
        .9148907_r8, .9956173_r8, .7504584_r8, .8131335_r8, .9401905_r8, .9999999_r8/)
      rsrpiza(13, :) = (/ &
        .9148907_r8, .9956173_r8, .7504584_r8, .8131335_r8, .9401905_r8, .9999999_r8/)
      rsrpiza(14, :) = (/ &
        .5230504_r8, .7868518_r8, .8531531_r8, .4048149_r8, .8748231_r8, .2355667_r8/)

      rsrasya( 1, :) = (/ &
        0.700610_r8, 0.818871_r8, 0.702399_r8, 0.689886_r8, .4629866_r8, .1907639_r8/)
      rsrasya( 2, :) = (/ &
        0.700610_r8, 0.818871_r8, 0.702399_r8, 0.689886_r8, .4629866_r8, .1907639_r8/)
      rsrasya( 3, :) = (/ &
        0.636342_r8, 0.802467_r8, 0.691305_r8, 0.627497_r8, .6105750_r8, .4760794_r8/)
      rsrasya( 4, :) = (/ &
        0.636342_r8, 0.802467_r8, 0.691305_r8, 0.627497_r8, .6105750_r8, .4760794_r8/)
      rsrasya( 5, :) = (/ &
        0.636342_r8, 0.802467_r8, 0.691305_r8, 0.627497_r8, .6105750_r8, .4760794_r8/)
      rsrasya( 6, :) = (/ &
        0.636342_r8, 0.802467_r8, 0.691305_r8, 0.627497_r8, .6105750_r8, .4760794_r8/)
      rsrasya( 7, :) = (/ &
        0.636342_r8, 0.802467_r8, 0.691305_r8, 0.627497_r8, .6105750_r8, .4760794_r8/)
      rsrasya( 8, :) = (/ &
        0.668431_r8, 0.788530_r8, 0.698682_r8, 0.657422_r8, .6735182_r8, .6519706_r8/)
      rsrasya( 9, :) = (/ &
        0.668431_r8, 0.788530_r8, 0.698682_r8, 0.657422_r8, .6735182_r8, .6519706_r8/)
      rsrasya(10, :) = (/ &
        0.729019_r8, 0.803129_r8, 0.784592_r8, 0.712208_r8, .7008249_r8, .7270548_r8/)
      rsrasya(11, :) = (/ &
        0.729019_r8, 0.803129_r8, 0.784592_r8, 0.712208_r8, .7008249_r8, .7270548_r8/)
      rsrasya(12, :) = (/ &
        0.729019_r8, 0.803129_r8, 0.784592_r8, 0.712208_r8, .7008249_r8, .7270548_r8/)
      rsrasya(13, :) = (/ &
        0.729019_r8, 0.803129_r8, 0.784592_r8, 0.712208_r8, .7008249_r8, .7270548_r8/)
      rsrasya(14, :) = (/ &
        0.700610_r8, 0.818871_r8, 0.702399_r8, 0.689886_r8, .4629866_r8, .1907639_r8/)

      end subroutine swaerpr
 
!***********************************************************************
      subroutine swcldpr
!***********************************************************************

! Purpose: Define cloud extinction coefficient, single scattering albedo
!          and asymmetry parameter data.
!

! ------- Modules -------

      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari

      save

!-----------------------------------------------------------------------
!
! Explanation of the method for each value of INFLAG.  A value of
!  0 for INFLAG do not distingish being liquid and ice clouds.
!  INFLAG = 2 does distinguish between liquid and ice clouds, and
!    requires further user input to specify the method to be used to 
!    compute the aborption due to each.
!  INFLAG = 0:  For each cloudy layer, the cloud fraction, the cloud optical
!    depth, the cloud single-scattering albedo, and the
!    moments of the phase function (0:NSTREAM).  Note
!    that these values are delta-m scaled within this
!    subroutine.

!  INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!    water path (g/m2), and cloud ice fraction are input.
!  ICEFLAG = 2:  The ice effective radius (microns) is input and the
!    optical properties due to ice clouds are computed from
!    the optical properties stored in the RT code, STREAMER v3.0 
!    (Reference: Key. J., Streamer User's Guide, Cooperative 
!    Institute for Meteorological Satellite Studies, 2001, 96 pp.).
!    Valid range of values for re are between 5.0 and
!    131.0 micron.
!    This version uses Ebert and Curry, JGR, (1992) method for 
!    ice particles larger than 131.0 microns. 
!  ICEFLAG = 3:  The ice generalized effective size (dge) is input
!    and the optical depths, single-scattering albedo,
!    and phase function moments are calculated as in
!    Q. Fu, J. Climate, (1996). Q. Fu provided high resolution
!    tables which were appropriately averaged for the
!    bands in RRTM_SW.  Linear interpolation is used to
!    get the coefficients from the stored tables.
!    Valid range of values for dge are between 5.0 and
!    140.0 micron. 
!    This version uses Ebert and Curry, JGR, (1992) method for 
!    ice particles larger than 140.0 microns. 
!  LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!    and the optical depths due to water clouds are computed 
!    as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!    The values for absorption coefficients appropriate for
!    the spectral bands in RRTM have been obtained for a 
!    range of effective radii by an averaging procedure 
!    based on the work of J. Pinto (private communication).
!    Linear interpolation is used to get the absorption 
!    coefficients for the input effective radius.
!
!     ------------------------------------------------------------------

! Everything below is for INFLAG = 2.

! Coefficients for Ebert and Curry method
      abari(:) = (/ &
        & 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8,3.448e-03_r8,3.448e-03_r8 /)
      bbari(:) = (/ &
        & 2.431e+00_r8,2.431e+00_r8,2.431e+00_r8,2.431e+00_r8,2.431e+00_r8 /)
      cbari(:) = (/ &
        & 1.000e-05_r8,1.100e-04_r8,1.240e-02_r8,3.779e-02_r8,4.666e-01_r8 /)
      dbari(:) = (/ &
        & 0.000e+00_r8,1.405e-05_r8,6.867e-04_r8,1.284e-03_r8,2.050e-05_r8 /)
      ebari(:) = (/ &
        & 7.661e-01_r8,7.730e-01_r8,7.865e-01_r8,8.172e-01_r8,9.595e-01_r8 /)
      fbari(:) = (/ &
        & 5.851e-04_r8,5.665e-04_r8,7.204e-04_r8,7.463e-04_r8,1.076e-04_r8 /)

! Extinction coefficient
      extliq1(:, 16) = (/ &
        & 8.981463e-01_r8,6.317895e-01_r8,4.557508e-01_r8,3.481624e-01_r8,2.797950e-01_r8,&
        & 2.342753e-01_r8,2.026934e-01_r8,1.800102e-01_r8,1.632408e-01_r8,1.505384e-01_r8,&
        & 1.354524e-01_r8,1.246520e-01_r8,1.154342e-01_r8,1.074756e-01_r8,1.005353e-01_r8,&
        & 9.442987e-02_r8,8.901760e-02_r8,8.418693e-02_r8,7.984904e-02_r8,7.593229e-02_r8,&
        & 7.237827e-02_r8,6.913887e-02_r8,6.617415e-02_r8,6.345061e-02_r8,6.094001e-02_r8,&
        & 5.861834e-02_r8,5.646506e-02_r8,5.446250e-02_r8,5.249596e-02_r8,5.081114e-02_r8,&
        & 4.922243e-02_r8,4.772189e-02_r8,4.630243e-02_r8,4.495766e-02_r8,4.368189e-02_r8,&
        & 4.246995e-02_r8,4.131720e-02_r8,4.021941e-02_r8,3.917276e-02_r8,3.817376e-02_r8,&
        & 3.721926e-02_r8,3.630635e-02_r8,3.543237e-02_r8,3.459491e-02_r8,3.379171e-02_r8,&
        & 3.302073e-02_r8,3.228007e-02_r8,3.156798e-02_r8,3.088284e-02_r8,3.022315e-02_r8,&
        & 2.958753e-02_r8,2.897468e-02_r8,2.838340e-02_r8,2.781258e-02_r8,2.726117e-02_r8,&
        & 2.672821e-02_r8,2.621278e-02_r8,2.5714e-02_r8 /)
      extliq1(:, 17) = (/ &
        & 8.293797e-01_r8,6.048371e-01_r8,4.465706e-01_r8,3.460387e-01_r8,2.800064e-01_r8,&
        & 2.346584e-01_r8,2.022399e-01_r8,1.782626e-01_r8,1.600153e-01_r8,1.457903e-01_r8,&
        & 1.334061e-01_r8,1.228548e-01_r8,1.138396e-01_r8,1.060486e-01_r8,9.924856e-02_r8,&
        & 9.326208e-02_r8,8.795158e-02_r8,8.320883e-02_r8,7.894750e-02_r8,7.509792e-02_r8,&
        & 7.160323e-02_r8,6.841653e-02_r8,6.549889e-02_r8,6.281763e-02_r8,6.034516e-02_r8,&
        & 5.805802e-02_r8,5.593615e-02_r8,5.396226e-02_r8,5.202302e-02_r8,5.036246e-02_r8,&
        & 4.879606e-02_r8,4.731610e-02_r8,4.591565e-02_r8,4.458852e-02_r8,4.332912e-02_r8,&
        & 4.213243e-02_r8,4.099390e-02_r8,3.990941e-02_r8,3.887522e-02_r8,3.788792e-02_r8,&
        & 3.694440e-02_r8,3.604183e-02_r8,3.517760e-02_r8,3.434934e-02_r8,3.355485e-02_r8,&
        & 3.279211e-02_r8,3.205925e-02_r8,3.135458e-02_r8,3.067648e-02_r8,3.002349e-02_r8,&
        & 2.939425e-02_r8,2.878748e-02_r8,2.820200e-02_r8,2.763673e-02_r8,2.709062e-02_r8,&
        & 2.656272e-02_r8,2.605214e-02_r8,2.5558e-02_r8 /)
      extliq1(:, 18) = (/ &
        & 9.193685e-01_r8,6.128292e-01_r8,4.344150e-01_r8,3.303048e-01_r8,2.659500e-01_r8,&
        & 2.239727e-01_r8,1.953457e-01_r8,1.751012e-01_r8,1.603515e-01_r8,1.493360e-01_r8,&
        & 1.323791e-01_r8,1.219335e-01_r8,1.130076e-01_r8,1.052926e-01_r8,9.855839e-02_r8,&
        & 9.262925e-02_r8,8.736918e-02_r8,8.267112e-02_r8,7.844965e-02_r8,7.463585e-02_r8,&
        & 7.117343e-02_r8,6.801601e-02_r8,6.512503e-02_r8,6.246815e-02_r8,6.001806e-02_r8,&
        & 5.775154e-02_r8,5.564872e-02_r8,5.369250e-02_r8,5.176284e-02_r8,5.011536e-02_r8,&
        & 4.856099e-02_r8,4.709211e-02_r8,4.570193e-02_r8,4.438430e-02_r8,4.313375e-02_r8,&
        & 4.194529e-02_r8,4.081443e-02_r8,3.973712e-02_r8,3.870966e-02_r8,3.772866e-02_r8,&
        & 3.679108e-02_r8,3.589409e-02_r8,3.503514e-02_r8,3.421185e-02_r8,3.342206e-02_r8,&
        & 3.266377e-02_r8,3.193513e-02_r8,3.123447e-02_r8,3.056018e-02_r8,2.991081e-02_r8,&
        & 2.928502e-02_r8,2.868154e-02_r8,2.809920e-02_r8,2.753692e-02_r8,2.699367e-02_r8,&
        & 2.646852e-02_r8,2.596057e-02_r8,2.5469e-02_r8 /)
      extliq1(:, 19) = (/ &
        & 9.136931e-01_r8,5.743244e-01_r8,4.080708e-01_r8,3.150572e-01_r8,2.577261e-01_r8,&
        & 2.197900e-01_r8,1.933037e-01_r8,1.740212e-01_r8,1.595056e-01_r8,1.482756e-01_r8,&
        & 1.312164e-01_r8,1.209246e-01_r8,1.121227e-01_r8,1.045095e-01_r8,9.785967e-02_r8,&
        & 9.200149e-02_r8,8.680170e-02_r8,8.215531e-02_r8,7.797850e-02_r8,7.420361e-02_r8,&
        & 7.077530e-02_r8,6.764798e-02_r8,6.478369e-02_r8,6.215063e-02_r8,5.972189e-02_r8,&
        & 5.747458e-02_r8,5.538913e-02_r8,5.344866e-02_r8,5.153216e-02_r8,4.989745e-02_r8,&
        & 4.835476e-02_r8,4.689661e-02_r8,4.551629e-02_r8,4.420777e-02_r8,4.296563e-02_r8,&
        & 4.178497e-02_r8,4.066137e-02_r8,3.959081e-02_r8,3.856963e-02_r8,3.759452e-02_r8,&
        & 3.666244e-02_r8,3.577061e-02_r8,3.491650e-02_r8,3.409777e-02_r8,3.331227e-02_r8,&
        & 3.255803e-02_r8,3.183322e-02_r8,3.113617e-02_r8,3.046530e-02_r8,2.981918e-02_r8,&
        & 2.919646e-02_r8,2.859591e-02_r8,2.801635e-02_r8,2.745671e-02_r8,2.691599e-02_r8,&
        & 2.639324e-02_r8,2.588759e-02_r8,2.5398e-02_r8 /)
      extliq1(:, 20) = (/ &
        & 8.447548e-01_r8,5.326840e-01_r8,3.921523e-01_r8,3.119082e-01_r8,2.597055e-01_r8,&
        & 2.228737e-01_r8,1.954157e-01_r8,1.741155e-01_r8,1.570881e-01_r8,1.431520e-01_r8,&
        & 1.302034e-01_r8,1.200491e-01_r8,1.113571e-01_r8,1.038330e-01_r8,9.725657e-02_r8,&
        & 9.145949e-02_r8,8.631112e-02_r8,8.170840e-02_r8,7.756901e-02_r8,7.382641e-02_r8,&
        & 7.042616e-02_r8,6.732338e-02_r8,6.448069e-02_r8,6.186672e-02_r8,5.945494e-02_r8,&
        & 5.722277e-02_r8,5.515089e-02_r8,5.322262e-02_r8,5.132153e-02_r8,4.969799e-02_r8,&
        & 4.816556e-02_r8,4.671686e-02_r8,4.534525e-02_r8,4.404480e-02_r8,4.281014e-02_r8,&
        & 4.163643e-02_r8,4.051930e-02_r8,3.945479e-02_r8,3.843927e-02_r8,3.746945e-02_r8,&
        & 3.654234e-02_r8,3.565518e-02_r8,3.480547e-02_r8,3.399088e-02_r8,3.320930e-02_r8,&
        & 3.245876e-02_r8,3.173745e-02_r8,3.104371e-02_r8,3.037600e-02_r8,2.973287e-02_r8,&
        & 2.911300e-02_r8,2.851516e-02_r8,2.793818e-02_r8,2.738101e-02_r8,2.684264e-02_r8,&
        & 2.632214e-02_r8,2.581863e-02_r8,2.5331e-02_r8 /)
      extliq1(:, 21) = (/ &
        & 7.727642e-01_r8,5.034865e-01_r8,3.808673e-01_r8,3.080333e-01_r8,2.586453e-01_r8,&
        & 2.224989e-01_r8,1.947060e-01_r8,1.725821e-01_r8,1.545096e-01_r8,1.394456e-01_r8,&
        & 1.288683e-01_r8,1.188852e-01_r8,1.103317e-01_r8,1.029214e-01_r8,9.643967e-02_r8,&
        & 9.072239e-02_r8,8.564194e-02_r8,8.109758e-02_r8,7.700875e-02_r8,7.331026e-02_r8,&
        & 6.994879e-02_r8,6.688028e-02_r8,6.406807e-02_r8,6.148133e-02_r8,5.909400e-02_r8,&
        & 5.688388e-02_r8,5.483197e-02_r8,5.292185e-02_r8,5.103763e-02_r8,4.942905e-02_r8,&
        & 4.791039e-02_r8,4.647438e-02_r8,4.511453e-02_r8,4.382497e-02_r8,4.260043e-02_r8,&
        & 4.143616e-02_r8,4.032784e-02_r8,3.927155e-02_r8,3.826375e-02_r8,3.730117e-02_r8,&
        & 3.638087e-02_r8,3.550013e-02_r8,3.465646e-02_r8,3.384759e-02_r8,3.307141e-02_r8,&
        & 3.232598e-02_r8,3.160953e-02_r8,3.092040e-02_r8,3.025706e-02_r8,2.961810e-02_r8,&
        & 2.900220e-02_r8,2.840814e-02_r8,2.783478e-02_r8,2.728106e-02_r8,2.674599e-02_r8,&
        & 2.622864e-02_r8,2.572816e-02_r8,2.5244e-02_r8 /)
      extliq1(:, 22) = (/ &
        & 7.416833e-01_r8,4.959591e-01_r8,3.775057e-01_r8,3.056353e-01_r8,2.565943e-01_r8,&
        & 2.206935e-01_r8,1.931479e-01_r8,1.712860e-01_r8,1.534837e-01_r8,1.386906e-01_r8,&
        & 1.281198e-01_r8,1.182344e-01_r8,1.097595e-01_r8,1.024137e-01_r8,9.598552e-02_r8,&
        & 9.031320e-02_r8,8.527093e-02_r8,8.075927e-02_r8,7.669869e-02_r8,7.302481e-02_r8,&
        & 6.968491e-02_r8,6.663542e-02_r8,6.384008e-02_r8,6.126838e-02_r8,5.889452e-02_r8,&
        & 5.669654e-02_r8,5.465558e-02_r8,5.275540e-02_r8,5.087937e-02_r8,4.927904e-02_r8,&
        & 4.776796e-02_r8,4.633895e-02_r8,4.498557e-02_r8,4.370202e-02_r8,4.248306e-02_r8,&
        & 4.132399e-02_r8,4.022052e-02_r8,3.916878e-02_r8,3.816523e-02_r8,3.720665e-02_r8,&
        & 3.629011e-02_r8,3.541290e-02_r8,3.457257e-02_r8,3.376685e-02_r8,3.299365e-02_r8,&
        & 3.225105e-02_r8,3.153728e-02_r8,3.085069e-02_r8,3.018977e-02_r8,2.955310e-02_r8,&
        & 2.893940e-02_r8,2.834742e-02_r8,2.777606e-02_r8,2.722424e-02_r8,2.669099e-02_r8,&
        & 2.617539e-02_r8,2.567658e-02_r8,2.5194e-02_r8 /)
      extliq1(:, 23) = (/ &
        & 7.058580e-01_r8,4.866573e-01_r8,3.712238e-01_r8,2.998638e-01_r8,2.513441e-01_r8,&
        & 2.161972e-01_r8,1.895576e-01_r8,1.686669e-01_r8,1.518437e-01_r8,1.380046e-01_r8,&
        & 1.267564e-01_r8,1.170399e-01_r8,1.087026e-01_r8,1.014704e-01_r8,9.513729e-02_r8,&
        & 8.954555e-02_r8,8.457221e-02_r8,8.012009e-02_r8,7.611136e-02_r8,7.248294e-02_r8,&
        & 6.918317e-02_r8,6.616934e-02_r8,6.340584e-02_r8,6.086273e-02_r8,5.851465e-02_r8,&
        & 5.634001e-02_r8,5.432027e-02_r8,5.243946e-02_r8,5.058070e-02_r8,4.899628e-02_r8,&
        & 4.749975e-02_r8,4.608411e-02_r8,4.474303e-02_r8,4.347082e-02_r8,4.226237e-02_r8,&
        & 4.111303e-02_r8,4.001861e-02_r8,3.897528e-02_r8,3.797959e-02_r8,3.702835e-02_r8,&
        & 3.611867e-02_r8,3.524791e-02_r8,3.441364e-02_r8,3.361360e-02_r8,3.284577e-02_r8,&
        & 3.210823e-02_r8,3.139923e-02_r8,3.071716e-02_r8,3.006052e-02_r8,2.942791e-02_r8,&
        & 2.881806e-02_r8,2.822974e-02_r8,2.766185e-02_r8,2.711335e-02_r8,2.658326e-02_r8,&
        & 2.607066e-02_r8,2.557473e-02_r8,2.5095e-02_r8 /)
      extliq1(:, 24) = (/ &
        & 6.822779e-01_r8,4.750373e-01_r8,3.634834e-01_r8,2.940726e-01_r8,2.468060e-01_r8,&
        & 2.125768e-01_r8,1.866586e-01_r8,1.663588e-01_r8,1.500326e-01_r8,1.366192e-01_r8,&
        & 1.253472e-01_r8,1.158052e-01_r8,1.076101e-01_r8,1.004954e-01_r8,9.426089e-02_r8,&
        & 8.875268e-02_r8,8.385090e-02_r8,7.946063e-02_r8,7.550578e-02_r8,7.192466e-02_r8,&
        & 6.866669e-02_r8,6.569001e-02_r8,6.295971e-02_r8,6.044642e-02_r8,5.812526e-02_r8,&
        & 5.597500e-02_r8,5.397746e-02_r8,5.211690e-02_r8,5.027505e-02_r8,4.870703e-02_r8,&
        & 4.722555e-02_r8,4.582373e-02_r8,4.449540e-02_r8,4.323497e-02_r8,4.203742e-02_r8,&
        & 4.089821e-02_r8,3.981321e-02_r8,3.877867e-02_r8,3.779118e-02_r8,3.684762e-02_r8,&
        & 3.594514e-02_r8,3.508114e-02_r8,3.425322e-02_r8,3.345917e-02_r8,3.269698e-02_r8,&
        & 3.196477e-02_r8,3.126082e-02_r8,3.058352e-02_r8,2.993141e-02_r8,2.930310e-02_r8,&
        & 2.869732e-02_r8,2.811289e-02_r8,2.754869e-02_r8,2.700371e-02_r8,2.647698e-02_r8,&
        & 2.596760e-02_r8,2.547473e-02_r8,2.4998e-02_r8 /)
      extliq1(:, 25) = (/ &
        & 6.666233e-01_r8,4.662044e-01_r8,3.579517e-01_r8,2.902984e-01_r8,2.440475e-01_r8,&
        & 2.104431e-01_r8,1.849277e-01_r8,1.648970e-01_r8,1.487555e-01_r8,1.354714e-01_r8,&
        & 1.244173e-01_r8,1.149913e-01_r8,1.068903e-01_r8,9.985323e-02_r8,9.368351e-02_r8,&
        & 8.823009e-02_r8,8.337507e-02_r8,7.902511e-02_r8,7.510529e-02_r8,7.155482e-02_r8,&
        & 6.832386e-02_r8,6.537113e-02_r8,6.266218e-02_r8,6.016802e-02_r8,5.786408e-02_r8,&
        & 5.572939e-02_r8,5.374598e-02_r8,5.189830e-02_r8,5.006825e-02_r8,4.851081e-02_r8,&
        & 4.703906e-02_r8,4.564623e-02_r8,4.432621e-02_r8,4.307349e-02_r8,4.188312e-02_r8,&
        & 4.075060e-02_r8,3.967183e-02_r8,3.864313e-02_r8,3.766111e-02_r8,3.672269e-02_r8,&
        & 3.582505e-02_r8,3.496559e-02_r8,3.414196e-02_r8,3.335198e-02_r8,3.259362e-02_r8,&
        & 3.186505e-02_r8,3.116454e-02_r8,3.049052e-02_r8,2.984152e-02_r8,2.921617e-02_r8,&
        & 2.861322e-02_r8,2.803148e-02_r8,2.746986e-02_r8,2.692733e-02_r8,2.640295e-02_r8,&
        & 2.589582e-02_r8,2.540510e-02_r8,2.4930e-02_r8 /)
      extliq1(:, 26) = (/ &
        & 6.535669e-01_r8,4.585865e-01_r8,3.529226e-01_r8,2.867245e-01_r8,2.413848e-01_r8,&
        & 2.083956e-01_r8,1.833191e-01_r8,1.636150e-01_r8,1.477247e-01_r8,1.346392e-01_r8,&
        & 1.236449e-01_r8,1.143095e-01_r8,1.062828e-01_r8,9.930773e-02_r8,9.319029e-02_r8,&
        & 8.778150e-02_r8,8.296497e-02_r8,7.864847e-02_r8,7.475799e-02_r8,7.123343e-02_r8,&
        & 6.802549e-02_r8,6.509332e-02_r8,6.240285e-02_r8,5.992538e-02_r8,5.763657e-02_r8,&
        & 5.551566e-02_r8,5.354483e-02_r8,5.170870e-02_r8,4.988866e-02_r8,4.834061e-02_r8,&
        & 4.687751e-02_r8,4.549264e-02_r8,4.417999e-02_r8,4.293410e-02_r8,4.175006e-02_r8,&
        & 4.062344e-02_r8,3.955019e-02_r8,3.852663e-02_r8,3.754943e-02_r8,3.661553e-02_r8,&
        & 3.572214e-02_r8,3.486669e-02_r8,3.404683e-02_r8,3.326040e-02_r8,3.250542e-02_r8,&
        & 3.178003e-02_r8,3.108254e-02_r8,3.041139e-02_r8,2.976511e-02_r8,2.914235e-02_r8,&
        & 2.854187e-02_r8,2.796247e-02_r8,2.740309e-02_r8,2.686271e-02_r8,2.634038e-02_r8,&
        & 2.583520e-02_r8,2.534636e-02_r8,2.4873e-02_r8 /)
      extliq1(:, 27) = (/ &
        & 6.448790e-01_r8,4.541425e-01_r8,3.503348e-01_r8,2.850494e-01_r8,2.401966e-01_r8,&
        & 2.074811e-01_r8,1.825631e-01_r8,1.629515e-01_r8,1.471142e-01_r8,1.340574e-01_r8,&
        & 1.231462e-01_r8,1.138628e-01_r8,1.058802e-01_r8,9.894286e-02_r8,9.285818e-02_r8,&
        & 8.747802e-02_r8,8.268676e-02_r8,7.839271e-02_r8,7.452230e-02_r8,7.101580e-02_r8,&
        & 6.782418e-02_r8,6.490685e-02_r8,6.222991e-02_r8,5.976484e-02_r8,5.748742e-02_r8,&
        & 5.537703e-02_r8,5.341593e-02_r8,5.158883e-02_r8,4.977355e-02_r8,4.823172e-02_r8,&
        & 4.677430e-02_r8,4.539465e-02_r8,4.408680e-02_r8,4.284533e-02_r8,4.166539e-02_r8,&
        & 4.054257e-02_r8,3.947283e-02_r8,3.845256e-02_r8,3.747842e-02_r8,3.654737e-02_r8,&
        & 3.565665e-02_r8,3.480370e-02_r8,3.398620e-02_r8,3.320198e-02_r8,3.244908e-02_r8,&
        & 3.172566e-02_r8,3.103002e-02_r8,3.036062e-02_r8,2.971600e-02_r8,2.909482e-02_r8,&
        & 2.849582e-02_r8,2.791785e-02_r8,2.735982e-02_r8,2.682072e-02_r8,2.629960e-02_r8,&
        & 2.579559e-02_r8,2.530786e-02_r8,2.4836e-02_r8 /)
      extliq1(:, 28) = (/ &
        & 6.422688e-01_r8,4.528453e-01_r8,3.497232e-01_r8,2.847724e-01_r8,2.400815e-01_r8,&
        & 2.074403e-01_r8,1.825502e-01_r8,1.629415e-01_r8,1.470934e-01_r8,1.340183e-01_r8,&
        & 1.230935e-01_r8,1.138049e-01_r8,1.058201e-01_r8,9.888245e-02_r8,9.279878e-02_r8,&
        & 8.742053e-02_r8,8.263175e-02_r8,7.834058e-02_r8,7.447327e-02_r8,7.097000e-02_r8,&
        & 6.778167e-02_r8,6.486765e-02_r8,6.219400e-02_r8,5.973215e-02_r8,5.745790e-02_r8,&
        & 5.535059e-02_r8,5.339250e-02_r8,5.156831e-02_r8,4.975308e-02_r8,4.821235e-02_r8,&
        & 4.675596e-02_r8,4.537727e-02_r8,4.407030e-02_r8,4.282968e-02_r8,4.165053e-02_r8,&
        & 4.052845e-02_r8,3.945941e-02_r8,3.843980e-02_r8,3.746628e-02_r8,3.653583e-02_r8,&
        & 3.564567e-02_r8,3.479326e-02_r8,3.397626e-02_r8,3.319253e-02_r8,3.244008e-02_r8,&
        & 3.171711e-02_r8,3.102189e-02_r8,3.035289e-02_r8,2.970866e-02_r8,2.908784e-02_r8,&
        & 2.848920e-02_r8,2.791156e-02_r8,2.735385e-02_r8,2.681507e-02_r8,2.629425e-02_r8,&
        & 2.579053e-02_r8,2.530308e-02_r8,2.4831e-02_r8 /)
      extliq1(:, 29) = (/ &
        & 4.614710e-01_r8,4.556116e-01_r8,4.056568e-01_r8,3.529833e-01_r8,3.060334e-01_r8,&
        & 2.658127e-01_r8,2.316095e-01_r8,2.024325e-01_r8,1.773749e-01_r8,1.556867e-01_r8,&
        & 1.455558e-01_r8,1.332882e-01_r8,1.229052e-01_r8,1.140067e-01_r8,1.062981e-01_r8,&
        & 9.955703e-02_r8,9.361333e-02_r8,8.833420e-02_r8,8.361467e-02_r8,7.937071e-02_r8,&
        & 7.553420e-02_r8,7.204942e-02_r8,6.887031e-02_r8,6.595851e-02_r8,6.328178e-02_r8,&
        & 6.081286e-02_r8,5.852854e-02_r8,5.640892e-02_r8,5.431269e-02_r8,5.252561e-02_r8,&
        & 5.084345e-02_r8,4.925727e-02_r8,4.775910e-02_r8,4.634182e-02_r8,4.499907e-02_r8,&
        & 4.372512e-02_r8,4.251484e-02_r8,4.136357e-02_r8,4.026710e-02_r8,3.922162e-02_r8,&
        & 3.822365e-02_r8,3.727004e-02_r8,3.635790e-02_r8,3.548457e-02_r8,3.464764e-02_r8,&
        & 3.384488e-02_r8,3.307424e-02_r8,3.233384e-02_r8,3.162192e-02_r8,3.093688e-02_r8,&
        & 3.027723e-02_r8,2.964158e-02_r8,2.902864e-02_r8,2.843722e-02_r8,2.786621e-02_r8,&
        & 2.731457e-02_r8,2.678133e-02_r8,2.6266e-02_r8 /)

! Single scattering albedo     
      ssaliq1(:, 16) = (/ &
        & 8.143821e-01_r8,7.836739e-01_r8,7.550722e-01_r8,7.306269e-01_r8,7.105612e-01_r8,&
        & 6.946649e-01_r8,6.825556e-01_r8,6.737762e-01_r8,6.678448e-01_r8,6.642830e-01_r8,&
        & 6.679741e-01_r8,6.584607e-01_r8,6.505598e-01_r8,6.440951e-01_r8,6.388901e-01_r8,&
        & 6.347689e-01_r8,6.315549e-01_r8,6.290718e-01_r8,6.271432e-01_r8,6.255928e-01_r8,&
        & 6.242441e-01_r8,6.229207e-01_r8,6.214464e-01_r8,6.196445e-01_r8,6.173388e-01_r8,&
        & 6.143527e-01_r8,6.105099e-01_r8,6.056339e-01_r8,6.108290e-01_r8,6.073939e-01_r8,&
        & 6.043073e-01_r8,6.015473e-01_r8,5.990913e-01_r8,5.969173e-01_r8,5.950028e-01_r8,&
        & 5.933257e-01_r8,5.918636e-01_r8,5.905944e-01_r8,5.894957e-01_r8,5.885453e-01_r8,&
        & 5.877209e-01_r8,5.870003e-01_r8,5.863611e-01_r8,5.857811e-01_r8,5.852381e-01_r8,&
        & 5.847098e-01_r8,5.841738e-01_r8,5.836081e-01_r8,5.829901e-01_r8,5.822979e-01_r8,&
        & 5.815089e-01_r8,5.806011e-01_r8,5.795521e-01_r8,5.783396e-01_r8,5.769413e-01_r8,&
        & 5.753351e-01_r8,5.734986e-01_r8,5.7141e-01_r8 /)
      ssaliq1(:, 17) = (/ &
        & 8.165821e-01_r8,8.002015e-01_r8,7.816921e-01_r8,7.634131e-01_r8,7.463721e-01_r8,&
        & 7.312469e-01_r8,7.185883e-01_r8,7.088975e-01_r8,7.026671e-01_r8,7.004020e-01_r8,&
        & 7.042138e-01_r8,6.960930e-01_r8,6.894243e-01_r8,6.840459e-01_r8,6.797957e-01_r8,&
        & 6.765119e-01_r8,6.740325e-01_r8,6.721955e-01_r8,6.708391e-01_r8,6.698013e-01_r8,&
        & 6.689201e-01_r8,6.680339e-01_r8,6.669805e-01_r8,6.655982e-01_r8,6.637250e-01_r8,&
        & 6.611992e-01_r8,6.578588e-01_r8,6.535420e-01_r8,6.584449e-01_r8,6.553992e-01_r8,&
        & 6.526547e-01_r8,6.501917e-01_r8,6.479905e-01_r8,6.460313e-01_r8,6.442945e-01_r8,&
        & 6.427605e-01_r8,6.414094e-01_r8,6.402217e-01_r8,6.391775e-01_r8,6.382573e-01_r8,&
        & 6.374413e-01_r8,6.367099e-01_r8,6.360433e-01_r8,6.354218e-01_r8,6.348257e-01_r8,&
        & 6.342355e-01_r8,6.336313e-01_r8,6.329935e-01_r8,6.323023e-01_r8,6.315383e-01_r8,&
        & 6.306814e-01_r8,6.297122e-01_r8,6.286110e-01_r8,6.273579e-01_r8,6.259333e-01_r8,&
        & 6.243176e-01_r8,6.224910e-01_r8,6.2043e-01_r8 /)
      ssaliq1(:, 18) = (/ &
        & 9.900163e-01_r8,9.854307e-01_r8,9.797730e-01_r8,9.733113e-01_r8,9.664245e-01_r8,&
        & 9.594976e-01_r8,9.529055e-01_r8,9.470112e-01_r8,9.421695e-01_r8,9.387304e-01_r8,&
        & 9.344918e-01_r8,9.305302e-01_r8,9.267048e-01_r8,9.230072e-01_r8,9.194289e-01_r8,&
        & 9.159616e-01_r8,9.125968e-01_r8,9.093260e-01_r8,9.061409e-01_r8,9.030330e-01_r8,&
        & 8.999940e-01_r8,8.970154e-01_r8,8.940888e-01_r8,8.912058e-01_r8,8.883579e-01_r8,&
        & 8.855368e-01_r8,8.827341e-01_r8,8.799413e-01_r8,8.777423e-01_r8,8.749566e-01_r8,&
        & 8.722298e-01_r8,8.695605e-01_r8,8.669469e-01_r8,8.643875e-01_r8,8.618806e-01_r8,&
        & 8.594246e-01_r8,8.570179e-01_r8,8.546589e-01_r8,8.523459e-01_r8,8.500773e-01_r8,&
        & 8.478516e-01_r8,8.456670e-01_r8,8.435219e-01_r8,8.414148e-01_r8,8.393439e-01_r8,&
        & 8.373078e-01_r8,8.353047e-01_r8,8.333330e-01_r8,8.313911e-01_r8,8.294774e-01_r8,&
        & 8.275904e-01_r8,8.257282e-01_r8,8.238893e-01_r8,8.220721e-01_r8,8.202751e-01_r8,&
        & 8.184965e-01_r8,8.167346e-01_r8,8.1499e-01_r8 /)
      ssaliq1(:, 19) = (/ &
        & 9.999916e-01_r8,9.987396e-01_r8,9.966900e-01_r8,9.950738e-01_r8,9.937531e-01_r8,&
        & 9.925912e-01_r8,9.914525e-01_r8,9.902018e-01_r8,9.887046e-01_r8,9.868263e-01_r8,&
        & 9.849039e-01_r8,9.832372e-01_r8,9.815265e-01_r8,9.797770e-01_r8,9.779940e-01_r8,&
        & 9.761827e-01_r8,9.743481e-01_r8,9.724955e-01_r8,9.706303e-01_r8,9.687575e-01_r8,&
        & 9.668823e-01_r8,9.650100e-01_r8,9.631457e-01_r8,9.612947e-01_r8,9.594622e-01_r8,&
        & 9.576534e-01_r8,9.558734e-01_r8,9.541275e-01_r8,9.522059e-01_r8,9.504258e-01_r8,&
        & 9.486459e-01_r8,9.468676e-01_r8,9.450921e-01_r8,9.433208e-01_r8,9.415548e-01_r8,&
        & 9.397955e-01_r8,9.380441e-01_r8,9.363022e-01_r8,9.345706e-01_r8,9.328510e-01_r8,&
        & 9.311445e-01_r8,9.294524e-01_r8,9.277761e-01_r8,9.261167e-01_r8,9.244755e-01_r8,&
        & 9.228540e-01_r8,9.212534e-01_r8,9.196748e-01_r8,9.181197e-01_r8,9.165894e-01_r8,&
        & 9.150851e-01_r8,9.136080e-01_r8,9.121596e-01_r8,9.107410e-01_r8,9.093536e-01_r8,&
        & 9.079987e-01_r8,9.066775e-01_r8,9.0539e-01_r8 /)
      ssaliq1(:, 20) = (/ &
        & 9.979493e-01_r8,9.964113e-01_r8,9.950014e-01_r8,9.937045e-01_r8,9.924964e-01_r8,&
        & 9.913546e-01_r8,9.902575e-01_r8,9.891843e-01_r8,9.881136e-01_r8,9.870238e-01_r8,&
        & 9.859934e-01_r8,9.849372e-01_r8,9.838873e-01_r8,9.828434e-01_r8,9.818052e-01_r8,&
        & 9.807725e-01_r8,9.797450e-01_r8,9.787225e-01_r8,9.777047e-01_r8,9.766914e-01_r8,&
        & 9.756823e-01_r8,9.746771e-01_r8,9.736756e-01_r8,9.726775e-01_r8,9.716827e-01_r8,&
        & 9.706907e-01_r8,9.697014e-01_r8,9.687145e-01_r8,9.678060e-01_r8,9.668108e-01_r8,&
        & 9.658218e-01_r8,9.648391e-01_r8,9.638629e-01_r8,9.628936e-01_r8,9.619313e-01_r8,&
        & 9.609763e-01_r8,9.600287e-01_r8,9.590888e-01_r8,9.581569e-01_r8,9.572330e-01_r8,&
        & 9.563176e-01_r8,9.554108e-01_r8,9.545128e-01_r8,9.536239e-01_r8,9.527443e-01_r8,&
        & 9.518741e-01_r8,9.510137e-01_r8,9.501633e-01_r8,9.493230e-01_r8,9.484931e-01_r8,&
        & 9.476740e-01_r8,9.468656e-01_r8,9.460683e-01_r8,9.452824e-01_r8,9.445080e-01_r8,&
        & 9.437454e-01_r8,9.429948e-01_r8,9.4226e-01_r8 /)
      ssaliq1(:, 21) = (/ &
        & 9.988742e-01_r8,9.982668e-01_r8,9.976935e-01_r8,9.971497e-01_r8,9.966314e-01_r8,&
        & 9.961344e-01_r8,9.956545e-01_r8,9.951873e-01_r8,9.947286e-01_r8,9.942741e-01_r8,&
        & 9.938457e-01_r8,9.933947e-01_r8,9.929473e-01_r8,9.925032e-01_r8,9.920621e-01_r8,&
        & 9.916237e-01_r8,9.911875e-01_r8,9.907534e-01_r8,9.903209e-01_r8,9.898898e-01_r8,&
        & 9.894597e-01_r8,9.890304e-01_r8,9.886015e-01_r8,9.881726e-01_r8,9.877435e-01_r8,&
        & 9.873138e-01_r8,9.868833e-01_r8,9.864516e-01_r8,9.860698e-01_r8,9.856317e-01_r8,&
        & 9.851957e-01_r8,9.847618e-01_r8,9.843302e-01_r8,9.839008e-01_r8,9.834739e-01_r8,&
        & 9.830494e-01_r8,9.826275e-01_r8,9.822083e-01_r8,9.817918e-01_r8,9.813782e-01_r8,&
        & 9.809675e-01_r8,9.805598e-01_r8,9.801552e-01_r8,9.797538e-01_r8,9.793556e-01_r8,&
        & 9.789608e-01_r8,9.785695e-01_r8,9.781817e-01_r8,9.777975e-01_r8,9.774171e-01_r8,&
        & 9.770404e-01_r8,9.766676e-01_r8,9.762988e-01_r8,9.759340e-01_r8,9.755733e-01_r8,&
        & 9.752169e-01_r8,9.748649e-01_r8,9.7452e-01_r8 /)
      ssaliq1(:, 22) = (/ &
        & 9.994441e-01_r8,9.991608e-01_r8,9.988949e-01_r8,9.986439e-01_r8,9.984054e-01_r8,&
        & 9.981768e-01_r8,9.979557e-01_r8,9.977396e-01_r8,9.975258e-01_r8,9.973120e-01_r8,&
        & 9.971011e-01_r8,9.968852e-01_r8,9.966708e-01_r8,9.964578e-01_r8,9.962462e-01_r8,&
        & 9.960357e-01_r8,9.958264e-01_r8,9.956181e-01_r8,9.954108e-01_r8,9.952043e-01_r8,&
        & 9.949987e-01_r8,9.947937e-01_r8,9.945892e-01_r8,9.943853e-01_r8,9.941818e-01_r8,&
        & 9.939786e-01_r8,9.937757e-01_r8,9.935728e-01_r8,9.933922e-01_r8,9.931825e-01_r8,&
        & 9.929739e-01_r8,9.927661e-01_r8,9.925592e-01_r8,9.923534e-01_r8,9.921485e-01_r8,&
        & 9.919447e-01_r8,9.917421e-01_r8,9.915406e-01_r8,9.913403e-01_r8,9.911412e-01_r8,&
        & 9.909435e-01_r8,9.907470e-01_r8,9.905519e-01_r8,9.903581e-01_r8,9.901659e-01_r8,&
        & 9.899751e-01_r8,9.897858e-01_r8,9.895981e-01_r8,9.894120e-01_r8,9.892276e-01_r8,&
        & 9.890447e-01_r8,9.888637e-01_r8,9.886845e-01_r8,9.885070e-01_r8,9.883314e-01_r8,&
        & 9.881576e-01_r8,9.879859e-01_r8,9.8782e-01_r8 /)
      ssaliq1(:, 23) = (/ &
        & 9.999138e-01_r8,9.998730e-01_r8,9.998338e-01_r8,9.997965e-01_r8,9.997609e-01_r8,&
        & 9.997270e-01_r8,9.996944e-01_r8,9.996629e-01_r8,9.996321e-01_r8,9.996016e-01_r8,&
        & 9.995690e-01_r8,9.995372e-01_r8,9.995057e-01_r8,9.994744e-01_r8,9.994433e-01_r8,&
        & 9.994124e-01_r8,9.993817e-01_r8,9.993510e-01_r8,9.993206e-01_r8,9.992903e-01_r8,&
        & 9.992600e-01_r8,9.992299e-01_r8,9.991998e-01_r8,9.991698e-01_r8,9.991398e-01_r8,&
        & 9.991098e-01_r8,9.990799e-01_r8,9.990499e-01_r8,9.990231e-01_r8,9.989920e-01_r8,&
        & 9.989611e-01_r8,9.989302e-01_r8,9.988996e-01_r8,9.988690e-01_r8,9.988386e-01_r8,&
        & 9.988084e-01_r8,9.987783e-01_r8,9.987485e-01_r8,9.987187e-01_r8,9.986891e-01_r8,&
        & 9.986598e-01_r8,9.986306e-01_r8,9.986017e-01_r8,9.985729e-01_r8,9.985443e-01_r8,&
        & 9.985160e-01_r8,9.984879e-01_r8,9.984600e-01_r8,9.984324e-01_r8,9.984050e-01_r8,&
        & 9.983778e-01_r8,9.983509e-01_r8,9.983243e-01_r8,9.982980e-01_r8,9.982719e-01_r8,&
        & 9.982461e-01_r8,9.982206e-01_r8,9.9820e-01_r8 /)
      ssaliq1(:, 24) = (/ &
        & 9.999985e-01_r8,9.999979e-01_r8,9.999972e-01_r8,9.999966e-01_r8,9.999961e-01_r8,&
        & 9.999955e-01_r8,9.999950e-01_r8,9.999944e-01_r8,9.999938e-01_r8,9.999933e-01_r8,&
        & 9.999927e-01_r8,9.999921e-01_r8,9.999915e-01_r8,9.999910e-01_r8,9.999904e-01_r8,&
        & 9.999899e-01_r8,9.999893e-01_r8,9.999888e-01_r8,9.999882e-01_r8,9.999877e-01_r8,&
        & 9.999871e-01_r8,9.999866e-01_r8,9.999861e-01_r8,9.999855e-01_r8,9.999850e-01_r8,&
        & 9.999844e-01_r8,9.999839e-01_r8,9.999833e-01_r8,9.999828e-01_r8,9.999823e-01_r8,&
        & 9.999817e-01_r8,9.999812e-01_r8,9.999807e-01_r8,9.999801e-01_r8,9.999796e-01_r8,&
        & 9.999791e-01_r8,9.999786e-01_r8,9.999781e-01_r8,9.999776e-01_r8,9.999770e-01_r8,&
        & 9.999765e-01_r8,9.999761e-01_r8,9.999756e-01_r8,9.999751e-01_r8,9.999746e-01_r8,&
        & 9.999741e-01_r8,9.999736e-01_r8,9.999732e-01_r8,9.999727e-01_r8,9.999722e-01_r8,&
        & 9.999718e-01_r8,9.999713e-01_r8,9.999709e-01_r8,9.999705e-01_r8,9.999701e-01_r8,&
        & 9.999697e-01_r8,9.999692e-01_r8,9.9997e-01_r8 /)
      ssaliq1(:, 25) = (/ &
        & 9.999999e-01_r8,9.999998e-01_r8,9.999997e-01_r8,9.999997e-01_r8,9.999997e-01_r8,&
        & 9.999996e-01_r8,9.999996e-01_r8,9.999995e-01_r8,9.999995e-01_r8,9.999994e-01_r8,&
        & 9.999994e-01_r8,9.999993e-01_r8,9.999993e-01_r8,9.999992e-01_r8,9.999992e-01_r8,&
        & 9.999991e-01_r8,9.999991e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999989e-01_r8,9.999989e-01_r8,9.999988e-01_r8,9.999988e-01_r8,9.999987e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999986e-01_r8,9.999985e-01_r8,9.999985e-01_r8,&
        & 9.999984e-01_r8,9.999984e-01_r8,9.999984e-01_r8,9.999983e-01_r8,9.999983e-01_r8,&
        & 9.999982e-01_r8,9.999982e-01_r8,9.999982e-01_r8,9.999981e-01_r8,9.999980e-01_r8,&
        & 9.999980e-01_r8,9.999980e-01_r8,9.999979e-01_r8,9.999979e-01_r8,9.999978e-01_r8,&
        & 9.999978e-01_r8,9.999977e-01_r8,9.999977e-01_r8,9.999977e-01_r8,9.999976e-01_r8,&
        & 9.999976e-01_r8,9.999975e-01_r8,9.999975e-01_r8,9.999974e-01_r8,9.999974e-01_r8,&
        & 9.999974e-01_r8,9.999973e-01_r8,1.0000e+00_r8 /)
      ssaliq1(:, 26) = (/ &
        & 9.999997e-01_r8,9.999995e-01_r8,9.999993e-01_r8,9.999992e-01_r8,9.999990e-01_r8,&
        & 9.999989e-01_r8,9.999988e-01_r8,9.999987e-01_r8,9.999986e-01_r8,9.999985e-01_r8,&
        & 9.999984e-01_r8,9.999983e-01_r8,9.999982e-01_r8,9.999981e-01_r8,9.999980e-01_r8,&
        & 9.999978e-01_r8,9.999977e-01_r8,9.999976e-01_r8,9.999975e-01_r8,9.999974e-01_r8,&
        & 9.999973e-01_r8,9.999972e-01_r8,9.999970e-01_r8,9.999969e-01_r8,9.999968e-01_r8,&
        & 9.999967e-01_r8,9.999966e-01_r8,9.999965e-01_r8,9.999964e-01_r8,9.999963e-01_r8,&
        & 9.999962e-01_r8,9.999961e-01_r8,9.999959e-01_r8,9.999958e-01_r8,9.999957e-01_r8,&
        & 9.999956e-01_r8,9.999955e-01_r8,9.999954e-01_r8,9.999953e-01_r8,9.999952e-01_r8,&
        & 9.999951e-01_r8,9.999949e-01_r8,9.999949e-01_r8,9.999947e-01_r8,9.999946e-01_r8,&
        & 9.999945e-01_r8,9.999944e-01_r8,9.999943e-01_r8,9.999942e-01_r8,9.999941e-01_r8,&
        & 9.999940e-01_r8,9.999939e-01_r8,9.999938e-01_r8,9.999937e-01_r8,9.999936e-01_r8,&
        & 9.999935e-01_r8,9.999934e-01_r8,9.9999e-01_r8 /)
      ssaliq1(:, 27) = (/ &
        & 9.999984e-01_r8,9.999976e-01_r8,9.999969e-01_r8,9.999962e-01_r8,9.999956e-01_r8,&
        & 9.999950e-01_r8,9.999945e-01_r8,9.999940e-01_r8,9.999935e-01_r8,9.999931e-01_r8,&
        & 9.999926e-01_r8,9.999920e-01_r8,9.999914e-01_r8,9.999908e-01_r8,9.999903e-01_r8,&
        & 9.999897e-01_r8,9.999891e-01_r8,9.999886e-01_r8,9.999880e-01_r8,9.999874e-01_r8,&
        & 9.999868e-01_r8,9.999863e-01_r8,9.999857e-01_r8,9.999851e-01_r8,9.999846e-01_r8,&
        & 9.999840e-01_r8,9.999835e-01_r8,9.999829e-01_r8,9.999824e-01_r8,9.999818e-01_r8,&
        & 9.999812e-01_r8,9.999806e-01_r8,9.999800e-01_r8,9.999795e-01_r8,9.999789e-01_r8,&
        & 9.999783e-01_r8,9.999778e-01_r8,9.999773e-01_r8,9.999767e-01_r8,9.999761e-01_r8,&
        & 9.999756e-01_r8,9.999750e-01_r8,9.999745e-01_r8,9.999739e-01_r8,9.999734e-01_r8,&
        & 9.999729e-01_r8,9.999723e-01_r8,9.999718e-01_r8,9.999713e-01_r8,9.999708e-01_r8,&
        & 9.999703e-01_r8,9.999697e-01_r8,9.999692e-01_r8,9.999687e-01_r8,9.999683e-01_r8,&
        & 9.999678e-01_r8,9.999673e-01_r8,9.9997e-01_r8 /)
      ssaliq1(:, 28) = (/ &
        & 9.999981e-01_r8,9.999973e-01_r8,9.999965e-01_r8,9.999958e-01_r8,9.999951e-01_r8,&
        & 9.999943e-01_r8,9.999937e-01_r8,9.999930e-01_r8,9.999924e-01_r8,9.999918e-01_r8,&
        & 9.999912e-01_r8,9.999905e-01_r8,9.999897e-01_r8,9.999890e-01_r8,9.999883e-01_r8,&
        & 9.999876e-01_r8,9.999869e-01_r8,9.999862e-01_r8,9.999855e-01_r8,9.999847e-01_r8,&
        & 9.999840e-01_r8,9.999834e-01_r8,9.999827e-01_r8,9.999819e-01_r8,9.999812e-01_r8,&
        & 9.999805e-01_r8,9.999799e-01_r8,9.999791e-01_r8,9.999785e-01_r8,9.999778e-01_r8,&
        & 9.999771e-01_r8,9.999764e-01_r8,9.999757e-01_r8,9.999750e-01_r8,9.999743e-01_r8,&
        & 9.999736e-01_r8,9.999729e-01_r8,9.999722e-01_r8,9.999715e-01_r8,9.999709e-01_r8,&
        & 9.999701e-01_r8,9.999695e-01_r8,9.999688e-01_r8,9.999682e-01_r8,9.999675e-01_r8,&
        & 9.999669e-01_r8,9.999662e-01_r8,9.999655e-01_r8,9.999649e-01_r8,9.999642e-01_r8,&
        & 9.999636e-01_r8,9.999630e-01_r8,9.999624e-01_r8,9.999618e-01_r8,9.999612e-01_r8,&
        & 9.999606e-01_r8,9.999600e-01_r8,9.9996e-01_r8 /)
      ssaliq1(:, 29) = (/ &
        & 8.505737e-01_r8,8.465102e-01_r8,8.394829e-01_r8,8.279508e-01_r8,8.110806e-01_r8,&
        & 7.900397e-01_r8,7.669615e-01_r8,7.444422e-01_r8,7.253055e-01_r8,7.124831e-01_r8,&
        & 7.016434e-01_r8,6.885485e-01_r8,6.767340e-01_r8,6.661029e-01_r8,6.565577e-01_r8,&
        & 6.480013e-01_r8,6.403373e-01_r8,6.334697e-01_r8,6.273034e-01_r8,6.217440e-01_r8,&
        & 6.166983e-01_r8,6.120740e-01_r8,6.077796e-01_r8,6.037249e-01_r8,5.998207e-01_r8,&
        & 5.959788e-01_r8,5.921123e-01_r8,5.881354e-01_r8,5.891285e-01_r8,5.851143e-01_r8,&
        & 5.814653e-01_r8,5.781606e-01_r8,5.751792e-01_r8,5.724998e-01_r8,5.701016e-01_r8,&
        & 5.679634e-01_r8,5.660642e-01_r8,5.643829e-01_r8,5.628984e-01_r8,5.615898e-01_r8,&
        & 5.604359e-01_r8,5.594158e-01_r8,5.585083e-01_r8,5.576924e-01_r8,5.569470e-01_r8,&
        & 5.562512e-01_r8,5.555838e-01_r8,5.549239e-01_r8,5.542503e-01_r8,5.535420e-01_r8,&
        & 5.527781e-01_r8,5.519374e-01_r8,5.509989e-01_r8,5.499417e-01_r8,5.487445e-01_r8,&
        & 5.473865e-01_r8,5.458466e-01_r8,5.4410e-01_r8 /)

! asymmetry parameter
      asyliq1(:, 16) = (/ &
        & 8.133297e-01_r8,8.133528e-01_r8,8.173865e-01_r8,8.243205e-01_r8,8.333063e-01_r8,&
        & 8.436317e-01_r8,8.546611e-01_r8,8.657934e-01_r8,8.764345e-01_r8,8.859837e-01_r8,&
        & 8.627394e-01_r8,8.824569e-01_r8,8.976887e-01_r8,9.089541e-01_r8,9.167699e-01_r8,&
        & 9.216517e-01_r8,9.241147e-01_r8,9.246743e-01_r8,9.238469e-01_r8,9.221504e-01_r8,&
        & 9.201045e-01_r8,9.182299e-01_r8,9.170491e-01_r8,9.170862e-01_r8,9.188653e-01_r8,&
        & 9.229111e-01_r8,9.297468e-01_r8,9.398950e-01_r8,9.203269e-01_r8,9.260693e-01_r8,&
        & 9.309373e-01_r8,9.349918e-01_r8,9.382935e-01_r8,9.409030e-01_r8,9.428809e-01_r8,&
        & 9.442881e-01_r8,9.451851e-01_r8,9.456331e-01_r8,9.456926e-01_r8,9.454247e-01_r8,&
        & 9.448902e-01_r8,9.441503e-01_r8,9.432661e-01_r8,9.422987e-01_r8,9.413094e-01_r8,&
        & 9.403594e-01_r8,9.395102e-01_r8,9.388230e-01_r8,9.383594e-01_r8,9.381810e-01_r8,&
        & 9.383489e-01_r8,9.389251e-01_r8,9.399707e-01_r8,9.415475e-01_r8,9.437167e-01_r8,&
        & 9.465399e-01_r8,9.500786e-01_r8,9.5439e-01_r8 /)
      asyliq1(:, 17) = (/ &
        & 8.794448e-01_r8,8.819306e-01_r8,8.837667e-01_r8,8.853832e-01_r8,8.871010e-01_r8,&
        & 8.892675e-01_r8,8.922584e-01_r8,8.964666e-01_r8,9.022940e-01_r8,9.101456e-01_r8,&
        & 8.839999e-01_r8,9.035610e-01_r8,9.184568e-01_r8,9.292315e-01_r8,9.364282e-01_r8,&
        & 9.405887e-01_r8,9.422554e-01_r8,9.419703e-01_r8,9.402759e-01_r8,9.377159e-01_r8,&
        & 9.348345e-01_r8,9.321769e-01_r8,9.302888e-01_r8,9.297166e-01_r8,9.310075e-01_r8,&
        & 9.347080e-01_r8,9.413643e-01_r8,9.515216e-01_r8,9.306286e-01_r8,9.361781e-01_r8,&
        & 9.408374e-01_r8,9.446692e-01_r8,9.477363e-01_r8,9.501013e-01_r8,9.518268e-01_r8,&
        & 9.529756e-01_r8,9.536105e-01_r8,9.537938e-01_r8,9.535886e-01_r8,9.530574e-01_r8,&
        & 9.522633e-01_r8,9.512688e-01_r8,9.501370e-01_r8,9.489306e-01_r8,9.477126e-01_r8,&
        & 9.465459e-01_r8,9.454934e-01_r8,9.446183e-01_r8,9.439833e-01_r8,9.436519e-01_r8,&
        & 9.436866e-01_r8,9.441508e-01_r8,9.451073e-01_r8,9.466195e-01_r8,9.487501e-01_r8,&
        & 9.515621e-01_r8,9.551185e-01_r8,9.5948e-01_r8 /)
      asyliq1(:, 18) = (/ &
        & 8.478817e-01_r8,8.269312e-01_r8,8.161352e-01_r8,8.135960e-01_r8,8.173586e-01_r8,&
        & 8.254167e-01_r8,8.357072e-01_r8,8.461167e-01_r8,8.544952e-01_r8,8.586776e-01_r8,&
        & 8.335562e-01_r8,8.524273e-01_r8,8.669052e-01_r8,8.775014e-01_r8,8.847277e-01_r8,&
        & 8.890958e-01_r8,8.911173e-01_r8,8.913038e-01_r8,8.901669e-01_r8,8.882182e-01_r8,&
        & 8.859692e-01_r8,8.839315e-01_r8,8.826164e-01_r8,8.825356e-01_r8,8.842004e-01_r8,&
        & 8.881223e-01_r8,8.948131e-01_r8,9.047837e-01_r8,8.855951e-01_r8,8.911796e-01_r8,&
        & 8.959229e-01_r8,8.998837e-01_r8,9.031209e-01_r8,9.056939e-01_r8,9.076609e-01_r8,&
        & 9.090812e-01_r8,9.100134e-01_r8,9.105167e-01_r8,9.106496e-01_r8,9.104712e-01_r8,&
        & 9.100404e-01_r8,9.094159e-01_r8,9.086568e-01_r8,9.078218e-01_r8,9.069697e-01_r8,&
        & 9.061595e-01_r8,9.054499e-01_r8,9.048999e-01_r8,9.045683e-01_r8,9.045142e-01_r8,&
        & 9.047962e-01_r8,9.054730e-01_r8,9.066037e-01_r8,9.082472e-01_r8,9.104623e-01_r8,&
        & 9.133079e-01_r8,9.168427e-01_r8,9.2113e-01_r8 /)
      asyliq1(:, 19) = (/ &
        & 8.216697e-01_r8,7.982871e-01_r8,7.891147e-01_r8,7.909083e-01_r8,8.003833e-01_r8,&
        & 8.142516e-01_r8,8.292290e-01_r8,8.420356e-01_r8,8.493945e-01_r8,8.480316e-01_r8,&
        & 8.212381e-01_r8,8.394984e-01_r8,8.534095e-01_r8,8.634813e-01_r8,8.702242e-01_r8,&
        & 8.741483e-01_r8,8.757638e-01_r8,8.755808e-01_r8,8.741095e-01_r8,8.718604e-01_r8,&
        & 8.693433e-01_r8,8.670686e-01_r8,8.655464e-01_r8,8.652872e-01_r8,8.668006e-01_r8,&
        & 8.705973e-01_r8,8.771874e-01_r8,8.870809e-01_r8,8.678284e-01_r8,8.732315e-01_r8,&
        & 8.778084e-01_r8,8.816166e-01_r8,8.847146e-01_r8,8.871603e-01_r8,8.890116e-01_r8,&
        & 8.903266e-01_r8,8.911632e-01_r8,8.915796e-01_r8,8.916337e-01_r8,8.913834e-01_r8,&
        & 8.908869e-01_r8,8.902022e-01_r8,8.893873e-01_r8,8.885001e-01_r8,8.875986e-01_r8,&
        & 8.867411e-01_r8,8.859852e-01_r8,8.853891e-01_r8,8.850111e-01_r8,8.849089e-01_r8,&
        & 8.851405e-01_r8,8.857639e-01_r8,8.868372e-01_r8,8.884185e-01_r8,8.905656e-01_r8,&
        & 8.933368e-01_r8,8.967899e-01_r8,9.0098e-01_r8 /)
      asyliq1(:, 20) = (/ &
        & 8.063610e-01_r8,7.938147e-01_r8,7.921304e-01_r8,7.985092e-01_r8,8.101339e-01_r8,&
        & 8.242175e-01_r8,8.379913e-01_r8,8.486920e-01_r8,8.535547e-01_r8,8.498083e-01_r8,&
        & 8.224849e-01_r8,8.405509e-01_r8,8.542436e-01_r8,8.640770e-01_r8,8.705653e-01_r8,&
        & 8.742227e-01_r8,8.755630e-01_r8,8.751004e-01_r8,8.733491e-01_r8,8.708231e-01_r8,&
        & 8.680365e-01_r8,8.655035e-01_r8,8.637381e-01_r8,8.632544e-01_r8,8.645665e-01_r8,&
        & 8.681885e-01_r8,8.746346e-01_r8,8.844188e-01_r8,8.648180e-01_r8,8.700563e-01_r8,&
        & 8.744672e-01_r8,8.781087e-01_r8,8.810393e-01_r8,8.833174e-01_r8,8.850011e-01_r8,&
        & 8.861485e-01_r8,8.868183e-01_r8,8.870687e-01_r8,8.869579e-01_r8,8.865441e-01_r8,&
        & 8.858857e-01_r8,8.850412e-01_r8,8.840686e-01_r8,8.830263e-01_r8,8.819726e-01_r8,&
        & 8.809658e-01_r8,8.800642e-01_r8,8.793260e-01_r8,8.788099e-01_r8,8.785737e-01_r8,&
        & 8.786758e-01_r8,8.791746e-01_r8,8.801283e-01_r8,8.815955e-01_r8,8.836340e-01_r8,&
        & 8.863024e-01_r8,8.896592e-01_r8,8.9376e-01_r8 /)
      asyliq1(:, 21) = (/ &
        & 7.885899e-01_r8,7.937172e-01_r8,8.020658e-01_r8,8.123971e-01_r8,8.235502e-01_r8,&
        & 8.343776e-01_r8,8.437336e-01_r8,8.504711e-01_r8,8.534421e-01_r8,8.514978e-01_r8,&
        & 8.238888e-01_r8,8.417463e-01_r8,8.552057e-01_r8,8.647853e-01_r8,8.710038e-01_r8,&
        & 8.743798e-01_r8,8.754319e-01_r8,8.746786e-01_r8,8.726386e-01_r8,8.698303e-01_r8,&
        & 8.667724e-01_r8,8.639836e-01_r8,8.619823e-01_r8,8.612870e-01_r8,8.624165e-01_r8,&
        & 8.658893e-01_r8,8.722241e-01_r8,8.819394e-01_r8,8.620216e-01_r8,8.671239e-01_r8,&
        & 8.713983e-01_r8,8.749032e-01_r8,8.776970e-01_r8,8.798385e-01_r8,8.813860e-01_r8,&
        & 8.823980e-01_r8,8.829332e-01_r8,8.830500e-01_r8,8.828068e-01_r8,8.822623e-01_r8,&
        & 8.814750e-01_r8,8.805031e-01_r8,8.794056e-01_r8,8.782407e-01_r8,8.770672e-01_r8,&
        & 8.759432e-01_r8,8.749275e-01_r8,8.740784e-01_r8,8.734547e-01_r8,8.731146e-01_r8,&
        & 8.731170e-01_r8,8.735199e-01_r8,8.743823e-01_r8,8.757625e-01_r8,8.777191e-01_r8,&
        & 8.803105e-01_r8,8.835953e-01_r8,8.8763e-01_r8 /)
      asyliq1(:, 22) = (/ &
        & 7.811516e-01_r8,7.962229e-01_r8,8.096199e-01_r8,8.212996e-01_r8,8.312212e-01_r8,&
        & 8.393430e-01_r8,8.456236e-01_r8,8.500214e-01_r8,8.524950e-01_r8,8.530031e-01_r8,&
        & 8.251485e-01_r8,8.429043e-01_r8,8.562461e-01_r8,8.656954e-01_r8,8.717737e-01_r8,&
        & 8.750020e-01_r8,8.759022e-01_r8,8.749953e-01_r8,8.728027e-01_r8,8.698461e-01_r8,&
        & 8.666466e-01_r8,8.637257e-01_r8,8.616047e-01_r8,8.608051e-01_r8,8.618483e-01_r8,&
        & 8.652557e-01_r8,8.715487e-01_r8,8.812485e-01_r8,8.611645e-01_r8,8.662052e-01_r8,&
        & 8.704173e-01_r8,8.738594e-01_r8,8.765901e-01_r8,8.786678e-01_r8,8.801517e-01_r8,&
        & 8.810999e-01_r8,8.815713e-01_r8,8.816246e-01_r8,8.813185e-01_r8,8.807114e-01_r8,&
        & 8.798621e-01_r8,8.788290e-01_r8,8.776713e-01_r8,8.764470e-01_r8,8.752152e-01_r8,&
        & 8.740343e-01_r8,8.729631e-01_r8,8.720602e-01_r8,8.713842e-01_r8,8.709936e-01_r8,&
        & 8.709475e-01_r8,8.713041e-01_r8,8.721221e-01_r8,8.734602e-01_r8,8.753774e-01_r8,&
        & 8.779319e-01_r8,8.811825e-01_r8,8.8519e-01_r8 /)
      asyliq1(:, 23) = (/ &
        & 7.865744e-01_r8,8.093340e-01_r8,8.257596e-01_r8,8.369940e-01_r8,8.441574e-01_r8,&
        & 8.483602e-01_r8,8.507096e-01_r8,8.523139e-01_r8,8.542834e-01_r8,8.577321e-01_r8,&
        & 8.288960e-01_r8,8.465308e-01_r8,8.597175e-01_r8,8.689830e-01_r8,8.748542e-01_r8,&
        & 8.778584e-01_r8,8.785222e-01_r8,8.773728e-01_r8,8.749370e-01_r8,8.717419e-01_r8,&
        & 8.683145e-01_r8,8.651816e-01_r8,8.628704e-01_r8,8.619077e-01_r8,8.628205e-01_r8,&
        & 8.661356e-01_r8,8.723803e-01_r8,8.820815e-01_r8,8.616715e-01_r8,8.666389e-01_r8,&
        & 8.707753e-01_r8,8.741398e-01_r8,8.767912e-01_r8,8.787885e-01_r8,8.801908e-01_r8,&
        & 8.810570e-01_r8,8.814460e-01_r8,8.814167e-01_r8,8.810283e-01_r8,8.803395e-01_r8,&
        & 8.794095e-01_r8,8.782971e-01_r8,8.770613e-01_r8,8.757610e-01_r8,8.744553e-01_r8,&
        & 8.732031e-01_r8,8.720634e-01_r8,8.710951e-01_r8,8.703572e-01_r8,8.699086e-01_r8,&
        & 8.698084e-01_r8,8.701155e-01_r8,8.708887e-01_r8,8.721872e-01_r8,8.740698e-01_r8,&
        & 8.765957e-01_r8,8.798235e-01_r8,8.8381e-01_r8 /)
      asyliq1(:, 24) = (/ &
        & 8.069513e-01_r8,8.262939e-01_r8,8.398241e-01_r8,8.486352e-01_r8,8.538213e-01_r8,&
        & 8.564743e-01_r8,8.576854e-01_r8,8.585455e-01_r8,8.601452e-01_r8,8.635755e-01_r8,&
        & 8.337383e-01_r8,8.512655e-01_r8,8.643049e-01_r8,8.733896e-01_r8,8.790535e-01_r8,&
        & 8.818295e-01_r8,8.822518e-01_r8,8.808533e-01_r8,8.781676e-01_r8,8.747284e-01_r8,&
        & 8.710690e-01_r8,8.677229e-01_r8,8.652236e-01_r8,8.641047e-01_r8,8.648993e-01_r8,&
        & 8.681413e-01_r8,8.743640e-01_r8,8.841007e-01_r8,8.633558e-01_r8,8.682719e-01_r8,&
        & 8.723543e-01_r8,8.756621e-01_r8,8.782547e-01_r8,8.801915e-01_r8,8.815318e-01_r8,&
        & 8.823347e-01_r8,8.826598e-01_r8,8.825663e-01_r8,8.821135e-01_r8,8.813608e-01_r8,&
        & 8.803674e-01_r8,8.791928e-01_r8,8.778960e-01_r8,8.765366e-01_r8,8.751738e-01_r8,&
        & 8.738670e-01_r8,8.726755e-01_r8,8.716585e-01_r8,8.708755e-01_r8,8.703856e-01_r8,&
        & 8.702483e-01_r8,8.705229e-01_r8,8.712687e-01_r8,8.725448e-01_r8,8.744109e-01_r8,&
        & 8.769260e-01_r8,8.801496e-01_r8,8.8414e-01_r8 /)
      asyliq1(:, 25) = (/ &
        & 8.252182e-01_r8,8.379244e-01_r8,8.471709e-01_r8,8.535760e-01_r8,8.577540e-01_r8,&
        & 8.603183e-01_r8,8.618820e-01_r8,8.630578e-01_r8,8.644587e-01_r8,8.666970e-01_r8,&
        & 8.362159e-01_r8,8.536817e-01_r8,8.666387e-01_r8,8.756240e-01_r8,8.811746e-01_r8,&
        & 8.838273e-01_r8,8.841191e-01_r8,8.825871e-01_r8,8.797681e-01_r8,8.761992e-01_r8,&
        & 8.724174e-01_r8,8.689593e-01_r8,8.663623e-01_r8,8.651632e-01_r8,8.658988e-01_r8,&
        & 8.691064e-01_r8,8.753226e-01_r8,8.850847e-01_r8,8.641620e-01_r8,8.690500e-01_r8,&
        & 8.731026e-01_r8,8.763795e-01_r8,8.789400e-01_r8,8.808438e-01_r8,8.821503e-01_r8,&
        & 8.829191e-01_r8,8.832095e-01_r8,8.830813e-01_r8,8.825938e-01_r8,8.818064e-01_r8,&
        & 8.807787e-01_r8,8.795704e-01_r8,8.782408e-01_r8,8.768493e-01_r8,8.754557e-01_r8,&
        & 8.741193e-01_r8,8.728995e-01_r8,8.718561e-01_r8,8.710484e-01_r8,8.705360e-01_r8,&
        & 8.703782e-01_r8,8.706347e-01_r8,8.713650e-01_r8,8.726285e-01_r8,8.744849e-01_r8,&
        & 8.769933e-01_r8,8.802136e-01_r8,8.8421e-01_r8 /)
      asyliq1(:, 26) = (/ &
        & 8.370583e-01_r8,8.467920e-01_r8,8.537769e-01_r8,8.585136e-01_r8,8.615034e-01_r8,&
        & 8.632474e-01_r8,8.642468e-01_r8,8.650026e-01_r8,8.660161e-01_r8,8.677882e-01_r8,&
        & 8.369760e-01_r8,8.543821e-01_r8,8.672699e-01_r8,8.761782e-01_r8,8.816454e-01_r8,&
        & 8.842103e-01_r8,8.844114e-01_r8,8.827872e-01_r8,8.798766e-01_r8,8.762179e-01_r8,&
        & 8.723500e-01_r8,8.688112e-01_r8,8.661403e-01_r8,8.648758e-01_r8,8.655563e-01_r8,&
        & 8.687206e-01_r8,8.749072e-01_r8,8.846546e-01_r8,8.636289e-01_r8,8.684849e-01_r8,&
        & 8.725054e-01_r8,8.757501e-01_r8,8.782785e-01_r8,8.801503e-01_r8,8.814249e-01_r8,&
        & 8.821620e-01_r8,8.824211e-01_r8,8.822620e-01_r8,8.817440e-01_r8,8.809268e-01_r8,&
        & 8.798699e-01_r8,8.786330e-01_r8,8.772756e-01_r8,8.758572e-01_r8,8.744374e-01_r8,&
        & 8.730760e-01_r8,8.718323e-01_r8,8.707660e-01_r8,8.699366e-01_r8,8.694039e-01_r8,&
        & 8.692271e-01_r8,8.694661e-01_r8,8.701803e-01_r8,8.714293e-01_r8,8.732727e-01_r8,&
        & 8.757702e-01_r8,8.789811e-01_r8,8.8297e-01_r8 /)
      asyliq1(:, 27) = (/ &
        & 8.430819e-01_r8,8.510060e-01_r8,8.567270e-01_r8,8.606533e-01_r8,8.631934e-01_r8,&
        & 8.647554e-01_r8,8.657471e-01_r8,8.665760e-01_r8,8.676496e-01_r8,8.693754e-01_r8,&
        & 8.384298e-01_r8,8.557913e-01_r8,8.686214e-01_r8,8.774605e-01_r8,8.828495e-01_r8,&
        & 8.853287e-01_r8,8.854393e-01_r8,8.837215e-01_r8,8.807161e-01_r8,8.769639e-01_r8,&
        & 8.730053e-01_r8,8.693812e-01_r8,8.666321e-01_r8,8.652988e-01_r8,8.659219e-01_r8,&
        & 8.690419e-01_r8,8.751999e-01_r8,8.849360e-01_r8,8.638013e-01_r8,8.686371e-01_r8,&
        & 8.726369e-01_r8,8.758605e-01_r8,8.783674e-01_r8,8.802176e-01_r8,8.814705e-01_r8,&
        & 8.821859e-01_r8,8.824234e-01_r8,8.822429e-01_r8,8.817038e-01_r8,8.808658e-01_r8,&
        & 8.797887e-01_r8,8.785323e-01_r8,8.771560e-01_r8,8.757196e-01_r8,8.742828e-01_r8,&
        & 8.729052e-01_r8,8.716467e-01_r8,8.705666e-01_r8,8.697250e-01_r8,8.691812e-01_r8,&
        & 8.689950e-01_r8,8.692264e-01_r8,8.699346e-01_r8,8.711795e-01_r8,8.730209e-01_r8,&
        & 8.755181e-01_r8,8.787312e-01_r8,8.8272e-01_r8 /)
      asyliq1(:, 28) = (/ &
        & 8.452284e-01_r8,8.522700e-01_r8,8.572973e-01_r8,8.607031e-01_r8,8.628802e-01_r8,&
        & 8.642215e-01_r8,8.651198e-01_r8,8.659679e-01_r8,8.671588e-01_r8,8.690853e-01_r8,&
        & 8.383803e-01_r8,8.557485e-01_r8,8.685851e-01_r8,8.774303e-01_r8,8.828245e-01_r8,&
        & 8.853077e-01_r8,8.854207e-01_r8,8.837034e-01_r8,8.806962e-01_r8,8.769398e-01_r8,&
        & 8.729740e-01_r8,8.693393e-01_r8,8.665761e-01_r8,8.652247e-01_r8,8.658253e-01_r8,&
        & 8.689182e-01_r8,8.750438e-01_r8,8.847424e-01_r8,8.636140e-01_r8,8.684449e-01_r8,&
        & 8.724400e-01_r8,8.756589e-01_r8,8.781613e-01_r8,8.800072e-01_r8,8.812559e-01_r8,&
        & 8.819671e-01_r8,8.822007e-01_r8,8.820165e-01_r8,8.814737e-01_r8,8.806322e-01_r8,&
        & 8.795518e-01_r8,8.782923e-01_r8,8.769129e-01_r8,8.754737e-01_r8,8.740342e-01_r8,&
        & 8.726542e-01_r8,8.713934e-01_r8,8.703111e-01_r8,8.694677e-01_r8,8.689222e-01_r8,&
        & 8.687344e-01_r8,8.689646e-01_r8,8.696715e-01_r8,8.709156e-01_r8,8.727563e-01_r8,&
        & 8.752531e-01_r8,8.784659e-01_r8,8.8245e-01_r8 /)
      asyliq1(:, 29) = (/ &
        & 7.800869e-01_r8,8.091120e-01_r8,8.325369e-01_r8,8.466266e-01_r8,8.515495e-01_r8,&
        & 8.499371e-01_r8,8.456203e-01_r8,8.430521e-01_r8,8.470286e-01_r8,8.625431e-01_r8,&
        & 8.402261e-01_r8,8.610822e-01_r8,8.776608e-01_r8,8.904485e-01_r8,8.999294e-01_r8,&
        & 9.065860e-01_r8,9.108995e-01_r8,9.133503e-01_r8,9.144187e-01_r8,9.145855e-01_r8,&
        & 9.143320e-01_r8,9.141402e-01_r8,9.144933e-01_r8,9.158754e-01_r8,9.187716e-01_r8,&
        & 9.236677e-01_r8,9.310503e-01_r8,9.414058e-01_r8,9.239108e-01_r8,9.300719e-01_r8,&
        & 9.353612e-01_r8,9.398378e-01_r8,9.435609e-01_r8,9.465895e-01_r8,9.489829e-01_r8,&
        & 9.508000e-01_r8,9.521002e-01_r8,9.529424e-01_r8,9.533860e-01_r8,9.534902e-01_r8,&
        & 9.533143e-01_r8,9.529177e-01_r8,9.523596e-01_r8,9.516997e-01_r8,9.509973e-01_r8,&
        & 9.503121e-01_r8,9.497037e-01_r8,9.492317e-01_r8,9.489558e-01_r8,9.489356e-01_r8,&
        & 9.492311e-01_r8,9.499019e-01_r8,9.510077e-01_r8,9.526084e-01_r8,9.547636e-01_r8,&
        & 9.575331e-01_r8,9.609766e-01_r8,9.6515e-01_r8 /)

! Spherical Ice Particle Parameterization
! extinction units (ext coef/iwc): [(m^-1)/(g m^-3)]
      extice2(:, 16) = (/ &
! band 16
        & 4.101824e-01_r8,2.435514e-01_r8,1.713697e-01_r8,1.314865e-01_r8,1.063406e-01_r8,&
        & 8.910701e-02_r8,7.659480e-02_r8,6.711784e-02_r8,5.970353e-02_r8,5.375249e-02_r8,&
        & 4.887577e-02_r8,4.481025e-02_r8,4.137171e-02_r8,3.842744e-02_r8,3.587948e-02_r8,&
        & 3.365396e-02_r8,3.169419e-02_r8,2.995593e-02_r8,2.840419e-02_r8,2.701091e-02_r8,&
        & 2.575336e-02_r8,2.461293e-02_r8,2.357423e-02_r8,2.262443e-02_r8,2.175276e-02_r8,&
        & 2.095012e-02_r8,2.020875e-02_r8,1.952199e-02_r8,1.888412e-02_r8,1.829018e-02_r8,&
        & 1.773586e-02_r8,1.721738e-02_r8,1.673144e-02_r8,1.627510e-02_r8,1.584579e-02_r8,&
        & 1.544122e-02_r8,1.505934e-02_r8,1.469833e-02_r8,1.435654e-02_r8,1.403251e-02_r8,&
        & 1.372492e-02_r8,1.343255e-02_r8,1.315433e-02_r8 /)
      extice2(:, 17) = (/ &
! band 17
        & 3.836650e-01_r8,2.304055e-01_r8,1.637265e-01_r8,1.266681e-01_r8,1.031602e-01_r8,&
        & 8.695191e-02_r8,7.511544e-02_r8,6.610009e-02_r8,5.900909e-02_r8,5.328833e-02_r8,&
        & 4.857728e-02_r8,4.463133e-02_r8,4.127880e-02_r8,3.839567e-02_r8,3.589013e-02_r8,&
        & 3.369280e-02_r8,3.175027e-02_r8,3.002079e-02_r8,2.847121e-02_r8,2.707493e-02_r8,&
        & 2.581031e-02_r8,2.465962e-02_r8,2.360815e-02_r8,2.264363e-02_r8,2.175571e-02_r8,&
        & 2.093563e-02_r8,2.017592e-02_r8,1.947015e-02_r8,1.881278e-02_r8,1.819901e-02_r8,&
        & 1.762463e-02_r8,1.708598e-02_r8,1.657982e-02_r8,1.610330e-02_r8,1.565390e-02_r8,&
        & 1.522937e-02_r8,1.482768e-02_r8,1.444706e-02_r8,1.408588e-02_r8,1.374270e-02_r8,&
        & 1.341619e-02_r8,1.310517e-02_r8,1.280857e-02_r8 /)
      extice2(:, 18) = (/ &
! band 18
        & 4.152673e-01_r8,2.436816e-01_r8,1.702243e-01_r8,1.299704e-01_r8,1.047528e-01_r8,&
        & 8.756039e-02_r8,7.513327e-02_r8,6.575690e-02_r8,5.844616e-02_r8,5.259609e-02_r8,&
        & 4.781531e-02_r8,4.383980e-02_r8,4.048517e-02_r8,3.761891e-02_r8,3.514342e-02_r8,&
        & 3.298525e-02_r8,3.108814e-02_r8,2.940825e-02_r8,2.791096e-02_r8,2.656858e-02_r8,&
        & 2.535869e-02_r8,2.426297e-02_r8,2.326627e-02_r8,2.235602e-02_r8,2.152164e-02_r8,&
        & 2.075420e-02_r8,2.004613e-02_r8,1.939091e-02_r8,1.878296e-02_r8,1.821744e-02_r8,&
        & 1.769015e-02_r8,1.719741e-02_r8,1.673600e-02_r8,1.630308e-02_r8,1.589615e-02_r8,&
        & 1.551298e-02_r8,1.515159e-02_r8,1.481021e-02_r8,1.448726e-02_r8,1.418131e-02_r8,&
        & 1.389109e-02_r8,1.361544e-02_r8,1.335330e-02_r8 /)
      extice2(:, 19) = (/ &
! band 19
        & 3.873250e-01_r8,2.331609e-01_r8,1.655002e-01_r8,1.277753e-01_r8,1.038247e-01_r8,&
        & 8.731780e-02_r8,7.527638e-02_r8,6.611873e-02_r8,5.892850e-02_r8,5.313885e-02_r8,&
        & 4.838068e-02_r8,4.440356e-02_r8,4.103167e-02_r8,3.813804e-02_r8,3.562870e-02_r8,&
        & 3.343269e-02_r8,3.149539e-02_r8,2.977414e-02_r8,2.823510e-02_r8,2.685112e-02_r8,&
        & 2.560015e-02_r8,2.446411e-02_r8,2.342805e-02_r8,2.247948e-02_r8,2.160789e-02_r8,&
        & 2.080438e-02_r8,2.006139e-02_r8,1.937238e-02_r8,1.873177e-02_r8,1.813469e-02_r8,&
        & 1.757689e-02_r8,1.705468e-02_r8,1.656479e-02_r8,1.610435e-02_r8,1.567081e-02_r8,&
        & 1.526192e-02_r8,1.487565e-02_r8,1.451020e-02_r8,1.416396e-02_r8,1.383546e-02_r8,&
        & 1.352339e-02_r8,1.322657e-02_r8,1.294392e-02_r8 /)
      extice2(:, 20) = (/ &
! band 20
        & 3.784280e-01_r8,2.291396e-01_r8,1.632551e-01_r8,1.263775e-01_r8,1.028944e-01_r8,&
        & 8.666975e-02_r8,7.480952e-02_r8,6.577335e-02_r8,5.866714e-02_r8,5.293694e-02_r8,&
        & 4.822153e-02_r8,4.427547e-02_r8,4.092626e-02_r8,3.804918e-02_r8,3.555184e-02_r8,&
        & 3.336440e-02_r8,3.143307e-02_r8,2.971577e-02_r8,2.817912e-02_r8,2.679632e-02_r8,&
        & 2.554558e-02_r8,2.440903e-02_r8,2.337187e-02_r8,2.242173e-02_r8,2.154821e-02_r8,&
        & 2.074249e-02_r8,1.999706e-02_r8,1.930546e-02_r8,1.866212e-02_r8,1.806221e-02_r8,&
        & 1.750152e-02_r8,1.697637e-02_r8,1.648352e-02_r8,1.602010e-02_r8,1.558358e-02_r8,&
        & 1.517172e-02_r8,1.478250e-02_r8,1.441413e-02_r8,1.406498e-02_r8,1.373362e-02_r8,&
        & 1.341872e-02_r8,1.311911e-02_r8,1.283371e-02_r8 /)
      extice2(:, 21) = (/ &
! band 21
        & 3.719909e-01_r8,2.259490e-01_r8,1.613144e-01_r8,1.250648e-01_r8,1.019462e-01_r8,&
        & 8.595358e-02_r8,7.425064e-02_r8,6.532618e-02_r8,5.830218e-02_r8,5.263421e-02_r8,&
        & 4.796697e-02_r8,4.405891e-02_r8,4.074013e-02_r8,3.788776e-02_r8,3.541071e-02_r8,&
        & 3.324008e-02_r8,3.132280e-02_r8,2.961733e-02_r8,2.809071e-02_r8,2.671645e-02_r8,&
        & 2.547302e-02_r8,2.434276e-02_r8,2.331102e-02_r8,2.236558e-02_r8,2.149614e-02_r8,&
        & 2.069397e-02_r8,1.995163e-02_r8,1.926272e-02_r8,1.862174e-02_r8,1.802389e-02_r8,&
        & 1.746500e-02_r8,1.694142e-02_r8,1.644994e-02_r8,1.598772e-02_r8,1.555225e-02_r8,&
        & 1.514129e-02_r8,1.475286e-02_r8,1.438515e-02_r8,1.403659e-02_r8,1.370572e-02_r8,&
        & 1.339124e-02_r8,1.309197e-02_r8,1.280685e-02_r8 /)
      extice2(:, 22) = (/ &
! band 22
        & 3.713158e-01_r8,2.253816e-01_r8,1.608461e-01_r8,1.246718e-01_r8,1.016109e-01_r8,&
        & 8.566332e-02_r8,7.399666e-02_r8,6.510199e-02_r8,5.810290e-02_r8,5.245608e-02_r8,&
        & 4.780702e-02_r8,4.391478e-02_r8,4.060989e-02_r8,3.776982e-02_r8,3.530374e-02_r8,&
        & 3.314296e-02_r8,3.123458e-02_r8,2.953719e-02_r8,2.801794e-02_r8,2.665043e-02_r8,&
        & 2.541321e-02_r8,2.428868e-02_r8,2.326224e-02_r8,2.232173e-02_r8,2.145688e-02_r8,&
        & 2.065899e-02_r8,1.992067e-02_r8,1.923552e-02_r8,1.859808e-02_r8,1.800356e-02_r8,&
        & 1.744782e-02_r8,1.692721e-02_r8,1.643855e-02_r8,1.597900e-02_r8,1.554606e-02_r8,&
        & 1.513751e-02_r8,1.475137e-02_r8,1.438586e-02_r8,1.403938e-02_r8,1.371050e-02_r8,&
        & 1.339793e-02_r8,1.310050e-02_r8,1.281713e-02_r8 /)
      extice2(:, 23) = (/ &
! band 23
        & 3.605883e-01_r8,2.204388e-01_r8,1.580431e-01_r8,1.229033e-01_r8,1.004203e-01_r8,&
        & 8.482616e-02_r8,7.338941e-02_r8,6.465105e-02_r8,5.776176e-02_r8,5.219398e-02_r8,&
        & 4.760288e-02_r8,4.375369e-02_r8,4.048111e-02_r8,3.766539e-02_r8,3.521771e-02_r8,&
        & 3.307079e-02_r8,3.117277e-02_r8,2.948303e-02_r8,2.796929e-02_r8,2.660560e-02_r8,&
        & 2.537086e-02_r8,2.424772e-02_r8,2.322182e-02_r8,2.228114e-02_r8,2.141556e-02_r8,&
        & 2.061649e-02_r8,1.987661e-02_r8,1.918962e-02_r8,1.855009e-02_r8,1.795330e-02_r8,&
        & 1.739514e-02_r8,1.687199e-02_r8,1.638069e-02_r8,1.591845e-02_r8,1.548276e-02_r8,&
        & 1.507143e-02_r8,1.468249e-02_r8,1.431416e-02_r8,1.396486e-02_r8,1.363318e-02_r8,&
        & 1.331781e-02_r8,1.301759e-02_r8,1.273147e-02_r8 /)
      extice2(:, 24) = (/ &
! band 24
        & 3.527890e-01_r8,2.168469e-01_r8,1.560090e-01_r8,1.216216e-01_r8,9.955787e-02_r8,&
        & 8.421942e-02_r8,7.294827e-02_r8,6.432192e-02_r8,5.751081e-02_r8,5.199888e-02_r8,&
        & 4.744835e-02_r8,4.362899e-02_r8,4.037847e-02_r8,3.757910e-02_r8,3.514351e-02_r8,&
        & 3.300546e-02_r8,3.111382e-02_r8,2.942853e-02_r8,2.791775e-02_r8,2.655584e-02_r8,&
        & 2.532195e-02_r8,2.419892e-02_r8,2.317255e-02_r8,2.223092e-02_r8,2.136402e-02_r8,&
        & 2.056334e-02_r8,1.982160e-02_r8,1.913258e-02_r8,1.849087e-02_r8,1.789178e-02_r8,&
        & 1.733124e-02_r8,1.680565e-02_r8,1.631187e-02_r8,1.584711e-02_r8,1.540889e-02_r8,&
        & 1.499502e-02_r8,1.460354e-02_r8,1.423269e-02_r8,1.388088e-02_r8,1.354670e-02_r8,&
        & 1.322887e-02_r8,1.292620e-02_r8,1.263767e-02_r8 /)
      extice2(:, 25) = (/ &
! band 25
        & 3.477874e-01_r8,2.143515e-01_r8,1.544887e-01_r8,1.205942e-01_r8,9.881779e-02_r8,&
        & 8.366261e-02_r8,7.251586e-02_r8,6.397790e-02_r8,5.723183e-02_r8,5.176908e-02_r8,&
        & 4.725658e-02_r8,4.346715e-02_r8,4.024055e-02_r8,3.746055e-02_r8,3.504080e-02_r8,&
        & 3.291583e-02_r8,3.103507e-02_r8,2.935891e-02_r8,2.785582e-02_r8,2.650042e-02_r8,&
        & 2.527206e-02_r8,2.415376e-02_r8,2.313142e-02_r8,2.219326e-02_r8,2.132934e-02_r8,&
        & 2.053122e-02_r8,1.979169e-02_r8,1.910456e-02_r8,1.846448e-02_r8,1.786680e-02_r8,&
        & 1.730745e-02_r8,1.678289e-02_r8,1.628998e-02_r8,1.582595e-02_r8,1.538835e-02_r8,&
        & 1.497499e-02_r8,1.458393e-02_r8,1.421341e-02_r8,1.386187e-02_r8,1.352788e-02_r8,&
        & 1.321019e-02_r8,1.290762e-02_r8,1.261913e-02_r8 /)
      extice2(:, 26) = (/ &
! band 26
        & 3.453721e-01_r8,2.130744e-01_r8,1.536698e-01_r8,1.200140e-01_r8,9.838078e-02_r8,&
        & 8.331940e-02_r8,7.223803e-02_r8,6.374775e-02_r8,5.703770e-02_r8,5.160290e-02_r8,&
        & 4.711259e-02_r8,4.334110e-02_r8,4.012923e-02_r8,3.736150e-02_r8,3.495208e-02_r8,&
        & 3.283589e-02_r8,3.096267e-02_r8,2.929302e-02_r8,2.779560e-02_r8,2.644517e-02_r8,&
        & 2.522119e-02_r8,2.410677e-02_r8,2.308788e-02_r8,2.215281e-02_r8,2.129165e-02_r8,&
        & 2.049602e-02_r8,1.975874e-02_r8,1.907365e-02_r8,1.843542e-02_r8,1.783943e-02_r8,&
        & 1.728162e-02_r8,1.675847e-02_r8,1.626685e-02_r8,1.580401e-02_r8,1.536750e-02_r8,&
        & 1.495515e-02_r8,1.456502e-02_r8,1.419537e-02_r8,1.384463e-02_r8,1.351139e-02_r8,&
        & 1.319438e-02_r8,1.289246e-02_r8,1.260456e-02_r8 /)
      extice2(:, 27) = (/ &
! band 27
        & 3.417883e-01_r8,2.113379e-01_r8,1.526395e-01_r8,1.193347e-01_r8,9.790253e-02_r8,&
        & 8.296715e-02_r8,7.196979e-02_r8,6.353806e-02_r8,5.687024e-02_r8,5.146670e-02_r8,&
        & 4.700001e-02_r8,4.324667e-02_r8,4.004894e-02_r8,3.729233e-02_r8,3.489172e-02_r8,&
        & 3.278257e-02_r8,3.091499e-02_r8,2.924987e-02_r8,2.775609e-02_r8,2.640859e-02_r8,&
        & 2.518695e-02_r8,2.407439e-02_r8,2.305697e-02_r8,2.212303e-02_r8,2.126273e-02_r8,&
        & 2.046774e-02_r8,1.973090e-02_r8,1.904610e-02_r8,1.840801e-02_r8,1.781204e-02_r8,&
        & 1.725417e-02_r8,1.673086e-02_r8,1.623902e-02_r8,1.577590e-02_r8,1.533906e-02_r8,&
        & 1.492634e-02_r8,1.453580e-02_r8,1.416571e-02_r8,1.381450e-02_r8,1.348078e-02_r8,&
        & 1.316327e-02_r8,1.286082e-02_r8,1.257240e-02_r8 /)
      extice2(:, 28) = (/ &
! band 28
        & 3.416111e-01_r8,2.114124e-01_r8,1.527734e-01_r8,1.194809e-01_r8,9.804612e-02_r8,&
        & 8.310287e-02_r8,7.209595e-02_r8,6.365442e-02_r8,5.697710e-02_r8,5.156460e-02_r8,&
        & 4.708957e-02_r8,4.332850e-02_r8,4.012361e-02_r8,3.736037e-02_r8,3.495364e-02_r8,&
        & 3.283879e-02_r8,3.096593e-02_r8,2.929589e-02_r8,2.779751e-02_r8,2.644571e-02_r8,&
        & 2.522004e-02_r8,2.410369e-02_r8,2.308271e-02_r8,2.214542e-02_r8,2.128195e-02_r8,&
        & 2.048396e-02_r8,1.974429e-02_r8,1.905679e-02_r8,1.841614e-02_r8,1.781774e-02_r8,&
        & 1.725754e-02_r8,1.673203e-02_r8,1.623807e-02_r8,1.577293e-02_r8,1.533416e-02_r8,&
        & 1.491958e-02_r8,1.452727e-02_r8,1.415547e-02_r8,1.380262e-02_r8,1.346732e-02_r8,&
        & 1.314830e-02_r8,1.284439e-02_r8,1.255456e-02_r8 /)
      extice2(:, 29) = (/ &
! band 29
        & 4.196611e-01_r8,2.493642e-01_r8,1.761261e-01_r8,1.357197e-01_r8,1.102161e-01_r8,&
        & 9.269376e-02_r8,7.992985e-02_r8,7.022538e-02_r8,6.260168e-02_r8,5.645603e-02_r8,&
        & 5.139732e-02_r8,4.716088e-02_r8,4.356133e-02_r8,4.046498e-02_r8,3.777303e-02_r8,&
        & 3.541094e-02_r8,3.332137e-02_r8,3.145954e-02_r8,2.978998e-02_r8,2.828419e-02_r8,&
        & 2.691905e-02_r8,2.567559e-02_r8,2.453811e-02_r8,2.349350e-02_r8,2.253072e-02_r8,&
        & 2.164042e-02_r8,2.081464e-02_r8,2.004652e-02_r8,1.933015e-02_r8,1.866041e-02_r8,&
        & 1.803283e-02_r8,1.744348e-02_r8,1.688894e-02_r8,1.636616e-02_r8,1.587244e-02_r8,&
        & 1.540539e-02_r8,1.496287e-02_r8,1.454295e-02_r8,1.414392e-02_r8,1.376423e-02_r8,&
        & 1.340247e-02_r8,1.305739e-02_r8,1.272784e-02_r8 /)

! single-scattering albedo: unitless
      ssaice2(:, 16) = (/ &
! band 16
        & 6.630615e-01_r8,6.451169e-01_r8,6.333696e-01_r8,6.246927e-01_r8,6.178420e-01_r8,&
        & 6.121976e-01_r8,6.074069e-01_r8,6.032505e-01_r8,5.995830e-01_r8,5.963030e-01_r8,&
        & 5.933372e-01_r8,5.906311e-01_r8,5.881427e-01_r8,5.858395e-01_r8,5.836955e-01_r8,&
        & 5.816896e-01_r8,5.798046e-01_r8,5.780264e-01_r8,5.763429e-01_r8,5.747441e-01_r8,&
        & 5.732213e-01_r8,5.717672e-01_r8,5.703754e-01_r8,5.690403e-01_r8,5.677571e-01_r8,&
        & 5.665215e-01_r8,5.653297e-01_r8,5.641782e-01_r8,5.630643e-01_r8,5.619850e-01_r8,&
        & 5.609381e-01_r8,5.599214e-01_r8,5.589328e-01_r8,5.579707e-01_r8,5.570333e-01_r8,&
        & 5.561193e-01_r8,5.552272e-01_r8,5.543558e-01_r8,5.535041e-01_r8,5.526708e-01_r8,&
        & 5.518551e-01_r8,5.510561e-01_r8,5.502729e-01_r8 /)
      ssaice2(:, 17) = (/ &
! band 17
        & 7.689749e-01_r8,7.398171e-01_r8,7.205819e-01_r8,7.065690e-01_r8,6.956928e-01_r8,&
        & 6.868989e-01_r8,6.795813e-01_r8,6.733606e-01_r8,6.679838e-01_r8,6.632742e-01_r8,&
        & 6.591036e-01_r8,6.553766e-01_r8,6.520197e-01_r8,6.489757e-01_r8,6.461991e-01_r8,&
        & 6.436531e-01_r8,6.413075e-01_r8,6.391375e-01_r8,6.371221e-01_r8,6.352438e-01_r8,&
        & 6.334876e-01_r8,6.318406e-01_r8,6.302918e-01_r8,6.288315e-01_r8,6.274512e-01_r8,&
        & 6.261436e-01_r8,6.249022e-01_r8,6.237211e-01_r8,6.225953e-01_r8,6.215201e-01_r8,&
        & 6.204914e-01_r8,6.195055e-01_r8,6.185592e-01_r8,6.176492e-01_r8,6.167730e-01_r8,&
        & 6.159280e-01_r8,6.151120e-01_r8,6.143228e-01_r8,6.135587e-01_r8,6.128177e-01_r8,&
        & 6.120984e-01_r8,6.113993e-01_r8,6.107189e-01_r8 /)
      ssaice2(:, 18) = (/ &
! band 18
        & 9.956167e-01_r8,9.814770e-01_r8,9.716104e-01_r8,9.639746e-01_r8,9.577179e-01_r8,&
        & 9.524010e-01_r8,9.477672e-01_r8,9.436527e-01_r8,9.399467e-01_r8,9.365708e-01_r8,&
        & 9.334672e-01_r8,9.305921e-01_r8,9.279118e-01_r8,9.253993e-01_r8,9.230330e-01_r8,&
        & 9.207954e-01_r8,9.186719e-01_r8,9.166501e-01_r8,9.147199e-01_r8,9.128722e-01_r8,&
        & 9.110997e-01_r8,9.093956e-01_r8,9.077544e-01_r8,9.061708e-01_r8,9.046406e-01_r8,&
        & 9.031598e-01_r8,9.017248e-01_r8,9.003326e-01_r8,8.989804e-01_r8,8.976655e-01_r8,&
        & 8.963857e-01_r8,8.951389e-01_r8,8.939233e-01_r8,8.927370e-01_r8,8.915785e-01_r8,&
        & 8.904464e-01_r8,8.893392e-01_r8,8.882559e-01_r8,8.871951e-01_r8,8.861559e-01_r8,&
        & 8.851373e-01_r8,8.841383e-01_r8,8.831581e-01_r8 /)
      ssaice2(:, 19) = (/ &
! band 19
        & 9.723177e-01_r8,9.452119e-01_r8,9.267592e-01_r8,9.127393e-01_r8,9.014238e-01_r8,&
        & 8.919334e-01_r8,8.837584e-01_r8,8.765773e-01_r8,8.701736e-01_r8,8.643950e-01_r8,&
        & 8.591299e-01_r8,8.542942e-01_r8,8.498230e-01_r8,8.456651e-01_r8,8.417794e-01_r8,&
        & 8.381324e-01_r8,8.346964e-01_r8,8.314484e-01_r8,8.283687e-01_r8,8.254408e-01_r8,&
        & 8.226505e-01_r8,8.199854e-01_r8,8.174348e-01_r8,8.149891e-01_r8,8.126403e-01_r8,&
        & 8.103808e-01_r8,8.082041e-01_r8,8.061044e-01_r8,8.040765e-01_r8,8.021156e-01_r8,&
        & 8.002174e-01_r8,7.983781e-01_r8,7.965941e-01_r8,7.948622e-01_r8,7.931795e-01_r8,&
        & 7.915432e-01_r8,7.899508e-01_r8,7.884002e-01_r8,7.868891e-01_r8,7.854156e-01_r8,&
        & 7.839779e-01_r8,7.825742e-01_r8,7.812031e-01_r8 /)
      ssaice2(:, 20) = (/ &
! band 20
        & 9.933294e-01_r8,9.860917e-01_r8,9.811564e-01_r8,9.774008e-01_r8,9.743652e-01_r8,&
        & 9.718155e-01_r8,9.696159e-01_r8,9.676810e-01_r8,9.659531e-01_r8,9.643915e-01_r8,&
        & 9.629667e-01_r8,9.616561e-01_r8,9.604426e-01_r8,9.593125e-01_r8,9.582548e-01_r8,&
        & 9.572607e-01_r8,9.563227e-01_r8,9.554347e-01_r8,9.545915e-01_r8,9.537888e-01_r8,&
        & 9.530226e-01_r8,9.522898e-01_r8,9.515874e-01_r8,9.509130e-01_r8,9.502643e-01_r8,&
        & 9.496394e-01_r8,9.490366e-01_r8,9.484542e-01_r8,9.478910e-01_r8,9.473456e-01_r8,&
        & 9.468169e-01_r8,9.463039e-01_r8,9.458056e-01_r8,9.453212e-01_r8,9.448499e-01_r8,&
        & 9.443910e-01_r8,9.439438e-01_r8,9.435077e-01_r8,9.430821e-01_r8,9.426666e-01_r8,&
        & 9.422607e-01_r8,9.418638e-01_r8,9.414756e-01_r8 /)
      ssaice2(:, 21) = (/ &
! band 21
        & 9.900787e-01_r8,9.828880e-01_r8,9.779258e-01_r8,9.741173e-01_r8,9.710184e-01_r8,&
        & 9.684012e-01_r8,9.661332e-01_r8,9.641301e-01_r8,9.623352e-01_r8,9.607083e-01_r8,&
        & 9.592198e-01_r8,9.578474e-01_r8,9.565739e-01_r8,9.553856e-01_r8,9.542715e-01_r8,&
        & 9.532226e-01_r8,9.522314e-01_r8,9.512919e-01_r8,9.503986e-01_r8,9.495472e-01_r8,&
        & 9.487337e-01_r8,9.479549e-01_r8,9.472077e-01_r8,9.464897e-01_r8,9.457985e-01_r8,&
        & 9.451322e-01_r8,9.444890e-01_r8,9.438673e-01_r8,9.432656e-01_r8,9.426826e-01_r8,&
        & 9.421173e-01_r8,9.415684e-01_r8,9.410351e-01_r8,9.405164e-01_r8,9.400115e-01_r8,&
        & 9.395198e-01_r8,9.390404e-01_r8,9.385728e-01_r8,9.381164e-01_r8,9.376707e-01_r8,&
        & 9.372350e-01_r8,9.368091e-01_r8,9.363923e-01_r8 /)
      ssaice2(:, 22) = (/ &
! band 22
        & 9.986793e-01_r8,9.985239e-01_r8,9.983911e-01_r8,9.982715e-01_r8,9.981606e-01_r8,&
        & 9.980562e-01_r8,9.979567e-01_r8,9.978613e-01_r8,9.977691e-01_r8,9.976798e-01_r8,&
        & 9.975929e-01_r8,9.975081e-01_r8,9.974251e-01_r8,9.973438e-01_r8,9.972640e-01_r8,&
        & 9.971855e-01_r8,9.971083e-01_r8,9.970322e-01_r8,9.969571e-01_r8,9.968830e-01_r8,&
        & 9.968099e-01_r8,9.967375e-01_r8,9.966660e-01_r8,9.965951e-01_r8,9.965250e-01_r8,&
        & 9.964555e-01_r8,9.963867e-01_r8,9.963185e-01_r8,9.962508e-01_r8,9.961836e-01_r8,&
        & 9.961170e-01_r8,9.960508e-01_r8,9.959851e-01_r8,9.959198e-01_r8,9.958550e-01_r8,&
        & 9.957906e-01_r8,9.957266e-01_r8,9.956629e-01_r8,9.955997e-01_r8,9.955367e-01_r8,&
        & 9.954742e-01_r8,9.954119e-01_r8,9.953500e-01_r8 /)
      ssaice2(:, 23) = (/ &
! band 23
        & 9.997944e-01_r8,9.997791e-01_r8,9.997664e-01_r8,9.997547e-01_r8,9.997436e-01_r8,&
        & 9.997327e-01_r8,9.997219e-01_r8,9.997110e-01_r8,9.996999e-01_r8,9.996886e-01_r8,&
        & 9.996771e-01_r8,9.996653e-01_r8,9.996533e-01_r8,9.996409e-01_r8,9.996282e-01_r8,&
        & 9.996152e-01_r8,9.996019e-01_r8,9.995883e-01_r8,9.995743e-01_r8,9.995599e-01_r8,&
        & 9.995453e-01_r8,9.995302e-01_r8,9.995149e-01_r8,9.994992e-01_r8,9.994831e-01_r8,&
        & 9.994667e-01_r8,9.994500e-01_r8,9.994329e-01_r8,9.994154e-01_r8,9.993976e-01_r8,&
        & 9.993795e-01_r8,9.993610e-01_r8,9.993422e-01_r8,9.993230e-01_r8,9.993035e-01_r8,&
        & 9.992837e-01_r8,9.992635e-01_r8,9.992429e-01_r8,9.992221e-01_r8,9.992008e-01_r8,&
        & 9.991793e-01_r8,9.991574e-01_r8,9.991352e-01_r8 /)
      ssaice2(:, 24) = (/ &
! band 24
        & 9.999949e-01_r8,9.999947e-01_r8,9.999943e-01_r8,9.999939e-01_r8,9.999934e-01_r8,&
        & 9.999927e-01_r8,9.999920e-01_r8,9.999913e-01_r8,9.999904e-01_r8,9.999895e-01_r8,&
        & 9.999885e-01_r8,9.999874e-01_r8,9.999863e-01_r8,9.999851e-01_r8,9.999838e-01_r8,&
        & 9.999824e-01_r8,9.999810e-01_r8,9.999795e-01_r8,9.999780e-01_r8,9.999764e-01_r8,&
        & 9.999747e-01_r8,9.999729e-01_r8,9.999711e-01_r8,9.999692e-01_r8,9.999673e-01_r8,&
        & 9.999653e-01_r8,9.999632e-01_r8,9.999611e-01_r8,9.999589e-01_r8,9.999566e-01_r8,&
        & 9.999543e-01_r8,9.999519e-01_r8,9.999495e-01_r8,9.999470e-01_r8,9.999444e-01_r8,&
        & 9.999418e-01_r8,9.999392e-01_r8,9.999364e-01_r8,9.999336e-01_r8,9.999308e-01_r8,&
        & 9.999279e-01_r8,9.999249e-01_r8,9.999219e-01_r8 /)
      ssaice2(:, 25) = (/ &
! band 25
        & 9.999997e-01_r8,9.999997e-01_r8,9.999997e-01_r8,9.999996e-01_r8,9.999996e-01_r8,&
        & 9.999995e-01_r8,9.999994e-01_r8,9.999993e-01_r8,9.999993e-01_r8,9.999992e-01_r8,&
        & 9.999991e-01_r8,9.999989e-01_r8,9.999988e-01_r8,9.999987e-01_r8,9.999986e-01_r8,&
        & 9.999984e-01_r8,9.999983e-01_r8,9.999981e-01_r8,9.999980e-01_r8,9.999978e-01_r8,&
        & 9.999976e-01_r8,9.999974e-01_r8,9.999972e-01_r8,9.999971e-01_r8,9.999969e-01_r8,&
        & 9.999966e-01_r8,9.999964e-01_r8,9.999962e-01_r8,9.999960e-01_r8,9.999957e-01_r8,&
        & 9.999955e-01_r8,9.999953e-01_r8,9.999950e-01_r8,9.999947e-01_r8,9.999945e-01_r8,&
        & 9.999942e-01_r8,9.999939e-01_r8,9.999936e-01_r8,9.999934e-01_r8,9.999931e-01_r8,&
        & 9.999928e-01_r8,9.999925e-01_r8,9.999921e-01_r8 /)
      ssaice2(:, 26) = (/ &
! band 26
        & 9.999997e-01_r8,9.999996e-01_r8,9.999996e-01_r8,9.999995e-01_r8,9.999994e-01_r8,&
        & 9.999993e-01_r8,9.999992e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,9.999982e-01_r8,9.999980e-01_r8,&
        & 9.999978e-01_r8,9.999976e-01_r8,9.999974e-01_r8,9.999972e-01_r8,9.999970e-01_r8,&
        & 9.999967e-01_r8,9.999965e-01_r8,9.999962e-01_r8,9.999959e-01_r8,9.999956e-01_r8,&
        & 9.999954e-01_r8,9.999951e-01_r8,9.999947e-01_r8,9.999944e-01_r8,9.999941e-01_r8,&
        & 9.999938e-01_r8,9.999934e-01_r8,9.999931e-01_r8,9.999927e-01_r8,9.999923e-01_r8,&
        & 9.999920e-01_r8,9.999916e-01_r8,9.999912e-01_r8,9.999908e-01_r8,9.999904e-01_r8,&
        & 9.999899e-01_r8,9.999895e-01_r8,9.999891e-01_r8 /)
      ssaice2(:, 27) = (/ &
! band 27
        & 9.999987e-01_r8,9.999987e-01_r8,9.999985e-01_r8,9.999984e-01_r8,9.999982e-01_r8,&
        & 9.999980e-01_r8,9.999978e-01_r8,9.999976e-01_r8,9.999973e-01_r8,9.999970e-01_r8,&
        & 9.999967e-01_r8,9.999964e-01_r8,9.999960e-01_r8,9.999956e-01_r8,9.999952e-01_r8,&
        & 9.999948e-01_r8,9.999944e-01_r8,9.999939e-01_r8,9.999934e-01_r8,9.999929e-01_r8,&
        & 9.999924e-01_r8,9.999918e-01_r8,9.999913e-01_r8,9.999907e-01_r8,9.999901e-01_r8,&
        & 9.999894e-01_r8,9.999888e-01_r8,9.999881e-01_r8,9.999874e-01_r8,9.999867e-01_r8,&
        & 9.999860e-01_r8,9.999853e-01_r8,9.999845e-01_r8,9.999837e-01_r8,9.999829e-01_r8,&
        & 9.999821e-01_r8,9.999813e-01_r8,9.999804e-01_r8,9.999796e-01_r8,9.999787e-01_r8,&
        & 9.999778e-01_r8,9.999768e-01_r8,9.999759e-01_r8 /)
      ssaice2(:, 28) = (/ &
! band 28
        & 9.999989e-01_r8,9.999989e-01_r8,9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,&
        & 9.999982e-01_r8,9.999980e-01_r8,9.999978e-01_r8,9.999975e-01_r8,9.999972e-01_r8,&
        & 9.999969e-01_r8,9.999966e-01_r8,9.999962e-01_r8,9.999958e-01_r8,9.999954e-01_r8,&
        & 9.999950e-01_r8,9.999945e-01_r8,9.999941e-01_r8,9.999936e-01_r8,9.999931e-01_r8,&
        & 9.999925e-01_r8,9.999920e-01_r8,9.999914e-01_r8,9.999908e-01_r8,9.999902e-01_r8,&
        & 9.999896e-01_r8,9.999889e-01_r8,9.999883e-01_r8,9.999876e-01_r8,9.999869e-01_r8,&
        & 9.999861e-01_r8,9.999854e-01_r8,9.999846e-01_r8,9.999838e-01_r8,9.999830e-01_r8,&
        & 9.999822e-01_r8,9.999814e-01_r8,9.999805e-01_r8,9.999796e-01_r8,9.999787e-01_r8,&
        & 9.999778e-01_r8,9.999769e-01_r8,9.999759e-01_r8 /)
      ssaice2(:, 29) = (/ &
! band 29
        & 7.042143e-01_r8,6.691161e-01_r8,6.463240e-01_r8,6.296590e-01_r8,6.166381e-01_r8,&
        & 6.060183e-01_r8,5.970908e-01_r8,5.894144e-01_r8,5.826968e-01_r8,5.767343e-01_r8,&
        & 5.713804e-01_r8,5.665256e-01_r8,5.620867e-01_r8,5.579987e-01_r8,5.542101e-01_r8,&
        & 5.506794e-01_r8,5.473727e-01_r8,5.442620e-01_r8,5.413239e-01_r8,5.385389e-01_r8,&
        & 5.358901e-01_r8,5.333633e-01_r8,5.309460e-01_r8,5.286277e-01_r8,5.263988e-01_r8,&
        & 5.242512e-01_r8,5.221777e-01_r8,5.201719e-01_r8,5.182280e-01_r8,5.163410e-01_r8,&
        & 5.145062e-01_r8,5.127197e-01_r8,5.109776e-01_r8,5.092766e-01_r8,5.076137e-01_r8,&
        & 5.059860e-01_r8,5.043911e-01_r8,5.028266e-01_r8,5.012904e-01_r8,4.997805e-01_r8,&
        & 4.982951e-01_r8,4.968326e-01_r8,4.953913e-01_r8 /)

! asymmetry factor: unitless
      asyice2(:, 16) = (/ &
! band 16
        & 7.946655e-01_r8,8.547685e-01_r8,8.806016e-01_r8,8.949880e-01_r8,9.041676e-01_r8,&
        & 9.105399e-01_r8,9.152249e-01_r8,9.188160e-01_r8,9.216573e-01_r8,9.239620e-01_r8,&
        & 9.258695e-01_r8,9.274745e-01_r8,9.288441e-01_r8,9.300267e-01_r8,9.310584e-01_r8,&
        & 9.319665e-01_r8,9.327721e-01_r8,9.334918e-01_r8,9.341387e-01_r8,9.347236e-01_r8,&
        & 9.352551e-01_r8,9.357402e-01_r8,9.361850e-01_r8,9.365942e-01_r8,9.369722e-01_r8,&
        & 9.373225e-01_r8,9.376481e-01_r8,9.379516e-01_r8,9.382352e-01_r8,9.385010e-01_r8,&
        & 9.387505e-01_r8,9.389854e-01_r8,9.392070e-01_r8,9.394163e-01_r8,9.396145e-01_r8,&
        & 9.398024e-01_r8,9.399809e-01_r8,9.401508e-01_r8,9.403126e-01_r8,9.404670e-01_r8,&
        & 9.406144e-01_r8,9.407555e-01_r8,9.408906e-01_r8 /)
      asyice2(:, 17) = (/ &
! band 17
        & 9.078091e-01_r8,9.195850e-01_r8,9.267250e-01_r8,9.317083e-01_r8,9.354632e-01_r8,&
        & 9.384323e-01_r8,9.408597e-01_r8,9.428935e-01_r8,9.446301e-01_r8,9.461351e-01_r8,&
        & 9.474555e-01_r8,9.486259e-01_r8,9.496722e-01_r8,9.506146e-01_r8,9.514688e-01_r8,&
        & 9.522476e-01_r8,9.529612e-01_r8,9.536181e-01_r8,9.542251e-01_r8,9.547883e-01_r8,&
        & 9.553124e-01_r8,9.558019e-01_r8,9.562601e-01_r8,9.566904e-01_r8,9.570953e-01_r8,&
        & 9.574773e-01_r8,9.578385e-01_r8,9.581806e-01_r8,9.585054e-01_r8,9.588142e-01_r8,&
        & 9.591083e-01_r8,9.593888e-01_r8,9.596569e-01_r8,9.599135e-01_r8,9.601593e-01_r8,&
        & 9.603952e-01_r8,9.606219e-01_r8,9.608399e-01_r8,9.610499e-01_r8,9.612523e-01_r8,&
        & 9.614477e-01_r8,9.616365e-01_r8,9.618192e-01_r8 /)
      asyice2(:, 18) = (/ &
! band 18
        & 8.322045e-01_r8,8.528693e-01_r8,8.648167e-01_r8,8.729163e-01_r8,8.789054e-01_r8,&
        & 8.835845e-01_r8,8.873819e-01_r8,8.905511e-01_r8,8.932532e-01_r8,8.955965e-01_r8,&
        & 8.976567e-01_r8,8.994887e-01_r8,9.011334e-01_r8,9.026221e-01_r8,9.039791e-01_r8,&
        & 9.052237e-01_r8,9.063715e-01_r8,9.074349e-01_r8,9.084245e-01_r8,9.093489e-01_r8,&
        & 9.102154e-01_r8,9.110303e-01_r8,9.117987e-01_r8,9.125253e-01_r8,9.132140e-01_r8,&
        & 9.138682e-01_r8,9.144910e-01_r8,9.150850e-01_r8,9.156524e-01_r8,9.161955e-01_r8,&
        & 9.167160e-01_r8,9.172157e-01_r8,9.176959e-01_r8,9.181581e-01_r8,9.186034e-01_r8,&
        & 9.190330e-01_r8,9.194478e-01_r8,9.198488e-01_r8,9.202368e-01_r8,9.206126e-01_r8,&
        & 9.209768e-01_r8,9.213301e-01_r8,9.216731e-01_r8 /)
      asyice2(:, 19) = (/ &
! band 19
        & 8.116560e-01_r8,8.488278e-01_r8,8.674331e-01_r8,8.788148e-01_r8,8.865810e-01_r8,&
        & 8.922595e-01_r8,8.966149e-01_r8,9.000747e-01_r8,9.028980e-01_r8,9.052513e-01_r8,&
        & 9.072468e-01_r8,9.089632e-01_r8,9.104574e-01_r8,9.117713e-01_r8,9.129371e-01_r8,&
        & 9.139793e-01_r8,9.149174e-01_r8,9.157668e-01_r8,9.165400e-01_r8,9.172473e-01_r8,&
        & 9.178970e-01_r8,9.184962e-01_r8,9.190508e-01_r8,9.195658e-01_r8,9.200455e-01_r8,&
        & 9.204935e-01_r8,9.209130e-01_r8,9.213067e-01_r8,9.216771e-01_r8,9.220262e-01_r8,&
        & 9.223560e-01_r8,9.226680e-01_r8,9.229636e-01_r8,9.232443e-01_r8,9.235112e-01_r8,&
        & 9.237652e-01_r8,9.240074e-01_r8,9.242385e-01_r8,9.244594e-01_r8,9.246708e-01_r8,&
        & 9.248733e-01_r8,9.250674e-01_r8,9.252536e-01_r8 /)
      asyice2(:, 20) = (/ &
! band 20
        & 8.047113e-01_r8,8.402864e-01_r8,8.570332e-01_r8,8.668455e-01_r8,8.733206e-01_r8,&
        & 8.779272e-01_r8,8.813796e-01_r8,8.840676e-01_r8,8.862225e-01_r8,8.879904e-01_r8,&
        & 8.894682e-01_r8,8.907228e-01_r8,8.918019e-01_r8,8.927404e-01_r8,8.935645e-01_r8,&
        & 8.942943e-01_r8,8.949452e-01_r8,8.955296e-01_r8,8.960574e-01_r8,8.965366e-01_r8,&
        & 8.969736e-01_r8,8.973740e-01_r8,8.977422e-01_r8,8.980820e-01_r8,8.983966e-01_r8,&
        & 8.986889e-01_r8,8.989611e-01_r8,8.992153e-01_r8,8.994533e-01_r8,8.996766e-01_r8,&
        & 8.998865e-01_r8,9.000843e-01_r8,9.002709e-01_r8,9.004474e-01_r8,9.006146e-01_r8,&
        & 9.007731e-01_r8,9.009237e-01_r8,9.010670e-01_r8,9.012034e-01_r8,9.013336e-01_r8,&
        & 9.014579e-01_r8,9.015767e-01_r8,9.016904e-01_r8 /)
      asyice2(:, 21) = (/ &
! band 21
        & 8.179122e-01_r8,8.480726e-01_r8,8.621945e-01_r8,8.704354e-01_r8,8.758555e-01_r8,&
        & 8.797007e-01_r8,8.825750e-01_r8,8.848078e-01_r8,8.865939e-01_r8,8.880564e-01_r8,&
        & 8.892765e-01_r8,8.903105e-01_r8,8.911982e-01_r8,8.919689e-01_r8,8.926446e-01_r8,&
        & 8.932419e-01_r8,8.937738e-01_r8,8.942506e-01_r8,8.946806e-01_r8,8.950702e-01_r8,&
        & 8.954251e-01_r8,8.957497e-01_r8,8.960477e-01_r8,8.963223e-01_r8,8.965762e-01_r8,&
        & 8.968116e-01_r8,8.970306e-01_r8,8.972347e-01_r8,8.974255e-01_r8,8.976042e-01_r8,&
        & 8.977720e-01_r8,8.979298e-01_r8,8.980784e-01_r8,8.982188e-01_r8,8.983515e-01_r8,&
        & 8.984771e-01_r8,8.985963e-01_r8,8.987095e-01_r8,8.988171e-01_r8,8.989195e-01_r8,&
        & 8.990172e-01_r8,8.991104e-01_r8,8.991994e-01_r8 /)
      asyice2(:, 22) = (/ &
! band 22
        & 8.169789e-01_r8,8.455024e-01_r8,8.586925e-01_r8,8.663283e-01_r8,8.713217e-01_r8,&
        & 8.748488e-01_r8,8.774765e-01_r8,8.795122e-01_r8,8.811370e-01_r8,8.824649e-01_r8,&
        & 8.835711e-01_r8,8.845073e-01_r8,8.853103e-01_r8,8.860068e-01_r8,8.866170e-01_r8,&
        & 8.871560e-01_r8,8.876358e-01_r8,8.880658e-01_r8,8.884533e-01_r8,8.888044e-01_r8,&
        & 8.891242e-01_r8,8.894166e-01_r8,8.896851e-01_r8,8.899324e-01_r8,8.901612e-01_r8,&
        & 8.903733e-01_r8,8.905706e-01_r8,8.907545e-01_r8,8.909265e-01_r8,8.910876e-01_r8,&
        & 8.912388e-01_r8,8.913812e-01_r8,8.915153e-01_r8,8.916419e-01_r8,8.917617e-01_r8,&
        & 8.918752e-01_r8,8.919829e-01_r8,8.920851e-01_r8,8.921824e-01_r8,8.922751e-01_r8,&
        & 8.923635e-01_r8,8.924478e-01_r8,8.925284e-01_r8 /)
      asyice2(:, 23) = (/ &
! band 23
        & 8.387642e-01_r8,8.569979e-01_r8,8.658630e-01_r8,8.711825e-01_r8,8.747605e-01_r8,&
        & 8.773472e-01_r8,8.793129e-01_r8,8.808621e-01_r8,8.821179e-01_r8,8.831583e-01_r8,&
        & 8.840361e-01_r8,8.847875e-01_r8,8.854388e-01_r8,8.860094e-01_r8,8.865138e-01_r8,&
        & 8.869634e-01_r8,8.873668e-01_r8,8.877310e-01_r8,8.880617e-01_r8,8.883635e-01_r8,&
        & 8.886401e-01_r8,8.888947e-01_r8,8.891298e-01_r8,8.893477e-01_r8,8.895504e-01_r8,&
        & 8.897393e-01_r8,8.899159e-01_r8,8.900815e-01_r8,8.902370e-01_r8,8.903833e-01_r8,&
        & 8.905214e-01_r8,8.906518e-01_r8,8.907753e-01_r8,8.908924e-01_r8,8.910036e-01_r8,&
        & 8.911094e-01_r8,8.912101e-01_r8,8.913062e-01_r8,8.913979e-01_r8,8.914856e-01_r8,&
        & 8.915695e-01_r8,8.916498e-01_r8,8.917269e-01_r8 /)
      asyice2(:, 24) = (/ &
! band 24
        & 8.522208e-01_r8,8.648132e-01_r8,8.711224e-01_r8,8.749901e-01_r8,8.776354e-01_r8,&
        & 8.795743e-01_r8,8.810649e-01_r8,8.822518e-01_r8,8.832225e-01_r8,8.840333e-01_r8,&
        & 8.847224e-01_r8,8.853162e-01_r8,8.858342e-01_r8,8.862906e-01_r8,8.866962e-01_r8,&
        & 8.870595e-01_r8,8.873871e-01_r8,8.876842e-01_r8,8.879551e-01_r8,8.882032e-01_r8,&
        & 8.884316e-01_r8,8.886425e-01_r8,8.888380e-01_r8,8.890199e-01_r8,8.891895e-01_r8,&
        & 8.893481e-01_r8,8.894968e-01_r8,8.896366e-01_r8,8.897683e-01_r8,8.898926e-01_r8,&
        & 8.900102e-01_r8,8.901215e-01_r8,8.902272e-01_r8,8.903276e-01_r8,8.904232e-01_r8,&
        & 8.905144e-01_r8,8.906014e-01_r8,8.906845e-01_r8,8.907640e-01_r8,8.908402e-01_r8,&
        & 8.909132e-01_r8,8.909834e-01_r8,8.910507e-01_r8 /)
      asyice2(:, 25) = (/ &
! band 25
        & 8.578202e-01_r8,8.683033e-01_r8,8.735431e-01_r8,8.767488e-01_r8,8.789378e-01_r8,&
        & 8.805399e-01_r8,8.817701e-01_r8,8.827485e-01_r8,8.835480e-01_r8,8.842152e-01_r8,&
        & 8.847817e-01_r8,8.852696e-01_r8,8.856949e-01_r8,8.860694e-01_r8,8.864020e-01_r8,&
        & 8.866997e-01_r8,8.869681e-01_r8,8.872113e-01_r8,8.874330e-01_r8,8.876360e-01_r8,&
        & 8.878227e-01_r8,8.879951e-01_r8,8.881548e-01_r8,8.883033e-01_r8,8.884418e-01_r8,&
        & 8.885712e-01_r8,8.886926e-01_r8,8.888066e-01_r8,8.889139e-01_r8,8.890152e-01_r8,&
        & 8.891110e-01_r8,8.892017e-01_r8,8.892877e-01_r8,8.893695e-01_r8,8.894473e-01_r8,&
        & 8.895214e-01_r8,8.895921e-01_r8,8.896597e-01_r8,8.897243e-01_r8,8.897862e-01_r8,&
        & 8.898456e-01_r8,8.899025e-01_r8,8.899572e-01_r8 /)
      asyice2(:, 26) = (/ &
! band 26
        & 8.625615e-01_r8,8.713831e-01_r8,8.755799e-01_r8,8.780560e-01_r8,8.796983e-01_r8,&
        & 8.808714e-01_r8,8.817534e-01_r8,8.824420e-01_r8,8.829953e-01_r8,8.834501e-01_r8,&
        & 8.838310e-01_r8,8.841549e-01_r8,8.844338e-01_r8,8.846767e-01_r8,8.848902e-01_r8,&
        & 8.850795e-01_r8,8.852484e-01_r8,8.854002e-01_r8,8.855374e-01_r8,8.856620e-01_r8,&
        & 8.857758e-01_r8,8.858800e-01_r8,8.859759e-01_r8,8.860644e-01_r8,8.861464e-01_r8,&
        & 8.862225e-01_r8,8.862935e-01_r8,8.863598e-01_r8,8.864218e-01_r8,8.864800e-01_r8,&
        & 8.865347e-01_r8,8.865863e-01_r8,8.866349e-01_r8,8.866809e-01_r8,8.867245e-01_r8,&
        & 8.867658e-01_r8,8.868050e-01_r8,8.868423e-01_r8,8.868778e-01_r8,8.869117e-01_r8,&
        & 8.869440e-01_r8,8.869749e-01_r8,8.870044e-01_r8 /)
      asyice2(:, 27) = (/ &
! band 27
        & 8.587495e-01_r8,8.684764e-01_r8,8.728189e-01_r8,8.752872e-01_r8,8.768846e-01_r8,&
        & 8.780060e-01_r8,8.788386e-01_r8,8.794824e-01_r8,8.799960e-01_r8,8.804159e-01_r8,&
        & 8.807660e-01_r8,8.810626e-01_r8,8.813175e-01_r8,8.815390e-01_r8,8.817335e-01_r8,&
        & 8.819057e-01_r8,8.820593e-01_r8,8.821973e-01_r8,8.823220e-01_r8,8.824353e-01_r8,&
        & 8.825387e-01_r8,8.826336e-01_r8,8.827209e-01_r8,8.828016e-01_r8,8.828764e-01_r8,&
        & 8.829459e-01_r8,8.830108e-01_r8,8.830715e-01_r8,8.831283e-01_r8,8.831817e-01_r8,&
        & 8.832320e-01_r8,8.832795e-01_r8,8.833244e-01_r8,8.833668e-01_r8,8.834071e-01_r8,&
        & 8.834454e-01_r8,8.834817e-01_r8,8.835164e-01_r8,8.835495e-01_r8,8.835811e-01_r8,&
        & 8.836113e-01_r8,8.836402e-01_r8,8.836679e-01_r8 /)
      asyice2(:, 28) = (/ &
! band 28
        & 8.561110e-01_r8,8.678583e-01_r8,8.727554e-01_r8,8.753892e-01_r8,8.770154e-01_r8,&
        & 8.781109e-01_r8,8.788949e-01_r8,8.794812e-01_r8,8.799348e-01_r8,8.802952e-01_r8,&
        & 8.805880e-01_r8,8.808300e-01_r8,8.810331e-01_r8,8.812058e-01_r8,8.813543e-01_r8,&
        & 8.814832e-01_r8,8.815960e-01_r8,8.816956e-01_r8,8.817839e-01_r8,8.818629e-01_r8,&
        & 8.819339e-01_r8,8.819979e-01_r8,8.820560e-01_r8,8.821089e-01_r8,8.821573e-01_r8,&
        & 8.822016e-01_r8,8.822425e-01_r8,8.822801e-01_r8,8.823150e-01_r8,8.823474e-01_r8,&
        & 8.823775e-01_r8,8.824056e-01_r8,8.824318e-01_r8,8.824564e-01_r8,8.824795e-01_r8,&
        & 8.825011e-01_r8,8.825215e-01_r8,8.825408e-01_r8,8.825589e-01_r8,8.825761e-01_r8,&
        & 8.825924e-01_r8,8.826078e-01_r8,8.826224e-01_r8 /)
      asyice2(:, 29) = (/ &
! band 29
        & 8.311124e-01_r8,8.688197e-01_r8,8.900274e-01_r8,9.040696e-01_r8,9.142334e-01_r8,&
        & 9.220181e-01_r8,9.282195e-01_r8,9.333048e-01_r8,9.375689e-01_r8,9.412085e-01_r8,&
        & 9.443604e-01_r8,9.471230e-01_r8,9.495694e-01_r8,9.517549e-01_r8,9.537224e-01_r8,&
        & 9.555057e-01_r8,9.571316e-01_r8,9.586222e-01_r8,9.599952e-01_r8,9.612656e-01_r8,&
        & 9.624458e-01_r8,9.635461e-01_r8,9.645756e-01_r8,9.655418e-01_r8,9.664513e-01_r8,&
        & 9.673098e-01_r8,9.681222e-01_r8,9.688928e-01_r8,9.696256e-01_r8,9.703237e-01_r8,&
        & 9.709903e-01_r8,9.716280e-01_r8,9.722391e-01_r8,9.728258e-01_r8,9.733901e-01_r8,&
        & 9.739336e-01_r8,9.744579e-01_r8,9.749645e-01_r8,9.754546e-01_r8,9.759294e-01_r8,&
        & 9.763901e-01_r8,9.768376e-01_r8,9.772727e-01_r8 /)

! Hexagonal Ice Particle Parameterization
! extinction units (ext coef/iwc): [(m^-1)/(g m^-3)]
      extice3(:, 16) = (/ &
! band 16
        & 5.194013e-01_r8,3.215089e-01_r8,2.327917e-01_r8,1.824424e-01_r8,1.499977e-01_r8,&
        & 1.273492e-01_r8,1.106421e-01_r8,9.780982e-02_r8,8.764435e-02_r8,7.939266e-02_r8,&
        & 7.256081e-02_r8,6.681137e-02_r8,6.190600e-02_r8,5.767154e-02_r8,5.397915e-02_r8,&
        & 5.073102e-02_r8,4.785151e-02_r8,4.528125e-02_r8,4.297296e-02_r8,4.088853e-02_r8,&
        & 3.899690e-02_r8,3.727251e-02_r8,3.569411e-02_r8,3.424393e-02_r8,3.290694e-02_r8,&
        & 3.167040e-02_r8,3.052340e-02_r8,2.945654e-02_r8,2.846172e-02_r8,2.753188e-02_r8,&
        & 2.666085e-02_r8,2.584322e-02_r8,2.507423e-02_r8,2.434967e-02_r8,2.366579e-02_r8,&
        & 2.301926e-02_r8,2.240711e-02_r8,2.182666e-02_r8,2.127551e-02_r8,2.075150e-02_r8,&
        & 2.025267e-02_r8,1.977725e-02_r8,1.932364e-02_r8,1.889035e-02_r8,1.847607e-02_r8,&
        & 1.807956e-02_r8 /)
      extice3(:, 17) = (/ &
! band 17
        & 4.901155e-01_r8,3.065286e-01_r8,2.230800e-01_r8,1.753951e-01_r8,1.445402e-01_r8,&
        & 1.229417e-01_r8,1.069777e-01_r8,9.469760e-02_r8,8.495824e-02_r8,7.704501e-02_r8,&
        & 7.048834e-02_r8,6.496693e-02_r8,6.025353e-02_r8,5.618286e-02_r8,5.263186e-02_r8,&
        & 4.950698e-02_r8,4.673585e-02_r8,4.426164e-02_r8,4.203904e-02_r8,4.003153e-02_r8,&
        & 3.820932e-02_r8,3.654790e-02_r8,3.502688e-02_r8,3.362919e-02_r8,3.234041e-02_r8,&
        & 3.114829e-02_r8,3.004234e-02_r8,2.901356e-02_r8,2.805413e-02_r8,2.715727e-02_r8,&
        & 2.631705e-02_r8,2.552828e-02_r8,2.478637e-02_r8,2.408725e-02_r8,2.342734e-02_r8,&
        & 2.280343e-02_r8,2.221264e-02_r8,2.165242e-02_r8,2.112043e-02_r8,2.061461e-02_r8,&
        & 2.013308e-02_r8,1.967411e-02_r8,1.923616e-02_r8,1.881783e-02_r8,1.841781e-02_r8,&
        & 1.803494e-02_r8 /)
      extice3(:, 18) = (/ &
! band 18
        & 5.056264e-01_r8,3.160261e-01_r8,2.298442e-01_r8,1.805973e-01_r8,1.487318e-01_r8,&
        & 1.264258e-01_r8,1.099389e-01_r8,9.725656e-02_r8,8.719819e-02_r8,7.902576e-02_r8,&
        & 7.225433e-02_r8,6.655206e-02_r8,6.168427e-02_r8,5.748028e-02_r8,5.381296e-02_r8,&
        & 5.058572e-02_r8,4.772383e-02_r8,4.516857e-02_r8,4.287317e-02_r8,4.079990e-02_r8,&
        & 3.891801e-02_r8,3.720217e-02_r8,3.563133e-02_r8,3.418786e-02_r8,3.285686e-02_r8,&
        & 3.162569e-02_r8,3.048352e-02_r8,2.942104e-02_r8,2.843018e-02_r8,2.750395e-02_r8,&
        & 2.663621e-02_r8,2.582160e-02_r8,2.505539e-02_r8,2.433337e-02_r8,2.365185e-02_r8,&
        & 2.300750e-02_r8,2.239736e-02_r8,2.181878e-02_r8,2.126937e-02_r8,2.074699e-02_r8,&
        & 2.024968e-02_r8,1.977567e-02_r8,1.932338e-02_r8,1.889134e-02_r8,1.847823e-02_r8,&
        & 1.808281e-02_r8 /)
      extice3(:, 19) = (/ &
! band 19
        & 4.881605e-01_r8,3.055237e-01_r8,2.225070e-01_r8,1.750688e-01_r8,1.443736e-01_r8,&
        & 1.228869e-01_r8,1.070054e-01_r8,9.478893e-02_r8,8.509997e-02_r8,7.722769e-02_r8,&
        & 7.070495e-02_r8,6.521211e-02_r8,6.052311e-02_r8,5.647351e-02_r8,5.294088e-02_r8,&
        & 4.983217e-02_r8,4.707539e-02_r8,4.461398e-02_r8,4.240288e-02_r8,4.040575e-02_r8,&
        & 3.859298e-02_r8,3.694016e-02_r8,3.542701e-02_r8,3.403655e-02_r8,3.275444e-02_r8,&
        & 3.156849e-02_r8,3.046827e-02_r8,2.944481e-02_r8,2.849034e-02_r8,2.759812e-02_r8,&
        & 2.676226e-02_r8,2.597757e-02_r8,2.523949e-02_r8,2.454400e-02_r8,2.388750e-02_r8,&
        & 2.326682e-02_r8,2.267909e-02_r8,2.212176e-02_r8,2.159253e-02_r8,2.108933e-02_r8,&
        & 2.061028e-02_r8,2.015369e-02_r8,1.971801e-02_r8,1.930184e-02_r8,1.890389e-02_r8,&
        & 1.852300e-02_r8 /)
      extice3(:, 20) = (/ &
! band 20
        & 5.103703e-01_r8,3.188144e-01_r8,2.317435e-01_r8,1.819887e-01_r8,1.497944e-01_r8,&
        & 1.272584e-01_r8,1.106013e-01_r8,9.778822e-02_r8,8.762610e-02_r8,7.936938e-02_r8,&
        & 7.252809e-02_r8,6.676701e-02_r8,6.184901e-02_r8,5.760165e-02_r8,5.389651e-02_r8,&
        & 5.063598e-02_r8,4.774457e-02_r8,4.516295e-02_r8,4.284387e-02_r8,4.074922e-02_r8,&
        & 3.884792e-02_r8,3.711438e-02_r8,3.552734e-02_r8,3.406898e-02_r8,3.272425e-02_r8,&
        & 3.148038e-02_r8,3.032643e-02_r8,2.925299e-02_r8,2.825191e-02_r8,2.731612e-02_r8,&
        & 2.643943e-02_r8,2.561642e-02_r8,2.484230e-02_r8,2.411284e-02_r8,2.342429e-02_r8,&
        & 2.277329e-02_r8,2.215686e-02_r8,2.157231e-02_r8,2.101724e-02_r8,2.048946e-02_r8,&
        & 1.998702e-02_r8,1.950813e-02_r8,1.905118e-02_r8,1.861468e-02_r8,1.819730e-02_r8,&
        & 1.779781e-02_r8 /)
      extice3(:, 21) = (/ &
! band 21
        & 5.031161e-01_r8,3.144511e-01_r8,2.286942e-01_r8,1.796903e-01_r8,1.479819e-01_r8,&
        & 1.257860e-01_r8,1.093803e-01_r8,9.676059e-02_r8,8.675183e-02_r8,7.861971e-02_r8,&
        & 7.188168e-02_r8,6.620754e-02_r8,6.136376e-02_r8,5.718050e-02_r8,5.353127e-02_r8,&
        & 5.031995e-02_r8,4.747218e-02_r8,4.492952e-02_r8,4.264544e-02_r8,4.058240e-02_r8,&
        & 3.870979e-02_r8,3.700242e-02_r8,3.543933e-02_r8,3.400297e-02_r8,3.267854e-02_r8,&
        & 3.145345e-02_r8,3.031691e-02_r8,2.925967e-02_r8,2.827370e-02_r8,2.735203e-02_r8,&
        & 2.648858e-02_r8,2.567798e-02_r8,2.491555e-02_r8,2.419710e-02_r8,2.351893e-02_r8,&
        & 2.287776e-02_r8,2.227063e-02_r8,2.169491e-02_r8,2.114821e-02_r8,2.062840e-02_r8,&
        & 2.013354e-02_r8,1.966188e-02_r8,1.921182e-02_r8,1.878191e-02_r8,1.837083e-02_r8,&
        & 1.797737e-02_r8 /)
      extice3(:, 22) = (/ &
! band 22
        & 4.949453e-01_r8,3.095918e-01_r8,2.253402e-01_r8,1.771964e-01_r8,1.460446e-01_r8,&
        & 1.242383e-01_r8,1.081206e-01_r8,9.572235e-02_r8,8.588928e-02_r8,7.789990e-02_r8,&
        & 7.128013e-02_r8,6.570559e-02_r8,6.094684e-02_r8,5.683701e-02_r8,5.325183e-02_r8,&
        & 5.009688e-02_r8,4.729909e-02_r8,4.480106e-02_r8,4.255708e-02_r8,4.053025e-02_r8,&
        & 3.869051e-02_r8,3.701310e-02_r8,3.547745e-02_r8,3.406631e-02_r8,3.276512e-02_r8,&
        & 3.156153e-02_r8,3.044494e-02_r8,2.940626e-02_r8,2.843759e-02_r8,2.753211e-02_r8,&
        & 2.668381e-02_r8,2.588744e-02_r8,2.513839e-02_r8,2.443255e-02_r8,2.376629e-02_r8,&
        & 2.313637e-02_r8,2.253990e-02_r8,2.197428e-02_r8,2.143718e-02_r8,2.092649e-02_r8,&
        & 2.044032e-02_r8,1.997694e-02_r8,1.953478e-02_r8,1.911241e-02_r8,1.870855e-02_r8,&
        & 1.832199e-02_r8 /)
      extice3(:, 23) = (/ &
! band 23
        & 5.052816e-01_r8,3.157665e-01_r8,2.296233e-01_r8,1.803986e-01_r8,1.485473e-01_r8,&
        & 1.262514e-01_r8,1.097718e-01_r8,9.709524e-02_r8,8.704139e-02_r8,7.887264e-02_r8,&
        & 7.210424e-02_r8,6.640454e-02_r8,6.153894e-02_r8,5.733683e-02_r8,5.367116e-02_r8,&
        & 5.044537e-02_r8,4.758477e-02_r8,4.503066e-02_r8,4.273629e-02_r8,4.066395e-02_r8,&
        & 3.878291e-02_r8,3.706784e-02_r8,3.549771e-02_r8,3.405488e-02_r8,3.272448e-02_r8,&
        & 3.149387e-02_r8,3.035221e-02_r8,2.929020e-02_r8,2.829979e-02_r8,2.737397e-02_r8,&
        & 2.650663e-02_r8,2.569238e-02_r8,2.492651e-02_r8,2.420482e-02_r8,2.352361e-02_r8,&
        & 2.287954e-02_r8,2.226968e-02_r8,2.169136e-02_r8,2.114220e-02_r8,2.062005e-02_r8,&
        & 2.012296e-02_r8,1.964917e-02_r8,1.919709e-02_r8,1.876524e-02_r8,1.835231e-02_r8,&
        & 1.795707e-02_r8 /)
      extice3(:, 24) = (/ &
! band 24
        & 5.042067e-01_r8,3.151195e-01_r8,2.291708e-01_r8,1.800573e-01_r8,1.482779e-01_r8,&
        & 1.260324e-01_r8,1.095900e-01_r8,9.694202e-02_r8,8.691087e-02_r8,7.876056e-02_r8,&
        & 7.200745e-02_r8,6.632062e-02_r8,6.146600e-02_r8,5.727338e-02_r8,5.361599e-02_r8,&
        & 5.039749e-02_r8,4.754334e-02_r8,4.499500e-02_r8,4.270580e-02_r8,4.063815e-02_r8,&
        & 3.876135e-02_r8,3.705016e-02_r8,3.548357e-02_r8,3.404400e-02_r8,3.271661e-02_r8,&
        & 3.148877e-02_r8,3.034969e-02_r8,2.929008e-02_r8,2.830191e-02_r8,2.737818e-02_r8,&
        & 2.651279e-02_r8,2.570039e-02_r8,2.493624e-02_r8,2.421618e-02_r8,2.353650e-02_r8,&
        & 2.289390e-02_r8,2.228541e-02_r8,2.170840e-02_r8,2.116048e-02_r8,2.063950e-02_r8,&
        & 2.014354e-02_r8,1.967082e-02_r8,1.921975e-02_r8,1.878888e-02_r8,1.837688e-02_r8,&
        & 1.798254e-02_r8 /)
      extice3(:, 25) = (/ &
! band 25
        & 5.022507e-01_r8,3.139246e-01_r8,2.283218e-01_r8,1.794059e-01_r8,1.477544e-01_r8,&
        & 1.255984e-01_r8,1.092222e-01_r8,9.662516e-02_r8,8.663439e-02_r8,7.851688e-02_r8,&
        & 7.179095e-02_r8,6.612700e-02_r8,6.129193e-02_r8,5.711618e-02_r8,5.347351e-02_r8,&
        & 5.026796e-02_r8,4.742530e-02_r8,4.488721e-02_r8,4.260724e-02_r8,4.054790e-02_r8,&
        & 3.867866e-02_r8,3.697435e-02_r8,3.541407e-02_r8,3.398029e-02_r8,3.265824e-02_r8,&
        & 3.143535e-02_r8,3.030085e-02_r8,2.924551e-02_r8,2.826131e-02_r8,2.734130e-02_r8,&
        & 2.647939e-02_r8,2.567026e-02_r8,2.490919e-02_r8,2.419203e-02_r8,2.351509e-02_r8,&
        & 2.287507e-02_r8,2.226903e-02_r8,2.169434e-02_r8,2.114862e-02_r8,2.062975e-02_r8,&
        & 2.013578e-02_r8,1.966496e-02_r8,1.921571e-02_r8,1.878658e-02_r8,1.837623e-02_r8,&
        & 1.798348e-02_r8 /)
      extice3(:, 26) = (/ &
! band 26
        & 5.068316e-01_r8,3.166869e-01_r8,2.302576e-01_r8,1.808693e-01_r8,1.489122e-01_r8,&
        & 1.265423e-01_r8,1.100080e-01_r8,9.728926e-02_r8,8.720201e-02_r8,7.900612e-02_r8,&
        & 7.221524e-02_r8,6.649660e-02_r8,6.161484e-02_r8,5.739877e-02_r8,5.372093e-02_r8,&
        & 5.048442e-02_r8,4.761431e-02_r8,4.505172e-02_r8,4.274972e-02_r8,4.067050e-02_r8,&
        & 3.878321e-02_r8,3.706244e-02_r8,3.548710e-02_r8,3.403948e-02_r8,3.270466e-02_r8,&
        & 3.146995e-02_r8,3.032450e-02_r8,2.925897e-02_r8,2.826527e-02_r8,2.733638e-02_r8,&
        & 2.646615e-02_r8,2.564920e-02_r8,2.488078e-02_r8,2.415670e-02_r8,2.347322e-02_r8,&
        & 2.282702e-02_r8,2.221513e-02_r8,2.163489e-02_r8,2.108390e-02_r8,2.056002e-02_r8,&
        & 2.006128e-02_r8,1.958591e-02_r8,1.913232e-02_r8,1.869904e-02_r8,1.828474e-02_r8,&
        & 1.788819e-02_r8 /)
      extice3(:, 27) = (/ &
! band 27
        & 5.077707e-01_r8,3.172636e-01_r8,2.306695e-01_r8,1.811871e-01_r8,1.491691e-01_r8,&
        & 1.267565e-01_r8,1.101907e-01_r8,9.744773e-02_r8,8.734125e-02_r8,7.912973e-02_r8,&
        & 7.232591e-02_r8,6.659637e-02_r8,6.170530e-02_r8,5.748120e-02_r8,5.379634e-02_r8,&
        & 5.055367e-02_r8,4.767809e-02_r8,4.511061e-02_r8,4.280423e-02_r8,4.072104e-02_r8,&
        & 3.883015e-02_r8,3.710611e-02_r8,3.552776e-02_r8,3.407738e-02_r8,3.274002e-02_r8,&
        & 3.150296e-02_r8,3.035532e-02_r8,2.928776e-02_r8,2.829216e-02_r8,2.736150e-02_r8,&
        & 2.648961e-02_r8,2.567111e-02_r8,2.490123e-02_r8,2.417576e-02_r8,2.349098e-02_r8,&
        & 2.284354e-02_r8,2.223049e-02_r8,2.164914e-02_r8,2.109711e-02_r8,2.057222e-02_r8,&
        & 2.007253e-02_r8,1.959626e-02_r8,1.914181e-02_r8,1.870770e-02_r8,1.829261e-02_r8,&
        & 1.789531e-02_r8 /)
      extice3(:, 28) = (/ &
! band 28
        & 5.062281e-01_r8,3.163402e-01_r8,2.300275e-01_r8,1.807060e-01_r8,1.487921e-01_r8,&
        & 1.264523e-01_r8,1.099403e-01_r8,9.723879e-02_r8,8.716516e-02_r8,7.898034e-02_r8,&
        & 7.219863e-02_r8,6.648771e-02_r8,6.161254e-02_r8,5.740217e-02_r8,5.372929e-02_r8,&
        & 5.049716e-02_r8,4.763092e-02_r8,4.507179e-02_r8,4.277290e-02_r8,4.069649e-02_r8,&
        & 3.881175e-02_r8,3.709331e-02_r8,3.552008e-02_r8,3.407442e-02_r8,3.274141e-02_r8,&
        & 3.150837e-02_r8,3.036447e-02_r8,2.930037e-02_r8,2.830801e-02_r8,2.738037e-02_r8,&
        & 2.651132e-02_r8,2.569547e-02_r8,2.492810e-02_r8,2.420499e-02_r8,2.352243e-02_r8,&
        & 2.287710e-02_r8,2.226604e-02_r8,2.168658e-02_r8,2.113634e-02_r8,2.061316e-02_r8,&
        & 2.011510e-02_r8,1.964038e-02_r8,1.918740e-02_r8,1.875471e-02_r8,1.834096e-02_r8,&
        & 1.794495e-02_r8 /)
      extice3(:, 29) = (/ &
! band 29
        & 1.338834e-01_r8,1.924912e-01_r8,1.755523e-01_r8,1.534793e-01_r8,1.343937e-01_r8,&
        & 1.187883e-01_r8,1.060654e-01_r8,9.559106e-02_r8,8.685880e-02_r8,7.948698e-02_r8,&
        & 7.319086e-02_r8,6.775669e-02_r8,6.302215e-02_r8,5.886236e-02_r8,5.517996e-02_r8,&
        & 5.189810e-02_r8,4.895539e-02_r8,4.630225e-02_r8,4.389823e-02_r8,4.171002e-02_r8,&
        & 3.970998e-02_r8,3.787493e-02_r8,3.618537e-02_r8,3.462471e-02_r8,3.317880e-02_r8,&
        & 3.183547e-02_r8,3.058421e-02_r8,2.941590e-02_r8,2.832256e-02_r8,2.729724e-02_r8,&
        & 2.633377e-02_r8,2.542675e-02_r8,2.457136e-02_r8,2.376332e-02_r8,2.299882e-02_r8,&
        & 2.227443e-02_r8,2.158707e-02_r8,2.093400e-02_r8,2.031270e-02_r8,1.972091e-02_r8,&
        & 1.915659e-02_r8,1.861787e-02_r8,1.810304e-02_r8,1.761055e-02_r8,1.713899e-02_r8,&
        & 1.668704e-02_r8 /)

! single-scattering albedo: unitless
      ssaice3(:, 16) = (/ &
! band 16
        & 6.749442e-01_r8,6.649947e-01_r8,6.565828e-01_r8,6.489928e-01_r8,6.420046e-01_r8,&
        & 6.355231e-01_r8,6.294964e-01_r8,6.238901e-01_r8,6.186783e-01_r8,6.138395e-01_r8,&
        & 6.093543e-01_r8,6.052049e-01_r8,6.013742e-01_r8,5.978457e-01_r8,5.946030e-01_r8,&
        & 5.916302e-01_r8,5.889115e-01_r8,5.864310e-01_r8,5.841731e-01_r8,5.821221e-01_r8,&
        & 5.802624e-01_r8,5.785785e-01_r8,5.770549e-01_r8,5.756759e-01_r8,5.744262e-01_r8,&
        & 5.732901e-01_r8,5.722524e-01_r8,5.712974e-01_r8,5.704097e-01_r8,5.695739e-01_r8,&
        & 5.687747e-01_r8,5.679964e-01_r8,5.672238e-01_r8,5.664415e-01_r8,5.656340e-01_r8,&
        & 5.647860e-01_r8,5.638821e-01_r8,5.629070e-01_r8,5.618452e-01_r8,5.606815e-01_r8,&
        & 5.594006e-01_r8,5.579870e-01_r8,5.564255e-01_r8,5.547008e-01_r8,5.527976e-01_r8,&
        & 5.507005e-01_r8 /)
      ssaice3(:, 17) = (/ &
! band 17
        & 7.628550e-01_r8,7.567297e-01_r8,7.508463e-01_r8,7.451972e-01_r8,7.397745e-01_r8,&
        & 7.345705e-01_r8,7.295775e-01_r8,7.247881e-01_r8,7.201945e-01_r8,7.157894e-01_r8,&
        & 7.115652e-01_r8,7.075145e-01_r8,7.036300e-01_r8,6.999044e-01_r8,6.963304e-01_r8,&
        & 6.929007e-01_r8,6.896083e-01_r8,6.864460e-01_r8,6.834067e-01_r8,6.804833e-01_r8,&
        & 6.776690e-01_r8,6.749567e-01_r8,6.723397e-01_r8,6.698109e-01_r8,6.673637e-01_r8,&
        & 6.649913e-01_r8,6.626870e-01_r8,6.604441e-01_r8,6.582561e-01_r8,6.561163e-01_r8,&
        & 6.540182e-01_r8,6.519554e-01_r8,6.499215e-01_r8,6.479099e-01_r8,6.459145e-01_r8,&
        & 6.439289e-01_r8,6.419468e-01_r8,6.399621e-01_r8,6.379686e-01_r8,6.359601e-01_r8,&
        & 6.339306e-01_r8,6.318740e-01_r8,6.297845e-01_r8,6.276559e-01_r8,6.254825e-01_r8,&
        & 6.232583e-01_r8 /)
      ssaice3(:, 18) = (/ &
! band 18
        & 9.924147e-01_r8,9.882792e-01_r8,9.842257e-01_r8,9.802522e-01_r8,9.763566e-01_r8,&
        & 9.725367e-01_r8,9.687905e-01_r8,9.651157e-01_r8,9.615104e-01_r8,9.579725e-01_r8,&
        & 9.544997e-01_r8,9.510901e-01_r8,9.477416e-01_r8,9.444520e-01_r8,9.412194e-01_r8,&
        & 9.380415e-01_r8,9.349165e-01_r8,9.318421e-01_r8,9.288164e-01_r8,9.258373e-01_r8,&
        & 9.229027e-01_r8,9.200106e-01_r8,9.171589e-01_r8,9.143457e-01_r8,9.115688e-01_r8,&
        & 9.088263e-01_r8,9.061161e-01_r8,9.034362e-01_r8,9.007846e-01_r8,8.981592e-01_r8,&
        & 8.955581e-01_r8,8.929792e-01_r8,8.904206e-01_r8,8.878803e-01_r8,8.853562e-01_r8,&
        & 8.828464e-01_r8,8.803488e-01_r8,8.778616e-01_r8,8.753827e-01_r8,8.729102e-01_r8,&
        & 8.704421e-01_r8,8.679764e-01_r8,8.655112e-01_r8,8.630445e-01_r8,8.605744e-01_r8,&
        & 8.580989e-01_r8 /)
      ssaice3(:, 19) = (/ &
! band 19
        & 9.629413e-01_r8,9.517182e-01_r8,9.409209e-01_r8,9.305366e-01_r8,9.205529e-01_r8,&
        & 9.109569e-01_r8,9.017362e-01_r8,8.928780e-01_r8,8.843699e-01_r8,8.761992e-01_r8,&
        & 8.683536e-01_r8,8.608204e-01_r8,8.535873e-01_r8,8.466417e-01_r8,8.399712e-01_r8,&
        & 8.335635e-01_r8,8.274062e-01_r8,8.214868e-01_r8,8.157932e-01_r8,8.103129e-01_r8,&
        & 8.050336e-01_r8,7.999432e-01_r8,7.950294e-01_r8,7.902798e-01_r8,7.856825e-01_r8,&
        & 7.812250e-01_r8,7.768954e-01_r8,7.726815e-01_r8,7.685711e-01_r8,7.645522e-01_r8,&
        & 7.606126e-01_r8,7.567404e-01_r8,7.529234e-01_r8,7.491498e-01_r8,7.454074e-01_r8,&
        & 7.416844e-01_r8,7.379688e-01_r8,7.342485e-01_r8,7.305118e-01_r8,7.267468e-01_r8,&
        & 7.229415e-01_r8,7.190841e-01_r8,7.151628e-01_r8,7.111657e-01_r8,7.070811e-01_r8,&
        & 7.028972e-01_r8 /)
      ssaice3(:, 20) = (/ &
! band 20
        & 9.942270e-01_r8,9.909206e-01_r8,9.876775e-01_r8,9.844960e-01_r8,9.813746e-01_r8,&
        & 9.783114e-01_r8,9.753049e-01_r8,9.723535e-01_r8,9.694553e-01_r8,9.666088e-01_r8,&
        & 9.638123e-01_r8,9.610641e-01_r8,9.583626e-01_r8,9.557060e-01_r8,9.530928e-01_r8,&
        & 9.505211e-01_r8,9.479895e-01_r8,9.454961e-01_r8,9.430393e-01_r8,9.406174e-01_r8,&
        & 9.382288e-01_r8,9.358717e-01_r8,9.335446e-01_r8,9.312456e-01_r8,9.289731e-01_r8,&
        & 9.267255e-01_r8,9.245010e-01_r8,9.222980e-01_r8,9.201147e-01_r8,9.179496e-01_r8,&
        & 9.158008e-01_r8,9.136667e-01_r8,9.115457e-01_r8,9.094359e-01_r8,9.073358e-01_r8,&
        & 9.052436e-01_r8,9.031577e-01_r8,9.010763e-01_r8,8.989977e-01_r8,8.969203e-01_r8,&
        & 8.948423e-01_r8,8.927620e-01_r8,8.906778e-01_r8,8.885879e-01_r8,8.864907e-01_r8,&
        & 8.843843e-01_r8 /)
      ssaice3(:, 21) = (/ &
! band 21
        & 9.934014e-01_r8,9.899331e-01_r8,9.865537e-01_r8,9.832610e-01_r8,9.800523e-01_r8,&
        & 9.769254e-01_r8,9.738777e-01_r8,9.709069e-01_r8,9.680106e-01_r8,9.651862e-01_r8,&
        & 9.624315e-01_r8,9.597439e-01_r8,9.571212e-01_r8,9.545608e-01_r8,9.520605e-01_r8,&
        & 9.496177e-01_r8,9.472301e-01_r8,9.448954e-01_r8,9.426111e-01_r8,9.403749e-01_r8,&
        & 9.381843e-01_r8,9.360370e-01_r8,9.339307e-01_r8,9.318629e-01_r8,9.298313e-01_r8,&
        & 9.278336e-01_r8,9.258673e-01_r8,9.239302e-01_r8,9.220198e-01_r8,9.201338e-01_r8,&
        & 9.182700e-01_r8,9.164258e-01_r8,9.145991e-01_r8,9.127874e-01_r8,9.109884e-01_r8,&
        & 9.091999e-01_r8,9.074194e-01_r8,9.056447e-01_r8,9.038735e-01_r8,9.021033e-01_r8,&
        & 9.003320e-01_r8,8.985572e-01_r8,8.967766e-01_r8,8.949879e-01_r8,8.931888e-01_r8,&
        & 8.913770e-01_r8 /)
      ssaice3(:, 22) = (/ &
! band 22
        & 9.994833e-01_r8,9.992055e-01_r8,9.989278e-01_r8,9.986500e-01_r8,9.983724e-01_r8,&
        & 9.980947e-01_r8,9.978172e-01_r8,9.975397e-01_r8,9.972623e-01_r8,9.969849e-01_r8,&
        & 9.967077e-01_r8,9.964305e-01_r8,9.961535e-01_r8,9.958765e-01_r8,9.955997e-01_r8,&
        & 9.953230e-01_r8,9.950464e-01_r8,9.947699e-01_r8,9.944936e-01_r8,9.942174e-01_r8,&
        & 9.939414e-01_r8,9.936656e-01_r8,9.933899e-01_r8,9.931144e-01_r8,9.928390e-01_r8,&
        & 9.925639e-01_r8,9.922889e-01_r8,9.920141e-01_r8,9.917396e-01_r8,9.914652e-01_r8,&
        & 9.911911e-01_r8,9.909171e-01_r8,9.906434e-01_r8,9.903700e-01_r8,9.900967e-01_r8,&
        & 9.898237e-01_r8,9.895510e-01_r8,9.892784e-01_r8,9.890062e-01_r8,9.887342e-01_r8,&
        & 9.884625e-01_r8,9.881911e-01_r8,9.879199e-01_r8,9.876490e-01_r8,9.873784e-01_r8,&
        & 9.871081e-01_r8 /)
      ssaice3(:, 23) = (/ &
! band 23
        & 9.999343e-01_r8,9.998917e-01_r8,9.998492e-01_r8,9.998067e-01_r8,9.997642e-01_r8,&
        & 9.997218e-01_r8,9.996795e-01_r8,9.996372e-01_r8,9.995949e-01_r8,9.995528e-01_r8,&
        & 9.995106e-01_r8,9.994686e-01_r8,9.994265e-01_r8,9.993845e-01_r8,9.993426e-01_r8,&
        & 9.993007e-01_r8,9.992589e-01_r8,9.992171e-01_r8,9.991754e-01_r8,9.991337e-01_r8,&
        & 9.990921e-01_r8,9.990505e-01_r8,9.990089e-01_r8,9.989674e-01_r8,9.989260e-01_r8,&
        & 9.988846e-01_r8,9.988432e-01_r8,9.988019e-01_r8,9.987606e-01_r8,9.987194e-01_r8,&
        & 9.986782e-01_r8,9.986370e-01_r8,9.985959e-01_r8,9.985549e-01_r8,9.985139e-01_r8,&
        & 9.984729e-01_r8,9.984319e-01_r8,9.983910e-01_r8,9.983502e-01_r8,9.983094e-01_r8,&
        & 9.982686e-01_r8,9.982279e-01_r8,9.981872e-01_r8,9.981465e-01_r8,9.981059e-01_r8,&
        & 9.980653e-01_r8 /)
      ssaice3(:, 24) = (/ &
! band 24
        & 9.999978e-01_r8,9.999965e-01_r8,9.999952e-01_r8,9.999939e-01_r8,9.999926e-01_r8,&
        & 9.999913e-01_r8,9.999900e-01_r8,9.999887e-01_r8,9.999873e-01_r8,9.999860e-01_r8,&
        & 9.999847e-01_r8,9.999834e-01_r8,9.999821e-01_r8,9.999808e-01_r8,9.999795e-01_r8,&
        & 9.999782e-01_r8,9.999769e-01_r8,9.999756e-01_r8,9.999743e-01_r8,9.999730e-01_r8,&
        & 9.999717e-01_r8,9.999704e-01_r8,9.999691e-01_r8,9.999678e-01_r8,9.999665e-01_r8,&
        & 9.999652e-01_r8,9.999639e-01_r8,9.999626e-01_r8,9.999613e-01_r8,9.999600e-01_r8,&
        & 9.999587e-01_r8,9.999574e-01_r8,9.999561e-01_r8,9.999548e-01_r8,9.999535e-01_r8,&
        & 9.999522e-01_r8,9.999509e-01_r8,9.999496e-01_r8,9.999483e-01_r8,9.999470e-01_r8,&
        & 9.999457e-01_r8,9.999444e-01_r8,9.999431e-01_r8,9.999418e-01_r8,9.999405e-01_r8,&
        & 9.999392e-01_r8 /)
      ssaice3(:, 25) = (/ &
! band 25
        & 9.999994e-01_r8,9.999993e-01_r8,9.999991e-01_r8,9.999990e-01_r8,9.999989e-01_r8,&
        & 9.999987e-01_r8,9.999986e-01_r8,9.999984e-01_r8,9.999983e-01_r8,9.999982e-01_r8,&
        & 9.999980e-01_r8,9.999979e-01_r8,9.999977e-01_r8,9.999976e-01_r8,9.999975e-01_r8,&
        & 9.999973e-01_r8,9.999972e-01_r8,9.999970e-01_r8,9.999969e-01_r8,9.999967e-01_r8,&
        & 9.999966e-01_r8,9.999965e-01_r8,9.999963e-01_r8,9.999962e-01_r8,9.999960e-01_r8,&
        & 9.999959e-01_r8,9.999957e-01_r8,9.999956e-01_r8,9.999954e-01_r8,9.999953e-01_r8,&
        & 9.999952e-01_r8,9.999950e-01_r8,9.999949e-01_r8,9.999947e-01_r8,9.999946e-01_r8,&
        & 9.999944e-01_r8,9.999943e-01_r8,9.999941e-01_r8,9.999940e-01_r8,9.999939e-01_r8,&
        & 9.999937e-01_r8,9.999936e-01_r8,9.999934e-01_r8,9.999933e-01_r8,9.999931e-01_r8,&
        & 9.999930e-01_r8 /)
      ssaice3(:, 26) = (/ &
! band 26
        & 9.999997e-01_r8,9.999995e-01_r8,9.999992e-01_r8,9.999990e-01_r8,9.999987e-01_r8,&
        & 9.999985e-01_r8,9.999983e-01_r8,9.999980e-01_r8,9.999978e-01_r8,9.999976e-01_r8,&
        & 9.999973e-01_r8,9.999971e-01_r8,9.999969e-01_r8,9.999967e-01_r8,9.999965e-01_r8,&
        & 9.999963e-01_r8,9.999960e-01_r8,9.999958e-01_r8,9.999956e-01_r8,9.999954e-01_r8,&
        & 9.999952e-01_r8,9.999950e-01_r8,9.999948e-01_r8,9.999946e-01_r8,9.999944e-01_r8,&
        & 9.999942e-01_r8,9.999939e-01_r8,9.999937e-01_r8,9.999935e-01_r8,9.999933e-01_r8,&
        & 9.999931e-01_r8,9.999929e-01_r8,9.999927e-01_r8,9.999925e-01_r8,9.999923e-01_r8,&
        & 9.999920e-01_r8,9.999918e-01_r8,9.999916e-01_r8,9.999914e-01_r8,9.999911e-01_r8,&
        & 9.999909e-01_r8,9.999907e-01_r8,9.999905e-01_r8,9.999902e-01_r8,9.999900e-01_r8,&
        & 9.999897e-01_r8 /)
      ssaice3(:, 27) = (/ &
! band 27
        & 9.999991e-01_r8,9.999985e-01_r8,9.999980e-01_r8,9.999974e-01_r8,9.999968e-01_r8,&
        & 9.999963e-01_r8,9.999957e-01_r8,9.999951e-01_r8,9.999946e-01_r8,9.999940e-01_r8,&
        & 9.999934e-01_r8,9.999929e-01_r8,9.999923e-01_r8,9.999918e-01_r8,9.999912e-01_r8,&
        & 9.999907e-01_r8,9.999901e-01_r8,9.999896e-01_r8,9.999891e-01_r8,9.999885e-01_r8,&
        & 9.999880e-01_r8,9.999874e-01_r8,9.999869e-01_r8,9.999863e-01_r8,9.999858e-01_r8,&
        & 9.999853e-01_r8,9.999847e-01_r8,9.999842e-01_r8,9.999836e-01_r8,9.999831e-01_r8,&
        & 9.999826e-01_r8,9.999820e-01_r8,9.999815e-01_r8,9.999809e-01_r8,9.999804e-01_r8,&
        & 9.999798e-01_r8,9.999793e-01_r8,9.999787e-01_r8,9.999782e-01_r8,9.999776e-01_r8,&
        & 9.999770e-01_r8,9.999765e-01_r8,9.999759e-01_r8,9.999754e-01_r8,9.999748e-01_r8,&
        & 9.999742e-01_r8 /)
      ssaice3(:, 28) = (/ &
! band 28
        & 9.999975e-01_r8,9.999961e-01_r8,9.999946e-01_r8,9.999931e-01_r8,9.999917e-01_r8,&
        & 9.999903e-01_r8,9.999888e-01_r8,9.999874e-01_r8,9.999859e-01_r8,9.999845e-01_r8,&
        & 9.999831e-01_r8,9.999816e-01_r8,9.999802e-01_r8,9.999788e-01_r8,9.999774e-01_r8,&
        & 9.999759e-01_r8,9.999745e-01_r8,9.999731e-01_r8,9.999717e-01_r8,9.999702e-01_r8,&
        & 9.999688e-01_r8,9.999674e-01_r8,9.999660e-01_r8,9.999646e-01_r8,9.999631e-01_r8,&
        & 9.999617e-01_r8,9.999603e-01_r8,9.999589e-01_r8,9.999574e-01_r8,9.999560e-01_r8,&
        & 9.999546e-01_r8,9.999532e-01_r8,9.999517e-01_r8,9.999503e-01_r8,9.999489e-01_r8,&
        & 9.999474e-01_r8,9.999460e-01_r8,9.999446e-01_r8,9.999431e-01_r8,9.999417e-01_r8,&
        & 9.999403e-01_r8,9.999388e-01_r8,9.999374e-01_r8,9.999359e-01_r8,9.999345e-01_r8,&
        & 9.999330e-01_r8 /)
      ssaice3(:, 29) = (/ &
! band 29
        & 4.526500e-01_r8,5.287890e-01_r8,5.410487e-01_r8,5.459865e-01_r8,5.485149e-01_r8,&
        & 5.498914e-01_r8,5.505895e-01_r8,5.508310e-01_r8,5.507364e-01_r8,5.503793e-01_r8,&
        & 5.498090e-01_r8,5.490612e-01_r8,5.481637e-01_r8,5.471395e-01_r8,5.460083e-01_r8,&
        & 5.447878e-01_r8,5.434946e-01_r8,5.421442e-01_r8,5.407514e-01_r8,5.393309e-01_r8,&
        & 5.378970e-01_r8,5.364641e-01_r8,5.350464e-01_r8,5.336582e-01_r8,5.323140e-01_r8,&
        & 5.310283e-01_r8,5.298158e-01_r8,5.286914e-01_r8,5.276704e-01_r8,5.267680e-01_r8,&
        & 5.260000e-01_r8,5.253823e-01_r8,5.249311e-01_r8,5.246629e-01_r8,5.245946e-01_r8,&
        & 5.247434e-01_r8,5.251268e-01_r8,5.257626e-01_r8,5.266693e-01_r8,5.278653e-01_r8,&
        & 5.293698e-01_r8,5.312022e-01_r8,5.333823e-01_r8,5.359305e-01_r8,5.388676e-01_r8,&
        & 5.422146e-01_r8 /)

! asymmetry factor: unitless
      asyice3(:, 16) = (/ &
! band 16
        & 8.340752e-01_r8,8.435170e-01_r8,8.517487e-01_r8,8.592064e-01_r8,8.660387e-01_r8,&
        & 8.723204e-01_r8,8.780997e-01_r8,8.834137e-01_r8,8.882934e-01_r8,8.927662e-01_r8,&
        & 8.968577e-01_r8,9.005914e-01_r8,9.039899e-01_r8,9.070745e-01_r8,9.098659e-01_r8,&
        & 9.123836e-01_r8,9.146466e-01_r8,9.166734e-01_r8,9.184817e-01_r8,9.200886e-01_r8,&
        & 9.215109e-01_r8,9.227648e-01_r8,9.238661e-01_r8,9.248304e-01_r8,9.256727e-01_r8,&
        & 9.264078e-01_r8,9.270505e-01_r8,9.276150e-01_r8,9.281156e-01_r8,9.285662e-01_r8,&
        & 9.289806e-01_r8,9.293726e-01_r8,9.297557e-01_r8,9.301435e-01_r8,9.305491e-01_r8,&
        & 9.309859e-01_r8,9.314671e-01_r8,9.320055e-01_r8,9.326140e-01_r8,9.333053e-01_r8,&
        & 9.340919e-01_r8,9.349861e-01_r8,9.360000e-01_r8,9.371451e-01_r8,9.384329e-01_r8,&
        & 9.398744e-01_r8 /)
      asyice3(:, 17) = (/ &
! band 17
        & 8.728160e-01_r8,8.777333e-01_r8,8.823754e-01_r8,8.867535e-01_r8,8.908785e-01_r8,&
        & 8.947611e-01_r8,8.984118e-01_r8,9.018408e-01_r8,9.050582e-01_r8,9.080739e-01_r8,&
        & 9.108976e-01_r8,9.135388e-01_r8,9.160068e-01_r8,9.183106e-01_r8,9.204595e-01_r8,&
        & 9.224620e-01_r8,9.243271e-01_r8,9.260632e-01_r8,9.276788e-01_r8,9.291822e-01_r8,&
        & 9.305817e-01_r8,9.318853e-01_r8,9.331012e-01_r8,9.342372e-01_r8,9.353013e-01_r8,&
        & 9.363013e-01_r8,9.372450e-01_r8,9.381400e-01_r8,9.389939e-01_r8,9.398145e-01_r8,&
        & 9.406092e-01_r8,9.413856e-01_r8,9.421511e-01_r8,9.429131e-01_r8,9.436790e-01_r8,&
        & 9.444561e-01_r8,9.452517e-01_r8,9.460729e-01_r8,9.469270e-01_r8,9.478209e-01_r8,&
        & 9.487617e-01_r8,9.497562e-01_r8,9.508112e-01_r8,9.519335e-01_r8,9.531294e-01_r8,&
        & 9.544055e-01_r8 /)
      asyice3(:, 18) = (/ &
! band 18
        & 7.897566e-01_r8,7.948704e-01_r8,7.998041e-01_r8,8.045623e-01_r8,8.091495e-01_r8,&
        & 8.135702e-01_r8,8.178290e-01_r8,8.219305e-01_r8,8.258790e-01_r8,8.296792e-01_r8,&
        & 8.333355e-01_r8,8.368524e-01_r8,8.402343e-01_r8,8.434856e-01_r8,8.466108e-01_r8,&
        & 8.496143e-01_r8,8.525004e-01_r8,8.552737e-01_r8,8.579384e-01_r8,8.604990e-01_r8,&
        & 8.629597e-01_r8,8.653250e-01_r8,8.675992e-01_r8,8.697867e-01_r8,8.718916e-01_r8,&
        & 8.739185e-01_r8,8.758715e-01_r8,8.777551e-01_r8,8.795734e-01_r8,8.813308e-01_r8,&
        & 8.830315e-01_r8,8.846799e-01_r8,8.862802e-01_r8,8.878366e-01_r8,8.893534e-01_r8,&
        & 8.908350e-01_r8,8.922854e-01_r8,8.937090e-01_r8,8.951099e-01_r8,8.964925e-01_r8,&
        & 8.978609e-01_r8,8.992192e-01_r8,9.005718e-01_r8,9.019229e-01_r8,9.032765e-01_r8,&
        & 9.046369e-01_r8 /)
      asyice3(:, 19) = (/ &
! band 19
        & 7.812615e-01_r8,7.887764e-01_r8,7.959664e-01_r8,8.028413e-01_r8,8.094109e-01_r8,&
        & 8.156849e-01_r8,8.216730e-01_r8,8.273846e-01_r8,8.328294e-01_r8,8.380166e-01_r8,&
        & 8.429556e-01_r8,8.476556e-01_r8,8.521258e-01_r8,8.563753e-01_r8,8.604131e-01_r8,&
        & 8.642481e-01_r8,8.678893e-01_r8,8.713455e-01_r8,8.746254e-01_r8,8.777378e-01_r8,&
        & 8.806914e-01_r8,8.834948e-01_r8,8.861566e-01_r8,8.886854e-01_r8,8.910897e-01_r8,&
        & 8.933779e-01_r8,8.955586e-01_r8,8.976402e-01_r8,8.996311e-01_r8,9.015398e-01_r8,&
        & 9.033745e-01_r8,9.051436e-01_r8,9.068555e-01_r8,9.085185e-01_r8,9.101410e-01_r8,&
        & 9.117311e-01_r8,9.132972e-01_r8,9.148476e-01_r8,9.163905e-01_r8,9.179340e-01_r8,&
        & 9.194864e-01_r8,9.210559e-01_r8,9.226505e-01_r8,9.242784e-01_r8,9.259476e-01_r8,&
        & 9.276661e-01_r8 /)
      asyice3(:, 20) = (/ &
! band 20
        & 7.640720e-01_r8,7.691119e-01_r8,7.739941e-01_r8,7.787222e-01_r8,7.832998e-01_r8,&
        & 7.877304e-01_r8,7.920177e-01_r8,7.961652e-01_r8,8.001765e-01_r8,8.040551e-01_r8,&
        & 8.078044e-01_r8,8.114280e-01_r8,8.149294e-01_r8,8.183119e-01_r8,8.215791e-01_r8,&
        & 8.247344e-01_r8,8.277812e-01_r8,8.307229e-01_r8,8.335629e-01_r8,8.363046e-01_r8,&
        & 8.389514e-01_r8,8.415067e-01_r8,8.439738e-01_r8,8.463560e-01_r8,8.486568e-01_r8,&
        & 8.508795e-01_r8,8.530274e-01_r8,8.551039e-01_r8,8.571122e-01_r8,8.590558e-01_r8,&
        & 8.609378e-01_r8,8.627618e-01_r8,8.645309e-01_r8,8.662485e-01_r8,8.679178e-01_r8,&
        & 8.695423e-01_r8,8.711251e-01_r8,8.726697e-01_r8,8.741792e-01_r8,8.756571e-01_r8,&
        & 8.771065e-01_r8,8.785307e-01_r8,8.799331e-01_r8,8.813169e-01_r8,8.826854e-01_r8,&
        & 8.840419e-01_r8 /)
      asyice3(:, 21) = (/ &
! band 21
        & 7.602598e-01_r8,7.651572e-01_r8,7.699014e-01_r8,7.744962e-01_r8,7.789452e-01_r8,&
        & 7.832522e-01_r8,7.874205e-01_r8,7.914538e-01_r8,7.953555e-01_r8,7.991290e-01_r8,&
        & 8.027777e-01_r8,8.063049e-01_r8,8.097140e-01_r8,8.130081e-01_r8,8.161906e-01_r8,&
        & 8.192645e-01_r8,8.222331e-01_r8,8.250993e-01_r8,8.278664e-01_r8,8.305374e-01_r8,&
        & 8.331153e-01_r8,8.356030e-01_r8,8.380037e-01_r8,8.403201e-01_r8,8.425553e-01_r8,&
        & 8.447121e-01_r8,8.467935e-01_r8,8.488022e-01_r8,8.507412e-01_r8,8.526132e-01_r8,&
        & 8.544210e-01_r8,8.561675e-01_r8,8.578554e-01_r8,8.594875e-01_r8,8.610665e-01_r8,&
        & 8.625951e-01_r8,8.640760e-01_r8,8.655119e-01_r8,8.669055e-01_r8,8.682594e-01_r8,&
        & 8.695763e-01_r8,8.708587e-01_r8,8.721094e-01_r8,8.733308e-01_r8,8.745255e-01_r8,&
        & 8.756961e-01_r8 /)
      asyice3(:, 22) = (/ &
! band 22
        & 7.568957e-01_r8,7.606995e-01_r8,7.644072e-01_r8,7.680204e-01_r8,7.715402e-01_r8,&
        & 7.749682e-01_r8,7.783057e-01_r8,7.815541e-01_r8,7.847148e-01_r8,7.877892e-01_r8,&
        & 7.907786e-01_r8,7.936846e-01_r8,7.965084e-01_r8,7.992515e-01_r8,8.019153e-01_r8,&
        & 8.045011e-01_r8,8.070103e-01_r8,8.094444e-01_r8,8.118048e-01_r8,8.140927e-01_r8,&
        & 8.163097e-01_r8,8.184571e-01_r8,8.205364e-01_r8,8.225488e-01_r8,8.244958e-01_r8,&
        & 8.263789e-01_r8,8.281993e-01_r8,8.299586e-01_r8,8.316580e-01_r8,8.332991e-01_r8,&
        & 8.348831e-01_r8,8.364115e-01_r8,8.378857e-01_r8,8.393071e-01_r8,8.406770e-01_r8,&
        & 8.419969e-01_r8,8.432682e-01_r8,8.444923e-01_r8,8.456706e-01_r8,8.468044e-01_r8,&
        & 8.478952e-01_r8,8.489444e-01_r8,8.499533e-01_r8,8.509234e-01_r8,8.518561e-01_r8,&
        & 8.527528e-01_r8 /)
      asyice3(:, 23) = (/ &
! band 23
        & 7.575066e-01_r8,7.606912e-01_r8,7.638236e-01_r8,7.669035e-01_r8,7.699306e-01_r8,&
        & 7.729046e-01_r8,7.758254e-01_r8,7.786926e-01_r8,7.815060e-01_r8,7.842654e-01_r8,&
        & 7.869705e-01_r8,7.896211e-01_r8,7.922168e-01_r8,7.947574e-01_r8,7.972428e-01_r8,&
        & 7.996726e-01_r8,8.020466e-01_r8,8.043646e-01_r8,8.066262e-01_r8,8.088313e-01_r8,&
        & 8.109796e-01_r8,8.130709e-01_r8,8.151049e-01_r8,8.170814e-01_r8,8.190001e-01_r8,&
        & 8.208608e-01_r8,8.226632e-01_r8,8.244071e-01_r8,8.260924e-01_r8,8.277186e-01_r8,&
        & 8.292856e-01_r8,8.307932e-01_r8,8.322411e-01_r8,8.336291e-01_r8,8.349570e-01_r8,&
        & 8.362244e-01_r8,8.374312e-01_r8,8.385772e-01_r8,8.396621e-01_r8,8.406856e-01_r8,&
        & 8.416476e-01_r8,8.425479e-01_r8,8.433861e-01_r8,8.441620e-01_r8,8.448755e-01_r8,&
        & 8.455263e-01_r8 /)
      asyice3(:, 24) = (/ &
! band 24
        & 7.568829e-01_r8,7.597947e-01_r8,7.626745e-01_r8,7.655212e-01_r8,7.683337e-01_r8,&
        & 7.711111e-01_r8,7.738523e-01_r8,7.765565e-01_r8,7.792225e-01_r8,7.818494e-01_r8,&
        & 7.844362e-01_r8,7.869819e-01_r8,7.894854e-01_r8,7.919459e-01_r8,7.943623e-01_r8,&
        & 7.967337e-01_r8,7.990590e-01_r8,8.013373e-01_r8,8.035676e-01_r8,8.057488e-01_r8,&
        & 8.078802e-01_r8,8.099605e-01_r8,8.119890e-01_r8,8.139645e-01_r8,8.158862e-01_r8,&
        & 8.177530e-01_r8,8.195641e-01_r8,8.213183e-01_r8,8.230149e-01_r8,8.246527e-01_r8,&
        & 8.262308e-01_r8,8.277483e-01_r8,8.292042e-01_r8,8.305976e-01_r8,8.319275e-01_r8,&
        & 8.331929e-01_r8,8.343929e-01_r8,8.355265e-01_r8,8.365928e-01_r8,8.375909e-01_r8,&
        & 8.385197e-01_r8,8.393784e-01_r8,8.401659e-01_r8,8.408815e-01_r8,8.415240e-01_r8,&
        & 8.420926e-01_r8 /)
      asyice3(:, 25) = (/ &
! band 25
        & 7.548616e-01_r8,7.575454e-01_r8,7.602153e-01_r8,7.628696e-01_r8,7.655067e-01_r8,&
        & 7.681249e-01_r8,7.707225e-01_r8,7.732978e-01_r8,7.758492e-01_r8,7.783750e-01_r8,&
        & 7.808735e-01_r8,7.833430e-01_r8,7.857819e-01_r8,7.881886e-01_r8,7.905612e-01_r8,&
        & 7.928983e-01_r8,7.951980e-01_r8,7.974588e-01_r8,7.996789e-01_r8,8.018567e-01_r8,&
        & 8.039905e-01_r8,8.060787e-01_r8,8.081196e-01_r8,8.101115e-01_r8,8.120527e-01_r8,&
        & 8.139416e-01_r8,8.157764e-01_r8,8.175557e-01_r8,8.192776e-01_r8,8.209405e-01_r8,&
        & 8.225427e-01_r8,8.240826e-01_r8,8.255585e-01_r8,8.269688e-01_r8,8.283117e-01_r8,&
        & 8.295856e-01_r8,8.307889e-01_r8,8.319198e-01_r8,8.329767e-01_r8,8.339579e-01_r8,&
        & 8.348619e-01_r8,8.356868e-01_r8,8.364311e-01_r8,8.370930e-01_r8,8.376710e-01_r8,&
        & 8.381633e-01_r8 /)
      asyice3(:, 26) = (/ &
! band 26
        & 7.491854e-01_r8,7.518523e-01_r8,7.545089e-01_r8,7.571534e-01_r8,7.597839e-01_r8,&
        & 7.623987e-01_r8,7.649959e-01_r8,7.675737e-01_r8,7.701303e-01_r8,7.726639e-01_r8,&
        & 7.751727e-01_r8,7.776548e-01_r8,7.801084e-01_r8,7.825318e-01_r8,7.849230e-01_r8,&
        & 7.872804e-01_r8,7.896020e-01_r8,7.918862e-01_r8,7.941309e-01_r8,7.963345e-01_r8,&
        & 7.984951e-01_r8,8.006109e-01_r8,8.026802e-01_r8,8.047009e-01_r8,8.066715e-01_r8,&
        & 8.085900e-01_r8,8.104546e-01_r8,8.122636e-01_r8,8.140150e-01_r8,8.157072e-01_r8,&
        & 8.173382e-01_r8,8.189063e-01_r8,8.204096e-01_r8,8.218464e-01_r8,8.232148e-01_r8,&
        & 8.245130e-01_r8,8.257391e-01_r8,8.268915e-01_r8,8.279682e-01_r8,8.289675e-01_r8,&
        & 8.298875e-01_r8,8.307264e-01_r8,8.314824e-01_r8,8.321537e-01_r8,8.327385e-01_r8,&
        & 8.332350e-01_r8 /)
      asyice3(:, 27) = (/ &
! band 27
        & 7.397086e-01_r8,7.424069e-01_r8,7.450955e-01_r8,7.477725e-01_r8,7.504362e-01_r8,&
        & 7.530846e-01_r8,7.557159e-01_r8,7.583283e-01_r8,7.609199e-01_r8,7.634888e-01_r8,&
        & 7.660332e-01_r8,7.685512e-01_r8,7.710411e-01_r8,7.735009e-01_r8,7.759288e-01_r8,&
        & 7.783229e-01_r8,7.806814e-01_r8,7.830024e-01_r8,7.852841e-01_r8,7.875246e-01_r8,&
        & 7.897221e-01_r8,7.918748e-01_r8,7.939807e-01_r8,7.960380e-01_r8,7.980449e-01_r8,&
        & 7.999995e-01_r8,8.019000e-01_r8,8.037445e-01_r8,8.055311e-01_r8,8.072581e-01_r8,&
        & 8.089235e-01_r8,8.105255e-01_r8,8.120623e-01_r8,8.135319e-01_r8,8.149326e-01_r8,&
        & 8.162626e-01_r8,8.175198e-01_r8,8.187025e-01_r8,8.198089e-01_r8,8.208371e-01_r8,&
        & 8.217852e-01_r8,8.226514e-01_r8,8.234338e-01_r8,8.241306e-01_r8,8.247399e-01_r8,&
        & 8.252599e-01_r8 /)
      asyice3(:, 28) = (/ &
! band 28
        & 7.224533e-01_r8,7.251681e-01_r8,7.278728e-01_r8,7.305654e-01_r8,7.332444e-01_r8,&
        & 7.359078e-01_r8,7.385539e-01_r8,7.411808e-01_r8,7.437869e-01_r8,7.463702e-01_r8,&
        & 7.489291e-01_r8,7.514616e-01_r8,7.539661e-01_r8,7.564408e-01_r8,7.588837e-01_r8,&
        & 7.612933e-01_r8,7.636676e-01_r8,7.660049e-01_r8,7.683034e-01_r8,7.705612e-01_r8,&
        & 7.727767e-01_r8,7.749480e-01_r8,7.770733e-01_r8,7.791509e-01_r8,7.811789e-01_r8,&
        & 7.831556e-01_r8,7.850791e-01_r8,7.869478e-01_r8,7.887597e-01_r8,7.905131e-01_r8,&
        & 7.922062e-01_r8,7.938372e-01_r8,7.954044e-01_r8,7.969059e-01_r8,7.983399e-01_r8,&
        & 7.997047e-01_r8,8.009985e-01_r8,8.022195e-01_r8,8.033658e-01_r8,8.044357e-01_r8,&
        & 8.054275e-01_r8,8.063392e-01_r8,8.071692e-01_r8,8.079157e-01_r8,8.085768e-01_r8,&
        & 8.091507e-01_r8 /)
      asyice3(:, 29) = (/ &
! band 29
        & 8.850026e-01_r8,9.005489e-01_r8,9.069242e-01_r8,9.121799e-01_r8,9.168987e-01_r8,&
        & 9.212259e-01_r8,9.252176e-01_r8,9.289028e-01_r8,9.323000e-01_r8,9.354235e-01_r8,&
        & 9.382858e-01_r8,9.408985e-01_r8,9.432734e-01_r8,9.454218e-01_r8,9.473557e-01_r8,&
        & 9.490871e-01_r8,9.506282e-01_r8,9.519917e-01_r8,9.531904e-01_r8,9.542374e-01_r8,&
        & 9.551461e-01_r8,9.559298e-01_r8,9.566023e-01_r8,9.571775e-01_r8,9.576692e-01_r8,&
        & 9.580916e-01_r8,9.584589e-01_r8,9.587853e-01_r8,9.590851e-01_r8,9.593729e-01_r8,&
        & 9.596632e-01_r8,9.599705e-01_r8,9.603096e-01_r8,9.606954e-01_r8,9.611427e-01_r8,&
        & 9.616667e-01_r8,9.622826e-01_r8,9.630060e-01_r8,9.638524e-01_r8,9.648379e-01_r8,&
        & 9.659788e-01_r8,9.672916e-01_r8,9.687933e-01_r8,9.705014e-01_r8,9.724337e-01_r8,&
        & 9.746084e-01_r8 /)

! fdelta: unitless
      fdlice3(:, 16) = (/ &
! band 16
        & 4.959277e-02_r8,4.685292e-02_r8,4.426104e-02_r8,4.181231e-02_r8,3.950191e-02_r8,&
        & 3.732500e-02_r8,3.527675e-02_r8,3.335235e-02_r8,3.154697e-02_r8,2.985578e-02_r8,&
        & 2.827395e-02_r8,2.679666e-02_r8,2.541909e-02_r8,2.413640e-02_r8,2.294378e-02_r8,&
        & 2.183639e-02_r8,2.080940e-02_r8,1.985801e-02_r8,1.897736e-02_r8,1.816265e-02_r8,&
        & 1.740905e-02_r8,1.671172e-02_r8,1.606585e-02_r8,1.546661e-02_r8,1.490917e-02_r8,&
        & 1.438870e-02_r8,1.390038e-02_r8,1.343939e-02_r8,1.300089e-02_r8,1.258006e-02_r8,&
        & 1.217208e-02_r8,1.177212e-02_r8,1.137536e-02_r8,1.097696e-02_r8,1.057210e-02_r8,&
        & 1.015596e-02_r8,9.723704e-03_r8,9.270516e-03_r8,8.791565e-03_r8,8.282026e-03_r8,&
        & 7.737072e-03_r8,7.151879e-03_r8,6.521619e-03_r8,5.841467e-03_r8,5.106597e-03_r8,&
        & 4.312183e-03_r8 /)
      fdlice3(:, 17) = (/ &
! band 17
        & 5.071224e-02_r8,5.000217e-02_r8,4.933872e-02_r8,4.871992e-02_r8,4.814380e-02_r8,&
        & 4.760839e-02_r8,4.711170e-02_r8,4.665177e-02_r8,4.622662e-02_r8,4.583426e-02_r8,&
        & 4.547274e-02_r8,4.514007e-02_r8,4.483428e-02_r8,4.455340e-02_r8,4.429544e-02_r8,&
        & 4.405844e-02_r8,4.384041e-02_r8,4.363939e-02_r8,4.345340e-02_r8,4.328047e-02_r8,&
        & 4.311861e-02_r8,4.296586e-02_r8,4.282024e-02_r8,4.267977e-02_r8,4.254248e-02_r8,&
        & 4.240640e-02_r8,4.226955e-02_r8,4.212995e-02_r8,4.198564e-02_r8,4.183462e-02_r8,&
        & 4.167494e-02_r8,4.150462e-02_r8,4.132167e-02_r8,4.112413e-02_r8,4.091003e-02_r8,&
        & 4.067737e-02_r8,4.042420e-02_r8,4.014854e-02_r8,3.984840e-02_r8,3.952183e-02_r8,&
        & 3.916683e-02_r8,3.878144e-02_r8,3.836368e-02_r8,3.791158e-02_r8,3.742316e-02_r8,&
        & 3.689645e-02_r8 /)
      fdlice3(:, 18) = (/ &
! band 18
        & 1.062938e-01_r8,1.065234e-01_r8,1.067822e-01_r8,1.070682e-01_r8,1.073793e-01_r8,&
        & 1.077137e-01_r8,1.080693e-01_r8,1.084442e-01_r8,1.088364e-01_r8,1.092439e-01_r8,&
        & 1.096647e-01_r8,1.100970e-01_r8,1.105387e-01_r8,1.109878e-01_r8,1.114423e-01_r8,&
        & 1.119004e-01_r8,1.123599e-01_r8,1.128190e-01_r8,1.132757e-01_r8,1.137279e-01_r8,&
        & 1.141738e-01_r8,1.146113e-01_r8,1.150385e-01_r8,1.154534e-01_r8,1.158540e-01_r8,&
        & 1.162383e-01_r8,1.166045e-01_r8,1.169504e-01_r8,1.172741e-01_r8,1.175738e-01_r8,&
        & 1.178472e-01_r8,1.180926e-01_r8,1.183080e-01_r8,1.184913e-01_r8,1.186405e-01_r8,&
        & 1.187538e-01_r8,1.188291e-01_r8,1.188645e-01_r8,1.188580e-01_r8,1.188076e-01_r8,&
        & 1.187113e-01_r8,1.185672e-01_r8,1.183733e-01_r8,1.181277e-01_r8,1.178282e-01_r8,&
        & 1.174731e-01_r8 /)
      fdlice3(:, 19) = (/ &
! band 19
        & 1.076195e-01_r8,1.065195e-01_r8,1.054696e-01_r8,1.044673e-01_r8,1.035099e-01_r8,&
        & 1.025951e-01_r8,1.017203e-01_r8,1.008831e-01_r8,1.000808e-01_r8,9.931116e-02_r8,&
        & 9.857151e-02_r8,9.785939e-02_r8,9.717230e-02_r8,9.650774e-02_r8,9.586322e-02_r8,&
        & 9.523623e-02_r8,9.462427e-02_r8,9.402484e-02_r8,9.343544e-02_r8,9.285358e-02_r8,&
        & 9.227675e-02_r8,9.170245e-02_r8,9.112818e-02_r8,9.055144e-02_r8,8.996974e-02_r8,&
        & 8.938056e-02_r8,8.878142e-02_r8,8.816981e-02_r8,8.754323e-02_r8,8.689919e-02_r8,&
        & 8.623517e-02_r8,8.554869e-02_r8,8.483724e-02_r8,8.409832e-02_r8,8.332943e-02_r8,&
        & 8.252807e-02_r8,8.169175e-02_r8,8.081795e-02_r8,7.990419e-02_r8,7.894796e-02_r8,&
        & 7.794676e-02_r8,7.689809e-02_r8,7.579945e-02_r8,7.464834e-02_r8,7.344227e-02_r8,&
        & 7.217872e-02_r8 /)
      fdlice3(:, 20) = (/ &
! band 20
        & 1.119014e-01_r8,1.122706e-01_r8,1.126690e-01_r8,1.130947e-01_r8,1.135456e-01_r8,&
        & 1.140199e-01_r8,1.145154e-01_r8,1.150302e-01_r8,1.155623e-01_r8,1.161096e-01_r8,&
        & 1.166703e-01_r8,1.172422e-01_r8,1.178233e-01_r8,1.184118e-01_r8,1.190055e-01_r8,&
        & 1.196025e-01_r8,1.202008e-01_r8,1.207983e-01_r8,1.213931e-01_r8,1.219832e-01_r8,&
        & 1.225665e-01_r8,1.231411e-01_r8,1.237050e-01_r8,1.242561e-01_r8,1.247926e-01_r8,&
        & 1.253122e-01_r8,1.258132e-01_r8,1.262934e-01_r8,1.267509e-01_r8,1.271836e-01_r8,&
        & 1.275896e-01_r8,1.279669e-01_r8,1.283134e-01_r8,1.286272e-01_r8,1.289063e-01_r8,&
        & 1.291486e-01_r8,1.293522e-01_r8,1.295150e-01_r8,1.296351e-01_r8,1.297104e-01_r8,&
        & 1.297390e-01_r8,1.297189e-01_r8,1.296480e-01_r8,1.295244e-01_r8,1.293460e-01_r8,&
        & 1.291109e-01_r8 /)
      fdlice3(:, 21) = (/ &
! band 21
        & 1.133298e-01_r8,1.136777e-01_r8,1.140556e-01_r8,1.144615e-01_r8,1.148934e-01_r8,&
        & 1.153492e-01_r8,1.158269e-01_r8,1.163243e-01_r8,1.168396e-01_r8,1.173706e-01_r8,&
        & 1.179152e-01_r8,1.184715e-01_r8,1.190374e-01_r8,1.196108e-01_r8,1.201897e-01_r8,&
        & 1.207720e-01_r8,1.213558e-01_r8,1.219389e-01_r8,1.225194e-01_r8,1.230951e-01_r8,&
        & 1.236640e-01_r8,1.242241e-01_r8,1.247733e-01_r8,1.253096e-01_r8,1.258309e-01_r8,&
        & 1.263352e-01_r8,1.268205e-01_r8,1.272847e-01_r8,1.277257e-01_r8,1.281415e-01_r8,&
        & 1.285300e-01_r8,1.288893e-01_r8,1.292173e-01_r8,1.295118e-01_r8,1.297710e-01_r8,&
        & 1.299927e-01_r8,1.301748e-01_r8,1.303154e-01_r8,1.304124e-01_r8,1.304637e-01_r8,&
        & 1.304673e-01_r8,1.304212e-01_r8,1.303233e-01_r8,1.301715e-01_r8,1.299638e-01_r8,&
        & 1.296983e-01_r8 /)
      fdlice3(:, 22) = (/ &
! band 22
        & 1.145360e-01_r8,1.153256e-01_r8,1.161453e-01_r8,1.169929e-01_r8,1.178666e-01_r8,&
        & 1.187641e-01_r8,1.196835e-01_r8,1.206227e-01_r8,1.215796e-01_r8,1.225522e-01_r8,&
        & 1.235383e-01_r8,1.245361e-01_r8,1.255433e-01_r8,1.265579e-01_r8,1.275779e-01_r8,&
        & 1.286011e-01_r8,1.296257e-01_r8,1.306494e-01_r8,1.316703e-01_r8,1.326862e-01_r8,&
        & 1.336951e-01_r8,1.346950e-01_r8,1.356838e-01_r8,1.366594e-01_r8,1.376198e-01_r8,&
        & 1.385629e-01_r8,1.394866e-01_r8,1.403889e-01_r8,1.412678e-01_r8,1.421212e-01_r8,&
        & 1.429469e-01_r8,1.437430e-01_r8,1.445074e-01_r8,1.452381e-01_r8,1.459329e-01_r8,&
        & 1.465899e-01_r8,1.472069e-01_r8,1.477819e-01_r8,1.483128e-01_r8,1.487976e-01_r8,&
        & 1.492343e-01_r8,1.496207e-01_r8,1.499548e-01_r8,1.502346e-01_r8,1.504579e-01_r8,&
        & 1.506227e-01_r8 /)
      fdlice3(:, 23) = (/ &
! band 23
        & 1.153263e-01_r8,1.161445e-01_r8,1.169932e-01_r8,1.178703e-01_r8,1.187738e-01_r8,&
        & 1.197016e-01_r8,1.206516e-01_r8,1.216217e-01_r8,1.226099e-01_r8,1.236141e-01_r8,&
        & 1.246322e-01_r8,1.256621e-01_r8,1.267017e-01_r8,1.277491e-01_r8,1.288020e-01_r8,&
        & 1.298584e-01_r8,1.309163e-01_r8,1.319736e-01_r8,1.330281e-01_r8,1.340778e-01_r8,&
        & 1.351207e-01_r8,1.361546e-01_r8,1.371775e-01_r8,1.381873e-01_r8,1.391820e-01_r8,&
        & 1.401593e-01_r8,1.411174e-01_r8,1.420540e-01_r8,1.429671e-01_r8,1.438547e-01_r8,&
        & 1.447146e-01_r8,1.455449e-01_r8,1.463433e-01_r8,1.471078e-01_r8,1.478364e-01_r8,&
        & 1.485270e-01_r8,1.491774e-01_r8,1.497857e-01_r8,1.503497e-01_r8,1.508674e-01_r8,&
        & 1.513367e-01_r8,1.517554e-01_r8,1.521216e-01_r8,1.524332e-01_r8,1.526880e-01_r8,&
        & 1.528840e-01_r8 /)
      fdlice3(:, 24) = (/ &
! band 24
        & 1.160842e-01_r8,1.169118e-01_r8,1.177697e-01_r8,1.186556e-01_r8,1.195676e-01_r8,&
        & 1.205036e-01_r8,1.214616e-01_r8,1.224394e-01_r8,1.234349e-01_r8,1.244463e-01_r8,&
        & 1.254712e-01_r8,1.265078e-01_r8,1.275539e-01_r8,1.286075e-01_r8,1.296664e-01_r8,&
        & 1.307287e-01_r8,1.317923e-01_r8,1.328550e-01_r8,1.339149e-01_r8,1.349699e-01_r8,&
        & 1.360179e-01_r8,1.370567e-01_r8,1.380845e-01_r8,1.390991e-01_r8,1.400984e-01_r8,&
        & 1.410803e-01_r8,1.420429e-01_r8,1.429840e-01_r8,1.439016e-01_r8,1.447936e-01_r8,&
        & 1.456579e-01_r8,1.464925e-01_r8,1.472953e-01_r8,1.480642e-01_r8,1.487972e-01_r8,&
        & 1.494923e-01_r8,1.501472e-01_r8,1.507601e-01_r8,1.513287e-01_r8,1.518511e-01_r8,&
        & 1.523252e-01_r8,1.527489e-01_r8,1.531201e-01_r8,1.534368e-01_r8,1.536969e-01_r8,&
        & 1.538984e-01_r8 /)
      fdlice3(:, 25) = (/ &
! band 25
        & 1.168725e-01_r8,1.177088e-01_r8,1.185747e-01_r8,1.194680e-01_r8,1.203867e-01_r8,&
        & 1.213288e-01_r8,1.222923e-01_r8,1.232750e-01_r8,1.242750e-01_r8,1.252903e-01_r8,&
        & 1.263187e-01_r8,1.273583e-01_r8,1.284069e-01_r8,1.294626e-01_r8,1.305233e-01_r8,&
        & 1.315870e-01_r8,1.326517e-01_r8,1.337152e-01_r8,1.347756e-01_r8,1.358308e-01_r8,&
        & 1.368788e-01_r8,1.379175e-01_r8,1.389449e-01_r8,1.399590e-01_r8,1.409577e-01_r8,&
        & 1.419389e-01_r8,1.429007e-01_r8,1.438410e-01_r8,1.447577e-01_r8,1.456488e-01_r8,&
        & 1.465123e-01_r8,1.473461e-01_r8,1.481483e-01_r8,1.489166e-01_r8,1.496492e-01_r8,&
        & 1.503439e-01_r8,1.509988e-01_r8,1.516118e-01_r8,1.521808e-01_r8,1.527038e-01_r8,&
        & 1.531788e-01_r8,1.536037e-01_r8,1.539764e-01_r8,1.542951e-01_r8,1.545575e-01_r8,&
        & 1.547617e-01_r8 /)
      fdlice3(:, 26) = (/ &
!band 26
        & 1.180509e-01_r8,1.189025e-01_r8,1.197820e-01_r8,1.206875e-01_r8,1.216171e-01_r8,&
        & 1.225687e-01_r8,1.235404e-01_r8,1.245303e-01_r8,1.255363e-01_r8,1.265564e-01_r8,&
        & 1.275888e-01_r8,1.286313e-01_r8,1.296821e-01_r8,1.307392e-01_r8,1.318006e-01_r8,&
        & 1.328643e-01_r8,1.339284e-01_r8,1.349908e-01_r8,1.360497e-01_r8,1.371029e-01_r8,&
        & 1.381486e-01_r8,1.391848e-01_r8,1.402095e-01_r8,1.412208e-01_r8,1.422165e-01_r8,&
        & 1.431949e-01_r8,1.441539e-01_r8,1.450915e-01_r8,1.460058e-01_r8,1.468947e-01_r8,&
        & 1.477564e-01_r8,1.485888e-01_r8,1.493900e-01_r8,1.501580e-01_r8,1.508907e-01_r8,&
        & 1.515864e-01_r8,1.522428e-01_r8,1.528582e-01_r8,1.534305e-01_r8,1.539578e-01_r8,&
        & 1.544380e-01_r8,1.548692e-01_r8,1.552494e-01_r8,1.555767e-01_r8,1.558490e-01_r8,&
        & 1.560645e-01_r8 /)
      fdlice3(:, 27) = (/ &
! band 27
        & 1.200480e-01_r8,1.209267e-01_r8,1.218304e-01_r8,1.227575e-01_r8,1.237059e-01_r8,&
        & 1.246739e-01_r8,1.256595e-01_r8,1.266610e-01_r8,1.276765e-01_r8,1.287041e-01_r8,&
        & 1.297420e-01_r8,1.307883e-01_r8,1.318412e-01_r8,1.328988e-01_r8,1.339593e-01_r8,&
        & 1.350207e-01_r8,1.360813e-01_r8,1.371393e-01_r8,1.381926e-01_r8,1.392396e-01_r8,&
        & 1.402783e-01_r8,1.413069e-01_r8,1.423235e-01_r8,1.433263e-01_r8,1.443134e-01_r8,&
        & 1.452830e-01_r8,1.462332e-01_r8,1.471622e-01_r8,1.480681e-01_r8,1.489490e-01_r8,&
        & 1.498032e-01_r8,1.506286e-01_r8,1.514236e-01_r8,1.521863e-01_r8,1.529147e-01_r8,&
        & 1.536070e-01_r8,1.542614e-01_r8,1.548761e-01_r8,1.554491e-01_r8,1.559787e-01_r8,&
        & 1.564629e-01_r8,1.568999e-01_r8,1.572879e-01_r8,1.576249e-01_r8,1.579093e-01_r8,&
        & 1.581390e-01_r8 /)
      fdlice3(:, 28) = (/ &
! band 28
        & 1.247813e-01_r8,1.256496e-01_r8,1.265417e-01_r8,1.274560e-01_r8,1.283905e-01_r8,&
        & 1.293436e-01_r8,1.303135e-01_r8,1.312983e-01_r8,1.322964e-01_r8,1.333060e-01_r8,&
        & 1.343252e-01_r8,1.353523e-01_r8,1.363855e-01_r8,1.374231e-01_r8,1.384632e-01_r8,&
        & 1.395042e-01_r8,1.405441e-01_r8,1.415813e-01_r8,1.426140e-01_r8,1.436404e-01_r8,&
        & 1.446587e-01_r8,1.456672e-01_r8,1.466640e-01_r8,1.476475e-01_r8,1.486157e-01_r8,&
        & 1.495671e-01_r8,1.504997e-01_r8,1.514117e-01_r8,1.523016e-01_r8,1.531673e-01_r8,&
        & 1.540073e-01_r8,1.548197e-01_r8,1.556026e-01_r8,1.563545e-01_r8,1.570734e-01_r8,&
        & 1.577576e-01_r8,1.584054e-01_r8,1.590149e-01_r8,1.595843e-01_r8,1.601120e-01_r8,&
        & 1.605962e-01_r8,1.610349e-01_r8,1.614266e-01_r8,1.617693e-01_r8,1.620614e-01_r8,&
        & 1.623011e-01_r8 /)
      fdlice3(:, 29) = (/ &
! band 29
        & 1.006055e-01_r8,9.549582e-02_r8,9.063960e-02_r8,8.602900e-02_r8,8.165612e-02_r8,&
        & 7.751308e-02_r8,7.359199e-02_r8,6.988496e-02_r8,6.638412e-02_r8,6.308156e-02_r8,&
        & 5.996942e-02_r8,5.703979e-02_r8,5.428481e-02_r8,5.169657e-02_r8,4.926719e-02_r8,&
        & 4.698880e-02_r8,4.485349e-02_r8,4.285339e-02_r8,4.098061e-02_r8,3.922727e-02_r8,&
        & 3.758547e-02_r8,3.604733e-02_r8,3.460497e-02_r8,3.325051e-02_r8,3.197604e-02_r8,&
        & 3.077369e-02_r8,2.963558e-02_r8,2.855381e-02_r8,2.752050e-02_r8,2.652776e-02_r8,&
        & 2.556772e-02_r8,2.463247e-02_r8,2.371415e-02_r8,2.280485e-02_r8,2.189670e-02_r8,&
        & 2.098180e-02_r8,2.005228e-02_r8,1.910024e-02_r8,1.811781e-02_r8,1.709709e-02_r8,&
        & 1.603020e-02_r8,1.490925e-02_r8,1.372635e-02_r8,1.247363e-02_r8,1.114319e-02_r8,&
        & 9.727157e-03_r8 /)

      end subroutine swcldpr

      end module rrtmg_sw_init


