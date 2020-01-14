!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_init.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/22 19:20:03 $
!
      module rrtmg_lw_init

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use rrlw_wvn
      use rrtmg_lw_setcoef, only: lwatmref, lwavplank

      implicit none

      contains

! **************************************************************************
      subroutine rrtmg_lw_ini
! **************************************************************************
!
!  Original version:       Michael J. Iacono; July, 1998
!  First revision for NCAR CCM:   September, 1998
!  Second revision for RRTM_V3.0:  September, 2002
!
!  This subroutine performs calculations necessary for the initialization
!  of the longwave model.  Lookup tables are computed for use in the LW
!  radiative transfer, and input absorption coefficient data for each
!  spectral band are reduced from 256 g-point intervals to 140.
! **************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_tbl, only: ntbl, tblint, pade, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrini, hnamini

! ------- Local -------

      integer :: itr, ibnd, igc, ig, ind, ipr 
      integer :: igcsm, iprsm

      real(kind=r8) :: wtsum, wtsm(mg)        !
      real(kind=r8) :: tfn                    !

! ------- Definitions -------
!     Arrays for 10000-point look-up tables:
!     TAU_TBL Clear-sky optical depth (used in cloudy radiative transfer)
!     EXP_TBL Exponential lookup table for ransmittance
!     TFN_TBL Tau transition function; i.e. the transition of the Planck
!             function from that for the mean layer temperature to that for
!             the layer boundary temperature as a function of optical depth.
!             The "linear in tau" method is used to make the table.
!     PADE    Pade approximation constant (= 0.278)
!     BPADE   Inverse of the Pade approximation constant
!

      hvrini = '$Revision: 1.2 $'

! Initialize model data
      call lwdatinit
      call lwcmbdat               ! g-point interval reduction data
      call lwcldpr                ! cloud optical properties
      call lwatmref               ! reference MLS profile
      call lwavplank              ! Planck function 
      call lw_kgb01               ! molecular absorption coefficients
      call lw_kgb02
      call lw_kgb03
      call lw_kgb04
      call lw_kgb05
      call lw_kgb06
      call lw_kgb07
      call lw_kgb08
      call lw_kgb09
      call lw_kgb10
      call lw_kgb11
      call lw_kgb12
      call lw_kgb13
      call lw_kgb14
      call lw_kgb15
      call lw_kgb16

! Compute lookup tables for transmittance, tau transition function,
! and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
! computed as a function of the tau transition function, transmittance 
! is calculated as a function of tau, and the tau transition function 
! is calculated using the linear in tau formulation at values of tau 
! above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
! are computed at intervals of 0.001.  The inverse of the constant used
! in the Pade approximation to the tau transition function is set to b.

      tau_tbl(0) = 0.0_r8
      tau_tbl(ntbl) = 1.e10_r8
      exp_tbl(0) = 1.0_r8
      exp_tbl(ntbl) = 0.0_r8
      tfn_tbl(0) = 0.0_r8
      tfn_tbl(ntbl) = 1.0_r8
      bpade = 1.0_r8 / pade
      do itr = 1, ntbl-1
         tfn = float(itr) / float(ntbl)
         tau_tbl(itr) = bpade * tfn / (1._r8 - tfn)
         exp_tbl(itr) = exp(-tau_tbl(itr))
         if (tau_tbl(itr) .lt. 0.06_r8) then
            tfn_tbl(itr) = tau_tbl(itr)/6._r8
         else
            tfn_tbl(itr) = 1._r8-2._r8*((1._r8/tau_tbl(itr))-(exp_tbl(itr)/(1.-exp_tbl(itr))))
         endif
      enddo

! Perform g-point reduction from 16 per band (256 total points) to
! a band dependant number (140 total points) for all absorption
! coefficient input data and Planck fraction input data.
! Compute relative weighting for new g-point combinations.

      igcsm = 0
      do ibnd = 1,nbndlw
         iprsm = 0
         if (ngc(ibnd).lt.mg) then
            do igc = 1,ngc(ibnd) 
               igcsm = igcsm + 1
               wtsum = 0._r8
               do ipr = 1, ngn(igcsm)
                  iprsm = iprsm + 1
                  wtsum = wtsum + wt(iprsm)
               enddo
               wtsm(igc) = wtsum
            enddo
            do ig = 1, ng(ibnd)
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = wt(ig)/wtsm(ngm(ind))
            enddo
         else
            do ig = 1, ng(ibnd)
               igcsm = igcsm + 1
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = 1.0_r8
            enddo
         endif
      enddo

! Reduce g-points for absorption coefficient data in each LW spectral band.

      call cmbgb1
      call cmbgb2
      call cmbgb3
      call cmbgb4
      call cmbgb5
      call cmbgb6
      call cmbgb7
      call cmbgb8
      call cmbgb9
      call cmbgb10
      call cmbgb11
      call cmbgb12
      call cmbgb13
      call cmbgb14
      call cmbgb15
      call cmbgb16

      end subroutine rrtmg_lw_ini

!***************************************************************************
      subroutine lwdatinit
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : maxxsec, maxinpx
      use rrlw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2 
      use shr_const_mod, only: shr_const_avogad
      use physconst,     only: cday, gravit, cpair
      use rrlw_vsn

      save 
 
! Longwave spectral band limits (wavenumbers)
      wavenum1(:) = (/ 10._r8, 350._r8, 500._r8, 630._r8, 700._r8, 820._r8, &
                      980._r8,1080._r8,1180._r8,1390._r8,1480._r8,1800._r8, &
                     2080._r8,2250._r8,2390._r8,2600._r8/)
      wavenum2(:) = (/350._r8, 500._r8, 630._r8, 700._r8, 820._r8, 980._r8, &
                     1080._r8,1180._r8,1390._r8,1480._r8,1800._r8,2080._r8, &
                     2250._r8,2390._r8,2600._r8,3250._r8/)
      delwave(:) =  (/340._r8, 150._r8, 130._r8,  70._r8, 120._r8, 160._r8, &
                      100._r8, 100._r8, 210._r8,  90._r8, 320._r8, 280._r8, &
                      170._r8, 130._r8, 220._r8, 650._r8/)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)
      nspb(:) = (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

! Use constants set in CAM for consistency
      grav = gravit
      avogad = shr_const_avogad * 1.e-3_r8

!     Heatfac is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
!      heatfac = 8.4391_r8

!     Modified values for consistency with CAM:
!        =  (9.80616)(86400)(1e-5)/(1.00464)
!      heatfac = 8.43339130434_r8

!     Calculate heatfac directly from CAM constants:
      heatfac = grav * cday * 1.e-5_r8 / (cpair * 1.e-3_r8)

!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
      nxmol = 4
      ixindx(1) = 1
      ixindx(2) = 2
      ixindx(3) = 3
      ixindx(4) = 4
      ixindx(5:maxinpx) = 0

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

      end subroutine lwdatinit

!***************************************************************************
      subroutine lwcmbdat
!***************************************************************************

      save
 
! ------- Definitions -------
!     Arrays for the g-point reduction from 256 to 140 for the 16 LW bands:
!     This mapping from 256 to 140 points has been carefully selected to 
!     minimize the effect on the resulting fluxes and cooling rates, and
!     caution should be used if the mapping is modified.  The full 256
!     g-point set can be restored with ngptlw=256, ngc=16*16, ngn=256*1., etc.
!     ngptlw  The total number of new g-points
!     ngc     The number of new g-points in each band
!     ngs     The cumulative sum of new g-points for each band
!     ngm     The index of each new g-point relative to the original
!             16 g-points for each band.  
!     ngn     The number of original g-points that are combined to make
!             each new g-point in each band.
!     ngb     The band index for each new g-point.
!     wt      RRTM weights for 16 g-points.

! ------- Data statements -------
      ngc(:) = (/10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
      ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
      ngm(:) = (/1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
                 1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
                 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
                 1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
                 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
                 1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
                 1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16
      ngn(:) = (/1,1,2,2,2,2,2,2,1,1, &                       ! band 1
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
                 1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
                 2,2,2,2,2,2,2,2, &                           ! band 6
                 2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
                 2,2,2,2,2,2,2,2, &                           ! band 8
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
                 2,2,2,2,4,4, &                               ! band 10
                 1,1,2,2,2,2,3,3, &                           ! band 11
                 1,1,1,1,2,2,4,4, &                           ! band 12
                 3,3,4,6, &                                   ! band 13
                 8,8, &                                       ! band 14
                 8,8, &                                       ! band 15
                 4,12/)                                       ! band 16
      ngb(:) = (/1,1,1,1,1,1,1,1,1,1, &                       ! band 1
                 2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
                 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
                 4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
                 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
                 6,6,6,6,6,6,6,6, &                           ! band 6
                 7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
                 8,8,8,8,8,8,8,8, &                           ! band 8
                 9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
                 10,10,10,10,10,10, &                         ! band 10
                 11,11,11,11,11,11,11,11, &                   ! band 11
                 12,12,12,12,12,12,12,12, &                   ! band 12
                 13,13,13,13, &                               ! band 13
                 14,14, &                                     ! band 14
                 15,15, &                                     ! band 15
                 16,16/)                                      ! band 16
      wt(:) = (/ 0.1527534276_r8, 0.1491729617_r8, 0.1420961469_r8, &
                 0.1316886544_r8, 0.1181945205_r8, 0.1019300893_r8, &
                 0.0832767040_r8, 0.0626720116_r8, 0.0424925000_r8, &
                 0.0046269894_r8, 0.0038279891_r8, 0.0030260086_r8, &
                 0.0022199750_r8, 0.0014140010_r8, 0.0005330000_r8, &
                 0.0000750000_r8/)

      end subroutine lwcmbdat

!***************************************************************************
      subroutine cmbgb1
!***************************************************************************
!
!  Original version:    MJIacono; July 1998
!  Revision for GCMs:   MJIacono; September 1998
!  Revision for RRTMG:  MJIacono, September 2002
!  Revision for F90 reformatting:  MJIacono, June 2006
!
!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.
!
!  band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                       (high key - h2o; high minor - n2)
!  note: previous versions of rrtm band 1: 
!        10-250 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng1
      use rrlw_kg01, only: fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
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

      do jt = 1,4
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

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(1)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
               sumk2 = sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
            enddo
            ka_mn2(jt,igc) = sumk1
            kb_mn2(jt,igc) = sumk2
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(1)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb1

!***************************************************************************
      subroutine cmbgb2
!***************************************************************************
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2: 
!           250 - 500 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng2
      use rrlw_kg02, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               kb(jt,jp,igc) = sumk
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

      iprsm = 0
      do igc = 1,ngc(2)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb2

!***************************************************************************
      subroutine cmbgb3
!***************************************************************************
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!
! old band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng3
      use rrlw_kg03, only: fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, &
                           selfref, forref

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
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(3)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               kb_mn2o(jn,jt,igc) = sumk
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

      do jt = 1,4
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
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb3

!***************************************************************************
      subroutine cmbgb4
!***************************************************************************
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!
! old band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng4
      use rrlw_kg04, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

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
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(4)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
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

      do jt = 1,4
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
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb4

!***************************************************************************
      subroutine cmbgb5
!***************************************************************************
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!
! old band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng5
      use rrlw_kg05, only: fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, &
                           selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(5)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+64)
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
               do igc = 1,ngc(5)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+64)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(5)
              sumk = 0.
               do ipr = 1, ngn(ngs(4)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mo3(jn,jt,iprsm)*rwgt(iprsm+64)
               enddo
               ka_mo3(jn,jt,igc) = sumk
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

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(5)
         sumk = 0.
         do ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + ccl4o(iprsm)*rwgt(iprsm+64)
         enddo
         ccl4(igc) = sumk
      enddo

      end subroutine cmbgb5

!***************************************************************************
      subroutine cmbgb6
!***************************************************************************
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!
! old band 6:  820-980 cm-1 (low - h2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng6
      use rrlw_kg06, only: fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, &
                           selfrefo, forrefo, &
                           fracrefa, ka, ka_mco2, cfc11adj, cfc12, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf, sumk1, sumk2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(6)
               sumk = 0.
               do ipr = 1, ngn(ngs(5)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+80)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kao_mco2(jt,iprsm)*rwgt(iprsm+80)
            enddo
            ka_mco2(jt,igc) = sumk
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

      iprsm = 0
      do igc = 1,ngc(6)
         sumf = 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefao(iprsm)
            sumk1= sumk1+ cfc11adjo(iprsm)*rwgt(iprsm+80)
            sumk2= sumk2+ cfc12o(iprsm)*rwgt(iprsm+80)
         enddo
         fracrefa(igc) = sumf
         cfc11adj(igc) = sumk1
         cfc12(igc) = sumk2
      enddo

      end subroutine cmbgb6

!***************************************************************************
      subroutine cmbgb7
!***************************************************************************
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!
! old band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng7
      use rrlw_kg07, only: fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, &
                           selfref, forref

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

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(7)
              sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+96)
               enddo
               ka_mco2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mco2(jt,iprsm)*rwgt(iprsm+96)
            enddo
            kb_mco2(jt,igc) = sumk
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

      do jt = 1,4
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
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(7)
         sumf = 0.
         do ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb7

!***************************************************************************
      subroutine cmbgb8
!***************************************************************************
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!
! old band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng8
      use rrlw_kg08, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                           kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                           cfc12o, cfc22adjo, &
                           fracrefa, fracrefb, ka, ka_mco2, ka_mn2o, &
                           ka_mo3, kb, kb_mco2, kb_mn2o, selfref, forref, &
                           cfc12, cfc22adj

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumk3, sumk4, sumk5, sumf1, sumf2


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
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(8)
               sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               kb(jt,jp,igc) = sumk
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

      do jt = 1,4
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

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(8)
            sumk1 = 0.
            sumk2 = 0.
            sumk3 = 0.
            sumk4 = 0.
            sumk5 = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk2 = sumk2 + kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk3 = sumk3 + kao_mo3(jt,iprsm)*rwgt(iprsm+112)
               sumk4 = sumk4 + kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
               sumk5 = sumk5 + kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
            enddo
            ka_mco2(jt,igc) = sumk1
            kb_mco2(jt,igc) = sumk2
            ka_mo3(jt,igc) = sumk3
            ka_mn2o(jt,igc) = sumk4
            kb_mn2o(jt,igc) = sumk5
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(8)
         sumf1= 0.
         sumf2= 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
            sumk1= sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
            sumk2= sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
         cfc12(igc) = sumk1
         cfc22adj(igc) = sumk2
      enddo

      end subroutine cmbgb8

!***************************************************************************
      subroutine cmbgb9
!***************************************************************************
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)!

! old band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng9
      use rrlw_kg09, only: fracrefao, fracrefbo, kao, kao_mn2o, &
                           kbo, kbo_mn2o, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mn2o, &
                           kb, kb_mn2o, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


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

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(9)
              sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+128)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mn2o(jt,iprsm)*rwgt(iprsm+128)
            enddo
            kb_mn2o(jt,igc) = sumk
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

      do jt = 1,4
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

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(9)
            sumf = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(9)
         sumf = 0.
         do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb9

!***************************************************************************
      subroutine cmbgb10
!***************************************************************************
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!
! old band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng10
      use rrlw_kg10, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


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

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(10)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb10

!***************************************************************************
      subroutine cmbgb11
!***************************************************************************
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!
! old band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng11
      use rrlw_kg11, only: fracrefao, fracrefbo, kao, kao_mo2, &
                           kbo, kbo_mo2, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mo2, &
                           kb, kb_mo2, selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(11)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mo2(jt,iprsm)*rwgt(iprsm+160)
               sumk2 = sumk2 + kbo_mo2(jt,iprsm)*rwgt(iprsm+160)
            enddo
            ka_mo2(jt,igc) = sumk1
            kb_mo2(jt,igc) = sumk2
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(11)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb11

!***************************************************************************
      subroutine cmbgb12
!***************************************************************************
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! old band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng12
      use rrlw_kg12, only: fracrefao, kao, selfrefo, forrefo, &
                           fracrefa, ka, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(12)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(11)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+176)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(12)
            sumf = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb12

!***************************************************************************
      subroutine cmbgb13
!***************************************************************************
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!
! old band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng13
      use rrlw_kg13, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                           kbo_mo3, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mco2, ka_mco, &
                           kb_mo3, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf


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

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(13)
              sumk1 = 0.
              sumk2 = 0.
               do ipr = 1, ngn(ngs(12)+igc)
                  iprsm = iprsm + 1
                  sumk1 = sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
                  sumk2 = sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
               enddo
               ka_mco2(jn,jt,igc) = sumk1
               ka_mco(jn,jt,igc) = sumk2
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
            enddo
            kb_mo3(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(13)
         sumf = 0.
         do ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(13)
            sumf = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb13

!***************************************************************************
      subroutine cmbgb14
!***************************************************************************
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! old band 14:  2250-2380 cm-1 (low - co2; high - co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng14
      use rrlw_kg14, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


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
      enddo

      do jt = 1,5
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
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb14

!***************************************************************************
      subroutine cmbgb15
!***************************************************************************
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!
! old band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng15
      use rrlw_kg15, only: fracrefao, kao, kao_mn2, selfrefo, forrefo, &
                           fracrefa, ka, ka_mn2, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(15)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(14)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+224)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(15)
              sumk = 0.
               do ipr = 1, ngn(ngs(14)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2(jn,jt,iprsm)*rwgt(iprsm+224)
               enddo
               ka_mn2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(15)
            sumf = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb15

!***************************************************************************
      subroutine cmbgb16
!***************************************************************************
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!
! old band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng16
      use rrlw_kg16, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(16)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(15)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+240)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(16)
               sumk = 0.
               do ipr = 1, ngn(ngs(15)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+240)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(16)
         sumf = 0.
         do ipr = 1, ngn(ngs(15)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(16)
            sumf = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb16

!***************************************************************************
      subroutine lwcldpr
!***************************************************************************

! --------- Modules ----------

      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3, &
                          absice4, extice4, ssaice4, asyice4 ! yihsuan@UMich 2018-08-28

      save

! ABSCLDn is the liquid water absorption coefficient (m2/g). 
! For INFLAG = 1.
      abscld1 = 0.0602410_r8
!  
! Everything below is for INFLAG = 2.

! ABSICEn(J,IB) are the parameters needed to compute the liquid water 
! absorption coefficient in spectral region IB for ICEFLAG=n.  The units
! of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units (microns (m2/g)).
! For ICEFLAG = 0.

      absice0(:)= (/0.005_r8,  1.0_r8/)

! For ICEFLAG = 1.
      absice1(1,:) = (/0.0036_r8, 0.0068_r8, 0.0003_r8, 0.0016_r8, 0.0020_r8/)
      absice1(2,:) = (/1.136_r8 , 0.600_r8 , 1.338_r8 , 1.166_r8 , 1.118_r8 /)

! For ICEFLAG = 2.  In each band, the absorption
! coefficients are listed for a range of effective radii from 5.0
! to 131.0 microns in increments of 3.0 microns.
! Spherical Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice2(:,1) = (/ &
! band 1
       7.798999e-02_r8,6.340479e-02_r8,5.417973e-02_r8,4.766245e-02_r8,4.272663e-02_r8, &
       3.880939e-02_r8,3.559544e-02_r8,3.289241e-02_r8,3.057511e-02_r8,2.855800e-02_r8, &
       2.678022e-02_r8,2.519712e-02_r8,2.377505e-02_r8,2.248806e-02_r8,2.131578e-02_r8, &
       2.024194e-02_r8,1.925337e-02_r8,1.833926e-02_r8,1.749067e-02_r8,1.670007e-02_r8, &
       1.596113e-02_r8,1.526845e-02_r8,1.461739e-02_r8,1.400394e-02_r8,1.342462e-02_r8, &
       1.287639e-02_r8,1.235656e-02_r8,1.186279e-02_r8,1.139297e-02_r8,1.094524e-02_r8, &
       1.051794e-02_r8,1.010956e-02_r8,9.718755e-03_r8,9.344316e-03_r8,8.985139e-03_r8, &
       8.640223e-03_r8,8.308656e-03_r8,7.989606e-03_r8,7.682312e-03_r8,7.386076e-03_r8, &
       7.100255e-03_r8,6.824258e-03_r8,6.557540e-03_r8/)
      absice2(:,2) = (/ &
! band 2
       2.784879e-02_r8,2.709863e-02_r8,2.619165e-02_r8,2.529230e-02_r8,2.443225e-02_r8, &
       2.361575e-02_r8,2.284021e-02_r8,2.210150e-02_r8,2.139548e-02_r8,2.071840e-02_r8, &
       2.006702e-02_r8,1.943856e-02_r8,1.883064e-02_r8,1.824120e-02_r8,1.766849e-02_r8, &
       1.711099e-02_r8,1.656737e-02_r8,1.603647e-02_r8,1.551727e-02_r8,1.500886e-02_r8, &
       1.451045e-02_r8,1.402132e-02_r8,1.354084e-02_r8,1.306842e-02_r8,1.260355e-02_r8, &
       1.214575e-02_r8,1.169460e-02_r8,1.124971e-02_r8,1.081072e-02_r8,1.037731e-02_r8, &
       9.949167e-03_r8,9.526021e-03_r8,9.107615e-03_r8,8.693714e-03_r8,8.284096e-03_r8, &
       7.878558e-03_r8,7.476910e-03_r8,7.078974e-03_r8,6.684586e-03_r8,6.293589e-03_r8, &
       5.905839e-03_r8,5.521200e-03_r8,5.139543e-03_r8/)
      absice2(:,3) = (/ &
! band 3
       1.065397e-01_r8,8.005726e-02_r8,6.546428e-02_r8,5.589131e-02_r8,4.898681e-02_r8, &
       4.369932e-02_r8,3.947901e-02_r8,3.600676e-02_r8,3.308299e-02_r8,3.057561e-02_r8, &
       2.839325e-02_r8,2.647040e-02_r8,2.475872e-02_r8,2.322164e-02_r8,2.183091e-02_r8, &
       2.056430e-02_r8,1.940407e-02_r8,1.833586e-02_r8,1.734787e-02_r8,1.643034e-02_r8, &
       1.557512e-02_r8,1.477530e-02_r8,1.402501e-02_r8,1.331924e-02_r8,1.265364e-02_r8, &
       1.202445e-02_r8,1.142838e-02_r8,1.086257e-02_r8,1.032445e-02_r8,9.811791e-03_r8, &
       9.322587e-03_r8,8.855053e-03_r8,8.407591e-03_r8,7.978763e-03_r8,7.567273e-03_r8, &
       7.171949e-03_r8,6.791728e-03_r8,6.425642e-03_r8,6.072809e-03_r8,5.732424e-03_r8, &
       5.403748e-03_r8,5.086103e-03_r8,4.778865e-03_r8/)
      absice2(:,4) = (/ &
! band 4
       1.804566e-01_r8,1.168987e-01_r8,8.680442e-02_r8,6.910060e-02_r8,5.738174e-02_r8, &
       4.902332e-02_r8,4.274585e-02_r8,3.784923e-02_r8,3.391734e-02_r8,3.068690e-02_r8, &
       2.798301e-02_r8,2.568480e-02_r8,2.370600e-02_r8,2.198337e-02_r8,2.046940e-02_r8, &
       1.912777e-02_r8,1.793016e-02_r8,1.685420e-02_r8,1.588193e-02_r8,1.499882e-02_r8, &
       1.419293e-02_r8,1.345440e-02_r8,1.277496e-02_r8,1.214769e-02_r8,1.156669e-02_r8, &
       1.102694e-02_r8,1.052412e-02_r8,1.005451e-02_r8,9.614854e-03_r8,9.202335e-03_r8, &
       8.814470e-03_r8,8.449077e-03_r8,8.104223e-03_r8,7.778195e-03_r8,7.469466e-03_r8, &
       7.176671e-03_r8,6.898588e-03_r8,6.634117e-03_r8,6.382264e-03_r8,6.142134e-03_r8, &
       5.912913e-03_r8,5.693862e-03_r8,5.484308e-03_r8/)
      absice2(:,5) = (/ &
! band 5
       2.131806e-01_r8,1.311372e-01_r8,9.407171e-02_r8,7.299442e-02_r8,5.941273e-02_r8, &
       4.994043e-02_r8,4.296242e-02_r8,3.761113e-02_r8,3.337910e-02_r8,2.994978e-02_r8, &
       2.711556e-02_r8,2.473461e-02_r8,2.270681e-02_r8,2.095943e-02_r8,1.943839e-02_r8, &
       1.810267e-02_r8,1.692057e-02_r8,1.586719e-02_r8,1.492275e-02_r8,1.407132e-02_r8, &
       1.329989e-02_r8,1.259780e-02_r8,1.195618e-02_r8,1.136761e-02_r8,1.082583e-02_r8, &
       1.032552e-02_r8,9.862158e-03_r8,9.431827e-03_r8,9.031157e-03_r8,8.657217e-03_r8, &
       8.307449e-03_r8,7.979609e-03_r8,7.671724e-03_r8,7.382048e-03_r8,7.109032e-03_r8, &
       6.851298e-03_r8,6.607615e-03_r8,6.376881e-03_r8,6.158105e-03_r8,5.950394e-03_r8, &
       5.752942e-03_r8,5.565019e-03_r8,5.385963e-03_r8/)
      absice2(:,6) = (/ &
! band 6
       1.546177e-01_r8,1.039251e-01_r8,7.910347e-02_r8,6.412429e-02_r8,5.399997e-02_r8, &
       4.664937e-02_r8,4.104237e-02_r8,3.660781e-02_r8,3.300218e-02_r8,3.000586e-02_r8, &
       2.747148e-02_r8,2.529633e-02_r8,2.340647e-02_r8,2.174723e-02_r8,2.027731e-02_r8, &
       1.896487e-02_r8,1.778492e-02_r8,1.671761e-02_r8,1.574692e-02_r8,1.485978e-02_r8, &
       1.404543e-02_r8,1.329489e-02_r8,1.260066e-02_r8,1.195636e-02_r8,1.135657e-02_r8, &
       1.079664e-02_r8,1.027257e-02_r8,9.780871e-03_r8,9.318505e-03_r8,8.882815e-03_r8, &
       8.471458e-03_r8,8.082364e-03_r8,7.713696e-03_r8,7.363817e-03_r8,7.031264e-03_r8, &
       6.714725e-03_r8,6.413021e-03_r8,6.125086e-03_r8,5.849958e-03_r8,5.586764e-03_r8, &
       5.334707e-03_r8,5.093066e-03_r8,4.861179e-03_r8/)
      absice2(:,7) = (/ &
! band 7
       7.583404e-02_r8,6.181558e-02_r8,5.312027e-02_r8,4.696039e-02_r8,4.225986e-02_r8, &
       3.849735e-02_r8,3.538340e-02_r8,3.274182e-02_r8,3.045798e-02_r8,2.845343e-02_r8, &
       2.667231e-02_r8,2.507353e-02_r8,2.362606e-02_r8,2.230595e-02_r8,2.109435e-02_r8, &
       1.997617e-02_r8,1.893916e-02_r8,1.797328e-02_r8,1.707016e-02_r8,1.622279e-02_r8, &
       1.542523e-02_r8,1.467241e-02_r8,1.395997e-02_r8,1.328414e-02_r8,1.264164e-02_r8, &
       1.202958e-02_r8,1.144544e-02_r8,1.088697e-02_r8,1.035218e-02_r8,9.839297e-03_r8, &
       9.346733e-03_r8,8.873057e-03_r8,8.416980e-03_r8,7.977335e-03_r8,7.553066e-03_r8, &
       7.143210e-03_r8,6.746888e-03_r8,6.363297e-03_r8,5.991700e-03_r8,5.631422e-03_r8, &
       5.281840e-03_r8,4.942378e-03_r8,4.612505e-03_r8/)
      absice2(:,8) = (/ &
! band 8
       9.022185e-02_r8,6.922700e-02_r8,5.710674e-02_r8,4.898377e-02_r8,4.305946e-02_r8, &
       3.849553e-02_r8,3.484183e-02_r8,3.183220e-02_r8,2.929794e-02_r8,2.712627e-02_r8, &
       2.523856e-02_r8,2.357810e-02_r8,2.210286e-02_r8,2.078089e-02_r8,1.958747e-02_r8, &
       1.850310e-02_r8,1.751218e-02_r8,1.660205e-02_r8,1.576232e-02_r8,1.498440e-02_r8, &
       1.426107e-02_r8,1.358624e-02_r8,1.295474e-02_r8,1.236212e-02_r8,1.180456e-02_r8, &
       1.127874e-02_r8,1.078175e-02_r8,1.031106e-02_r8,9.864433e-03_r8,9.439878e-03_r8, &
       9.035637e-03_r8,8.650140e-03_r8,8.281981e-03_r8,7.929895e-03_r8,7.592746e-03_r8, &
       7.269505e-03_r8,6.959238e-03_r8,6.661100e-03_r8,6.374317e-03_r8,6.098185e-03_r8, &
       5.832059e-03_r8,5.575347e-03_r8,5.327504e-03_r8/)
      absice2(:,9) = (/ &
! band 9
       1.294087e-01_r8,8.788217e-02_r8,6.728288e-02_r8,5.479720e-02_r8,4.635049e-02_r8, &
       4.022253e-02_r8,3.555576e-02_r8,3.187259e-02_r8,2.888498e-02_r8,2.640843e-02_r8, &
       2.431904e-02_r8,2.253038e-02_r8,2.098024e-02_r8,1.962267e-02_r8,1.842293e-02_r8, &
       1.735426e-02_r8,1.639571e-02_r8,1.553060e-02_r8,1.474552e-02_r8,1.402953e-02_r8, &
       1.337363e-02_r8,1.277033e-02_r8,1.221336e-02_r8,1.169741e-02_r8,1.121797e-02_r8, &
       1.077117e-02_r8,1.035369e-02_r8,9.962643e-03_r8,9.595509e-03_r8,9.250088e-03_r8, &
       8.924447e-03_r8,8.616876e-03_r8,8.325862e-03_r8,8.050057e-03_r8,7.788258e-03_r8, &
       7.539388e-03_r8,7.302478e-03_r8,7.076656e-03_r8,6.861134e-03_r8,6.655197e-03_r8, &
       6.458197e-03_r8,6.269543e-03_r8,6.088697e-03_r8/)
      absice2(:,10) = (/ &
! band 10
       1.593628e-01_r8,1.014552e-01_r8,7.458955e-02_r8,5.903571e-02_r8,4.887582e-02_r8, &
       4.171159e-02_r8,3.638480e-02_r8,3.226692e-02_r8,2.898717e-02_r8,2.631256e-02_r8, &
       2.408925e-02_r8,2.221156e-02_r8,2.060448e-02_r8,1.921325e-02_r8,1.799699e-02_r8, &
       1.692456e-02_r8,1.597177e-02_r8,1.511961e-02_r8,1.435289e-02_r8,1.365933e-02_r8, &
       1.302890e-02_r8,1.245334e-02_r8,1.192576e-02_r8,1.144037e-02_r8,1.099230e-02_r8, &
       1.057739e-02_r8,1.019208e-02_r8,9.833302e-03_r8,9.498395e-03_r8,9.185047e-03_r8, &
       8.891237e-03_r8,8.615185e-03_r8,8.355325e-03_r8,8.110267e-03_r8,7.878778e-03_r8, &
       7.659759e-03_r8,7.452224e-03_r8,7.255291e-03_r8,7.068166e-03_r8,6.890130e-03_r8, &
       6.720536e-03_r8,6.558794e-03_r8,6.404371e-03_r8/)
      absice2(:,11) = (/ &
! band 11
       1.656227e-01_r8,1.032129e-01_r8,7.487359e-02_r8,5.871431e-02_r8,4.828355e-02_r8, &
       4.099989e-02_r8,3.562924e-02_r8,3.150755e-02_r8,2.824593e-02_r8,2.560156e-02_r8, &
       2.341503e-02_r8,2.157740e-02_r8,2.001169e-02_r8,1.866199e-02_r8,1.748669e-02_r8, &
       1.645421e-02_r8,1.554015e-02_r8,1.472535e-02_r8,1.399457e-02_r8,1.333553e-02_r8, &
       1.273821e-02_r8,1.219440e-02_r8,1.169725e-02_r8,1.124104e-02_r8,1.082096e-02_r8, &
       1.043290e-02_r8,1.007336e-02_r8,9.739338e-03_r8,9.428223e-03_r8,9.137756e-03_r8, &
       8.865964e-03_r8,8.611115e-03_r8,8.371686e-03_r8,8.146330e-03_r8,7.933852e-03_r8, &
       7.733187e-03_r8,7.543386e-03_r8,7.363597e-03_r8,7.193056e-03_r8,7.031072e-03_r8, &
       6.877024e-03_r8,6.730348e-03_r8,6.590531e-03_r8/)
      absice2(:,12) = (/ &
! band 12
       9.194591e-02_r8,6.446867e-02_r8,4.962034e-02_r8,4.042061e-02_r8,3.418456e-02_r8, &
       2.968856e-02_r8,2.629900e-02_r8,2.365572e-02_r8,2.153915e-02_r8,1.980791e-02_r8, &
       1.836689e-02_r8,1.714979e-02_r8,1.610900e-02_r8,1.520946e-02_r8,1.442476e-02_r8, &
       1.373468e-02_r8,1.312345e-02_r8,1.257858e-02_r8,1.209010e-02_r8,1.164990e-02_r8, &
       1.125136e-02_r8,1.088901e-02_r8,1.055827e-02_r8,1.025531e-02_r8,9.976896e-03_r8, &
       9.720255e-03_r8,9.483022e-03_r8,9.263160e-03_r8,9.058902e-03_r8,8.868710e-03_r8, &
       8.691240e-03_r8,8.525312e-03_r8,8.369886e-03_r8,8.224042e-03_r8,8.086961e-03_r8, &
       7.957917e-03_r8,7.836258e-03_r8,7.721400e-03_r8,7.612821e-03_r8,7.510045e-03_r8, &
       7.412648e-03_r8,7.320242e-03_r8,7.232476e-03_r8/)
      absice2(:,13) = (/ &
! band 13
       1.437021e-01_r8,8.872535e-02_r8,6.392420e-02_r8,4.991833e-02_r8,4.096790e-02_r8, &
       3.477881e-02_r8,3.025782e-02_r8,2.681909e-02_r8,2.412102e-02_r8,2.195132e-02_r8, &
       2.017124e-02_r8,1.868641e-02_r8,1.743044e-02_r8,1.635529e-02_r8,1.542540e-02_r8, &
       1.461388e-02_r8,1.390003e-02_r8,1.326766e-02_r8,1.270395e-02_r8,1.219860e-02_r8, &
       1.174326e-02_r8,1.133107e-02_r8,1.095637e-02_r8,1.061442e-02_r8,1.030126e-02_r8, &
       1.001352e-02_r8,9.748340e-03_r8,9.503256e-03_r8,9.276155e-03_r8,9.065205e-03_r8, &
       8.868808e-03_r8,8.685571e-03_r8,8.514268e-03_r8,8.353820e-03_r8,8.203272e-03_r8, &
       8.061776e-03_r8,7.928578e-03_r8,7.803001e-03_r8,7.684443e-03_r8,7.572358e-03_r8, &
       7.466258e-03_r8,7.365701e-03_r8,7.270286e-03_r8/)
      absice2(:,14) = (/ &
! band 14
       1.288870e-01_r8,8.160295e-02_r8,5.964745e-02_r8,4.703790e-02_r8,3.888637e-02_r8, &
       3.320115e-02_r8,2.902017e-02_r8,2.582259e-02_r8,2.330224e-02_r8,2.126754e-02_r8, &
       1.959258e-02_r8,1.819130e-02_r8,1.700289e-02_r8,1.598320e-02_r8,1.509942e-02_r8, &
       1.432666e-02_r8,1.364572e-02_r8,1.304156e-02_r8,1.250220e-02_r8,1.201803e-02_r8, &
       1.158123e-02_r8,1.118537e-02_r8,1.082513e-02_r8,1.049605e-02_r8,1.019440e-02_r8, &
       9.916989e-03_r8,9.661116e-03_r8,9.424457e-03_r8,9.205005e-03_r8,9.001022e-03_r8, &
       8.810992e-03_r8,8.633588e-03_r8,8.467646e-03_r8,8.312137e-03_r8,8.166151e-03_r8, &
       8.028878e-03_r8,7.899597e-03_r8,7.777663e-03_r8,7.662498e-03_r8,7.553581e-03_r8, &
       7.450444e-03_r8,7.352662e-03_r8,7.259851e-03_r8/)
      absice2(:,15) = (/ &
! band 15
       8.254229e-02_r8,5.808787e-02_r8,4.492166e-02_r8,3.675028e-02_r8,3.119623e-02_r8, &
       2.718045e-02_r8,2.414450e-02_r8,2.177073e-02_r8,1.986526e-02_r8,1.830306e-02_r8, &
       1.699991e-02_r8,1.589698e-02_r8,1.495199e-02_r8,1.413374e-02_r8,1.341870e-02_r8, &
       1.278883e-02_r8,1.223002e-02_r8,1.173114e-02_r8,1.128322e-02_r8,1.087900e-02_r8, &
       1.051254e-02_r8,1.017890e-02_r8,9.873991e-03_r8,9.594347e-03_r8,9.337044e-03_r8, &
       9.099589e-03_r8,8.879842e-03_r8,8.675960e-03_r8,8.486341e-03_r8,8.309594e-03_r8, &
       8.144500e-03_r8,7.989986e-03_r8,7.845109e-03_r8,7.709031e-03_r8,7.581007e-03_r8, &
       7.460376e-03_r8,7.346544e-03_r8,7.238978e-03_r8,7.137201e-03_r8,7.040780e-03_r8, &
       6.949325e-03_r8,6.862483e-03_r8,6.779931e-03_r8/)
      absice2(:,16) = (/ &
! band 16
       1.382062e-01_r8,8.643227e-02_r8,6.282935e-02_r8,4.934783e-02_r8,4.063891e-02_r8, &
       3.455591e-02_r8,3.007059e-02_r8,2.662897e-02_r8,2.390631e-02_r8,2.169972e-02_r8, &
       1.987596e-02_r8,1.834393e-02_r8,1.703924e-02_r8,1.591513e-02_r8,1.493679e-02_r8, &
       1.407780e-02_r8,1.331775e-02_r8,1.264061e-02_r8,1.203364e-02_r8,1.148655e-02_r8, &
       1.099099e-02_r8,1.054006e-02_r8,1.012807e-02_r8,9.750215e-03_r8,9.402477e-03_r8, &
       9.081428e-03_r8,8.784143e-03_r8,8.508107e-03_r8,8.251146e-03_r8,8.011373e-03_r8, &
       7.787140e-03_r8,7.577002e-03_r8,7.379687e-03_r8,7.194071e-03_r8,7.019158e-03_r8, &
       6.854061e-03_r8,6.697986e-03_r8,6.550224e-03_r8,6.410138e-03_r8,6.277153e-03_r8, &
       6.150751e-03_r8,6.030462e-03_r8,5.915860e-03_r8/)

! ICEFLAG = 3; Fu parameterization. Particle size 5 - 140 micron in 
! increments of 3 microns.
! units = m2/g
! Hexagonal Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice3(:,1) = (/ &
! band 1
       3.110649e-03_r8,4.666352e-02_r8,6.606447e-02_r8,6.531678e-02_r8,6.012598e-02_r8, &
       5.437494e-02_r8,4.906411e-02_r8,4.441146e-02_r8,4.040585e-02_r8,3.697334e-02_r8, &
       3.403027e-02_r8,3.149979e-02_r8,2.931596e-02_r8,2.742365e-02_r8,2.577721e-02_r8, &
       2.433888e-02_r8,2.307732e-02_r8,2.196644e-02_r8,2.098437e-02_r8,2.011264e-02_r8, &
       1.933561e-02_r8,1.863992e-02_r8,1.801407e-02_r8,1.744812e-02_r8,1.693346e-02_r8, &
       1.646252e-02_r8,1.602866e-02_r8,1.562600e-02_r8,1.524933e-02_r8,1.489399e-02_r8, &
       1.455580e-02_r8,1.423098e-02_r8,1.391612e-02_r8,1.360812e-02_r8,1.330413e-02_r8, &
       1.300156e-02_r8,1.269801e-02_r8,1.239127e-02_r8,1.207928e-02_r8,1.176014e-02_r8, &
       1.143204e-02_r8,1.109334e-02_r8,1.074243e-02_r8,1.037786e-02_r8,9.998198e-03_r8, &
       9.602126e-03_r8/)
      absice3(:,2) = (/ &
! band 2
       3.984966e-04_r8,1.681097e-02_r8,2.627680e-02_r8,2.767465e-02_r8,2.700722e-02_r8, &
       2.579180e-02_r8,2.448677e-02_r8,2.323890e-02_r8,2.209096e-02_r8,2.104882e-02_r8, &
       2.010547e-02_r8,1.925003e-02_r8,1.847128e-02_r8,1.775883e-02_r8,1.710358e-02_r8, &
       1.649769e-02_r8,1.593449e-02_r8,1.540829e-02_r8,1.491429e-02_r8,1.444837e-02_r8, &
       1.400704e-02_r8,1.358729e-02_r8,1.318654e-02_r8,1.280258e-02_r8,1.243346e-02_r8, &
       1.207750e-02_r8,1.173325e-02_r8,1.139941e-02_r8,1.107487e-02_r8,1.075861e-02_r8, &
       1.044975e-02_r8,1.014753e-02_r8,9.851229e-03_r8,9.560240e-03_r8,9.274003e-03_r8, &
       8.992020e-03_r8,8.713845e-03_r8,8.439074e-03_r8,8.167346e-03_r8,7.898331e-03_r8, &
       7.631734e-03_r8,7.367286e-03_r8,7.104742e-03_r8,6.843882e-03_r8,6.584504e-03_r8, &
       6.326424e-03_r8/)
      absice3(:,3) = (/ &
! band 3
       6.933163e-02_r8,8.540475e-02_r8,7.701816e-02_r8,6.771158e-02_r8,5.986953e-02_r8, &
       5.348120e-02_r8,4.824962e-02_r8,4.390563e-02_r8,4.024411e-02_r8,3.711404e-02_r8, &
       3.440426e-02_r8,3.203200e-02_r8,2.993478e-02_r8,2.806474e-02_r8,2.638464e-02_r8, &
       2.486516e-02_r8,2.348288e-02_r8,2.221890e-02_r8,2.105780e-02_r8,1.998687e-02_r8, &
       1.899552e-02_r8,1.807490e-02_r8,1.721750e-02_r8,1.641693e-02_r8,1.566773e-02_r8, &
       1.496515e-02_r8,1.430509e-02_r8,1.368398e-02_r8,1.309865e-02_r8,1.254634e-02_r8, &
       1.202456e-02_r8,1.153114e-02_r8,1.106409e-02_r8,1.062166e-02_r8,1.020224e-02_r8, &
       9.804381e-03_r8,9.426771e-03_r8,9.068205e-03_r8,8.727578e-03_r8,8.403876e-03_r8, &
       8.096160e-03_r8,7.803564e-03_r8,7.525281e-03_r8,7.260560e-03_r8,7.008697e-03_r8, &
       6.769036e-03_r8/)
      absice3(:,4) = (/ &
! band 4
       1.765735e-01_r8,1.382700e-01_r8,1.095129e-01_r8,8.987475e-02_r8,7.591185e-02_r8, &
       6.554169e-02_r8,5.755500e-02_r8,5.122083e-02_r8,4.607610e-02_r8,4.181475e-02_r8, &
       3.822697e-02_r8,3.516432e-02_r8,3.251897e-02_r8,3.021073e-02_r8,2.817876e-02_r8, &
       2.637607e-02_r8,2.476582e-02_r8,2.331871e-02_r8,2.201113e-02_r8,2.082388e-02_r8, &
       1.974115e-02_r8,1.874983e-02_r8,1.783894e-02_r8,1.699922e-02_r8,1.622280e-02_r8, &
       1.550296e-02_r8,1.483390e-02_r8,1.421064e-02_r8,1.362880e-02_r8,1.308460e-02_r8, &
       1.257468e-02_r8,1.209611e-02_r8,1.164628e-02_r8,1.122287e-02_r8,1.082381e-02_r8, &
       1.044725e-02_r8,1.009154e-02_r8,9.755166e-03_r8,9.436783e-03_r8,9.135163e-03_r8, &
       8.849193e-03_r8,8.577856e-03_r8,8.320225e-03_r8,8.075451e-03_r8,7.842755e-03_r8, &
       7.621418e-03_r8/)
      absice3(:,5) = (/ &
! band 5
       2.339673e-01_r8,1.692124e-01_r8,1.291656e-01_r8,1.033837e-01_r8,8.562949e-02_r8, &
       7.273526e-02_r8,6.298262e-02_r8,5.537015e-02_r8,4.927787e-02_r8,4.430246e-02_r8, &
       4.017061e-02_r8,3.669072e-02_r8,3.372455e-02_r8,3.116995e-02_r8,2.894977e-02_r8, &
       2.700471e-02_r8,2.528842e-02_r8,2.376420e-02_r8,2.240256e-02_r8,2.117959e-02_r8, &
       2.007567e-02_r8,1.907456e-02_r8,1.816271e-02_r8,1.732874e-02_r8,1.656300e-02_r8, &
       1.585725e-02_r8,1.520445e-02_r8,1.459852e-02_r8,1.403419e-02_r8,1.350689e-02_r8, &
       1.301260e-02_r8,1.254781e-02_r8,1.210941e-02_r8,1.169468e-02_r8,1.130118e-02_r8, &
       1.092675e-02_r8,1.056945e-02_r8,1.022757e-02_r8,9.899560e-03_r8,9.584021e-03_r8, &
       9.279705e-03_r8,8.985479e-03_r8,8.700322e-03_r8,8.423306e-03_r8,8.153590e-03_r8, &
       7.890412e-03_r8/)
      absice3(:,6) = (/ &
! band 6
       1.145369e-01_r8,1.174566e-01_r8,9.917866e-02_r8,8.332990e-02_r8,7.104263e-02_r8, &
       6.153370e-02_r8,5.405472e-02_r8,4.806281e-02_r8,4.317918e-02_r8,3.913795e-02_r8, &
       3.574916e-02_r8,3.287437e-02_r8,3.041067e-02_r8,2.828017e-02_r8,2.642292e-02_r8, &
       2.479206e-02_r8,2.335051e-02_r8,2.206851e-02_r8,2.092195e-02_r8,1.989108e-02_r8, &
       1.895958e-02_r8,1.811385e-02_r8,1.734245e-02_r8,1.663573e-02_r8,1.598545e-02_r8, &
       1.538456e-02_r8,1.482700e-02_r8,1.430750e-02_r8,1.382150e-02_r8,1.336499e-02_r8, &
       1.293447e-02_r8,1.252685e-02_r8,1.213939e-02_r8,1.176968e-02_r8,1.141555e-02_r8, &
       1.107508e-02_r8,1.074655e-02_r8,1.042839e-02_r8,1.011923e-02_r8,9.817799e-03_r8, &
       9.522962e-03_r8,9.233688e-03_r8,8.949041e-03_r8,8.668171e-03_r8,8.390301e-03_r8, &
       8.114723e-03_r8/)
      absice3(:,7) = (/ &
! band 7
       1.222345e-02_r8,5.344230e-02_r8,5.523465e-02_r8,5.128759e-02_r8,4.676925e-02_r8, &
       4.266150e-02_r8,3.910561e-02_r8,3.605479e-02_r8,3.342843e-02_r8,3.115052e-02_r8, &
       2.915776e-02_r8,2.739935e-02_r8,2.583499e-02_r8,2.443266e-02_r8,2.316681e-02_r8, &
       2.201687e-02_r8,2.096619e-02_r8,2.000112e-02_r8,1.911044e-02_r8,1.828481e-02_r8, &
       1.751641e-02_r8,1.679866e-02_r8,1.612598e-02_r8,1.549360e-02_r8,1.489742e-02_r8, &
       1.433392e-02_r8,1.380002e-02_r8,1.329305e-02_r8,1.281068e-02_r8,1.235084e-02_r8, &
       1.191172e-02_r8,1.149171e-02_r8,1.108936e-02_r8,1.070341e-02_r8,1.033271e-02_r8, &
       9.976220e-03_r8,9.633021e-03_r8,9.302273e-03_r8,8.983216e-03_r8,8.675161e-03_r8, &
       8.377478e-03_r8,8.089595e-03_r8,7.810986e-03_r8,7.541170e-03_r8,7.279706e-03_r8, &
       7.026186e-03_r8/)
      absice3(:,8) = (/ &
! band 8
       6.711058e-02_r8,6.918198e-02_r8,6.127484e-02_r8,5.411944e-02_r8,4.836902e-02_r8, &
       4.375293e-02_r8,3.998077e-02_r8,3.683587e-02_r8,3.416508e-02_r8,3.186003e-02_r8, &
       2.984290e-02_r8,2.805671e-02_r8,2.645895e-02_r8,2.501733e-02_r8,2.370689e-02_r8, &
       2.250808e-02_r8,2.140532e-02_r8,2.038609e-02_r8,1.944018e-02_r8,1.855918e-02_r8, &
       1.773609e-02_r8,1.696504e-02_r8,1.624106e-02_r8,1.555990e-02_r8,1.491793e-02_r8, &
       1.431197e-02_r8,1.373928e-02_r8,1.319743e-02_r8,1.268430e-02_r8,1.219799e-02_r8, &
       1.173682e-02_r8,1.129925e-02_r8,1.088393e-02_r8,1.048961e-02_r8,1.011516e-02_r8, &
       9.759543e-03_r8,9.421813e-03_r8,9.101089e-03_r8,8.796559e-03_r8,8.507464e-03_r8, &
       8.233098e-03_r8,7.972798e-03_r8,7.725942e-03_r8,7.491940e-03_r8,7.270238e-03_r8, &
       7.060305e-03_r8/)
      absice3(:,9) = (/ &
! band 9
       1.236780e-01_r8,9.222386e-02_r8,7.383997e-02_r8,6.204072e-02_r8,5.381029e-02_r8, &
       4.770678e-02_r8,4.296928e-02_r8,3.916131e-02_r8,3.601540e-02_r8,3.335878e-02_r8, &
       3.107493e-02_r8,2.908247e-02_r8,2.732282e-02_r8,2.575276e-02_r8,2.433968e-02_r8, &
       2.305852e-02_r8,2.188966e-02_r8,2.081757e-02_r8,1.982974e-02_r8,1.891599e-02_r8, &
       1.806794e-02_r8,1.727865e-02_r8,1.654227e-02_r8,1.585387e-02_r8,1.520924e-02_r8, &
       1.460476e-02_r8,1.403730e-02_r8,1.350416e-02_r8,1.300293e-02_r8,1.253153e-02_r8, &
       1.208808e-02_r8,1.167094e-02_r8,1.127862e-02_r8,1.090979e-02_r8,1.056323e-02_r8, &
       1.023786e-02_r8,9.932665e-03_r8,9.646744e-03_r8,9.379250e-03_r8,9.129409e-03_r8, &
       8.896500e-03_r8,8.679856e-03_r8,8.478852e-03_r8,8.292904e-03_r8,8.121463e-03_r8, &
       7.964013e-03_r8/)
      absice3(:,10) = (/ &
! band 10
       1.655966e-01_r8,1.134205e-01_r8,8.714344e-02_r8,7.129241e-02_r8,6.063739e-02_r8, &
       5.294203e-02_r8,4.709309e-02_r8,4.247476e-02_r8,3.871892e-02_r8,3.559206e-02_r8, &
       3.293893e-02_r8,3.065226e-02_r8,2.865558e-02_r8,2.689288e-02_r8,2.532221e-02_r8, &
       2.391150e-02_r8,2.263582e-02_r8,2.147549e-02_r8,2.041476e-02_r8,1.944089e-02_r8, &
       1.854342e-02_r8,1.771371e-02_r8,1.694456e-02_r8,1.622989e-02_r8,1.556456e-02_r8, &
       1.494415e-02_r8,1.436491e-02_r8,1.382354e-02_r8,1.331719e-02_r8,1.284339e-02_r8, &
       1.239992e-02_r8,1.198486e-02_r8,1.159647e-02_r8,1.123323e-02_r8,1.089375e-02_r8, &
       1.057679e-02_r8,1.028124e-02_r8,1.000607e-02_r8,9.750376e-03_r8,9.513303e-03_r8, &
       9.294082e-03_r8,9.092003e-03_r8,8.906412e-03_r8,8.736702e-03_r8,8.582314e-03_r8, &
       8.442725e-03_r8/)
      absice3(:,11) = (/ &
! band 11
       1.775615e-01_r8,1.180046e-01_r8,8.929607e-02_r8,7.233500e-02_r8,6.108333e-02_r8, &
       5.303642e-02_r8,4.696927e-02_r8,4.221206e-02_r8,3.836768e-02_r8,3.518576e-02_r8, &
       3.250063e-02_r8,3.019825e-02_r8,2.819758e-02_r8,2.643943e-02_r8,2.487953e-02_r8, &
       2.348414e-02_r8,2.222705e-02_r8,2.108762e-02_r8,2.004936e-02_r8,1.909892e-02_r8, &
       1.822539e-02_r8,1.741975e-02_r8,1.667449e-02_r8,1.598330e-02_r8,1.534084e-02_r8, &
       1.474253e-02_r8,1.418446e-02_r8,1.366325e-02_r8,1.317597e-02_r8,1.272004e-02_r8, &
       1.229321e-02_r8,1.189350e-02_r8,1.151915e-02_r8,1.116859e-02_r8,1.084042e-02_r8, &
       1.053338e-02_r8,1.024636e-02_r8,9.978326e-03_r8,9.728357e-03_r8,9.495613e-03_r8, &
       9.279327e-03_r8,9.078798e-03_r8,8.893383e-03_r8,8.722488e-03_r8,8.565568e-03_r8, &
       8.422115e-03_r8/)
      absice3(:,12) = (/ &
! band 12
       9.465447e-02_r8,6.432047e-02_r8,5.060973e-02_r8,4.267283e-02_r8,3.741843e-02_r8, &
       3.363096e-02_r8,3.073531e-02_r8,2.842405e-02_r8,2.651789e-02_r8,2.490518e-02_r8, &
       2.351273e-02_r8,2.229056e-02_r8,2.120335e-02_r8,2.022541e-02_r8,1.933763e-02_r8, &
       1.852546e-02_r8,1.777763e-02_r8,1.708528e-02_r8,1.644134e-02_r8,1.584009e-02_r8, &
       1.527684e-02_r8,1.474774e-02_r8,1.424955e-02_r8,1.377957e-02_r8,1.333549e-02_r8, &
       1.291534e-02_r8,1.251743e-02_r8,1.214029e-02_r8,1.178265e-02_r8,1.144337e-02_r8, &
       1.112148e-02_r8,1.081609e-02_r8,1.052642e-02_r8,1.025178e-02_r8,9.991540e-03_r8, &
       9.745130e-03_r8,9.512038e-03_r8,9.291797e-03_r8,9.083980e-03_r8,8.888195e-03_r8, &
       8.704081e-03_r8,8.531306e-03_r8,8.369560e-03_r8,8.218558e-03_r8,8.078032e-03_r8, &
       7.947730e-03_r8/)
      absice3(:,13) = (/ &
! band 13
       1.560311e-01_r8,9.961097e-02_r8,7.502949e-02_r8,6.115022e-02_r8,5.214952e-02_r8, &
       4.578149e-02_r8,4.099731e-02_r8,3.724174e-02_r8,3.419343e-02_r8,3.165356e-02_r8, &
       2.949251e-02_r8,2.762222e-02_r8,2.598073e-02_r8,2.452322e-02_r8,2.321642e-02_r8, &
       2.203516e-02_r8,2.096002e-02_r8,1.997579e-02_r8,1.907036e-02_r8,1.823401e-02_r8, &
       1.745879e-02_r8,1.673819e-02_r8,1.606678e-02_r8,1.544003e-02_r8,1.485411e-02_r8, &
       1.430574e-02_r8,1.379215e-02_r8,1.331092e-02_r8,1.285996e-02_r8,1.243746e-02_r8, &
       1.204183e-02_r8,1.167164e-02_r8,1.132567e-02_r8,1.100281e-02_r8,1.070207e-02_r8, &
       1.042258e-02_r8,1.016352e-02_r8,9.924197e-03_r8,9.703953e-03_r8,9.502199e-03_r8, &
       9.318400e-03_r8,9.152066e-03_r8,9.002749e-03_r8,8.870038e-03_r8,8.753555e-03_r8, &
       8.652951e-03_r8/)
      absice3(:,14) = (/ &
! band 14
       1.559547e-01_r8,9.896700e-02_r8,7.441231e-02_r8,6.061469e-02_r8,5.168730e-02_r8, &
       4.537821e-02_r8,4.064106e-02_r8,3.692367e-02_r8,3.390714e-02_r8,3.139438e-02_r8, &
       2.925702e-02_r8,2.740783e-02_r8,2.578547e-02_r8,2.434552e-02_r8,2.305506e-02_r8, &
       2.188910e-02_r8,2.082842e-02_r8,1.985789e-02_r8,1.896553e-02_r8,1.814165e-02_r8, &
       1.737839e-02_r8,1.666927e-02_r8,1.600891e-02_r8,1.539279e-02_r8,1.481712e-02_r8, &
       1.427865e-02_r8,1.377463e-02_r8,1.330266e-02_r8,1.286068e-02_r8,1.244689e-02_r8, &
       1.205973e-02_r8,1.169780e-02_r8,1.135989e-02_r8,1.104492e-02_r8,1.075192e-02_r8, &
       1.048004e-02_r8,1.022850e-02_r8,9.996611e-03_r8,9.783753e-03_r8,9.589361e-03_r8, &
       9.412924e-03_r8,9.253977e-03_r8,9.112098e-03_r8,8.986903e-03_r8,8.878039e-03_r8, &
       8.785184e-03_r8/)
      absice3(:,15) = (/ &
! band 15
       1.102926e-01_r8,7.176622e-02_r8,5.530316e-02_r8,4.606056e-02_r8,4.006116e-02_r8, &
       3.579628e-02_r8,3.256909e-02_r8,3.001360e-02_r8,2.791920e-02_r8,2.615617e-02_r8, &
       2.464023e-02_r8,2.331426e-02_r8,2.213817e-02_r8,2.108301e-02_r8,2.012733e-02_r8, &
       1.925493e-02_r8,1.845331e-02_r8,1.771269e-02_r8,1.702531e-02_r8,1.638493e-02_r8, &
       1.578648e-02_r8,1.522579e-02_r8,1.469940e-02_r8,1.420442e-02_r8,1.373841e-02_r8, &
       1.329931e-02_r8,1.288535e-02_r8,1.249502e-02_r8,1.212700e-02_r8,1.178015e-02_r8, &
       1.145348e-02_r8,1.114612e-02_r8,1.085730e-02_r8,1.058633e-02_r8,1.033263e-02_r8, &
       1.009564e-02_r8,9.874895e-03_r8,9.669960e-03_r8,9.480449e-03_r8,9.306014e-03_r8, &
       9.146339e-03_r8,9.001138e-03_r8,8.870154e-03_r8,8.753148e-03_r8,8.649907e-03_r8, &
       8.560232e-03_r8/)
      absice3(:,16) = (/ &
! band 16
       1.688344e-01_r8,1.077072e-01_r8,7.994467e-02_r8,6.403862e-02_r8,5.369850e-02_r8, &
       4.641582e-02_r8,4.099331e-02_r8,3.678724e-02_r8,3.342069e-02_r8,3.065831e-02_r8, &
       2.834557e-02_r8,2.637680e-02_r8,2.467733e-02_r8,2.319286e-02_r8,2.188299e-02_r8, &
       2.071701e-02_r8,1.967121e-02_r8,1.872692e-02_r8,1.786931e-02_r8,1.708641e-02_r8, &
       1.636846e-02_r8,1.570743e-02_r8,1.509665e-02_r8,1.453052e-02_r8,1.400433e-02_r8, &
       1.351407e-02_r8,1.305631e-02_r8,1.262810e-02_r8,1.222688e-02_r8,1.185044e-02_r8, &
       1.149683e-02_r8,1.116436e-02_r8,1.085153e-02_r8,1.055701e-02_r8,1.027961e-02_r8, &
       1.001831e-02_r8,9.772141e-03_r8,9.540280e-03_r8,9.321966e-03_r8,9.116517e-03_r8, &
       8.923315e-03_r8,8.741803e-03_r8,8.571472e-03_r8,8.411860e-03_r8,8.262543e-03_r8, &
       8.123136e-03_r8/)

!>>> yihsuan@UMich 2018-08-28 add start >>>

! CPKuo@TAMU
! ICEFLAG = 4; MODIS collection 6 parameterization. 
! Effective diameter 3 - 370 micron.
! Roughened Aggregated Hexagonal Ice Particle shape.
! Data listed below are regression coefficients for polynomial fittings.
! For example, at ib spectral band, absice4(1,1:4,ib) are regression 
! coefficients of 3rd order of polynomial from 3rd order to 0th order 
! for effective diameter larger than 25 micron; absice4(2,1:8,ib) are 
! regression coefficients of 7th order of polynomial from 7th order to 
! 0th order for effective diameter smaller than 25 micron.

! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
! band  1
       absice4(1,:, 1) = (/  1.7768618851557e+02_r8, -5.0216289089028e+00_r8,  6.7638021981379e-01_r8,  4.4342605109832e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 1) = (/  3.4281290994587e+03_r8, -3.6178735111899e+03_r8,  1.6324931636039e+03_r8, -4.2559837121445e+02_r8, &
                             7.3296294563148e+01_r8, -8.7424474183646e+00_r8,  6.7638021981379e-01_r8,  4.4342605109832e-02_r8/)
! band  2
       absice4(1,:, 2) = (/  2.5956566494407e+02_r8,  1.4982891537084e+00_r8,  1.9473840863157e-01_r8,  2.2602355458432e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 2) = (/  1.7615520672443e+03_r8, -2.2501270067212e+03_r8,  1.1858399489772e+03_r8, -3.3320262234622e+02_r8, &
                             5.3643243729664e+01_r8, -4.8512284380817e+00_r8,  1.9473840863157e-01_r8,  2.2602355458432e-02_r8/)
! band  3
       absice4(1,:, 3) = (/  8.4208902545781e+01_r8, -5.6566668059479e+00_r8,  8.0881691866581e-01_r8,  4.5857932828471e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 3) = (/ -8.3927931421904e+03_r8,  8.6997560212212e+03_r8, -3.4494416605605e+03_r8,  6.1061436275872e+02_r8, &
                            -2.9188904779367e+01_r8, -5.5769656293431e+00_r8,  8.0881691866581e-01_r8,  4.5857932828471e-02_r8/)
! band  4
       absice4(1,:, 4) = (/ -1.8260195146217e+02_r8, -1.2363701672791e+01_r8,  1.2514656521124e+00_r8,  5.8638648700055e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 4) = (/ -8.8886840315571e+03_r8,  1.1484918815122e+04_r8, -5.7598106387407e+03_r8,  1.3956678090072e+03_r8, &
                            -1.5371557463940e+02_r8,  1.4251659366685e+00_r8,  1.2514656521124e+00_r8,  5.8638648700055e-02_r8/)
! band  5
       absice4(1,:, 5) = (/ -2.4717687565809e+02_r8, -1.0457756233011e+01_r8,  1.5467967235533e+00_r8,  6.3966828570920e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 5) = (/  3.9298527004489e+03_r8, -4.6797542148577e+02_r8, -1.5552846579777e+03_r8,  7.5451263737238e+02_r8, &
                            -1.3140133851088e+02_r8,  5.1141047474870e+00_r8,  1.5467967235533e+00_r8,  6.3966828570920e-02_r8/)
! band  6
       absice4(1,:, 6) = (/ -9.9462543587775e+01_r8, -6.9983524396956e+00_r8,  1.2816560841087e+00_r8,  5.6439999763648e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 6) = (/  3.2927144113559e+03_r8, -1.1671694213311e+03_r8, -5.5738967772723e+02_r8,  3.5373374160984e+02_r8, &
                            -6.1115681300347e+01_r8,  8.6236469369408e-01_r8,  1.2816560841087e+00_r8,  5.6439999763648e-02_r8/)
! band  7
       absice4(1,:, 7) = (/  1.7764385945777e+02_r8, -3.4786098081541e+00_r8,  5.5206735868209e-01_r8,  3.8286689343153e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 7) = (/  8.4619929490515e+02_r8, -5.9142238743104e+02_r8,  1.8334619099061e+02_r8, -6.2265768019186e+01_r8, &
                             2.2683541579886e+01_r8, -4.9653450636526e+00_r8,  5.5206735868209e-01_r8,  3.8286689343153e-02_r8/)
! band  8
       absice4(1,:, 8) = (/  1.7822197325529e+02_r8, -3.4845635563255e+00_r8,  5.3471915439578e-01_r8,  3.7701223902836e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 8) = (/  4.5781427812777e+02_r8, -1.8568527813523e+02_r8,  1.2291126531986e+01_r8, -2.5387268227944e+01_r8, &
                             1.8536932055404e+01_r8, -4.7496202611862e+00_r8,  5.3471915439578e-01_r8,  3.7701223902836e-02_r8/)
! band  9
       absice4(1,:, 9) = (/  1.2741914019179e+02_r8, -4.9519138638016e+00_r8,  6.6316044435397e-01_r8,  4.1859938578194e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:, 9) = (/ -2.2930824958501e+02_r8,  7.6738607072577e+02_r8, -4.9061257862572e+02_r8,  1.0171856463125e+02_r8, &
                             3.8581236508460e+00_r8, -4.5702574604381e+00_r8,  6.6316044435397e-01_r8,  4.1859938578194e-02_r8/)
! band 10
       absice4(1,:,10) = (/  3.3000397113931e+01_r8, -7.1398181236536e+00_r8,  8.9226328926866e-01_r8,  4.8634486647221e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,10) = (/ -1.0588319780351e+03_r8,  2.1413034318833e+03_r8, -1.3044576920366e+03_r8,  3.2921161154962e+02_r8, &
                            -2.6316689409232e+01_r8, -3.5983838123489e+00_r8,  8.9226328926866e-01_r8,  4.8634486647221e-02_r8/)
! band 11
       absice4(1,:,11) = (/  3.9243834245466e+00_r8, -7.5670694327206e+00_r8,  9.6366406281792e-01_r8,  5.0480007951557e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,11) = (/ -1.0359618956746e+03_r8,  2.3223585667469e+03_r8, -1.4826909869869e+03_r8,  3.9386970917859e+02_r8, &
                            -3.7172062438955e+01_r8, -2.9593669123465e+00_r8,  9.6366406281792e-01_r8,  5.0480007951557e-02_r8/)
! band 12
       absice4(1,:,12) = (/  2.0202201811418e+02_r8, -2.1384050572965e+00_r8,  4.2162329274198e-01_r8,  3.3075401292158e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,12) = (/  9.1782189430575e+02_r8, -6.7532410430445e+02_r8,  2.0562196328329e+02_r8, -5.5277314306493e+01_r8, &
                             1.7654501178747e+01_r8, -3.8498727935874e+00_r8,  4.2162329274198e-01_r8,  3.3075401292158e-02_r8/)
! band 13
       absice4(1,:,13) = (/  1.1655227371902e+02_r8, -4.9489096936260e+00_r8,  6.9472738357294e-01_r8,  4.2505428972267e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,13) = (/  2.6844840668137e+02_r8,  3.5429450138627e+02_r8, -3.7318325539933e+02_r8,  9.2886690671731e+01_r8, &
                             2.0108876027623e+00_r8, -4.2695481754025e+00_r8,  6.9472738357294e-01_r8,  4.2505428972267e-02_r8/)
! band 14
       absice4(1,:,14) = (/  1.0126681998680e+02_r8, -5.4433453077032e+00_r8,  7.3772246716454e-01_r8,  4.3986891764158e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,14) = (/  9.3245555431133e+01_r8,  6.1129492298628e+02_r8, -5.1384231104850e+02_r8,  1.2941198165182e+02_r8, &
                            -2.2865477294036e+00_r8, -4.2313649195526e+00_r8,  7.3772246716454e-01_r8,  4.3986891764158e-02_r8/)
! band 15
       absice4(1,:,15) = (/  1.7750967013494e+02_r8, -2.9728258428388e+00_r8,  5.0516550558990e-01_r8,  3.5997106210398e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,15) = (/  1.0355146375894e+03_r8, -6.7009498531165e+02_r8,  1.4559383373911e+02_r8, -3.0349676759069e+01_r8, &
                             1.4171442431211e+01_r8, -3.9488589643285e+00_r8,  5.0516550558990e-01_r8,  3.5997106210398e-02_r8/)
! band 16
       absice4(1,:,16) = (/  1.2774387352402e+02_r8, -2.5427891390490e+00_r8,  6.2719327733184e-01_r8,  3.7617781382362e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       absice4(2,:,16) = (/  1.8195871704194e+03_r8, -1.1191679917223e+03_r8,  1.4250877187512e+02_r8,  2.5663335916348e+01_r8, &
                            -2.0799472615670e+00_r8, -2.1225391366802e+00_r8,  6.2719327733184e-01_r8,  3.7617781382362e-02_r8/)

! extinction units (ext coef/iwc): [(m^-1)/(g m^-3)]
! band  1
       extice4(1,:, 1) = (/  7.0642297430674e+02_r8, -8.2251951415666e+00_r8,  8.3009825258186e-01_r8,  8.3326921971799e-02_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 1) = (/  1.9689976784551e+04_r8, -2.1881302456339e+04_r8,  9.9699025205862e+03_r8, -2.4018719542133e+03_r8, &
                             3.2801184451734e+02_r8, -2.4850363794509e+01_r8,  8.3009825258186e-01_r8,  8.3326921971799e-02_r8/)
! band  2
       extice4(1,:, 2) = (/  3.2076540791610e+02_r8, -1.6446108095976e+01_r8,  2.1086128916677e+00_r8,  1.2049222282429e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 2) = (/  9.8985398214186e+03_r8, -1.6976913936236e+04_r8,  1.0993161359046e+04_r8, -3.5951135796165e+03_r8, &
                             6.4487847779714e+02_r8, -6.1235483586471e+01_r8,  2.1086128916677e+00_r8,  1.2049222282429e-01_r8/)
! band  3
       extice4(1,:, 3) = (/ -3.3940015768006e+02_r8, -1.7644227631264e+01_r8,  3.2569481306562e+00_r8,  1.3674260906854e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 3) = (/ -6.1054298467816e+04_r8,  6.2728003155118e+04_r8, -2.4998839919915e+04_r8,  4.6274226135266e+03_r8, &
                            -3.0842320937379e+02_r8, -1.9830458525643e+01_r8,  3.2569481306562e+00_r8,  1.3674260906854e-01_r8/)
! band  4
       extice4(1,:, 4) = (/ -7.1575311189384e+02_r8, -3.1951270989680e+01_r8,  3.1522136016808e+00_r8,  1.3575431828245e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 4) = (/ -5.0207419959417e+04_r8,  5.8265471164018e+04_r8, -2.7008060051767e+04_r8,  6.2230154559326e+03_r8, &
                            -6.8127889683723e+02_r8,  1.4227122500662e+01_r8,  3.1522136016808e+00_r8,  1.3575431828245e-01_r8/)
! band  5
       extice4(1,:, 5) = (/ -4.9243253320764e+02_r8, -2.5003867653960e+01_r8,  2.9264042558399e+00_r8,  1.2840213460662e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 5) = (/ -7.7637671502239e+03_r8,  1.4581470515451e+04_r8, -9.2236586804328e+03_r8,  2.6881539704723e+03_r8, &
                            -3.6260705829930e+02_r8,  9.8636113517763e+00_r8,  2.9264042558399e+00_r8,  1.2840213460662e-01_r8/)
! band  6
       extice4(1,:, 6) = (/  1.8669724176797e+01_r8, -1.3929319309167e+01_r8,  2.0637171099576e+00_r8,  1.0507539706177e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 6) = (/  2.3078585240734e+03_r8, -1.8727776247408e+02_r8, -8.3002639050863e+02_r8,  3.1653618717193e+02_r8, &
                            -2.3519127177817e+01_r8, -7.8893486859270e+00_r8,  2.0637171099576e+00_r8,  1.0507539706177e-01_r8/)
! band  7
       extice4(1,:, 7) = (/ -2.7753978662197e+01_r8, -1.7855028234597e+01_r8,  2.3482033954918e+00_r8,  1.1721821489540e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 7) = (/ -1.0770560543943e+04_r8,  9.8246849219152e+03_r8, -2.9179964567477e+03_r8,  7.5956847454895e+01_r8, &
                             1.3926345295330e+02_r8, -3.1180815136536e+01_r8,  2.3482033954918e+00_r8,  1.1721821489540e-01_r8/)
! band  8
       extice4(1,:, 8) = (/ -2.7615130153785e+02_r8, -1.3926152712954e+01_r8,  3.1651444401631e+00_r8,  1.3162689595213e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 8) = (/ -4.9516172731536e+04_r8,  5.1477332245212e+04_r8, -2.0870466272835e+04_r8,  3.9749838456408e+03_r8, &
                            -2.8433704577698e+02_r8, -1.6069007018782e+01_r8,  3.1651444401631e+00_r8,  1.3162689595213e-01_r8/)
! band  9
       extice4(1,:, 9) = (/ -6.5342744608070e+02_r8, -2.6977073652688e+01_r8,  3.2374882477429e+00_r8,  1.3524022739846e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:, 9) = (/ -6.6507534864387e+04_r8,  7.3709390217482e+04_r8, -3.2618068046728e+04_r8,  7.1751966487619e+03_r8, &
                            -7.4976368961917e+02_r8,  1.4438977615802e+01_r8,  3.2374882477429e+00_r8,  1.3524022739846e-01_r8/)
! band 10
       extice4(1,:,10) = (/ -8.9390398570828e+02_r8, -3.8712015918337e+01_r8,  3.0636687574413e+00_r8,  1.3382266903802e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,10) = (/ -4.7953386795103e+04_r8,  5.8512635330783e+04_r8, -2.8570809222664e+04_r8,  7.0119723683464e+03_r8, &
                            -8.5345860323169e+02_r8,  3.0570039392808e+01_r8,  3.0636687574413e+00_r8,  1.3382266903802e-01_r8/)
! band 11
       extice4(1,:,11) = (/ -9.2560147670497e+02_r8, -4.0724421174880e+01_r8,  3.0044184616767e+00_r8,  1.3291936646152e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,11) = (/ -3.1766896666758e+04_r8,  4.3095334922533e+04_r8, -2.2947325284572e+04_r8,  6.0746262186972e+03_r8, &
                            -7.9775698474361e+02_r8,  3.2680108252255e+01_r8,  3.0044184616767e+00_r8,  1.3291936646152e-01_r8/)
! band 12
       extice4(1,:,12) = (/ -1.2567575876476e+03_r8, -5.5899659684856e+01_r8,  2.8893601852964e+00_r8,  1.3407344974440e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,12) = (/  1.1740577562475e+04_r8,  4.6284078015019e+03_r8, -1.0614658137762e+04_r8,  4.5455819207981e+03_r8, &
                            -8.0740562663262e+02_r8,  4.7442501326831e+01_r8,  2.8893601852964e+00_r8,  1.3407344974440e-01_r8/)
! band 13
       extice4(1,:,13) = (/ -8.3211126372767e+02_r8, -3.6034286534515e+01_r8,  3.0810709482932e+00_r8,  1.3373539554940e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,13) = (/  7.6729397037912e+04_r8, -6.3034831216861e+04_r8,  1.7450689836424e+04_r8, -1.2541140582288e+03_r8, &
                            -2.2052522410649e+02_r8,  2.6862176417705e+01_r8,  3.0810709482932e+00_r8,  1.3373539554940e-01_r8/)
! band 14
       extice4(1,:,14) = (/ -6.1991602935710e+02_r8, -2.5842543967416e+01_r8,  3.1775658885890e+00_r8,  1.3327721604313e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,14) = (/  8.9208900954035e+04_r8, -7.7524728752411e+04_r8,  2.4235042811446e+04_r8, -2.8597861564663e+03_r8, &
                            -3.0026046446311e+01_r8,  1.8388618646956e+01_r8,  3.1775658885890e+00_r8,  1.3327721604313e-01_r8/)
! band 15
       extice4(1,:,15) = (/ -2.4064622908476e+02_r8, -7.2306654976433e+00_r8,  3.4019926268719e+00_r8,  1.3382075237363e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,15) = (/  1.0358420064354e+05_r8, -9.5872417830220e+04_r8,  3.3583743925545e+04_r8, -5.2464100070897e+03_r8, &
                             2.7416494599033e+02_r8,  3.4469825641715e+00_r8,  3.4019926268719e+00_r8,  1.3382075237363e-01_r8/)
! band 16
       extice4(1,:,16) = (/  2.2183982824014e+02_r8,  1.5632298186052e+01_r8,  3.6691748500883e+00_r8,  1.3402580896768e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       extice4(2,:,16) = (/  6.5639587978015e+04_r8, -6.4625774828595e+04_r8,  2.5190775778991e+04_r8, -4.7851129445107e+03_r8, &
                             4.0595458359608e+02_r8, -1.0508121448097e+01_r8,  3.6691748500883e+00_r8,  1.3402580896768e-01_r8/)

! single-scattering albedo units : unitless
! band  1
       ssaice4(1,:, 1) = (/  2.7992888915477e+03_r8,  7.1281977694363e+01_r8, -2.6638617279438e+00_r8,  4.6751929266116e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 1) = (/  2.5102455060640e+04_r8, -3.5778375263444e+04_r8,  1.9639317859526e+04_r8, -5.3190625180687e+03_r8, &
                             7.2390634474944e+02_r8, -3.5621551943233e+01_r8, -2.6638617279438e+00_r8,  4.6751929266116e-01_r8/)
! band  2
       ssaice4(1,:, 2) = (/  3.1815610798476e+03_r8, -2.5613047978413e+01_r8,  1.2950397190275e+00_r8,  8.1297437201429e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 2) = (/ -3.1186671943573e+02_r8, -5.4749035635971e+03_r8,  5.1284138264248e+03_r8, -1.8749344274072e+03_r8, &
                             3.5816523616435e+02_r8, -3.9948465188978e+01_r8,  1.2950397190275e+00_r8,  8.1297437201429e-01_r8/)
! band  3
       ssaice4(1,:, 3) = (/ -1.5340146907669e+03_r8, -9.6998632132407e+01_r8,  1.8860146810504e+00_r8,  6.6492361585150e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 3) = (/  8.3342035990044e+03_r8, -1.1629318208977e+04_r8,  6.3791527531413e+03_r8, -1.8475831852622e+03_r8, &
                             3.2621302551884e+02_r8, -3.8652119791689e+01_r8,  1.8860146810504e+00_r8,  6.6492361585150e-01_r8/)
! band  4
       ssaice4(1,:, 4) = (/ -1.2917039254030e+03_r8, -4.2078600630886e+01_r8,  1.0344694040207e+00_r8,  5.6763342411629e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 4) = (/ -2.3129429671113e+04_r8,  2.1865069456760e+04_r8, -8.0121006412372e+03_r8,  1.3612241577579e+03_r8, &
                            -7.4092773347133e+01_r8, -1.0386977886490e+01_r8,  1.0344694040207e+00_r8,  5.6763342411629e-01_r8/)
! band  5
       ssaice4(1,:, 5) = (/ -1.4202753877598e+02_r8,  5.0642410031809e+00_r8, -6.0306768603553e-01_r8,  5.0165116208358e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 5) = (/ -1.7683757478756e+04_r8,  1.7866909970112e+04_r8, -7.1863937255739e+03_r8,  1.4385609518488e+03_r8, &
                            -1.3926603379945e+02_r8,  3.9287021968509e+00_r8, -6.0306768603553e-01_r8,  5.0165116208358e-01_r8/)
! band  6
       ssaice4(1,:, 6) = (/  9.4102736143094e+02_r8,  3.3472163487379e+01_r8, -1.7818201871850e+00_r8,  4.6311805303715e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 6) = (/  3.7820082683986e+03_r8, -5.8518630183046e+03_r8,  3.3514665585601e+03_r8, -9.2755448384446e+02_r8, &
                             1.2521838737597e+02_r8, -3.8485288076389e+00_r8, -1.7818201871850e+00_r8,  4.6311805303715e-01_r8/)
! band  7
       ssaice4(1,:, 7) = (/ -2.4018760766913e+02_r8, -8.1419427424145e+01_r8,  1.5175133216955e+00_r8,  6.7386099311082e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 7) = (/  1.5468349963839e+04_r8, -1.9876549184054e+04_r8,  1.0351343898490e+04_r8, -2.8607475891012e+03_r8, &
                             4.6468286684128e+02_r8, -4.5908384160813e+01_r8,  1.5175133216955e+00_r8,  6.7386099311082e-01_r8/)
! band  8
       ssaice4(1,:, 8) = (/ -1.4790113686459e+03_r8, -1.1712394826964e+02_r8,  2.6473193012760e+00_r8,  7.1384501974421e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 8) = (/  1.0579600114800e+04_r8, -1.4533964769145e+04_r8,  8.1445876425581e+03_r8, -2.4494137203721e+03_r8, &
                             4.3822066389628e+02_r8, -4.8997063197557e+01_r8,  2.6473193012760e+00_r8,  7.1384501974421e-01_r8/)
! band  9
       ssaice4(1,:, 9) = (/ -2.0561950999904e+03_r8, -1.1267356037207e+02_r8,  2.8520403865255e+00_r8,  6.8982247044383e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:, 9) = (/ -1.4290331250727e+04_r8,  1.2062768361957e+04_r8, -3.2362834782400e+03_r8,  2.1457296072509e+01_r8, &
                             1.5654937502289e+02_r8, -3.4390112011518e+01_r8,  2.8520403865255e+00_r8,  6.8982247044383e-01_r8/)
! band 10
       ssaice4(1,:,10) = (/ -2.1822257654115e+03_r8, -8.9910637602706e+01_r8,  2.4350127160217e+00_r8,  6.3506939908042e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,10) = (/ -3.9590782432739e+04_r8,  4.0029460129649e+04_r8, -1.5748653850137e+04_r8,  2.9162231426452e+03_r8, &
                            -2.0796307033026e+02_r8, -1.1046625229906e+01_r8,  2.4350127160217e+00_r8,  6.3506939908042e-01_r8/)
! band 11
       ssaice4(1,:,11) = (/ -2.0480331157500e+03_r8, -8.0701448418985e+01_r8,  2.1944164892573e+00_r8,  6.1853116529105e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,11) = (/ -4.2088501965548e+04_r8,  4.3824245988998e+04_r8, -1.7963094005464e+04_r8,  3.5670101954564e+03_r8, &
                            -3.1153069026420e+02_r8, -2.4215893735763e+00_r8,  2.1944164892573e+00_r8,  6.1853116529105e-01_r8/)
! band 12
       ssaice4(1,:,12) = (/ -1.2685736792681e+03_r8, -1.1672538018870e+02_r8,  3.1283075073336e+00_r8,  7.5129919478109e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,12) = (/ -1.3724227450025e+04_r8,  1.4759302602725e+04_r8, -5.9172004154424e+03_r8,  9.6777110413189e+02_r8, &
                             1.5282196820819e+00_r8, -2.2994186504011e+01_r8,  3.1283075073336e+00_r8,  7.5129919478109e-01_r8/)
! band 13
       ssaice4(1,:,13) = (/ -1.9833927654392e+03_r8, -1.0448516623767e+02_r8,  2.7632073757908e+00_r8,  6.8068368855964e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,13) = (/ -9.0825394491614e+03_r8,  1.3173085189049e+04_r8, -7.0766124497007e+03_r8,  1.7578325938517e+03_r8, &
                            -1.7819578275700e+02_r8, -5.5646270136109e+00_r8,  2.7632073757908e+00_r8,  6.8068368855964e-01_r8/)
! band 14
       ssaice4(1,:,14) = (/ -2.0154003443675e+03_r8, -9.9551109358248e+01_r8,  2.7139296724850e+00_r8,  6.6898125836682e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,14) = (/ -1.2667287827518e+03_r8,  5.4526820815200e+03_r8, -4.1529799508696e+03_r8,  1.2471658834147e+03_r8, &
                            -1.4322221755763e+02_r8, -5.1120412063899e+00_r8,  2.7139296724850e+00_r8,  6.6898125836682e-01_r8/)
! band 15
       ssaice4(1,:,15) = (/ -1.1779274424157e+03_r8, -1.0297037420037e+02_r8,  2.9725368910068e+00_r8,  7.3101289187912e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,15) = (/  1.8055666925760e+04_r8, -1.6422717039954e+04_r8,  5.8524375650952e+03_r8, -1.1055347012350e+03_r8, &
                             1.5318875081302e+02_r8, -2.3085494309606e+01_r8,  2.9725368910068e+00_r8,  7.3101289187912e-01_r8/)
! band 16
       ssaice4(1,:,16) = (/ -1.0853435758793e+02_r8, -6.5433019039226e+01_r8,  2.4540684129138e+00_r8,  7.2037469732487e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       ssaice4(2,:,16) = (/  3.4747300576651e+04_r8, -3.5751260857524e+04_r8,  1.4674041570864e+04_r8, -3.0750273417689e+03_r8, &
                             3.6064935644176e+02_r8, -2.8825293609070e+01_r8,  2.4540684129138e+00_r8,  7.2037469732487e-01_r8/)

! asymmetry factor units : unitless
! band  1
       asyice4(1,:, 1) = (/  1.1785375647732e+03_r8,  9.4575590883586e+01_r8, -2.3375325303111e+00_r8,  7.8250719666904e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 1) = (/ -7.0414320186460e+04_r8,  5.0579630552958e+04_r8, -1.0095776138067e+04_r8, -6.2071297853346e+02_r8, &
                             4.4763251873279e+02_r8, -4.4382654665922e+01_r8, -2.3375325303111e+00_r8,  7.8250719666904e-01_r8/)
! band  2
       asyice4(1,:, 2) = (/ -2.8388224121426e+03_r8,  5.4034598369159e+00_r8, -8.0352337143753e-01_r8,  7.6767644559685e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 2) = (/ -2.5606844231934e+05_r8,  2.4814807080467e+05_r8, -9.4347503216990e+04_r8,  1.7616306675951e+04_r8, &
                            -1.6084883429486e+03_r8,  5.3546376018265e+01_r8, -8.0352337143753e-01_r8,  7.6767644559685e-01_r8/)
! band  3
       asyice4(1,:, 3) = (/  5.0315739578391e+02_r8,  5.3142455048027e+01_r8, -1.9532145816670e+00_r8,  8.4039781434454e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 3) = (/ -6.3892678838576e+04_r8,  6.7438827986629e+04_r8, -2.8975609792060e+04_r8,  6.4470601642739e+03_r8, &
                            -7.6590483118622e+02_r8,  4.2889358923489e+01_r8, -1.9532145816670e+00_r8,  8.4039781434454e-01_r8/)
! band  4
       asyice4(1,:, 4) = (/  1.0305544278180e+03_r8,  4.3516394193686e+01_r8, -1.6199938175342e+00_r8,  8.8184033846652e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 4) = (/ -2.1264759705150e+04_r8,  2.0792565416058e+04_r8, -8.4912874769063e+03_r8,  1.8842635847086e+03_r8, &
                            -2.3861815709776e+02_r8,  1.5700809740304e+01_r8, -1.6199938175342e+00_r8,  8.8184033846652e-01_r8/)
! band  5
       asyice4(1,:, 5) = (/  4.6575529976098e+02_r8,  1.6276114457219e+01_r8, -1.1861732072678e+00_r8,  9.0293034948575e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 5) = (/ -1.5909461759446e+04_r8,  1.5290276721336e+04_r8, -6.0143968545773e+03_r8,  1.2371097935377e+03_r8, &
                            -1.3499269861958e+02_r8,  6.3616081593622e+00_r8, -1.1861732072678e+00_r8,  9.0293034948575e-01_r8/)
! band  6
       asyice4(1,:, 6) = (/  4.0966230084914e+02_r8,  1.4983222311175e+01_r8, -1.1010924581378e+00_r8,  9.2977755425733e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 6) = (/  8.1041852865586e+03_r8, -8.9404779089119e+03_r8,  3.5693793823053e+03_r8, -6.1899652864751e+02_r8, &
                             4.1954928353490e+01_r8, -5.3261916156777e-01_r8, -1.1010924581378e+00_r8,  9.2977755425733e-01_r8/)
! band  7
       asyice4(1,:, 7) = (/  8.2615766337341e+02_r8,  5.6544671806908e+01_r8, -1.3770461453160e+00_r8,  8.9728984758841e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 7) = (/  9.8942317662933e+04_r8, -9.8464229691926e+04_r8,  3.7313362473678e+04_r8, -6.4941707643104e+03_r8, &
                             4.5969927449852e+02_r8, -1.5874119703178e+00_r8, -1.3770461453160e+00_r8,  8.9728984758841e-01_r8/)
! band  8
       asyice4(1,:, 8) = (/  1.0368635192446e+03_r8,  7.5772771950089e+01_r8, -1.3182114883294e+00_r8,  8.8011198646507e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 8) = (/  7.7009895423753e+04_r8, -7.5213256098665e+04_r8,  2.7578199537678e+04_r8, -4.4676268236025e+03_r8, &
                             2.4420757967439e+02_r8,  8.7093867121748e+00_r8, -1.3182114883294e+00_r8,  8.8011198646507e-01_r8/)
! band  9
       asyice4(1,:, 9) = (/  1.2497257110604e+03_r8,  7.5047572362840e+01_r8, -1.3194005013220e+00_r8,  8.8929339619440e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:, 9) = (/  3.9910123117240e+04_r8, -3.8136807664452e+04_r8,  1.3255654375476e+04_r8, -1.8294051334532e+03_r8, &
                             1.7543593957863e+01_r8,  1.5967775640937e+01_r8, -1.3194005013220e+00_r8,  8.8929339619440e-01_r8/)
! band 10
       asyice4(1,:,10) = (/  1.1546912281299e+03_r8,  5.4822604099547e+01_r8, -1.2254187560289e+00_r8,  9.1436004079812e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,10) = (/  1.5305288259280e+04_r8, -1.4302595860626e+04_r8,  4.5745802051572e+03_r8, -4.2443984337943e+02_r8, &
                            -6.3091882319723e+01_r8,  1.4473355359629e+01_r8, -1.2254187560289e+00_r8,  9.1436004079812e-01_r8/)
! band 11
       asyice4(1,:,11) = (/  9.4700112027916e+02_r8,  4.3342297560352e+01_r8, -1.1228144838071e+00_r8,  9.2488966978705e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,11) = (/  9.7178530326009e+03_r8, -8.4852828169086e+03_r8,  2.3060048111363e+03_r8, -3.3471599634809e+01_r8, &
                            -8.6389930525525e+01_r8,  1.3806011557560e+01_r8, -1.1228144838071e+00_r8,  9.2488966978705e-01_r8/)
! band 12
       asyice4(1,:,12) = (/  4.0322059567152e+02_r8,  6.2417165690909e+01_r8, -1.3315114510178e+00_r8,  8.6204688036126e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,12) = (/  1.2522620941561e+04_r8, -7.7332693023977e+03_r8, -7.4374029413642e+00_r8,  1.0345159451180e+03_r8, &
                            -2.9895570468508e+02_r8,  3.2712813875463e+01_r8, -1.3315114510178e+00_r8,  8.6204688036126e-01_r8/)
! band 13
       asyice4(1,:,13) = (/  1.0872207412157e+03_r8,  6.1898348982546e+01_r8, -1.4369969293467e+00_r8,  8.9241984555694e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,13) = (/  1.7190457243086e+04_r8, -1.4521059802134e+04_r8,  3.9524691278949e+03_r8, -1.3437352458353e+02_r8, &
                            -1.2177506584726e+02_r8,  2.1644971221772e+01_r8, -1.4369969293467e+00_r8,  8.9241984555694e-01_r8/)
! band 14
       asyice4(1,:,14) = (/  1.2173050632803e+03_r8,  6.2134915415708e+01_r8, -1.4412970715043e+00_r8,  8.9853882311284e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,14) = (/  1.9681032921650e+04_r8, -1.7862239295561e+04_r8,  5.7524510511504e+03_r8, -6.2574989514910e+02_r8, &
                            -5.2656294140074e+01_r8,  1.7567971477352e+01_r8, -1.4412970715043e+00_r8,  8.9853882311284e-01_r8/)
! band 15
       asyice4(1,:,15) = (/  9.8450201984772e+02_r8,  7.7970630060755e+01_r8, -1.5492951241268e+00_r8,  8.5799122520656e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,15) = (/  1.9873268478438e+04_r8, -1.8311276665287e+04_r8,  5.9621136457063e+03_r8, -6.2876451202883e+02_r8, &
                            -7.0191140450781e+01_r8,  2.1346139739238e+01_r8, -1.5492951241268e+00_r8,  8.5799122520656e-01_r8/)
! band 16
       asyice4(1,:,16) = (/  6.5835756830549e+02_r8,  7.5078554758456e+01_r8, -1.4733055979538e+00_r8,  8.3799379857791e-01_r8, &
                             0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8,  0.0000000000000e+00_r8/) 
       asyice4(2,:,16) = (/  1.7476549279007e+04_r8, -1.7139525163433e+04_r8,  6.1834727475563e+03_r8, -8.6730112459330e+02_r8, &
                            -1.7400717666123e+01_r8,  1.7222668006228e+01_r8, -1.4733055979538e+00_r8,  8.3799379857791e-01_r8/)
!<<< yihsuan@UMich 2018-08-28 add end <<<


! For LIQFLAG = 0.
      absliq0 = 0.0903614_r8

! For LIQFLAG = 1.  In each band, the absorption
! coefficients are listed for a range of effective radii from 2.5
! to 59.5 microns in increments of 1.0 micron.
      absliq1(:, 1) = (/ &
! band  1
       1.64047e-03_r8, 6.90533e-02_r8, 7.72017e-02_r8, 7.78054e-02_r8, 7.69523e-02_r8, &
       7.58058e-02_r8, 7.46400e-02_r8, 7.35123e-02_r8, 7.24162e-02_r8, 7.13225e-02_r8, &
       6.99145e-02_r8, 6.66409e-02_r8, 6.36582e-02_r8, 6.09425e-02_r8, 5.84593e-02_r8, &
       5.61743e-02_r8, 5.40571e-02_r8, 5.20812e-02_r8, 5.02245e-02_r8, 4.84680e-02_r8, &
       4.67959e-02_r8, 4.51944e-02_r8, 4.36516e-02_r8, 4.21570e-02_r8, 4.07015e-02_r8, &
       3.92766e-02_r8, 3.78747e-02_r8, 3.64886e-02_r8, 3.53632e-02_r8, 3.41992e-02_r8, &
       3.31016e-02_r8, 3.20643e-02_r8, 3.10817e-02_r8, 3.01490e-02_r8, 2.92620e-02_r8, &
       2.84171e-02_r8, 2.76108e-02_r8, 2.68404e-02_r8, 2.61031e-02_r8, 2.53966e-02_r8, &
       2.47189e-02_r8, 2.40678e-02_r8, 2.34418e-02_r8, 2.28392e-02_r8, 2.22586e-02_r8, &
       2.16986e-02_r8, 2.11580e-02_r8, 2.06356e-02_r8, 2.01305e-02_r8, 1.96417e-02_r8, &
       1.91682e-02_r8, 1.87094e-02_r8, 1.82643e-02_r8, 1.78324e-02_r8, 1.74129e-02_r8, &
       1.70052e-02_r8, 1.66088e-02_r8, 1.62231e-02_r8/)
      absliq1(:, 2) = (/ &
! band  2
       2.19486e-01_r8, 1.80687e-01_r8, 1.59150e-01_r8, 1.44731e-01_r8, 1.33703e-01_r8, &
       1.24355e-01_r8, 1.15756e-01_r8, 1.07318e-01_r8, 9.86119e-02_r8, 8.92739e-02_r8, &
       8.34911e-02_r8, 7.70773e-02_r8, 7.15240e-02_r8, 6.66615e-02_r8, 6.23641e-02_r8, &
       5.85359e-02_r8, 5.51020e-02_r8, 5.20032e-02_r8, 4.91916e-02_r8, 4.66283e-02_r8, &
       4.42813e-02_r8, 4.21236e-02_r8, 4.01330e-02_r8, 3.82905e-02_r8, 3.65797e-02_r8, &
       3.49869e-02_r8, 3.35002e-02_r8, 3.21090e-02_r8, 3.08957e-02_r8, 2.97601e-02_r8, &
       2.86966e-02_r8, 2.76984e-02_r8, 2.67599e-02_r8, 2.58758e-02_r8, 2.50416e-02_r8, &
       2.42532e-02_r8, 2.35070e-02_r8, 2.27997e-02_r8, 2.21284e-02_r8, 2.14904e-02_r8, &
       2.08834e-02_r8, 2.03051e-02_r8, 1.97536e-02_r8, 1.92271e-02_r8, 1.87239e-02_r8, &
       1.82425e-02_r8, 1.77816e-02_r8, 1.73399e-02_r8, 1.69162e-02_r8, 1.65094e-02_r8, &
       1.61187e-02_r8, 1.57430e-02_r8, 1.53815e-02_r8, 1.50334e-02_r8, 1.46981e-02_r8, &
       1.43748e-02_r8, 1.40628e-02_r8, 1.37617e-02_r8/)
      absliq1(:, 3) = (/ &
! band  3
       2.95174e-01_r8, 2.34765e-01_r8, 1.98038e-01_r8, 1.72114e-01_r8, 1.52083e-01_r8, &
       1.35654e-01_r8, 1.21613e-01_r8, 1.09252e-01_r8, 9.81263e-02_r8, 8.79448e-02_r8, &
       8.12566e-02_r8, 7.44563e-02_r8, 6.86374e-02_r8, 6.36042e-02_r8, 5.92094e-02_r8, &
       5.53402e-02_r8, 5.19087e-02_r8, 4.88455e-02_r8, 4.60951e-02_r8, 4.36124e-02_r8, &
       4.13607e-02_r8, 3.93096e-02_r8, 3.74338e-02_r8, 3.57119e-02_r8, 3.41261e-02_r8, &
       3.26610e-02_r8, 3.13036e-02_r8, 3.00425e-02_r8, 2.88497e-02_r8, 2.78077e-02_r8, &
       2.68317e-02_r8, 2.59158e-02_r8, 2.50545e-02_r8, 2.42430e-02_r8, 2.34772e-02_r8, &
       2.27533e-02_r8, 2.20679e-02_r8, 2.14181e-02_r8, 2.08011e-02_r8, 2.02145e-02_r8, &
       1.96561e-02_r8, 1.91239e-02_r8, 1.86161e-02_r8, 1.81311e-02_r8, 1.76673e-02_r8, &
       1.72234e-02_r8, 1.67981e-02_r8, 1.63903e-02_r8, 1.59989e-02_r8, 1.56230e-02_r8, &
       1.52615e-02_r8, 1.49138e-02_r8, 1.45791e-02_r8, 1.42565e-02_r8, 1.39455e-02_r8, &
       1.36455e-02_r8, 1.33559e-02_r8, 1.30761e-02_r8/)
      absliq1(:, 4) = (/ &
! band  4
       3.00925e-01_r8, 2.36949e-01_r8, 1.96947e-01_r8, 1.68692e-01_r8, 1.47190e-01_r8, &
       1.29986e-01_r8, 1.15719e-01_r8, 1.03568e-01_r8, 9.30028e-02_r8, 8.36658e-02_r8, &
       7.71075e-02_r8, 7.07002e-02_r8, 6.52284e-02_r8, 6.05024e-02_r8, 5.63801e-02_r8, &
       5.27534e-02_r8, 4.95384e-02_r8, 4.66690e-02_r8, 4.40925e-02_r8, 4.17664e-02_r8, &
       3.96559e-02_r8, 3.77326e-02_r8, 3.59727e-02_r8, 3.43561e-02_r8, 3.28662e-02_r8, &
       3.14885e-02_r8, 3.02110e-02_r8, 2.90231e-02_r8, 2.78948e-02_r8, 2.69109e-02_r8, &
       2.59884e-02_r8, 2.51217e-02_r8, 2.43058e-02_r8, 2.35364e-02_r8, 2.28096e-02_r8, &
       2.21218e-02_r8, 2.14700e-02_r8, 2.08515e-02_r8, 2.02636e-02_r8, 1.97041e-02_r8, &
       1.91711e-02_r8, 1.86625e-02_r8, 1.81769e-02_r8, 1.77126e-02_r8, 1.72683e-02_r8, &
       1.68426e-02_r8, 1.64344e-02_r8, 1.60427e-02_r8, 1.56664e-02_r8, 1.53046e-02_r8, &
       1.49565e-02_r8, 1.46214e-02_r8, 1.42985e-02_r8, 1.39871e-02_r8, 1.36866e-02_r8, &
       1.33965e-02_r8, 1.31162e-02_r8, 1.28453e-02_r8/)
      absliq1(:, 5) = (/ &
! band  5
       2.64691e-01_r8, 2.12018e-01_r8, 1.78009e-01_r8, 1.53539e-01_r8, 1.34721e-01_r8, &
       1.19580e-01_r8, 1.06996e-01_r8, 9.62772e-02_r8, 8.69710e-02_r8, 7.87670e-02_r8, &
       7.29272e-02_r8, 6.70920e-02_r8, 6.20977e-02_r8, 5.77732e-02_r8, 5.39910e-02_r8, &
       5.06538e-02_r8, 4.76866e-02_r8, 4.50301e-02_r8, 4.26374e-02_r8, 4.04704e-02_r8, &
       3.84981e-02_r8, 3.66948e-02_r8, 3.50394e-02_r8, 3.35141e-02_r8, 3.21038e-02_r8, &
       3.07957e-02_r8, 2.95788e-02_r8, 2.84438e-02_r8, 2.73790e-02_r8, 2.64390e-02_r8, &
       2.55565e-02_r8, 2.47263e-02_r8, 2.39437e-02_r8, 2.32047e-02_r8, 2.25056e-02_r8, &
       2.18433e-02_r8, 2.12149e-02_r8, 2.06177e-02_r8, 2.00495e-02_r8, 1.95081e-02_r8, &
       1.89917e-02_r8, 1.84984e-02_r8, 1.80269e-02_r8, 1.75755e-02_r8, 1.71431e-02_r8, &
       1.67283e-02_r8, 1.63303e-02_r8, 1.59478e-02_r8, 1.55801e-02_r8, 1.52262e-02_r8, &
       1.48853e-02_r8, 1.45568e-02_r8, 1.42400e-02_r8, 1.39342e-02_r8, 1.36388e-02_r8, &
       1.33533e-02_r8, 1.30773e-02_r8, 1.28102e-02_r8/)
      absliq1(:, 6) = (/ &
! band  6
       8.81182e-02_r8, 1.06745e-01_r8, 9.79753e-02_r8, 8.99625e-02_r8, 8.35200e-02_r8, &
       7.81899e-02_r8, 7.35939e-02_r8, 6.94696e-02_r8, 6.56266e-02_r8, 6.19148e-02_r8, &
       5.83355e-02_r8, 5.49306e-02_r8, 5.19642e-02_r8, 4.93325e-02_r8, 4.69659e-02_r8, &
       4.48148e-02_r8, 4.28431e-02_r8, 4.10231e-02_r8, 3.93332e-02_r8, 3.77563e-02_r8, &
       3.62785e-02_r8, 3.48882e-02_r8, 3.35758e-02_r8, 3.23333e-02_r8, 3.11536e-02_r8, &
       3.00310e-02_r8, 2.89601e-02_r8, 2.79365e-02_r8, 2.70502e-02_r8, 2.62618e-02_r8, &
       2.55025e-02_r8, 2.47728e-02_r8, 2.40726e-02_r8, 2.34013e-02_r8, 2.27583e-02_r8, &
       2.21422e-02_r8, 2.15522e-02_r8, 2.09869e-02_r8, 2.04453e-02_r8, 1.99260e-02_r8, &
       1.94280e-02_r8, 1.89501e-02_r8, 1.84913e-02_r8, 1.80506e-02_r8, 1.76270e-02_r8, &
       1.72196e-02_r8, 1.68276e-02_r8, 1.64500e-02_r8, 1.60863e-02_r8, 1.57357e-02_r8, &
       1.53975e-02_r8, 1.50710e-02_r8, 1.47558e-02_r8, 1.44511e-02_r8, 1.41566e-02_r8, &
       1.38717e-02_r8, 1.35960e-02_r8, 1.33290e-02_r8/)
      absliq1(:, 7) = (/ &
! band  7
       4.32174e-02_r8, 7.36078e-02_r8, 6.98340e-02_r8, 6.65231e-02_r8, 6.41948e-02_r8, &
       6.23551e-02_r8, 6.06638e-02_r8, 5.88680e-02_r8, 5.67124e-02_r8, 5.38629e-02_r8, &
       4.99579e-02_r8, 4.86289e-02_r8, 4.70120e-02_r8, 4.52854e-02_r8, 4.35466e-02_r8, &
       4.18480e-02_r8, 4.02169e-02_r8, 3.86658e-02_r8, 3.71992e-02_r8, 3.58168e-02_r8, &
       3.45155e-02_r8, 3.32912e-02_r8, 3.21390e-02_r8, 3.10538e-02_r8, 3.00307e-02_r8, &
       2.90651e-02_r8, 2.81524e-02_r8, 2.72885e-02_r8, 2.62821e-02_r8, 2.55744e-02_r8, &
       2.48799e-02_r8, 2.42029e-02_r8, 2.35460e-02_r8, 2.29108e-02_r8, 2.22981e-02_r8, &
       2.17079e-02_r8, 2.11402e-02_r8, 2.05945e-02_r8, 2.00701e-02_r8, 1.95663e-02_r8, &
       1.90824e-02_r8, 1.86174e-02_r8, 1.81706e-02_r8, 1.77411e-02_r8, 1.73281e-02_r8, &
       1.69307e-02_r8, 1.65483e-02_r8, 1.61801e-02_r8, 1.58254e-02_r8, 1.54835e-02_r8, &
       1.51538e-02_r8, 1.48358e-02_r8, 1.45288e-02_r8, 1.42322e-02_r8, 1.39457e-02_r8, &
       1.36687e-02_r8, 1.34008e-02_r8, 1.31416e-02_r8/)
      absliq1(:, 8) = (/ &
! band  8
       1.41881e-01_r8, 7.15419e-02_r8, 6.30335e-02_r8, 6.11132e-02_r8, 6.01931e-02_r8, &
       5.92420e-02_r8, 5.78968e-02_r8, 5.58876e-02_r8, 5.28923e-02_r8, 4.84462e-02_r8, &
       4.60839e-02_r8, 4.56013e-02_r8, 4.45410e-02_r8, 4.31866e-02_r8, 4.17026e-02_r8, &
       4.01850e-02_r8, 3.86892e-02_r8, 3.72461e-02_r8, 3.58722e-02_r8, 3.45749e-02_r8, &
       3.33564e-02_r8, 3.22155e-02_r8, 3.11494e-02_r8, 3.01541e-02_r8, 2.92253e-02_r8, &
       2.83584e-02_r8, 2.75488e-02_r8, 2.67925e-02_r8, 2.57692e-02_r8, 2.50704e-02_r8, &
       2.43918e-02_r8, 2.37350e-02_r8, 2.31005e-02_r8, 2.24888e-02_r8, 2.18996e-02_r8, &
       2.13325e-02_r8, 2.07870e-02_r8, 2.02623e-02_r8, 1.97577e-02_r8, 1.92724e-02_r8, &
       1.88056e-02_r8, 1.83564e-02_r8, 1.79241e-02_r8, 1.75079e-02_r8, 1.71070e-02_r8, &
       1.67207e-02_r8, 1.63482e-02_r8, 1.59890e-02_r8, 1.56424e-02_r8, 1.53077e-02_r8, &
       1.49845e-02_r8, 1.46722e-02_r8, 1.43702e-02_r8, 1.40782e-02_r8, 1.37955e-02_r8, &
       1.35219e-02_r8, 1.32569e-02_r8, 1.30000e-02_r8/)
      absliq1(:, 9) = (/ &
! band  9
       6.72726e-02_r8, 6.61013e-02_r8, 6.47866e-02_r8, 6.33780e-02_r8, 6.18985e-02_r8, &
       6.03335e-02_r8, 5.86136e-02_r8, 5.65876e-02_r8, 5.39839e-02_r8, 5.03536e-02_r8, &
       4.71608e-02_r8, 4.63630e-02_r8, 4.50313e-02_r8, 4.34526e-02_r8, 4.17876e-02_r8, &
       4.01261e-02_r8, 3.85171e-02_r8, 3.69860e-02_r8, 3.55442e-02_r8, 3.41954e-02_r8, &
       3.29384e-02_r8, 3.17693e-02_r8, 3.06832e-02_r8, 2.96745e-02_r8, 2.87374e-02_r8, &
       2.78662e-02_r8, 2.70557e-02_r8, 2.63008e-02_r8, 2.52450e-02_r8, 2.45424e-02_r8, &
       2.38656e-02_r8, 2.32144e-02_r8, 2.25885e-02_r8, 2.19873e-02_r8, 2.14099e-02_r8, &
       2.08554e-02_r8, 2.03230e-02_r8, 1.98116e-02_r8, 1.93203e-02_r8, 1.88482e-02_r8, &
       1.83944e-02_r8, 1.79578e-02_r8, 1.75378e-02_r8, 1.71335e-02_r8, 1.67440e-02_r8, &
       1.63687e-02_r8, 1.60069e-02_r8, 1.56579e-02_r8, 1.53210e-02_r8, 1.49958e-02_r8, &
       1.46815e-02_r8, 1.43778e-02_r8, 1.40841e-02_r8, 1.37999e-02_r8, 1.35249e-02_r8, &
       1.32585e-02_r8, 1.30004e-02_r8, 1.27502e-02_r8/)
      absliq1(:,10) = (/ &
! band 10
       7.97040e-02_r8, 7.63844e-02_r8, 7.36499e-02_r8, 7.13525e-02_r8, 6.93043e-02_r8, &
       6.72807e-02_r8, 6.50227e-02_r8, 6.22395e-02_r8, 5.86093e-02_r8, 5.37815e-02_r8, &
       5.14682e-02_r8, 4.97214e-02_r8, 4.77392e-02_r8, 4.56961e-02_r8, 4.36858e-02_r8, &
       4.17569e-02_r8, 3.99328e-02_r8, 3.82224e-02_r8, 3.66265e-02_r8, 3.51416e-02_r8, &
       3.37617e-02_r8, 3.24798e-02_r8, 3.12887e-02_r8, 3.01812e-02_r8, 2.91505e-02_r8, &
       2.81900e-02_r8, 2.72939e-02_r8, 2.64568e-02_r8, 2.54165e-02_r8, 2.46832e-02_r8, &
       2.39783e-02_r8, 2.33017e-02_r8, 2.26531e-02_r8, 2.20314e-02_r8, 2.14359e-02_r8, &
       2.08653e-02_r8, 2.03187e-02_r8, 1.97947e-02_r8, 1.92924e-02_r8, 1.88106e-02_r8, &
       1.83483e-02_r8, 1.79043e-02_r8, 1.74778e-02_r8, 1.70678e-02_r8, 1.66735e-02_r8, &
       1.62941e-02_r8, 1.59286e-02_r8, 1.55766e-02_r8, 1.52371e-02_r8, 1.49097e-02_r8, &
       1.45937e-02_r8, 1.42885e-02_r8, 1.39936e-02_r8, 1.37085e-02_r8, 1.34327e-02_r8, &
       1.31659e-02_r8, 1.29075e-02_r8, 1.26571e-02_r8/)
      absliq1(:,11) = (/ &
! band 11
       1.49438e-01_r8, 1.33535e-01_r8, 1.21542e-01_r8, 1.11743e-01_r8, 1.03263e-01_r8, &
       9.55774e-02_r8, 8.83382e-02_r8, 8.12943e-02_r8, 7.42533e-02_r8, 6.70609e-02_r8, &
       6.38761e-02_r8, 5.97788e-02_r8, 5.59841e-02_r8, 5.25318e-02_r8, 4.94132e-02_r8, &
       4.66014e-02_r8, 4.40644e-02_r8, 4.17706e-02_r8, 3.96910e-02_r8, 3.77998e-02_r8, &
       3.60742e-02_r8, 3.44947e-02_r8, 3.30442e-02_r8, 3.17079e-02_r8, 3.04730e-02_r8, &
       2.93283e-02_r8, 2.82642e-02_r8, 2.72720e-02_r8, 2.61789e-02_r8, 2.53277e-02_r8, &
       2.45237e-02_r8, 2.37635e-02_r8, 2.30438e-02_r8, 2.23615e-02_r8, 2.17140e-02_r8, &
       2.10987e-02_r8, 2.05133e-02_r8, 1.99557e-02_r8, 1.94241e-02_r8, 1.89166e-02_r8, &
       1.84317e-02_r8, 1.79679e-02_r8, 1.75238e-02_r8, 1.70983e-02_r8, 1.66901e-02_r8, &
       1.62983e-02_r8, 1.59219e-02_r8, 1.55599e-02_r8, 1.52115e-02_r8, 1.48761e-02_r8, &
       1.45528e-02_r8, 1.42411e-02_r8, 1.39402e-02_r8, 1.36497e-02_r8, 1.33690e-02_r8, &
       1.30976e-02_r8, 1.28351e-02_r8, 1.25810e-02_r8/)
      absliq1(:,12) = (/ &
! band 12
       3.71985e-02_r8, 3.88586e-02_r8, 3.99070e-02_r8, 4.04351e-02_r8, 4.04610e-02_r8, &
       3.99834e-02_r8, 3.89953e-02_r8, 3.74886e-02_r8, 3.54551e-02_r8, 3.28870e-02_r8, &
       3.32576e-02_r8, 3.22444e-02_r8, 3.12384e-02_r8, 3.02584e-02_r8, 2.93146e-02_r8, &
       2.84120e-02_r8, 2.75525e-02_r8, 2.67361e-02_r8, 2.59618e-02_r8, 2.52280e-02_r8, &
       2.45327e-02_r8, 2.38736e-02_r8, 2.32487e-02_r8, 2.26558e-02_r8, 2.20929e-02_r8, &
       2.15579e-02_r8, 2.10491e-02_r8, 2.05648e-02_r8, 1.99749e-02_r8, 1.95704e-02_r8, &
       1.91731e-02_r8, 1.87839e-02_r8, 1.84032e-02_r8, 1.80315e-02_r8, 1.76689e-02_r8, &
       1.73155e-02_r8, 1.69712e-02_r8, 1.66362e-02_r8, 1.63101e-02_r8, 1.59928e-02_r8, &
       1.56842e-02_r8, 1.53840e-02_r8, 1.50920e-02_r8, 1.48080e-02_r8, 1.45318e-02_r8, &
       1.42631e-02_r8, 1.40016e-02_r8, 1.37472e-02_r8, 1.34996e-02_r8, 1.32586e-02_r8, &
       1.30239e-02_r8, 1.27954e-02_r8, 1.25728e-02_r8, 1.23559e-02_r8, 1.21445e-02_r8, &
       1.19385e-02_r8, 1.17376e-02_r8, 1.15417e-02_r8/)
      absliq1(:,13) = (/ &
! band 13
       3.11868e-02_r8, 4.48357e-02_r8, 4.90224e-02_r8, 4.96406e-02_r8, 4.86806e-02_r8, &
       4.69610e-02_r8, 4.48630e-02_r8, 4.25795e-02_r8, 4.02138e-02_r8, 3.78236e-02_r8, &
       3.74266e-02_r8, 3.60384e-02_r8, 3.47074e-02_r8, 3.34434e-02_r8, 3.22499e-02_r8, &
       3.11264e-02_r8, 3.00704e-02_r8, 2.90784e-02_r8, 2.81463e-02_r8, 2.72702e-02_r8, &
       2.64460e-02_r8, 2.56698e-02_r8, 2.49381e-02_r8, 2.42475e-02_r8, 2.35948e-02_r8, &
       2.29774e-02_r8, 2.23925e-02_r8, 2.18379e-02_r8, 2.11793e-02_r8, 2.07076e-02_r8, &
       2.02470e-02_r8, 1.97981e-02_r8, 1.93613e-02_r8, 1.89367e-02_r8, 1.85243e-02_r8, &
       1.81240e-02_r8, 1.77356e-02_r8, 1.73588e-02_r8, 1.69935e-02_r8, 1.66392e-02_r8, &
       1.62956e-02_r8, 1.59624e-02_r8, 1.56393e-02_r8, 1.53259e-02_r8, 1.50219e-02_r8, &
       1.47268e-02_r8, 1.44404e-02_r8, 1.41624e-02_r8, 1.38925e-02_r8, 1.36302e-02_r8, &
       1.33755e-02_r8, 1.31278e-02_r8, 1.28871e-02_r8, 1.26530e-02_r8, 1.24253e-02_r8, &
       1.22038e-02_r8, 1.19881e-02_r8, 1.17782e-02_r8/)
      absliq1(:,14) = (/ &
! band 14
       1.58988e-02_r8, 3.50652e-02_r8, 4.00851e-02_r8, 4.07270e-02_r8, 3.98101e-02_r8, &
       3.83306e-02_r8, 3.66829e-02_r8, 3.50327e-02_r8, 3.34497e-02_r8, 3.19609e-02_r8, &
       3.13712e-02_r8, 3.03348e-02_r8, 2.93415e-02_r8, 2.83973e-02_r8, 2.75037e-02_r8, &
       2.66604e-02_r8, 2.58654e-02_r8, 2.51161e-02_r8, 2.44100e-02_r8, 2.37440e-02_r8, &
       2.31154e-02_r8, 2.25215e-02_r8, 2.19599e-02_r8, 2.14282e-02_r8, 2.09242e-02_r8, &
       2.04459e-02_r8, 1.99915e-02_r8, 1.95594e-02_r8, 1.90254e-02_r8, 1.86598e-02_r8, &
       1.82996e-02_r8, 1.79455e-02_r8, 1.75983e-02_r8, 1.72584e-02_r8, 1.69260e-02_r8, &
       1.66013e-02_r8, 1.62843e-02_r8, 1.59752e-02_r8, 1.56737e-02_r8, 1.53799e-02_r8, &
       1.50936e-02_r8, 1.48146e-02_r8, 1.45429e-02_r8, 1.42782e-02_r8, 1.40203e-02_r8, &
       1.37691e-02_r8, 1.35243e-02_r8, 1.32858e-02_r8, 1.30534e-02_r8, 1.28270e-02_r8, &
       1.26062e-02_r8, 1.23909e-02_r8, 1.21810e-02_r8, 1.19763e-02_r8, 1.17766e-02_r8, &
       1.15817e-02_r8, 1.13915e-02_r8, 1.12058e-02_r8/)
      absliq1(:,15) = (/ &
! band 15
       5.02079e-03_r8, 2.17615e-02_r8, 2.55449e-02_r8, 2.59484e-02_r8, 2.53650e-02_r8, &
       2.45281e-02_r8, 2.36843e-02_r8, 2.29159e-02_r8, 2.22451e-02_r8, 2.16716e-02_r8, &
       2.11451e-02_r8, 2.05817e-02_r8, 2.00454e-02_r8, 1.95372e-02_r8, 1.90567e-02_r8, &
       1.86028e-02_r8, 1.81742e-02_r8, 1.77693e-02_r8, 1.73866e-02_r8, 1.70244e-02_r8, &
       1.66815e-02_r8, 1.63563e-02_r8, 1.60477e-02_r8, 1.57544e-02_r8, 1.54755e-02_r8, &
       1.52097e-02_r8, 1.49564e-02_r8, 1.47146e-02_r8, 1.43684e-02_r8, 1.41728e-02_r8, &
       1.39762e-02_r8, 1.37797e-02_r8, 1.35838e-02_r8, 1.33891e-02_r8, 1.31961e-02_r8, &
       1.30051e-02_r8, 1.28164e-02_r8, 1.26302e-02_r8, 1.24466e-02_r8, 1.22659e-02_r8, &
       1.20881e-02_r8, 1.19131e-02_r8, 1.17412e-02_r8, 1.15723e-02_r8, 1.14063e-02_r8, &
       1.12434e-02_r8, 1.10834e-02_r8, 1.09264e-02_r8, 1.07722e-02_r8, 1.06210e-02_r8, &
       1.04725e-02_r8, 1.03269e-02_r8, 1.01839e-02_r8, 1.00436e-02_r8, 9.90593e-03_r8, &
       9.77080e-03_r8, 9.63818e-03_r8, 9.50800e-03_r8/)
      absliq1(:,16) = (/ &
! band 16
       5.64971e-02_r8, 9.04736e-02_r8, 8.11726e-02_r8, 7.05450e-02_r8, 6.20052e-02_r8, &
       5.54286e-02_r8, 5.03503e-02_r8, 4.63791e-02_r8, 4.32290e-02_r8, 4.06959e-02_r8, &
       3.74690e-02_r8, 3.52964e-02_r8, 3.33799e-02_r8, 3.16774e-02_r8, 3.01550e-02_r8, &
       2.87856e-02_r8, 2.75474e-02_r8, 2.64223e-02_r8, 2.53953e-02_r8, 2.44542e-02_r8, &
       2.35885e-02_r8, 2.27894e-02_r8, 2.20494e-02_r8, 2.13622e-02_r8, 2.07222e-02_r8, &
       2.01246e-02_r8, 1.95654e-02_r8, 1.90408e-02_r8, 1.84398e-02_r8, 1.80021e-02_r8, &
       1.75816e-02_r8, 1.71775e-02_r8, 1.67889e-02_r8, 1.64152e-02_r8, 1.60554e-02_r8, &
       1.57089e-02_r8, 1.53751e-02_r8, 1.50531e-02_r8, 1.47426e-02_r8, 1.44428e-02_r8, &
       1.41532e-02_r8, 1.38734e-02_r8, 1.36028e-02_r8, 1.33410e-02_r8, 1.30875e-02_r8, &
       1.28420e-02_r8, 1.26041e-02_r8, 1.23735e-02_r8, 1.21497e-02_r8, 1.19325e-02_r8, &
       1.17216e-02_r8, 1.15168e-02_r8, 1.13177e-02_r8, 1.11241e-02_r8, 1.09358e-02_r8, &
       1.07525e-02_r8, 1.05741e-02_r8, 1.04003e-02_r8/)

      end subroutine lwcldpr

      end module rrtmg_lw_init

