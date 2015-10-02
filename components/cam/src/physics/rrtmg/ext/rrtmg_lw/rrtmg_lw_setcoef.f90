!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_setcoef.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/22 19:20:05 $
!
      module rrtmg_lw_setcoef

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
      use parrrtm, only : nbndlw, mg, maxxsec, mxmol
      use rrlw_wvn, only: totplnk, totplk16
      use rrlw_ref
      use rrlw_vsn, only: hvrset, hnamset

      implicit none

      contains

!----------------------------------------------------------------------------
      subroutine setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, &
                         coldry, wkl, wbroad, &
                         laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                         colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                         colbrd, fac00, fac01, fac10, fac11, &
                         rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                         rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                         rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                         selffac, selffrac, indself, forfac, forfrac, indfor, &
                         minorfrac, scaleminor, scaleminorn2, indminor)
!----------------------------------------------------------------------------
!
!  Purpose:  For a given atmosphere, calculate the indices and
!  fractions related to the pressure and temperature interpolations.
!  Also calculate the values of the integrated Planck functions 
!  for each band at the level and layer temperatures.

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: nlayers         ! total number of layers
      integer, intent(in) :: istart          ! beginning band of calculation

      real(kind=r8), intent(in) :: pavel(:)           ! layer pressures (mb) 
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: tavel(:)           ! layer temperatures (K)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: tz(0:)             ! level (interface) temperatures (K)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tbound             ! surface temperature (K)
      real(kind=r8), intent(in) :: coldry(:)          ! dry air column density (mol/cm2)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: wbroad(:)          ! broadening gas column density (mol/cm2)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: wkl(:,:)           ! molecular amounts (mol/cm-2)
                                                        !    Dimensions: (mxmol,nlayers)
      real(kind=r8), intent(in) :: semiss(:)          ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)

! ----- Output -----
      integer, intent(out) :: laytrop        ! tropopause layer index
      integer, intent(out) :: jp(:)          ! 
                                                        !    Dimensions: (nlayers)
      integer, intent(out) :: jt(:)          !
                                                        !    Dimensions: (nlayers)
      integer, intent(out) :: jt1(:)         !
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: planklay(:,:)     ! 
                                                        !    Dimensions: (nlayers,nbndlw)
      real(kind=r8), intent(out) :: planklev(0:,:)    ! 
                                                        !    Dimensions: (0:nlayers,nbndlw)
      real(kind=r8), intent(out) :: plankbnd(:)       ! 
                                                        !    Dimensions: (nbndlw)

      real(kind=r8), intent(out) :: colh2o(:)         ! column amount (h2o)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colco2(:)         ! column amount (co2)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colo3(:)          ! column amount (o3)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: coln2o(:)         ! column amount (n2o)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colco(:)          ! column amount (co)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colch4(:)         ! column amount (ch4)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colo2(:)          ! column amount (o2)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: colbrd(:)         ! column amount (broadening gases)
                                                        !    Dimensions: (nlayers)

      integer, intent(out) :: indself(:)
                                                        !    Dimensions: (nlayers)
      integer, intent(out) :: indfor(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: selffac(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: selffrac(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: forfac(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: forfrac(:)
                                                        !    Dimensions: (nlayers)

      integer, intent(out) :: indminor(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: minorfrac(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: scaleminor(:)
                                                        !    Dimensions: (nlayers)
      real(kind=r8), intent(out) :: scaleminorn2(:)
                                                        !    Dimensions: (nlayers)

      real(kind=r8), intent(out) :: &                 !
                         fac00(:), fac01(:), &          !    Dimensions: (nlayers)
                         fac10(:), fac11(:) 
                                                        
      real(kind=r8), intent(out) :: &                 !
                         rat_h2oco2(:),rat_h2oco2_1(:), &
                         rat_h2oo3(:),rat_h2oo3_1(:), & !    Dimensions: (nlayers)
                         rat_h2on2o(:),rat_h2on2o_1(:), &
                         rat_h2och4(:),rat_h2och4_1(:), &
                         rat_n2oco2(:),rat_n2oco2_1(:), &
                         rat_o3co2(:),rat_o3co2_1(:)
                                                        

! ----- Local -----
      integer :: indbound, indlev0
      integer :: lay, indlay, indlev, iband
      integer :: jp1
      real(kind=r8) :: stpfac, tbndfrac, t0frac, tlayfrac, tlevfrac
      real(kind=r8) :: dbdtlev, dbdtlay
      real(kind=r8) :: plog, fp, ft, ft1, water, scalefac, factor, compfp


      hvrset = '$Revision: 1.2 $'

      stpfac = 296._r8/1013._r8

      indbound = tbound - 159._r8
      if (indbound .lt. 1) then
         indbound = 1
      elseif (indbound .gt. 180) then
         indbound = 180
      endif
      tbndfrac = tbound - 159._r8 - float(indbound)
      indlev0 = tz(0) - 159._r8
      if (indlev0 .lt. 1) then
         indlev0 = 1
      elseif (indlev0 .gt. 180) then
         indlev0 = 180
      endif
      t0frac = tz(0) - 159._r8 - float(indlev0)
      laytrop = 0

! Begin layer loop 
!  Calculate the integrated Planck functions for each band at the
!  surface, level, and layer temperatures.
      do lay = 1, nlayers
         indlay = tavel(lay) - 159._r8
         if (indlay .lt. 1) then
            indlay = 1
         elseif (indlay .gt. 180) then
            indlay = 180
         endif
         tlayfrac = tavel(lay) - 159._r8 - float(indlay)
         indlev = tz(lay) - 159._r8
         if (indlev .lt. 1) then
            indlev = 1
         elseif (indlev .gt. 180) then
            indlev = 180
         endif
         tlevfrac = tz(lay) - 159._r8 - float(indlev)

! Begin spectral band loop 
         do iband = 1, 15
            if (lay.eq.1) then
               dbdtlev = totplnk(indbound+1,iband) - totplnk(indbound,iband)
               plankbnd(iband) = semiss(iband) * &
                   (totplnk(indbound,iband) + tbndfrac * dbdtlev)
               dbdtlev = totplnk(indlev0+1,iband)-totplnk(indlev0,iband)
               planklev(0,iband) = totplnk(indlev0,iband) + t0frac * dbdtlev
            endif
            dbdtlev = totplnk(indlev+1,iband) - totplnk(indlev,iband)
            dbdtlay = totplnk(indlay+1,iband) - totplnk(indlay,iband)
            planklay(lay,iband) = totplnk(indlay,iband) + tlayfrac * dbdtlay
            planklev(lay,iband) = totplnk(indlev,iband) + tlevfrac * dbdtlev
         enddo

!  For band 16, if radiative transfer will be performed on just
!  this band, use integrated Planck values up to 3250 cm-1.  
!  If radiative transfer will be performed across all 16 bands,
!  then include in the integrated Planck values for this band
!  contributions from 2600 cm-1 to infinity.
         iband = 16
         if (istart .eq. 16) then
            if (lay.eq.1) then
               dbdtlev = totplk16(indbound+1) - totplk16(indbound)
               plankbnd(iband) = semiss(iband) * &
                    (totplk16(indbound) + tbndfrac * dbdtlev)
               dbdtlev = totplnk(indlev0+1,iband)-totplnk(indlev0,iband)
               planklev(0,iband) = totplk16(indlev0) + &
                    t0frac * dbdtlev
            endif
            dbdtlev = totplk16(indlev+1) - totplk16(indlev)
            dbdtlay = totplk16(indlay+1) - totplk16(indlay)
            planklay(lay,iband) = totplk16(indlay) + tlayfrac * dbdtlay
            planklev(lay,iband) = totplk16(indlev) + tlevfrac * dbdtlev
         else
            if (lay.eq.1) then
               dbdtlev = totplnk(indbound+1,iband) - totplnk(indbound,iband)
               plankbnd(iband) = semiss(iband) * &
                    (totplnk(indbound,iband) + tbndfrac * dbdtlev)
               dbdtlev = totplnk(indlev0+1,iband)-totplnk(indlev0,iband)
               planklev(0,iband) = totplnk(indlev0,iband) + t0frac * dbdtlev
            endif
            dbdtlev = totplnk(indlev+1,iband) - totplnk(indlev,iband)
            dbdtlay = totplnk(indlay+1,iband) - totplnk(indlay,iband)
            planklay(lay,iband) = totplnk(indlay,iband) + tlayfrac * dbdtlay
            planklev(lay,iband) = totplnk(indlev,iband) + tlevfrac * dbdtlev
         endif

!  Find the two reference pressures on either side of the
!  layer pressure.  Store them in JP and JP1.  Store in FP the
!  fraction of the difference (in ln(pressure)) between these
!  two values that the layer pressure lies.
!         plog = alog(pavel(lay))
         plog = dlog(pavel(lay))
         jp(lay) = int(36._r8 - 5*(plog+0.04_r8))
         if (jp(lay) .lt. 1) then
            jp(lay) = 1
         elseif (jp(lay) .gt. 58) then
            jp(lay) = 58
         endif
         jp1 = jp(lay) + 1
         fp = 5._r8 *(preflog(jp(lay)) - plog)

!  Determine, for each reference pressure (JP and JP1), which
!  reference temperature (these are different for each  
!  reference pressure) is nearest the layer temperature but does
!  not exceed it.  Store these indices in JT and JT1, resp.
!  Store in FT (resp. FT1) the fraction of the way between JT
!  (JT1) and the next highest reference temperature that the 
!  layer temperature falls.
         jt(lay) = int(3._r8 + (tavel(lay)-tref(jp(lay)))/15._r8)
         if (jt(lay) .lt. 1) then
            jt(lay) = 1
         elseif (jt(lay) .gt. 4) then
            jt(lay) = 4
         endif
         ft = ((tavel(lay)-tref(jp(lay)))/15._r8) - float(jt(lay)-3)
         jt1(lay) = int(3._r8 + (tavel(lay)-tref(jp1))/15._r8)
         if (jt1(lay) .lt. 1) then
            jt1(lay) = 1
         elseif (jt1(lay) .gt. 4) then
            jt1(lay) = 4
         endif
         ft1 = ((tavel(lay)-tref(jp1))/15._r8) - float(jt1(lay)-3)
         water = wkl(1,lay)/coldry(lay)
         scalefac = pavel(lay) * stpfac / tavel(lay)

!  If the pressure is less than ~100mb, perform a different
!  set of species interpolations.
         if (plog .le. 4.56_r8) go to 5300
         laytrop =  laytrop + 1

         forfac(lay) = scalefac / (1.+water)
         factor = (332.0_r8-tavel(lay))/36.0_r8
         indfor(lay) = min(2, max(1, int(factor)))
         forfrac(lay) = factor - float(indfor(lay))

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(lay) = water * forfac(lay)
         factor = (tavel(lay)-188.0_r8)/7.2_r8
         indself(lay) = min(9, max(1, int(factor)-7))
         selffrac(lay) = factor - float(indself(lay) + 7)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(lay) = pavel(lay)/tavel(lay)
         scaleminorn2(lay) = (pavel(lay)/tavel(lay)) &
             *(wbroad(lay)/(coldry(lay)+wkl(1,lay)))
         factor = (tavel(lay)-180.8_r8)/7.2_r8
         indminor(lay) = min(18, max(1, int(factor)))
         minorfrac(lay) = factor - float(indminor(lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in lower atmosphere.
         rat_h2oco2(lay)=chi_mls(1,jp(lay))/chi_mls(2,jp(lay))
         rat_h2oco2_1(lay)=chi_mls(1,jp(lay)+1)/chi_mls(2,jp(lay)+1)

         rat_h2oo3(lay)=chi_mls(1,jp(lay))/chi_mls(3,jp(lay))
         rat_h2oo3_1(lay)=chi_mls(1,jp(lay)+1)/chi_mls(3,jp(lay)+1)

         rat_h2on2o(lay)=chi_mls(1,jp(lay))/chi_mls(4,jp(lay))
         rat_h2on2o_1(lay)=chi_mls(1,jp(lay)+1)/chi_mls(4,jp(lay)+1)

         rat_h2och4(lay)=chi_mls(1,jp(lay))/chi_mls(6,jp(lay))
         rat_h2och4_1(lay)=chi_mls(1,jp(lay)+1)/chi_mls(6,jp(lay)+1)

         rat_n2oco2(lay)=chi_mls(4,jp(lay))/chi_mls(2,jp(lay))
         rat_n2oco2_1(lay)=chi_mls(4,jp(lay)+1)/chi_mls(2,jp(lay)+1)

!  Calculate needed column amounts.
         colh2o(lay) = 1.e-20_r8 * wkl(1,lay)
         colco2(lay) = 1.e-20_r8 * wkl(2,lay)
         colo3(lay) = 1.e-20_r8 * wkl(3,lay)
         coln2o(lay) = 1.e-20_r8 * wkl(4,lay)
         colco(lay) = 1.e-20_r8 * wkl(5,lay)
         colch4(lay) = 1.e-20_r8 * wkl(6,lay)
         colo2(lay) = 1.e-20_r8 * wkl(7,lay)
         if (colco2(lay) .eq. 0._r8) colco2(lay) = 1.e-32_r8 * coldry(lay)
         if (colo3(lay) .eq. 0._r8) colo3(lay) = 1.e-32_r8 * coldry(lay)
         if (coln2o(lay) .eq. 0._r8) coln2o(lay) = 1.e-32_r8 * coldry(lay)
         if (colco(lay) .eq. 0._r8) colco(lay) = 1.e-32_r8 * coldry(lay)
         if (colch4(lay) .eq. 0._r8) colch4(lay) = 1.e-32_r8 * coldry(lay)
         colbrd(lay) = 1.e-20_r8 * wbroad(lay)
         go to 5400

!  Above laytrop.
 5300    continue

         forfac(lay) = scalefac / (1.+water)
         factor = (tavel(lay)-188.0_r8)/36.0_r8
         indfor(lay) = 3
         forfrac(lay) = factor - 1.0_r8

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(lay) = water * forfac(lay)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(lay) = pavel(lay)/tavel(lay)         
         scaleminorn2(lay) = (pavel(lay)/tavel(lay)) &
             * (wbroad(lay)/(coldry(lay)+wkl(1,lay)))
         factor = (tavel(lay)-180.8_r8)/7.2_r8
         indminor(lay) = min(18, max(1, int(factor)))
         minorfrac(lay) = factor - float(indminor(lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in upper atmosphere.
         rat_h2oco2(lay)=chi_mls(1,jp(lay))/chi_mls(2,jp(lay))
         rat_h2oco2_1(lay)=chi_mls(1,jp(lay)+1)/chi_mls(2,jp(lay)+1)         

         rat_o3co2(lay)=chi_mls(3,jp(lay))/chi_mls(2,jp(lay))
         rat_o3co2_1(lay)=chi_mls(3,jp(lay)+1)/chi_mls(2,jp(lay)+1)         

!  Calculate needed column amounts.
         colh2o(lay) = 1.e-20_r8 * wkl(1,lay)
         colco2(lay) = 1.e-20_r8 * wkl(2,lay)
         colo3(lay) = 1.e-20_r8 * wkl(3,lay)
         coln2o(lay) = 1.e-20_r8 * wkl(4,lay)
         colco(lay) = 1.e-20_r8 * wkl(5,lay)
         colch4(lay) = 1.e-20_r8 * wkl(6,lay)
         colo2(lay) = 1.e-20_r8 * wkl(7,lay)
         if (colco2(lay) .eq. 0._r8) colco2(lay) = 1.e-32_r8 * coldry(lay)
         if (colo3(lay) .eq. 0._r8) colo3(lay) = 1.e-32_r8 * coldry(lay)
         if (coln2o(lay) .eq. 0._r8) coln2o(lay) = 1.e-32_r8 * coldry(lay)
         if (colco(lay)  .eq. 0._r8) colco(lay) = 1.e-32_r8 * coldry(lay)
         if (colch4(lay) .eq. 0._r8) colch4(lay) = 1.e-32_r8 * coldry(lay)
         colbrd(lay) = 1.e-20_r8 * wbroad(lay)
 5400    continue

!  We have now isolated the layer ln pressure and temperature,
!  between two reference pressures and two reference temperatures 
!  (for each reference pressure).  We multiply the pressure 
!  fraction FP with the appropriate temperature fractions to get 
!  the factors that will be needed for the interpolation that yields
!  the optical depths (performed in routines TAUGBn for band n).`

         compfp = 1. - fp
         fac10(lay) = compfp * ft
         fac00(lay) = compfp * (1._r8 - ft)
         fac11(lay) = fp * ft1
         fac01(lay) = fp * (1._r8 - ft1)

!  Rescale selffac and forfac for use in taumol
         selffac(lay) = colh2o(lay)*selffac(lay)
         forfac(lay) = colh2o(lay)*forfac(lay)

! End layer loop
      enddo

      end subroutine setcoef

!***************************************************************************
      subroutine lwatmref
!***************************************************************************

      save
 
! These pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
! each subsequent ln(pressure) differs from the previous one by 0.2.

      pref(:) = (/ &
          1.05363e+03_r8,8.62642e+02_r8,7.06272e+02_r8,5.78246e+02_r8,4.73428e+02_r8, &
          3.87610e+02_r8,3.17348e+02_r8,2.59823e+02_r8,2.12725e+02_r8,1.74164e+02_r8, &
          1.42594e+02_r8,1.16746e+02_r8,9.55835e+01_r8,7.82571e+01_r8,6.40715e+01_r8, &
          5.24573e+01_r8,4.29484e+01_r8,3.51632e+01_r8,2.87892e+01_r8,2.35706e+01_r8, &
          1.92980e+01_r8,1.57998e+01_r8,1.29358e+01_r8,1.05910e+01_r8,8.67114e+00_r8, &
          7.09933e+00_r8,5.81244e+00_r8,4.75882e+00_r8,3.89619e+00_r8,3.18993e+00_r8, &
          2.61170e+00_r8,2.13828e+00_r8,1.75067e+00_r8,1.43333e+00_r8,1.17351e+00_r8, &
          9.60789e-01_r8,7.86628e-01_r8,6.44036e-01_r8,5.27292e-01_r8,4.31710e-01_r8, &
          3.53455e-01_r8,2.89384e-01_r8,2.36928e-01_r8,1.93980e-01_r8,1.58817e-01_r8, &
          1.30029e-01_r8,1.06458e-01_r8,8.71608e-02_r8,7.13612e-02_r8,5.84256e-02_r8, &
          4.78349e-02_r8,3.91639e-02_r8,3.20647e-02_r8,2.62523e-02_r8,2.14936e-02_r8, &
          1.75975e-02_r8,1.44076e-02_r8,1.17959e-02_r8,9.65769e-03_r8/)

      preflog(:) = (/ &
           6.9600e+00_r8, 6.7600e+00_r8, 6.5600e+00_r8, 6.3600e+00_r8, 6.1600e+00_r8, &
           5.9600e+00_r8, 5.7600e+00_r8, 5.5600e+00_r8, 5.3600e+00_r8, 5.1600e+00_r8, &
           4.9600e+00_r8, 4.7600e+00_r8, 4.5600e+00_r8, 4.3600e+00_r8, 4.1600e+00_r8, &
           3.9600e+00_r8, 3.7600e+00_r8, 3.5600e+00_r8, 3.3600e+00_r8, 3.1600e+00_r8, &
           2.9600e+00_r8, 2.7600e+00_r8, 2.5600e+00_r8, 2.3600e+00_r8, 2.1600e+00_r8, &
           1.9600e+00_r8, 1.7600e+00_r8, 1.5600e+00_r8, 1.3600e+00_r8, 1.1600e+00_r8, &
           9.6000e-01_r8, 7.6000e-01_r8, 5.6000e-01_r8, 3.6000e-01_r8, 1.6000e-01_r8, &
          -4.0000e-02_r8,-2.4000e-01_r8,-4.4000e-01_r8,-6.4000e-01_r8,-8.4000e-01_r8, &
          -1.0400e+00_r8,-1.2400e+00_r8,-1.4400e+00_r8,-1.6400e+00_r8,-1.8400e+00_r8, &
          -2.0400e+00_r8,-2.2400e+00_r8,-2.4400e+00_r8,-2.6400e+00_r8,-2.8400e+00_r8, &
          -3.0400e+00_r8,-3.2400e+00_r8,-3.4400e+00_r8,-3.6400e+00_r8,-3.8400e+00_r8, &
          -4.0400e+00_r8,-4.2400e+00_r8,-4.4400e+00_r8,-4.6400e+00_r8/)

! These are the temperatures associated with the respective 
! pressures for the mls standard atmosphere. 

      tref(:) = (/ &
           2.9420e+02_r8, 2.8799e+02_r8, 2.7894e+02_r8, 2.6925e+02_r8, 2.5983e+02_r8, &
           2.5017e+02_r8, 2.4077e+02_r8, 2.3179e+02_r8, 2.2306e+02_r8, 2.1578e+02_r8, &
           2.1570e+02_r8, 2.1570e+02_r8, 2.1570e+02_r8, 2.1706e+02_r8, 2.1858e+02_r8, &
           2.2018e+02_r8, 2.2174e+02_r8, 2.2328e+02_r8, 2.2479e+02_r8, 2.2655e+02_r8, &
           2.2834e+02_r8, 2.3113e+02_r8, 2.3401e+02_r8, 2.3703e+02_r8, 2.4022e+02_r8, &
           2.4371e+02_r8, 2.4726e+02_r8, 2.5085e+02_r8, 2.5457e+02_r8, 2.5832e+02_r8, &
           2.6216e+02_r8, 2.6606e+02_r8, 2.6999e+02_r8, 2.7340e+02_r8, 2.7536e+02_r8, &
           2.7568e+02_r8, 2.7372e+02_r8, 2.7163e+02_r8, 2.6955e+02_r8, 2.6593e+02_r8, &
           2.6211e+02_r8, 2.5828e+02_r8, 2.5360e+02_r8, 2.4854e+02_r8, 2.4348e+02_r8, &
           2.3809e+02_r8, 2.3206e+02_r8, 2.2603e+02_r8, 2.2000e+02_r8, 2.1435e+02_r8, &
           2.0887e+02_r8, 2.0340e+02_r8, 1.9792e+02_r8, 1.9290e+02_r8, 1.8809e+02_r8, &
           1.8329e+02_r8, 1.7849e+02_r8, 1.7394e+02_r8, 1.7212e+02_r8/)

       chi_mls(1,1:12) = (/ &
        1.8760e-02_r8, 1.2223e-02_r8, 5.8909e-03_r8, 2.7675e-03_r8, 1.4065e-03_r8, &
        7.5970e-04_r8, 3.8876e-04_r8, 1.6542e-04_r8, 3.7190e-05_r8, 7.4765e-06_r8, &
        4.3082e-06_r8, 3.3319e-06_r8/)
       chi_mls(1,13:59) = (/ &
        3.2039e-06_r8,  3.1619e-06_r8,  3.2524e-06_r8,  3.4226e-06_r8,  3.6288e-06_r8, &
        3.9148e-06_r8,  4.1488e-06_r8,  4.3081e-06_r8,  4.4420e-06_r8,  4.5778e-06_r8, &
        4.7087e-06_r8,  4.7943e-06_r8,  4.8697e-06_r8,  4.9260e-06_r8,  4.9669e-06_r8, &
        4.9963e-06_r8,  5.0527e-06_r8,  5.1266e-06_r8,  5.2503e-06_r8,  5.3571e-06_r8, &
        5.4509e-06_r8,  5.4830e-06_r8,  5.5000e-06_r8,  5.5000e-06_r8,  5.4536e-06_r8, &
        5.4047e-06_r8,  5.3558e-06_r8,  5.2533e-06_r8,  5.1436e-06_r8,  5.0340e-06_r8, &
        4.8766e-06_r8,  4.6979e-06_r8,  4.5191e-06_r8,  4.3360e-06_r8,  4.1442e-06_r8, &
        3.9523e-06_r8,  3.7605e-06_r8,  3.5722e-06_r8,  3.3855e-06_r8,  3.1988e-06_r8, &
        3.0121e-06_r8,  2.8262e-06_r8,  2.6407e-06_r8,  2.4552e-06_r8,  2.2696e-06_r8, &
        4.3360e-06_r8,  4.1442e-06_r8/)
       chi_mls(2,1:12) = (/ &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8/)
       chi_mls(2,13:59) = (/ &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8,  3.5500e-04_r8, &
        3.5500e-04_r8,  3.5471e-04_r8,  3.5427e-04_r8,  3.5384e-04_r8,  3.5340e-04_r8, &
        3.5500e-04_r8,  3.5500e-04_r8/)
       chi_mls(3,1:12) = (/ &
        3.0170e-08_r8,  3.4725e-08_r8,  4.2477e-08_r8,  5.2759e-08_r8,  6.6944e-08_r8, &
        8.7130e-08_r8,  1.1391e-07_r8,  1.5677e-07_r8,  2.1788e-07_r8,  3.2443e-07_r8, &
        4.6594e-07_r8,  5.6806e-07_r8/)
       chi_mls(3,13:59) = (/ &
        6.9607e-07_r8,  1.1186e-06_r8,  1.7618e-06_r8,  2.3269e-06_r8,  2.9577e-06_r8, &
        3.6593e-06_r8,  4.5950e-06_r8,  5.3189e-06_r8,  5.9618e-06_r8,  6.5113e-06_r8, &
        7.0635e-06_r8,  7.6917e-06_r8,  8.2577e-06_r8,  8.7082e-06_r8,  8.8325e-06_r8, &
        8.7149e-06_r8,  8.0943e-06_r8,  7.3307e-06_r8,  6.3101e-06_r8,  5.3672e-06_r8, &
        4.4829e-06_r8,  3.8391e-06_r8,  3.2827e-06_r8,  2.8235e-06_r8,  2.4906e-06_r8, &
        2.1645e-06_r8,  1.8385e-06_r8,  1.6618e-06_r8,  1.5052e-06_r8,  1.3485e-06_r8, &
        1.1972e-06_r8,  1.0482e-06_r8,  8.9926e-07_r8,  7.6343e-07_r8,  6.5381e-07_r8, &
        5.4419e-07_r8,  4.3456e-07_r8,  3.6421e-07_r8,  3.1194e-07_r8,  2.5967e-07_r8, &
        2.0740e-07_r8,  1.9146e-07_r8,  1.9364e-07_r8,  1.9582e-07_r8,  1.9800e-07_r8, &
        7.6343e-07_r8,  6.5381e-07_r8/)
       chi_mls(4,1:12) = (/ &
        3.2000e-07_r8,  3.2000e-07_r8,  3.2000e-07_r8,  3.2000e-07_r8,  3.2000e-07_r8, &
        3.1965e-07_r8,  3.1532e-07_r8,  3.0383e-07_r8,  2.9422e-07_r8,  2.8495e-07_r8, &
        2.7671e-07_r8,  2.6471e-07_r8/)
       chi_mls(4,13:59) = (/ &
        2.4285e-07_r8,  2.0955e-07_r8,  1.7195e-07_r8,  1.3749e-07_r8,  1.1332e-07_r8, &
        1.0035e-07_r8,  9.1281e-08_r8,  8.5463e-08_r8,  8.0363e-08_r8,  7.3372e-08_r8, &
        6.5975e-08_r8,  5.6039e-08_r8,  4.7090e-08_r8,  3.9977e-08_r8,  3.2979e-08_r8, &
        2.6064e-08_r8,  2.1066e-08_r8,  1.6592e-08_r8,  1.3017e-08_r8,  1.0090e-08_r8, &
        7.6249e-09_r8,  6.1159e-09_r8,  4.6672e-09_r8,  3.2857e-09_r8,  2.8484e-09_r8, &
        2.4620e-09_r8,  2.0756e-09_r8,  1.8551e-09_r8,  1.6568e-09_r8,  1.4584e-09_r8, &
        1.3195e-09_r8,  1.2072e-09_r8,  1.0948e-09_r8,  9.9780e-10_r8,  9.3126e-10_r8, &
        8.6472e-10_r8,  7.9818e-10_r8,  7.5138e-10_r8,  7.1367e-10_r8,  6.7596e-10_r8, &
        6.3825e-10_r8,  6.0981e-10_r8,  5.8600e-10_r8,  5.6218e-10_r8,  5.3837e-10_r8, &
        9.9780e-10_r8,  9.3126e-10_r8/)
       chi_mls(5,1:12) = (/ &
        1.5000e-07_r8,  1.4306e-07_r8,  1.3474e-07_r8,  1.3061e-07_r8,  1.2793e-07_r8, &
        1.2038e-07_r8,  1.0798e-07_r8,  9.4238e-08_r8,  7.9488e-08_r8,  6.1386e-08_r8, &
        4.5563e-08_r8,  3.3475e-08_r8/)
       chi_mls(5,13:59) = (/ &
        2.5118e-08_r8,  1.8671e-08_r8,  1.4349e-08_r8,  1.2501e-08_r8,  1.2407e-08_r8, &
        1.3472e-08_r8,  1.4900e-08_r8,  1.6079e-08_r8,  1.7156e-08_r8,  1.8616e-08_r8, &
        2.0106e-08_r8,  2.1654e-08_r8,  2.3096e-08_r8,  2.4340e-08_r8,  2.5643e-08_r8, &
        2.6990e-08_r8,  2.8456e-08_r8,  2.9854e-08_r8,  3.0943e-08_r8,  3.2023e-08_r8, &
        3.3101e-08_r8,  3.4260e-08_r8,  3.5360e-08_r8,  3.6397e-08_r8,  3.7310e-08_r8, &
        3.8217e-08_r8,  3.9123e-08_r8,  4.1303e-08_r8,  4.3652e-08_r8,  4.6002e-08_r8, &
        5.0289e-08_r8,  5.5446e-08_r8,  6.0603e-08_r8,  6.8946e-08_r8,  8.3652e-08_r8, &
        9.8357e-08_r8,  1.1306e-07_r8,  1.4766e-07_r8,  1.9142e-07_r8,  2.3518e-07_r8, &
        2.7894e-07_r8,  3.5001e-07_r8,  4.3469e-07_r8,  5.1938e-07_r8,  6.0407e-07_r8, &
        6.8946e-08_r8,  8.3652e-08_r8/)
       chi_mls(6,1:12) = (/ &
        1.7000e-06_r8,  1.7000e-06_r8,  1.6999e-06_r8,  1.6904e-06_r8,  1.6671e-06_r8, &
        1.6351e-06_r8,  1.6098e-06_r8,  1.5590e-06_r8,  1.5120e-06_r8,  1.4741e-06_r8, &
        1.4385e-06_r8,  1.4002e-06_r8/)
       chi_mls(6,13:59) = (/ &
        1.3573e-06_r8,  1.3130e-06_r8,  1.2512e-06_r8,  1.1668e-06_r8,  1.0553e-06_r8, &
        9.3281e-07_r8,  8.1217e-07_r8,  7.5239e-07_r8,  7.0728e-07_r8,  6.6722e-07_r8, &
        6.2733e-07_r8,  5.8604e-07_r8,  5.4769e-07_r8,  5.1480e-07_r8,  4.8206e-07_r8, &
        4.4943e-07_r8,  4.1702e-07_r8,  3.8460e-07_r8,  3.5200e-07_r8,  3.1926e-07_r8, &
        2.8646e-07_r8,  2.5498e-07_r8,  2.2474e-07_r8,  1.9588e-07_r8,  1.8295e-07_r8, &
        1.7089e-07_r8,  1.5882e-07_r8,  1.5536e-07_r8,  1.5304e-07_r8,  1.5072e-07_r8, &
        1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8, &
        1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8, &
        1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8,  1.5000e-07_r8, &
        1.5000e-07_r8,  1.5000e-07_r8/)
       chi_mls(7,1:12) = (/ &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8/)
       chi_mls(7,13:59) = (/ &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8,  0.2090_r8, &
        0.2090_r8,  0.2090_r8/)

      end subroutine lwatmref

!***************************************************************************
      subroutine lwavplank
!***************************************************************************

      save
 
      totplnk(1:50,  1) = (/ &
      0.14783e-05_r8,0.15006e-05_r8,0.15230e-05_r8,0.15455e-05_r8,0.15681e-05_r8, &
      0.15908e-05_r8,0.16136e-05_r8,0.16365e-05_r8,0.16595e-05_r8,0.16826e-05_r8, &
      0.17059e-05_r8,0.17292e-05_r8,0.17526e-05_r8,0.17762e-05_r8,0.17998e-05_r8, &
      0.18235e-05_r8,0.18473e-05_r8,0.18712e-05_r8,0.18953e-05_r8,0.19194e-05_r8, &
      0.19435e-05_r8,0.19678e-05_r8,0.19922e-05_r8,0.20166e-05_r8,0.20412e-05_r8, &
      0.20658e-05_r8,0.20905e-05_r8,0.21153e-05_r8,0.21402e-05_r8,0.21652e-05_r8, &
      0.21902e-05_r8,0.22154e-05_r8,0.22406e-05_r8,0.22659e-05_r8,0.22912e-05_r8, &
      0.23167e-05_r8,0.23422e-05_r8,0.23678e-05_r8,0.23934e-05_r8,0.24192e-05_r8, &
      0.24450e-05_r8,0.24709e-05_r8,0.24968e-05_r8,0.25229e-05_r8,0.25490e-05_r8, &
      0.25751e-05_r8,0.26014e-05_r8,0.26277e-05_r8,0.26540e-05_r8,0.26805e-05_r8/)
      totplnk(51:100,  1) = (/ &
      0.27070e-05_r8,0.27335e-05_r8,0.27602e-05_r8,0.27869e-05_r8,0.28136e-05_r8, &
      0.28404e-05_r8,0.28673e-05_r8,0.28943e-05_r8,0.29213e-05_r8,0.29483e-05_r8, &
      0.29754e-05_r8,0.30026e-05_r8,0.30298e-05_r8,0.30571e-05_r8,0.30845e-05_r8, &
      0.31119e-05_r8,0.31393e-05_r8,0.31669e-05_r8,0.31944e-05_r8,0.32220e-05_r8, &
      0.32497e-05_r8,0.32774e-05_r8,0.33052e-05_r8,0.33330e-05_r8,0.33609e-05_r8, &
      0.33888e-05_r8,0.34168e-05_r8,0.34448e-05_r8,0.34729e-05_r8,0.35010e-05_r8, &
      0.35292e-05_r8,0.35574e-05_r8,0.35857e-05_r8,0.36140e-05_r8,0.36424e-05_r8, &
      0.36708e-05_r8,0.36992e-05_r8,0.37277e-05_r8,0.37563e-05_r8,0.37848e-05_r8, &
      0.38135e-05_r8,0.38421e-05_r8,0.38708e-05_r8,0.38996e-05_r8,0.39284e-05_r8, &
      0.39572e-05_r8,0.39861e-05_r8,0.40150e-05_r8,0.40440e-05_r8,0.40730e-05_r8/)
      totplnk(101:150,  1) = (/ &
      0.41020e-05_r8,0.41311e-05_r8,0.41602e-05_r8,0.41893e-05_r8,0.42185e-05_r8, &
      0.42477e-05_r8,0.42770e-05_r8,0.43063e-05_r8,0.43356e-05_r8,0.43650e-05_r8, &
      0.43944e-05_r8,0.44238e-05_r8,0.44533e-05_r8,0.44828e-05_r8,0.45124e-05_r8, &
      0.45419e-05_r8,0.45715e-05_r8,0.46012e-05_r8,0.46309e-05_r8,0.46606e-05_r8, &
      0.46903e-05_r8,0.47201e-05_r8,0.47499e-05_r8,0.47797e-05_r8,0.48096e-05_r8, &
      0.48395e-05_r8,0.48695e-05_r8,0.48994e-05_r8,0.49294e-05_r8,0.49594e-05_r8, &
      0.49895e-05_r8,0.50196e-05_r8,0.50497e-05_r8,0.50798e-05_r8,0.51100e-05_r8, &
      0.51402e-05_r8,0.51704e-05_r8,0.52007e-05_r8,0.52309e-05_r8,0.52612e-05_r8, &
      0.52916e-05_r8,0.53219e-05_r8,0.53523e-05_r8,0.53827e-05_r8,0.54132e-05_r8, &
      0.54436e-05_r8,0.54741e-05_r8,0.55047e-05_r8,0.55352e-05_r8,0.55658e-05_r8/)
      totplnk(151:181,  1) = (/ &
      0.55964e-05_r8,0.56270e-05_r8,0.56576e-05_r8,0.56883e-05_r8,0.57190e-05_r8, &
      0.57497e-05_r8,0.57804e-05_r8,0.58112e-05_r8,0.58420e-05_r8,0.58728e-05_r8, &
      0.59036e-05_r8,0.59345e-05_r8,0.59653e-05_r8,0.59962e-05_r8,0.60272e-05_r8, &
      0.60581e-05_r8,0.60891e-05_r8,0.61201e-05_r8,0.61511e-05_r8,0.61821e-05_r8, &
      0.62131e-05_r8,0.62442e-05_r8,0.62753e-05_r8,0.63064e-05_r8,0.63376e-05_r8, &
      0.63687e-05_r8,0.63998e-05_r8,0.64310e-05_r8,0.64622e-05_r8,0.64935e-05_r8, &
      0.65247e-05_r8/)
      totplnk(1:50,  2) = (/ &
      0.20262e-05_r8,0.20757e-05_r8,0.21257e-05_r8,0.21763e-05_r8,0.22276e-05_r8, &
      0.22794e-05_r8,0.23319e-05_r8,0.23849e-05_r8,0.24386e-05_r8,0.24928e-05_r8, &
      0.25477e-05_r8,0.26031e-05_r8,0.26591e-05_r8,0.27157e-05_r8,0.27728e-05_r8, &
      0.28306e-05_r8,0.28889e-05_r8,0.29478e-05_r8,0.30073e-05_r8,0.30673e-05_r8, &
      0.31279e-05_r8,0.31890e-05_r8,0.32507e-05_r8,0.33129e-05_r8,0.33757e-05_r8, &
      0.34391e-05_r8,0.35029e-05_r8,0.35674e-05_r8,0.36323e-05_r8,0.36978e-05_r8, &
      0.37638e-05_r8,0.38304e-05_r8,0.38974e-05_r8,0.39650e-05_r8,0.40331e-05_r8, &
      0.41017e-05_r8,0.41708e-05_r8,0.42405e-05_r8,0.43106e-05_r8,0.43812e-05_r8, &
      0.44524e-05_r8,0.45240e-05_r8,0.45961e-05_r8,0.46687e-05_r8,0.47418e-05_r8, &
      0.48153e-05_r8,0.48894e-05_r8,0.49639e-05_r8,0.50389e-05_r8,0.51143e-05_r8/)
      totplnk(51:100,  2) = (/ &
      0.51902e-05_r8,0.52666e-05_r8,0.53434e-05_r8,0.54207e-05_r8,0.54985e-05_r8, &
      0.55767e-05_r8,0.56553e-05_r8,0.57343e-05_r8,0.58139e-05_r8,0.58938e-05_r8, &
      0.59742e-05_r8,0.60550e-05_r8,0.61362e-05_r8,0.62179e-05_r8,0.63000e-05_r8, &
      0.63825e-05_r8,0.64654e-05_r8,0.65487e-05_r8,0.66324e-05_r8,0.67166e-05_r8, &
      0.68011e-05_r8,0.68860e-05_r8,0.69714e-05_r8,0.70571e-05_r8,0.71432e-05_r8, &
      0.72297e-05_r8,0.73166e-05_r8,0.74039e-05_r8,0.74915e-05_r8,0.75796e-05_r8, &
      0.76680e-05_r8,0.77567e-05_r8,0.78459e-05_r8,0.79354e-05_r8,0.80252e-05_r8, &
      0.81155e-05_r8,0.82061e-05_r8,0.82970e-05_r8,0.83883e-05_r8,0.84799e-05_r8, &
      0.85719e-05_r8,0.86643e-05_r8,0.87569e-05_r8,0.88499e-05_r8,0.89433e-05_r8, &
      0.90370e-05_r8,0.91310e-05_r8,0.92254e-05_r8,0.93200e-05_r8,0.94150e-05_r8/)
      totplnk(101:150,  2) = (/ &
      0.95104e-05_r8,0.96060e-05_r8,0.97020e-05_r8,0.97982e-05_r8,0.98948e-05_r8, &
      0.99917e-05_r8,0.10089e-04_r8,0.10186e-04_r8,0.10284e-04_r8,0.10382e-04_r8, &
      0.10481e-04_r8,0.10580e-04_r8,0.10679e-04_r8,0.10778e-04_r8,0.10877e-04_r8, &
      0.10977e-04_r8,0.11077e-04_r8,0.11178e-04_r8,0.11279e-04_r8,0.11380e-04_r8, &
      0.11481e-04_r8,0.11583e-04_r8,0.11684e-04_r8,0.11786e-04_r8,0.11889e-04_r8, &
      0.11992e-04_r8,0.12094e-04_r8,0.12198e-04_r8,0.12301e-04_r8,0.12405e-04_r8, &
      0.12509e-04_r8,0.12613e-04_r8,0.12717e-04_r8,0.12822e-04_r8,0.12927e-04_r8, &
      0.13032e-04_r8,0.13138e-04_r8,0.13244e-04_r8,0.13349e-04_r8,0.13456e-04_r8, &
      0.13562e-04_r8,0.13669e-04_r8,0.13776e-04_r8,0.13883e-04_r8,0.13990e-04_r8, &
      0.14098e-04_r8,0.14206e-04_r8,0.14314e-04_r8,0.14422e-04_r8,0.14531e-04_r8/)
      totplnk(151:181,  2) = (/ &
      0.14639e-04_r8,0.14748e-04_r8,0.14857e-04_r8,0.14967e-04_r8,0.15076e-04_r8, &
      0.15186e-04_r8,0.15296e-04_r8,0.15407e-04_r8,0.15517e-04_r8,0.15628e-04_r8, &
      0.15739e-04_r8,0.15850e-04_r8,0.15961e-04_r8,0.16072e-04_r8,0.16184e-04_r8, &
      0.16296e-04_r8,0.16408e-04_r8,0.16521e-04_r8,0.16633e-04_r8,0.16746e-04_r8, &
      0.16859e-04_r8,0.16972e-04_r8,0.17085e-04_r8,0.17198e-04_r8,0.17312e-04_r8, &
      0.17426e-04_r8,0.17540e-04_r8,0.17654e-04_r8,0.17769e-04_r8,0.17883e-04_r8, &
      0.17998e-04_r8/)
      totplnk(1:50, 3) = (/ &
      1.34822e-06_r8,1.39134e-06_r8,1.43530e-06_r8,1.48010e-06_r8,1.52574e-06_r8, &
      1.57222e-06_r8,1.61956e-06_r8,1.66774e-06_r8,1.71678e-06_r8,1.76666e-06_r8, &
      1.81741e-06_r8,1.86901e-06_r8,1.92147e-06_r8,1.97479e-06_r8,2.02898e-06_r8, &
      2.08402e-06_r8,2.13993e-06_r8,2.19671e-06_r8,2.25435e-06_r8,2.31285e-06_r8, &
      2.37222e-06_r8,2.43246e-06_r8,2.49356e-06_r8,2.55553e-06_r8,2.61837e-06_r8, &
      2.68207e-06_r8,2.74664e-06_r8,2.81207e-06_r8,2.87837e-06_r8,2.94554e-06_r8, &
      3.01356e-06_r8,3.08245e-06_r8,3.15221e-06_r8,3.22282e-06_r8,3.29429e-06_r8, &
      3.36662e-06_r8,3.43982e-06_r8,3.51386e-06_r8,3.58876e-06_r8,3.66451e-06_r8, &
      3.74112e-06_r8,3.81857e-06_r8,3.89688e-06_r8,3.97602e-06_r8,4.05601e-06_r8, &
      4.13685e-06_r8,4.21852e-06_r8,4.30104e-06_r8,4.38438e-06_r8,4.46857e-06_r8/)
      totplnk(51:100, 3) = (/ &
      4.55358e-06_r8,4.63943e-06_r8,4.72610e-06_r8,4.81359e-06_r8,4.90191e-06_r8, &
      4.99105e-06_r8,5.08100e-06_r8,5.17176e-06_r8,5.26335e-06_r8,5.35573e-06_r8, &
      5.44892e-06_r8,5.54292e-06_r8,5.63772e-06_r8,5.73331e-06_r8,5.82970e-06_r8, &
      5.92688e-06_r8,6.02485e-06_r8,6.12360e-06_r8,6.22314e-06_r8,6.32346e-06_r8, &
      6.42455e-06_r8,6.52641e-06_r8,6.62906e-06_r8,6.73247e-06_r8,6.83664e-06_r8, &
      6.94156e-06_r8,7.04725e-06_r8,7.15370e-06_r8,7.26089e-06_r8,7.36883e-06_r8, &
      7.47752e-06_r8,7.58695e-06_r8,7.69712e-06_r8,7.80801e-06_r8,7.91965e-06_r8, &
      8.03201e-06_r8,8.14510e-06_r8,8.25891e-06_r8,8.37343e-06_r8,8.48867e-06_r8, &
      8.60463e-06_r8,8.72128e-06_r8,8.83865e-06_r8,8.95672e-06_r8,9.07548e-06_r8, &
      9.19495e-06_r8,9.31510e-06_r8,9.43594e-06_r8,9.55745e-06_r8,9.67966e-06_r8/)
      totplnk(101:150, 3) = (/ &
      9.80254e-06_r8,9.92609e-06_r8,1.00503e-05_r8,1.01752e-05_r8,1.03008e-05_r8, &
      1.04270e-05_r8,1.05539e-05_r8,1.06814e-05_r8,1.08096e-05_r8,1.09384e-05_r8, &
      1.10679e-05_r8,1.11980e-05_r8,1.13288e-05_r8,1.14601e-05_r8,1.15922e-05_r8, &
      1.17248e-05_r8,1.18581e-05_r8,1.19920e-05_r8,1.21265e-05_r8,1.22616e-05_r8, &
      1.23973e-05_r8,1.25337e-05_r8,1.26706e-05_r8,1.28081e-05_r8,1.29463e-05_r8, &
      1.30850e-05_r8,1.32243e-05_r8,1.33642e-05_r8,1.35047e-05_r8,1.36458e-05_r8, &
      1.37875e-05_r8,1.39297e-05_r8,1.40725e-05_r8,1.42159e-05_r8,1.43598e-05_r8, &
      1.45044e-05_r8,1.46494e-05_r8,1.47950e-05_r8,1.49412e-05_r8,1.50879e-05_r8, &
      1.52352e-05_r8,1.53830e-05_r8,1.55314e-05_r8,1.56803e-05_r8,1.58297e-05_r8, &
      1.59797e-05_r8,1.61302e-05_r8,1.62812e-05_r8,1.64327e-05_r8,1.65848e-05_r8/)
      totplnk(151:181, 3) = (/ &
      1.67374e-05_r8,1.68904e-05_r8,1.70441e-05_r8,1.71982e-05_r8,1.73528e-05_r8, &
      1.75079e-05_r8,1.76635e-05_r8,1.78197e-05_r8,1.79763e-05_r8,1.81334e-05_r8, &
      1.82910e-05_r8,1.84491e-05_r8,1.86076e-05_r8,1.87667e-05_r8,1.89262e-05_r8, &
      1.90862e-05_r8,1.92467e-05_r8,1.94076e-05_r8,1.95690e-05_r8,1.97309e-05_r8, &
      1.98932e-05_r8,2.00560e-05_r8,2.02193e-05_r8,2.03830e-05_r8,2.05472e-05_r8, &
      2.07118e-05_r8,2.08768e-05_r8,2.10423e-05_r8,2.12083e-05_r8,2.13747e-05_r8, &
      2.15414e-05_r8/)
      totplnk(1:50, 4) = (/ &
      8.90528e-07_r8,9.24222e-07_r8,9.58757e-07_r8,9.94141e-07_r8,1.03038e-06_r8, &
      1.06748e-06_r8,1.10545e-06_r8,1.14430e-06_r8,1.18403e-06_r8,1.22465e-06_r8, &
      1.26618e-06_r8,1.30860e-06_r8,1.35193e-06_r8,1.39619e-06_r8,1.44136e-06_r8, &
      1.48746e-06_r8,1.53449e-06_r8,1.58246e-06_r8,1.63138e-06_r8,1.68124e-06_r8, &
      1.73206e-06_r8,1.78383e-06_r8,1.83657e-06_r8,1.89028e-06_r8,1.94495e-06_r8, &
      2.00060e-06_r8,2.05724e-06_r8,2.11485e-06_r8,2.17344e-06_r8,2.23303e-06_r8, &
      2.29361e-06_r8,2.35519e-06_r8,2.41777e-06_r8,2.48134e-06_r8,2.54592e-06_r8, &
      2.61151e-06_r8,2.67810e-06_r8,2.74571e-06_r8,2.81433e-06_r8,2.88396e-06_r8, &
      2.95461e-06_r8,3.02628e-06_r8,3.09896e-06_r8,3.17267e-06_r8,3.24741e-06_r8, &
      3.32316e-06_r8,3.39994e-06_r8,3.47774e-06_r8,3.55657e-06_r8,3.63642e-06_r8/)
      totplnk(51:100, 4) = (/ &
      3.71731e-06_r8,3.79922e-06_r8,3.88216e-06_r8,3.96612e-06_r8,4.05112e-06_r8, &
      4.13714e-06_r8,4.22419e-06_r8,4.31227e-06_r8,4.40137e-06_r8,4.49151e-06_r8, &
      4.58266e-06_r8,4.67485e-06_r8,4.76806e-06_r8,4.86229e-06_r8,4.95754e-06_r8, &
      5.05383e-06_r8,5.15113e-06_r8,5.24946e-06_r8,5.34879e-06_r8,5.44916e-06_r8, &
      5.55053e-06_r8,5.65292e-06_r8,5.75632e-06_r8,5.86073e-06_r8,5.96616e-06_r8, &
      6.07260e-06_r8,6.18003e-06_r8,6.28848e-06_r8,6.39794e-06_r8,6.50838e-06_r8, &
      6.61983e-06_r8,6.73229e-06_r8,6.84573e-06_r8,6.96016e-06_r8,7.07559e-06_r8, &
      7.19200e-06_r8,7.30940e-06_r8,7.42779e-06_r8,7.54715e-06_r8,7.66749e-06_r8, &
      7.78882e-06_r8,7.91110e-06_r8,8.03436e-06_r8,8.15859e-06_r8,8.28379e-06_r8, &
      8.40994e-06_r8,8.53706e-06_r8,8.66515e-06_r8,8.79418e-06_r8,8.92416e-06_r8/)
      totplnk(101:150, 4) = (/ &
      9.05510e-06_r8,9.18697e-06_r8,9.31979e-06_r8,9.45356e-06_r8,9.58826e-06_r8, &
      9.72389e-06_r8,9.86046e-06_r8,9.99793e-06_r8,1.01364e-05_r8,1.02757e-05_r8, &
      1.04159e-05_r8,1.05571e-05_r8,1.06992e-05_r8,1.08422e-05_r8,1.09861e-05_r8, &
      1.11309e-05_r8,1.12766e-05_r8,1.14232e-05_r8,1.15707e-05_r8,1.17190e-05_r8, &
      1.18683e-05_r8,1.20184e-05_r8,1.21695e-05_r8,1.23214e-05_r8,1.24741e-05_r8, &
      1.26277e-05_r8,1.27822e-05_r8,1.29376e-05_r8,1.30939e-05_r8,1.32509e-05_r8, &
      1.34088e-05_r8,1.35676e-05_r8,1.37273e-05_r8,1.38877e-05_r8,1.40490e-05_r8, &
      1.42112e-05_r8,1.43742e-05_r8,1.45380e-05_r8,1.47026e-05_r8,1.48680e-05_r8, &
      1.50343e-05_r8,1.52014e-05_r8,1.53692e-05_r8,1.55379e-05_r8,1.57074e-05_r8, &
      1.58778e-05_r8,1.60488e-05_r8,1.62207e-05_r8,1.63934e-05_r8,1.65669e-05_r8/)
      totplnk(151:181, 4) = (/ &
      1.67411e-05_r8,1.69162e-05_r8,1.70920e-05_r8,1.72685e-05_r8,1.74459e-05_r8, &
      1.76240e-05_r8,1.78029e-05_r8,1.79825e-05_r8,1.81629e-05_r8,1.83440e-05_r8, &
      1.85259e-05_r8,1.87086e-05_r8,1.88919e-05_r8,1.90760e-05_r8,1.92609e-05_r8, &
      1.94465e-05_r8,1.96327e-05_r8,1.98199e-05_r8,2.00076e-05_r8,2.01961e-05_r8, &
      2.03853e-05_r8,2.05752e-05_r8,2.07658e-05_r8,2.09571e-05_r8,2.11491e-05_r8, &
      2.13418e-05_r8,2.15352e-05_r8,2.17294e-05_r8,2.19241e-05_r8,2.21196e-05_r8, &
      2.23158e-05_r8/)
      totplnk(1:50, 5) = (/ &
      5.70230e-07_r8,5.94788e-07_r8,6.20085e-07_r8,6.46130e-07_r8,6.72936e-07_r8, &
      7.00512e-07_r8,7.28869e-07_r8,7.58019e-07_r8,7.87971e-07_r8,8.18734e-07_r8, &
      8.50320e-07_r8,8.82738e-07_r8,9.15999e-07_r8,9.50110e-07_r8,9.85084e-07_r8, &
      1.02093e-06_r8,1.05765e-06_r8,1.09527e-06_r8,1.13378e-06_r8,1.17320e-06_r8, &
      1.21353e-06_r8,1.25479e-06_r8,1.29698e-06_r8,1.34011e-06_r8,1.38419e-06_r8, &
      1.42923e-06_r8,1.47523e-06_r8,1.52221e-06_r8,1.57016e-06_r8,1.61910e-06_r8, &
      1.66904e-06_r8,1.71997e-06_r8,1.77192e-06_r8,1.82488e-06_r8,1.87886e-06_r8, &
      1.93387e-06_r8,1.98991e-06_r8,2.04699e-06_r8,2.10512e-06_r8,2.16430e-06_r8, &
      2.22454e-06_r8,2.28584e-06_r8,2.34821e-06_r8,2.41166e-06_r8,2.47618e-06_r8, &
      2.54178e-06_r8,2.60847e-06_r8,2.67626e-06_r8,2.74514e-06_r8,2.81512e-06_r8/)
      totplnk(51:100, 5) = (/ &
      2.88621e-06_r8,2.95841e-06_r8,3.03172e-06_r8,3.10615e-06_r8,3.18170e-06_r8, &
      3.25838e-06_r8,3.33618e-06_r8,3.41511e-06_r8,3.49518e-06_r8,3.57639e-06_r8, &
      3.65873e-06_r8,3.74221e-06_r8,3.82684e-06_r8,3.91262e-06_r8,3.99955e-06_r8, &
      4.08763e-06_r8,4.17686e-06_r8,4.26725e-06_r8,4.35880e-06_r8,4.45150e-06_r8, &
      4.54537e-06_r8,4.64039e-06_r8,4.73659e-06_r8,4.83394e-06_r8,4.93246e-06_r8, &
      5.03215e-06_r8,5.13301e-06_r8,5.23504e-06_r8,5.33823e-06_r8,5.44260e-06_r8, &
      5.54814e-06_r8,5.65484e-06_r8,5.76272e-06_r8,5.87177e-06_r8,5.98199e-06_r8, &
      6.09339e-06_r8,6.20596e-06_r8,6.31969e-06_r8,6.43460e-06_r8,6.55068e-06_r8, &
      6.66793e-06_r8,6.78636e-06_r8,6.90595e-06_r8,7.02670e-06_r8,7.14863e-06_r8, &
      7.27173e-06_r8,7.39599e-06_r8,7.52142e-06_r8,7.64802e-06_r8,7.77577e-06_r8/)
      totplnk(101:150, 5) = (/ &
      7.90469e-06_r8,8.03477e-06_r8,8.16601e-06_r8,8.29841e-06_r8,8.43198e-06_r8, &
      8.56669e-06_r8,8.70256e-06_r8,8.83957e-06_r8,8.97775e-06_r8,9.11706e-06_r8, &
      9.25753e-06_r8,9.39915e-06_r8,9.54190e-06_r8,9.68580e-06_r8,9.83085e-06_r8, &
      9.97704e-06_r8,1.01243e-05_r8,1.02728e-05_r8,1.04224e-05_r8,1.05731e-05_r8, &
      1.07249e-05_r8,1.08779e-05_r8,1.10320e-05_r8,1.11872e-05_r8,1.13435e-05_r8, &
      1.15009e-05_r8,1.16595e-05_r8,1.18191e-05_r8,1.19799e-05_r8,1.21418e-05_r8, &
      1.23048e-05_r8,1.24688e-05_r8,1.26340e-05_r8,1.28003e-05_r8,1.29676e-05_r8, &
      1.31361e-05_r8,1.33056e-05_r8,1.34762e-05_r8,1.36479e-05_r8,1.38207e-05_r8, &
      1.39945e-05_r8,1.41694e-05_r8,1.43454e-05_r8,1.45225e-05_r8,1.47006e-05_r8, &
      1.48797e-05_r8,1.50600e-05_r8,1.52413e-05_r8,1.54236e-05_r8,1.56070e-05_r8/)
      totplnk(151:181, 5) = (/ &
      1.57914e-05_r8,1.59768e-05_r8,1.61633e-05_r8,1.63509e-05_r8,1.65394e-05_r8, &
      1.67290e-05_r8,1.69197e-05_r8,1.71113e-05_r8,1.73040e-05_r8,1.74976e-05_r8, &
      1.76923e-05_r8,1.78880e-05_r8,1.80847e-05_r8,1.82824e-05_r8,1.84811e-05_r8, &
      1.86808e-05_r8,1.88814e-05_r8,1.90831e-05_r8,1.92857e-05_r8,1.94894e-05_r8, &
      1.96940e-05_r8,1.98996e-05_r8,2.01061e-05_r8,2.03136e-05_r8,2.05221e-05_r8, &
      2.07316e-05_r8,2.09420e-05_r8,2.11533e-05_r8,2.13657e-05_r8,2.15789e-05_r8, &
      2.17931e-05_r8/)
      totplnk(1:50, 6) = (/ &
      2.73493e-07_r8,2.87408e-07_r8,3.01848e-07_r8,3.16825e-07_r8,3.32352e-07_r8, &
      3.48439e-07_r8,3.65100e-07_r8,3.82346e-07_r8,4.00189e-07_r8,4.18641e-07_r8, &
      4.37715e-07_r8,4.57422e-07_r8,4.77774e-07_r8,4.98784e-07_r8,5.20464e-07_r8, &
      5.42824e-07_r8,5.65879e-07_r8,5.89638e-07_r8,6.14115e-07_r8,6.39320e-07_r8, &
      6.65266e-07_r8,6.91965e-07_r8,7.19427e-07_r8,7.47666e-07_r8,7.76691e-07_r8, &
      8.06516e-07_r8,8.37151e-07_r8,8.68607e-07_r8,9.00896e-07_r8,9.34029e-07_r8, &
      9.68018e-07_r8,1.00287e-06_r8,1.03860e-06_r8,1.07522e-06_r8,1.11274e-06_r8, &
      1.15117e-06_r8,1.19052e-06_r8,1.23079e-06_r8,1.27201e-06_r8,1.31418e-06_r8, &
      1.35731e-06_r8,1.40141e-06_r8,1.44650e-06_r8,1.49257e-06_r8,1.53965e-06_r8, &
      1.58773e-06_r8,1.63684e-06_r8,1.68697e-06_r8,1.73815e-06_r8,1.79037e-06_r8/)
      totplnk(51:100, 6) = (/ &
      1.84365e-06_r8,1.89799e-06_r8,1.95341e-06_r8,2.00991e-06_r8,2.06750e-06_r8, &
      2.12619e-06_r8,2.18599e-06_r8,2.24691e-06_r8,2.30895e-06_r8,2.37212e-06_r8, &
      2.43643e-06_r8,2.50189e-06_r8,2.56851e-06_r8,2.63628e-06_r8,2.70523e-06_r8, &
      2.77536e-06_r8,2.84666e-06_r8,2.91916e-06_r8,2.99286e-06_r8,3.06776e-06_r8, &
      3.14387e-06_r8,3.22120e-06_r8,3.29975e-06_r8,3.37953e-06_r8,3.46054e-06_r8, &
      3.54280e-06_r8,3.62630e-06_r8,3.71105e-06_r8,3.79707e-06_r8,3.88434e-06_r8, &
      3.97288e-06_r8,4.06270e-06_r8,4.15380e-06_r8,4.24617e-06_r8,4.33984e-06_r8, &
      4.43479e-06_r8,4.53104e-06_r8,4.62860e-06_r8,4.72746e-06_r8,4.82763e-06_r8, &
      4.92911e-06_r8,5.03191e-06_r8,5.13603e-06_r8,5.24147e-06_r8,5.34824e-06_r8, &
      5.45634e-06_r8,5.56578e-06_r8,5.67656e-06_r8,5.78867e-06_r8,5.90213e-06_r8/)
      totplnk(101:150, 6) = (/ &
      6.01694e-06_r8,6.13309e-06_r8,6.25060e-06_r8,6.36947e-06_r8,6.48968e-06_r8, &
      6.61126e-06_r8,6.73420e-06_r8,6.85850e-06_r8,6.98417e-06_r8,7.11120e-06_r8, &
      7.23961e-06_r8,7.36938e-06_r8,7.50053e-06_r8,7.63305e-06_r8,7.76694e-06_r8, &
      7.90221e-06_r8,8.03887e-06_r8,8.17690e-06_r8,8.31632e-06_r8,8.45710e-06_r8, &
      8.59928e-06_r8,8.74282e-06_r8,8.88776e-06_r8,9.03409e-06_r8,9.18179e-06_r8, &
      9.33088e-06_r8,9.48136e-06_r8,9.63323e-06_r8,9.78648e-06_r8,9.94111e-06_r8, &
      1.00971e-05_r8,1.02545e-05_r8,1.04133e-05_r8,1.05735e-05_r8,1.07351e-05_r8, &
      1.08980e-05_r8,1.10624e-05_r8,1.12281e-05_r8,1.13952e-05_r8,1.15637e-05_r8, &
      1.17335e-05_r8,1.19048e-05_r8,1.20774e-05_r8,1.22514e-05_r8,1.24268e-05_r8, &
      1.26036e-05_r8,1.27817e-05_r8,1.29612e-05_r8,1.31421e-05_r8,1.33244e-05_r8/)
      totplnk(151:181, 6) = (/ &
      1.35080e-05_r8,1.36930e-05_r8,1.38794e-05_r8,1.40672e-05_r8,1.42563e-05_r8, &
      1.44468e-05_r8,1.46386e-05_r8,1.48318e-05_r8,1.50264e-05_r8,1.52223e-05_r8, &
      1.54196e-05_r8,1.56182e-05_r8,1.58182e-05_r8,1.60196e-05_r8,1.62223e-05_r8, &
      1.64263e-05_r8,1.66317e-05_r8,1.68384e-05_r8,1.70465e-05_r8,1.72559e-05_r8, &
      1.74666e-05_r8,1.76787e-05_r8,1.78921e-05_r8,1.81069e-05_r8,1.83230e-05_r8, &
      1.85404e-05_r8,1.87591e-05_r8,1.89791e-05_r8,1.92005e-05_r8,1.94232e-05_r8, &
      1.96471e-05_r8/)
      totplnk(1:50, 7) = (/ &
      1.25349e-07_r8,1.32735e-07_r8,1.40458e-07_r8,1.48527e-07_r8,1.56954e-07_r8, &
      1.65748e-07_r8,1.74920e-07_r8,1.84481e-07_r8,1.94443e-07_r8,2.04814e-07_r8, &
      2.15608e-07_r8,2.26835e-07_r8,2.38507e-07_r8,2.50634e-07_r8,2.63229e-07_r8, &
      2.76301e-07_r8,2.89864e-07_r8,3.03930e-07_r8,3.18508e-07_r8,3.33612e-07_r8, &
      3.49253e-07_r8,3.65443e-07_r8,3.82195e-07_r8,3.99519e-07_r8,4.17428e-07_r8, &
      4.35934e-07_r8,4.55050e-07_r8,4.74785e-07_r8,4.95155e-07_r8,5.16170e-07_r8, &
      5.37844e-07_r8,5.60186e-07_r8,5.83211e-07_r8,6.06929e-07_r8,6.31355e-07_r8, &
      6.56498e-07_r8,6.82373e-07_r8,7.08990e-07_r8,7.36362e-07_r8,7.64501e-07_r8, &
      7.93420e-07_r8,8.23130e-07_r8,8.53643e-07_r8,8.84971e-07_r8,9.17128e-07_r8, &
      9.50123e-07_r8,9.83969e-07_r8,1.01868e-06_r8,1.05426e-06_r8,1.09073e-06_r8/)
      totplnk(51:100, 7) = (/ &
      1.12810e-06_r8,1.16638e-06_r8,1.20558e-06_r8,1.24572e-06_r8,1.28680e-06_r8, &
      1.32883e-06_r8,1.37183e-06_r8,1.41581e-06_r8,1.46078e-06_r8,1.50675e-06_r8, &
      1.55374e-06_r8,1.60174e-06_r8,1.65078e-06_r8,1.70087e-06_r8,1.75200e-06_r8, &
      1.80421e-06_r8,1.85749e-06_r8,1.91186e-06_r8,1.96732e-06_r8,2.02389e-06_r8, &
      2.08159e-06_r8,2.14040e-06_r8,2.20035e-06_r8,2.26146e-06_r8,2.32372e-06_r8, &
      2.38714e-06_r8,2.45174e-06_r8,2.51753e-06_r8,2.58451e-06_r8,2.65270e-06_r8, &
      2.72210e-06_r8,2.79272e-06_r8,2.86457e-06_r8,2.93767e-06_r8,3.01201e-06_r8, &
      3.08761e-06_r8,3.16448e-06_r8,3.24261e-06_r8,3.32204e-06_r8,3.40275e-06_r8, &
      3.48476e-06_r8,3.56808e-06_r8,3.65271e-06_r8,3.73866e-06_r8,3.82595e-06_r8, &
      3.91456e-06_r8,4.00453e-06_r8,4.09584e-06_r8,4.18851e-06_r8,4.28254e-06_r8/)
      totplnk(101:150, 7) = (/ &
      4.37796e-06_r8,4.47475e-06_r8,4.57293e-06_r8,4.67249e-06_r8,4.77346e-06_r8, &
      4.87583e-06_r8,4.97961e-06_r8,5.08481e-06_r8,5.19143e-06_r8,5.29948e-06_r8, &
      5.40896e-06_r8,5.51989e-06_r8,5.63226e-06_r8,5.74608e-06_r8,5.86136e-06_r8, &
      5.97810e-06_r8,6.09631e-06_r8,6.21597e-06_r8,6.33713e-06_r8,6.45976e-06_r8, &
      6.58388e-06_r8,6.70950e-06_r8,6.83661e-06_r8,6.96521e-06_r8,7.09531e-06_r8, &
      7.22692e-06_r8,7.36005e-06_r8,7.49468e-06_r8,7.63084e-06_r8,7.76851e-06_r8, &
      7.90773e-06_r8,8.04846e-06_r8,8.19072e-06_r8,8.33452e-06_r8,8.47985e-06_r8, &
      8.62674e-06_r8,8.77517e-06_r8,8.92514e-06_r8,9.07666e-06_r8,9.22975e-06_r8, &
      9.38437e-06_r8,9.54057e-06_r8,9.69832e-06_r8,9.85762e-06_r8,1.00185e-05_r8, &
      1.01810e-05_r8,1.03450e-05_r8,1.05106e-05_r8,1.06777e-05_r8,1.08465e-05_r8/)
      totplnk(151:181, 7) = (/ &
      1.10168e-05_r8,1.11887e-05_r8,1.13621e-05_r8,1.15372e-05_r8,1.17138e-05_r8, &
      1.18920e-05_r8,1.20718e-05_r8,1.22532e-05_r8,1.24362e-05_r8,1.26207e-05_r8, &
      1.28069e-05_r8,1.29946e-05_r8,1.31839e-05_r8,1.33749e-05_r8,1.35674e-05_r8, &
      1.37615e-05_r8,1.39572e-05_r8,1.41544e-05_r8,1.43533e-05_r8,1.45538e-05_r8, &
      1.47558e-05_r8,1.49595e-05_r8,1.51647e-05_r8,1.53716e-05_r8,1.55800e-05_r8, &
      1.57900e-05_r8,1.60017e-05_r8,1.62149e-05_r8,1.64296e-05_r8,1.66460e-05_r8, &
      1.68640e-05_r8/)
      totplnk(1:50, 8) = (/ &
      6.74445e-08_r8,7.18176e-08_r8,7.64153e-08_r8,8.12456e-08_r8,8.63170e-08_r8, &
      9.16378e-08_r8,9.72168e-08_r8,1.03063e-07_r8,1.09184e-07_r8,1.15591e-07_r8, &
      1.22292e-07_r8,1.29296e-07_r8,1.36613e-07_r8,1.44253e-07_r8,1.52226e-07_r8, &
      1.60540e-07_r8,1.69207e-07_r8,1.78236e-07_r8,1.87637e-07_r8,1.97421e-07_r8, &
      2.07599e-07_r8,2.18181e-07_r8,2.29177e-07_r8,2.40598e-07_r8,2.52456e-07_r8, &
      2.64761e-07_r8,2.77523e-07_r8,2.90755e-07_r8,3.04468e-07_r8,3.18673e-07_r8, &
      3.33381e-07_r8,3.48603e-07_r8,3.64352e-07_r8,3.80638e-07_r8,3.97474e-07_r8, &
      4.14871e-07_r8,4.32841e-07_r8,4.51395e-07_r8,4.70547e-07_r8,4.90306e-07_r8, &
      5.10687e-07_r8,5.31699e-07_r8,5.53357e-07_r8,5.75670e-07_r8,5.98652e-07_r8, &
      6.22315e-07_r8,6.46672e-07_r8,6.71731e-07_r8,6.97511e-07_r8,7.24018e-07_r8/)
      totplnk(51:100, 8) = (/ &
      7.51266e-07_r8,7.79269e-07_r8,8.08038e-07_r8,8.37584e-07_r8,8.67922e-07_r8, &
      8.99061e-07_r8,9.31016e-07_r8,9.63797e-07_r8,9.97417e-07_r8,1.03189e-06_r8, &
      1.06722e-06_r8,1.10343e-06_r8,1.14053e-06_r8,1.17853e-06_r8,1.21743e-06_r8, &
      1.25726e-06_r8,1.29803e-06_r8,1.33974e-06_r8,1.38241e-06_r8,1.42606e-06_r8, &
      1.47068e-06_r8,1.51630e-06_r8,1.56293e-06_r8,1.61056e-06_r8,1.65924e-06_r8, &
      1.70894e-06_r8,1.75971e-06_r8,1.81153e-06_r8,1.86443e-06_r8,1.91841e-06_r8, &
      1.97350e-06_r8,2.02968e-06_r8,2.08699e-06_r8,2.14543e-06_r8,2.20500e-06_r8, &
      2.26573e-06_r8,2.32762e-06_r8,2.39068e-06_r8,2.45492e-06_r8,2.52036e-06_r8, &
      2.58700e-06_r8,2.65485e-06_r8,2.72393e-06_r8,2.79424e-06_r8,2.86580e-06_r8, &
      2.93861e-06_r8,3.01269e-06_r8,3.08803e-06_r8,3.16467e-06_r8,3.24259e-06_r8/)
      totplnk(101:150, 8) = (/ &
      3.32181e-06_r8,3.40235e-06_r8,3.48420e-06_r8,3.56739e-06_r8,3.65192e-06_r8, &
      3.73779e-06_r8,3.82502e-06_r8,3.91362e-06_r8,4.00359e-06_r8,4.09494e-06_r8, &
      4.18768e-06_r8,4.28182e-06_r8,4.37737e-06_r8,4.47434e-06_r8,4.57273e-06_r8, &
      4.67254e-06_r8,4.77380e-06_r8,4.87651e-06_r8,4.98067e-06_r8,5.08630e-06_r8, &
      5.19339e-06_r8,5.30196e-06_r8,5.41201e-06_r8,5.52356e-06_r8,5.63660e-06_r8, &
      5.75116e-06_r8,5.86722e-06_r8,5.98479e-06_r8,6.10390e-06_r8,6.22453e-06_r8, &
      6.34669e-06_r8,6.47042e-06_r8,6.59569e-06_r8,6.72252e-06_r8,6.85090e-06_r8, &
      6.98085e-06_r8,7.11238e-06_r8,7.24549e-06_r8,7.38019e-06_r8,7.51646e-06_r8, &
      7.65434e-06_r8,7.79382e-06_r8,7.93490e-06_r8,8.07760e-06_r8,8.22192e-06_r8, &
      8.36784e-06_r8,8.51540e-06_r8,8.66459e-06_r8,8.81542e-06_r8,8.96786e-06_r8/)
      totplnk(151:181, 8) = (/ &
      9.12197e-06_r8,9.27772e-06_r8,9.43513e-06_r8,9.59419e-06_r8,9.75490e-06_r8, &
      9.91728e-06_r8,1.00813e-05_r8,1.02471e-05_r8,1.04144e-05_r8,1.05835e-05_r8, &
      1.07543e-05_r8,1.09267e-05_r8,1.11008e-05_r8,1.12766e-05_r8,1.14541e-05_r8, &
      1.16333e-05_r8,1.18142e-05_r8,1.19969e-05_r8,1.21812e-05_r8,1.23672e-05_r8, &
      1.25549e-05_r8,1.27443e-05_r8,1.29355e-05_r8,1.31284e-05_r8,1.33229e-05_r8, &
      1.35193e-05_r8,1.37173e-05_r8,1.39170e-05_r8,1.41185e-05_r8,1.43217e-05_r8, &
      1.45267e-05_r8/)
      totplnk(1:50, 9) = (/ &
      2.61522e-08_r8,2.80613e-08_r8,3.00838e-08_r8,3.22250e-08_r8,3.44899e-08_r8, &
      3.68841e-08_r8,3.94129e-08_r8,4.20820e-08_r8,4.48973e-08_r8,4.78646e-08_r8, &
      5.09901e-08_r8,5.42799e-08_r8,5.77405e-08_r8,6.13784e-08_r8,6.52001e-08_r8, &
      6.92126e-08_r8,7.34227e-08_r8,7.78375e-08_r8,8.24643e-08_r8,8.73103e-08_r8, &
      9.23832e-08_r8,9.76905e-08_r8,1.03240e-07_r8,1.09039e-07_r8,1.15097e-07_r8, &
      1.21421e-07_r8,1.28020e-07_r8,1.34902e-07_r8,1.42075e-07_r8,1.49548e-07_r8, &
      1.57331e-07_r8,1.65432e-07_r8,1.73860e-07_r8,1.82624e-07_r8,1.91734e-07_r8, &
      2.01198e-07_r8,2.11028e-07_r8,2.21231e-07_r8,2.31818e-07_r8,2.42799e-07_r8, &
      2.54184e-07_r8,2.65983e-07_r8,2.78205e-07_r8,2.90862e-07_r8,3.03963e-07_r8, &
      3.17519e-07_r8,3.31541e-07_r8,3.46039e-07_r8,3.61024e-07_r8,3.76507e-07_r8/)
      totplnk(51:100, 9) = (/ &
      3.92498e-07_r8,4.09008e-07_r8,4.26050e-07_r8,4.43633e-07_r8,4.61769e-07_r8, &
      4.80469e-07_r8,4.99744e-07_r8,5.19606e-07_r8,5.40067e-07_r8,5.61136e-07_r8, &
      5.82828e-07_r8,6.05152e-07_r8,6.28120e-07_r8,6.51745e-07_r8,6.76038e-07_r8, &
      7.01010e-07_r8,7.26674e-07_r8,7.53041e-07_r8,7.80124e-07_r8,8.07933e-07_r8, &
      8.36482e-07_r8,8.65781e-07_r8,8.95845e-07_r8,9.26683e-07_r8,9.58308e-07_r8, &
      9.90732e-07_r8,1.02397e-06_r8,1.05803e-06_r8,1.09292e-06_r8,1.12866e-06_r8, &
      1.16526e-06_r8,1.20274e-06_r8,1.24109e-06_r8,1.28034e-06_r8,1.32050e-06_r8, &
      1.36158e-06_r8,1.40359e-06_r8,1.44655e-06_r8,1.49046e-06_r8,1.53534e-06_r8, &
      1.58120e-06_r8,1.62805e-06_r8,1.67591e-06_r8,1.72478e-06_r8,1.77468e-06_r8, &
      1.82561e-06_r8,1.87760e-06_r8,1.93066e-06_r8,1.98479e-06_r8,2.04000e-06_r8/)
      totplnk(101:150, 9) = (/ &
      2.09631e-06_r8,2.15373e-06_r8,2.21228e-06_r8,2.27196e-06_r8,2.33278e-06_r8, &
      2.39475e-06_r8,2.45790e-06_r8,2.52222e-06_r8,2.58773e-06_r8,2.65445e-06_r8, &
      2.72238e-06_r8,2.79152e-06_r8,2.86191e-06_r8,2.93354e-06_r8,3.00643e-06_r8, &
      3.08058e-06_r8,3.15601e-06_r8,3.23273e-06_r8,3.31075e-06_r8,3.39009e-06_r8, &
      3.47074e-06_r8,3.55272e-06_r8,3.63605e-06_r8,3.72072e-06_r8,3.80676e-06_r8, &
      3.89417e-06_r8,3.98297e-06_r8,4.07315e-06_r8,4.16474e-06_r8,4.25774e-06_r8, &
      4.35217e-06_r8,4.44802e-06_r8,4.54532e-06_r8,4.64406e-06_r8,4.74428e-06_r8, &
      4.84595e-06_r8,4.94911e-06_r8,5.05376e-06_r8,5.15990e-06_r8,5.26755e-06_r8, &
      5.37671e-06_r8,5.48741e-06_r8,5.59963e-06_r8,5.71340e-06_r8,5.82871e-06_r8, &
      5.94559e-06_r8,6.06403e-06_r8,6.18404e-06_r8,6.30565e-06_r8,6.42885e-06_r8/)
      totplnk(151:181, 9) = (/ &
      6.55364e-06_r8,6.68004e-06_r8,6.80806e-06_r8,6.93771e-06_r8,7.06898e-06_r8, &
      7.20190e-06_r8,7.33646e-06_r8,7.47267e-06_r8,7.61056e-06_r8,7.75010e-06_r8, &
      7.89133e-06_r8,8.03423e-06_r8,8.17884e-06_r8,8.32514e-06_r8,8.47314e-06_r8, &
      8.62284e-06_r8,8.77427e-06_r8,8.92743e-06_r8,9.08231e-06_r8,9.23893e-06_r8, &
      9.39729e-06_r8,9.55741e-06_r8,9.71927e-06_r8,9.88291e-06_r8,1.00483e-05_r8, &
      1.02155e-05_r8,1.03844e-05_r8,1.05552e-05_r8,1.07277e-05_r8,1.09020e-05_r8, &
      1.10781e-05_r8/)
      totplnk(1:50,10) = (/ &
      8.89300e-09_r8,9.63263e-09_r8,1.04235e-08_r8,1.12685e-08_r8,1.21703e-08_r8, &
      1.31321e-08_r8,1.41570e-08_r8,1.52482e-08_r8,1.64090e-08_r8,1.76428e-08_r8, &
      1.89533e-08_r8,2.03441e-08_r8,2.18190e-08_r8,2.33820e-08_r8,2.50370e-08_r8, &
      2.67884e-08_r8,2.86402e-08_r8,3.05969e-08_r8,3.26632e-08_r8,3.48436e-08_r8, &
      3.71429e-08_r8,3.95660e-08_r8,4.21179e-08_r8,4.48040e-08_r8,4.76294e-08_r8, &
      5.05996e-08_r8,5.37201e-08_r8,5.69966e-08_r8,6.04349e-08_r8,6.40411e-08_r8, &
      6.78211e-08_r8,7.17812e-08_r8,7.59276e-08_r8,8.02670e-08_r8,8.48059e-08_r8, &
      8.95508e-08_r8,9.45090e-08_r8,9.96873e-08_r8,1.05093e-07_r8,1.10733e-07_r8, &
      1.16614e-07_r8,1.22745e-07_r8,1.29133e-07_r8,1.35786e-07_r8,1.42711e-07_r8, &
      1.49916e-07_r8,1.57410e-07_r8,1.65202e-07_r8,1.73298e-07_r8,1.81709e-07_r8/)
      totplnk(51:100,10) = (/ &
      1.90441e-07_r8,1.99505e-07_r8,2.08908e-07_r8,2.18660e-07_r8,2.28770e-07_r8, &
      2.39247e-07_r8,2.50101e-07_r8,2.61340e-07_r8,2.72974e-07_r8,2.85013e-07_r8, &
      2.97467e-07_r8,3.10345e-07_r8,3.23657e-07_r8,3.37413e-07_r8,3.51623e-07_r8, &
      3.66298e-07_r8,3.81448e-07_r8,3.97082e-07_r8,4.13212e-07_r8,4.29848e-07_r8, &
      4.47000e-07_r8,4.64680e-07_r8,4.82898e-07_r8,5.01664e-07_r8,5.20991e-07_r8, &
      5.40888e-07_r8,5.61369e-07_r8,5.82440e-07_r8,6.04118e-07_r8,6.26410e-07_r8, &
      6.49329e-07_r8,6.72887e-07_r8,6.97095e-07_r8,7.21964e-07_r8,7.47506e-07_r8, &
      7.73732e-07_r8,8.00655e-07_r8,8.28287e-07_r8,8.56635e-07_r8,8.85717e-07_r8, &
      9.15542e-07_r8,9.46122e-07_r8,9.77469e-07_r8,1.00960e-06_r8,1.04251e-06_r8, &
      1.07623e-06_r8,1.11077e-06_r8,1.14613e-06_r8,1.18233e-06_r8,1.21939e-06_r8/)
      totplnk(101:150,10) = (/ &
      1.25730e-06_r8,1.29610e-06_r8,1.33578e-06_r8,1.37636e-06_r8,1.41785e-06_r8, &
      1.46027e-06_r8,1.50362e-06_r8,1.54792e-06_r8,1.59319e-06_r8,1.63942e-06_r8, &
      1.68665e-06_r8,1.73487e-06_r8,1.78410e-06_r8,1.83435e-06_r8,1.88564e-06_r8, &
      1.93797e-06_r8,1.99136e-06_r8,2.04582e-06_r8,2.10137e-06_r8,2.15801e-06_r8, &
      2.21576e-06_r8,2.27463e-06_r8,2.33462e-06_r8,2.39577e-06_r8,2.45806e-06_r8, &
      2.52153e-06_r8,2.58617e-06_r8,2.65201e-06_r8,2.71905e-06_r8,2.78730e-06_r8, &
      2.85678e-06_r8,2.92749e-06_r8,2.99946e-06_r8,3.07269e-06_r8,3.14720e-06_r8, &
      3.22299e-06_r8,3.30007e-06_r8,3.37847e-06_r8,3.45818e-06_r8,3.53923e-06_r8, &
      3.62161e-06_r8,3.70535e-06_r8,3.79046e-06_r8,3.87695e-06_r8,3.96481e-06_r8, &
      4.05409e-06_r8,4.14477e-06_r8,4.23687e-06_r8,4.33040e-06_r8,4.42538e-06_r8/)
      totplnk(151:181,10) = (/ &
      4.52180e-06_r8,4.61969e-06_r8,4.71905e-06_r8,4.81991e-06_r8,4.92226e-06_r8, &
      5.02611e-06_r8,5.13148e-06_r8,5.23839e-06_r8,5.34681e-06_r8,5.45681e-06_r8, &
      5.56835e-06_r8,5.68146e-06_r8,5.79614e-06_r8,5.91242e-06_r8,6.03030e-06_r8, &
      6.14978e-06_r8,6.27088e-06_r8,6.39360e-06_r8,6.51798e-06_r8,6.64398e-06_r8, &
      6.77165e-06_r8,6.90099e-06_r8,7.03198e-06_r8,7.16468e-06_r8,7.29906e-06_r8, &
      7.43514e-06_r8,7.57294e-06_r8,7.71244e-06_r8,7.85369e-06_r8,7.99666e-06_r8, &
      8.14138e-06_r8/)
      totplnk(1:50,11) = (/ &
      2.53767e-09_r8,2.77242e-09_r8,3.02564e-09_r8,3.29851e-09_r8,3.59228e-09_r8, &
      3.90825e-09_r8,4.24777e-09_r8,4.61227e-09_r8,5.00322e-09_r8,5.42219e-09_r8, &
      5.87080e-09_r8,6.35072e-09_r8,6.86370e-09_r8,7.41159e-09_r8,7.99628e-09_r8, &
      8.61974e-09_r8,9.28404e-09_r8,9.99130e-09_r8,1.07437e-08_r8,1.15436e-08_r8, &
      1.23933e-08_r8,1.32953e-08_r8,1.42522e-08_r8,1.52665e-08_r8,1.63410e-08_r8, &
      1.74786e-08_r8,1.86820e-08_r8,1.99542e-08_r8,2.12985e-08_r8,2.27179e-08_r8, &
      2.42158e-08_r8,2.57954e-08_r8,2.74604e-08_r8,2.92141e-08_r8,3.10604e-08_r8, &
      3.30029e-08_r8,3.50457e-08_r8,3.71925e-08_r8,3.94476e-08_r8,4.18149e-08_r8, &
      4.42991e-08_r8,4.69043e-08_r8,4.96352e-08_r8,5.24961e-08_r8,5.54921e-08_r8, &
      5.86277e-08_r8,6.19081e-08_r8,6.53381e-08_r8,6.89231e-08_r8,7.26681e-08_r8/)
      totplnk(51:100,11) = (/ &
      7.65788e-08_r8,8.06604e-08_r8,8.49187e-08_r8,8.93591e-08_r8,9.39879e-08_r8, &
      9.88106e-08_r8,1.03834e-07_r8,1.09063e-07_r8,1.14504e-07_r8,1.20165e-07_r8, &
      1.26051e-07_r8,1.32169e-07_r8,1.38525e-07_r8,1.45128e-07_r8,1.51982e-07_r8, &
      1.59096e-07_r8,1.66477e-07_r8,1.74132e-07_r8,1.82068e-07_r8,1.90292e-07_r8, &
      1.98813e-07_r8,2.07638e-07_r8,2.16775e-07_r8,2.26231e-07_r8,2.36015e-07_r8, &
      2.46135e-07_r8,2.56599e-07_r8,2.67415e-07_r8,2.78592e-07_r8,2.90137e-07_r8, &
      3.02061e-07_r8,3.14371e-07_r8,3.27077e-07_r8,3.40186e-07_r8,3.53710e-07_r8, &
      3.67655e-07_r8,3.82031e-07_r8,3.96848e-07_r8,4.12116e-07_r8,4.27842e-07_r8, &
      4.44039e-07_r8,4.60713e-07_r8,4.77876e-07_r8,4.95537e-07_r8,5.13706e-07_r8, &
      5.32392e-07_r8,5.51608e-07_r8,5.71360e-07_r8,5.91662e-07_r8,6.12521e-07_r8/)
      totplnk(101:150,11) = (/ &
      6.33950e-07_r8,6.55958e-07_r8,6.78556e-07_r8,7.01753e-07_r8,7.25562e-07_r8, &
      7.49992e-07_r8,7.75055e-07_r8,8.00760e-07_r8,8.27120e-07_r8,8.54145e-07_r8, &
      8.81845e-07_r8,9.10233e-07_r8,9.39318e-07_r8,9.69113e-07_r8,9.99627e-07_r8, &
      1.03087e-06_r8,1.06286e-06_r8,1.09561e-06_r8,1.12912e-06_r8,1.16340e-06_r8, &
      1.19848e-06_r8,1.23435e-06_r8,1.27104e-06_r8,1.30855e-06_r8,1.34690e-06_r8, &
      1.38609e-06_r8,1.42614e-06_r8,1.46706e-06_r8,1.50886e-06_r8,1.55155e-06_r8, &
      1.59515e-06_r8,1.63967e-06_r8,1.68512e-06_r8,1.73150e-06_r8,1.77884e-06_r8, &
      1.82715e-06_r8,1.87643e-06_r8,1.92670e-06_r8,1.97797e-06_r8,2.03026e-06_r8, &
      2.08356e-06_r8,2.13791e-06_r8,2.19330e-06_r8,2.24975e-06_r8,2.30728e-06_r8, &
      2.36589e-06_r8,2.42560e-06_r8,2.48641e-06_r8,2.54835e-06_r8,2.61142e-06_r8/)
      totplnk(151:181,11) = (/ &
      2.67563e-06_r8,2.74100e-06_r8,2.80754e-06_r8,2.87526e-06_r8,2.94417e-06_r8, &
      3.01429e-06_r8,3.08562e-06_r8,3.15819e-06_r8,3.23199e-06_r8,3.30704e-06_r8, &
      3.38336e-06_r8,3.46096e-06_r8,3.53984e-06_r8,3.62002e-06_r8,3.70151e-06_r8, &
      3.78433e-06_r8,3.86848e-06_r8,3.95399e-06_r8,4.04084e-06_r8,4.12907e-06_r8, &
      4.21868e-06_r8,4.30968e-06_r8,4.40209e-06_r8,4.49592e-06_r8,4.59117e-06_r8, &
      4.68786e-06_r8,4.78600e-06_r8,4.88561e-06_r8,4.98669e-06_r8,5.08926e-06_r8, &
      5.19332e-06_r8/)
      totplnk(1:50,12) = (/ &
      2.73921e-10_r8,3.04500e-10_r8,3.38056e-10_r8,3.74835e-10_r8,4.15099e-10_r8, &
      4.59126e-10_r8,5.07214e-10_r8,5.59679e-10_r8,6.16857e-10_r8,6.79103e-10_r8, &
      7.46796e-10_r8,8.20335e-10_r8,9.00144e-10_r8,9.86671e-10_r8,1.08039e-09_r8, &
      1.18180e-09_r8,1.29142e-09_r8,1.40982e-09_r8,1.53757e-09_r8,1.67529e-09_r8, &
      1.82363e-09_r8,1.98327e-09_r8,2.15492e-09_r8,2.33932e-09_r8,2.53726e-09_r8, &
      2.74957e-09_r8,2.97710e-09_r8,3.22075e-09_r8,3.48145e-09_r8,3.76020e-09_r8, &
      4.05801e-09_r8,4.37595e-09_r8,4.71513e-09_r8,5.07672e-09_r8,5.46193e-09_r8, &
      5.87201e-09_r8,6.30827e-09_r8,6.77205e-09_r8,7.26480e-09_r8,7.78794e-09_r8, &
      8.34304e-09_r8,8.93163e-09_r8,9.55537e-09_r8,1.02159e-08_r8,1.09151e-08_r8, &
      1.16547e-08_r8,1.24365e-08_r8,1.32625e-08_r8,1.41348e-08_r8,1.50554e-08_r8/)
      totplnk(51:100,12) = (/ &
      1.60264e-08_r8,1.70500e-08_r8,1.81285e-08_r8,1.92642e-08_r8,2.04596e-08_r8, &
      2.17171e-08_r8,2.30394e-08_r8,2.44289e-08_r8,2.58885e-08_r8,2.74209e-08_r8, &
      2.90290e-08_r8,3.07157e-08_r8,3.24841e-08_r8,3.43371e-08_r8,3.62782e-08_r8, &
      3.83103e-08_r8,4.04371e-08_r8,4.26617e-08_r8,4.49878e-08_r8,4.74190e-08_r8, &
      4.99589e-08_r8,5.26113e-08_r8,5.53801e-08_r8,5.82692e-08_r8,6.12826e-08_r8, &
      6.44245e-08_r8,6.76991e-08_r8,7.11105e-08_r8,7.46634e-08_r8,7.83621e-08_r8, &
      8.22112e-08_r8,8.62154e-08_r8,9.03795e-08_r8,9.47081e-08_r8,9.92066e-08_r8, &
      1.03879e-07_r8,1.08732e-07_r8,1.13770e-07_r8,1.18998e-07_r8,1.24422e-07_r8, &
      1.30048e-07_r8,1.35880e-07_r8,1.41924e-07_r8,1.48187e-07_r8,1.54675e-07_r8, &
      1.61392e-07_r8,1.68346e-07_r8,1.75543e-07_r8,1.82988e-07_r8,1.90688e-07_r8/)
      totplnk(101:150,12) = (/ &
      1.98650e-07_r8,2.06880e-07_r8,2.15385e-07_r8,2.24172e-07_r8,2.33247e-07_r8, &
      2.42617e-07_r8,2.52289e-07_r8,2.62272e-07_r8,2.72571e-07_r8,2.83193e-07_r8, &
      2.94147e-07_r8,3.05440e-07_r8,3.17080e-07_r8,3.29074e-07_r8,3.41430e-07_r8, &
      3.54155e-07_r8,3.67259e-07_r8,3.80747e-07_r8,3.94631e-07_r8,4.08916e-07_r8, &
      4.23611e-07_r8,4.38725e-07_r8,4.54267e-07_r8,4.70245e-07_r8,4.86666e-07_r8, &
      5.03541e-07_r8,5.20879e-07_r8,5.38687e-07_r8,5.56975e-07_r8,5.75751e-07_r8, &
      5.95026e-07_r8,6.14808e-07_r8,6.35107e-07_r8,6.55932e-07_r8,6.77293e-07_r8, &
      6.99197e-07_r8,7.21656e-07_r8,7.44681e-07_r8,7.68278e-07_r8,7.92460e-07_r8, &
      8.17235e-07_r8,8.42614e-07_r8,8.68606e-07_r8,8.95223e-07_r8,9.22473e-07_r8, &
      9.50366e-07_r8,9.78915e-07_r8,1.00813e-06_r8,1.03802e-06_r8,1.06859e-06_r8/)
      totplnk(151:181,12) = (/ &
      1.09986e-06_r8,1.13184e-06_r8,1.16453e-06_r8,1.19796e-06_r8,1.23212e-06_r8, &
      1.26703e-06_r8,1.30270e-06_r8,1.33915e-06_r8,1.37637e-06_r8,1.41440e-06_r8, &
      1.45322e-06_r8,1.49286e-06_r8,1.53333e-06_r8,1.57464e-06_r8,1.61679e-06_r8, &
      1.65981e-06_r8,1.70370e-06_r8,1.74847e-06_r8,1.79414e-06_r8,1.84071e-06_r8, &
      1.88821e-06_r8,1.93663e-06_r8,1.98599e-06_r8,2.03631e-06_r8,2.08759e-06_r8, &
      2.13985e-06_r8,2.19310e-06_r8,2.24734e-06_r8,2.30260e-06_r8,2.35888e-06_r8, &
      2.41619e-06_r8/)
      totplnk(1:50,13) = (/ &
      4.53634e-11_r8,5.11435e-11_r8,5.75754e-11_r8,6.47222e-11_r8,7.26531e-11_r8, &
      8.14420e-11_r8,9.11690e-11_r8,1.01921e-10_r8,1.13790e-10_r8,1.26877e-10_r8, &
      1.41288e-10_r8,1.57140e-10_r8,1.74555e-10_r8,1.93665e-10_r8,2.14613e-10_r8, &
      2.37548e-10_r8,2.62633e-10_r8,2.90039e-10_r8,3.19948e-10_r8,3.52558e-10_r8, &
      3.88073e-10_r8,4.26716e-10_r8,4.68719e-10_r8,5.14331e-10_r8,5.63815e-10_r8, &
      6.17448e-10_r8,6.75526e-10_r8,7.38358e-10_r8,8.06277e-10_r8,8.79625e-10_r8, &
      9.58770e-10_r8,1.04410e-09_r8,1.13602e-09_r8,1.23495e-09_r8,1.34135e-09_r8, &
      1.45568e-09_r8,1.57845e-09_r8,1.71017e-09_r8,1.85139e-09_r8,2.00268e-09_r8, &
      2.16464e-09_r8,2.33789e-09_r8,2.52309e-09_r8,2.72093e-09_r8,2.93212e-09_r8, &
      3.15740e-09_r8,3.39757e-09_r8,3.65341e-09_r8,3.92579e-09_r8,4.21559e-09_r8/)
      totplnk(51:100,13) = (/ &
      4.52372e-09_r8,4.85115e-09_r8,5.19886e-09_r8,5.56788e-09_r8,5.95928e-09_r8, &
      6.37419e-09_r8,6.81375e-09_r8,7.27917e-09_r8,7.77168e-09_r8,8.29256e-09_r8, &
      8.84317e-09_r8,9.42487e-09_r8,1.00391e-08_r8,1.06873e-08_r8,1.13710e-08_r8, &
      1.20919e-08_r8,1.28515e-08_r8,1.36514e-08_r8,1.44935e-08_r8,1.53796e-08_r8, &
      1.63114e-08_r8,1.72909e-08_r8,1.83201e-08_r8,1.94008e-08_r8,2.05354e-08_r8, &
      2.17258e-08_r8,2.29742e-08_r8,2.42830e-08_r8,2.56545e-08_r8,2.70910e-08_r8, &
      2.85950e-08_r8,3.01689e-08_r8,3.18155e-08_r8,3.35373e-08_r8,3.53372e-08_r8, &
      3.72177e-08_r8,3.91818e-08_r8,4.12325e-08_r8,4.33727e-08_r8,4.56056e-08_r8, &
      4.79342e-08_r8,5.03617e-08_r8,5.28915e-08_r8,5.55270e-08_r8,5.82715e-08_r8, &
      6.11286e-08_r8,6.41019e-08_r8,6.71951e-08_r8,7.04119e-08_r8,7.37560e-08_r8/)
      totplnk(101:150,13) = (/ &
      7.72315e-08_r8,8.08424e-08_r8,8.45927e-08_r8,8.84866e-08_r8,9.25281e-08_r8, &
      9.67218e-08_r8,1.01072e-07_r8,1.05583e-07_r8,1.10260e-07_r8,1.15107e-07_r8, &
      1.20128e-07_r8,1.25330e-07_r8,1.30716e-07_r8,1.36291e-07_r8,1.42061e-07_r8, &
      1.48031e-07_r8,1.54206e-07_r8,1.60592e-07_r8,1.67192e-07_r8,1.74015e-07_r8, &
      1.81064e-07_r8,1.88345e-07_r8,1.95865e-07_r8,2.03628e-07_r8,2.11643e-07_r8, &
      2.19912e-07_r8,2.28443e-07_r8,2.37244e-07_r8,2.46318e-07_r8,2.55673e-07_r8, &
      2.65316e-07_r8,2.75252e-07_r8,2.85489e-07_r8,2.96033e-07_r8,3.06891e-07_r8, &
      3.18070e-07_r8,3.29576e-07_r8,3.41417e-07_r8,3.53600e-07_r8,3.66133e-07_r8, &
      3.79021e-07_r8,3.92274e-07_r8,4.05897e-07_r8,4.19899e-07_r8,4.34288e-07_r8, &
      4.49071e-07_r8,4.64255e-07_r8,4.79850e-07_r8,4.95863e-07_r8,5.12300e-07_r8/)
      totplnk(151:181,13) = (/ &
      5.29172e-07_r8,5.46486e-07_r8,5.64250e-07_r8,5.82473e-07_r8,6.01164e-07_r8, &
      6.20329e-07_r8,6.39979e-07_r8,6.60122e-07_r8,6.80767e-07_r8,7.01922e-07_r8, &
      7.23596e-07_r8,7.45800e-07_r8,7.68539e-07_r8,7.91826e-07_r8,8.15669e-07_r8, &
      8.40076e-07_r8,8.65058e-07_r8,8.90623e-07_r8,9.16783e-07_r8,9.43544e-07_r8, &
      9.70917e-07_r8,9.98912e-07_r8,1.02754e-06_r8,1.05681e-06_r8,1.08673e-06_r8, &
      1.11731e-06_r8,1.14856e-06_r8,1.18050e-06_r8,1.21312e-06_r8,1.24645e-06_r8, &
      1.28049e-06_r8/)
      totplnk(1:50,14) = (/ &
      1.40113e-11_r8,1.59358e-11_r8,1.80960e-11_r8,2.05171e-11_r8,2.32266e-11_r8, &
      2.62546e-11_r8,2.96335e-11_r8,3.33990e-11_r8,3.75896e-11_r8,4.22469e-11_r8, &
      4.74164e-11_r8,5.31466e-11_r8,5.94905e-11_r8,6.65054e-11_r8,7.42522e-11_r8, &
      8.27975e-11_r8,9.22122e-11_r8,1.02573e-10_r8,1.13961e-10_r8,1.26466e-10_r8, &
      1.40181e-10_r8,1.55206e-10_r8,1.71651e-10_r8,1.89630e-10_r8,2.09265e-10_r8, &
      2.30689e-10_r8,2.54040e-10_r8,2.79467e-10_r8,3.07128e-10_r8,3.37190e-10_r8, &
      3.69833e-10_r8,4.05243e-10_r8,4.43623e-10_r8,4.85183e-10_r8,5.30149e-10_r8, &
      5.78755e-10_r8,6.31255e-10_r8,6.87910e-10_r8,7.49002e-10_r8,8.14824e-10_r8, &
      8.85687e-10_r8,9.61914e-10_r8,1.04385e-09_r8,1.13186e-09_r8,1.22631e-09_r8, &
      1.32761e-09_r8,1.43617e-09_r8,1.55243e-09_r8,1.67686e-09_r8,1.80992e-09_r8/)
      totplnk(51:100,14) = (/ &
      1.95212e-09_r8,2.10399e-09_r8,2.26607e-09_r8,2.43895e-09_r8,2.62321e-09_r8, &
      2.81949e-09_r8,3.02844e-09_r8,3.25073e-09_r8,3.48707e-09_r8,3.73820e-09_r8, &
      4.00490e-09_r8,4.28794e-09_r8,4.58819e-09_r8,4.90647e-09_r8,5.24371e-09_r8, &
      5.60081e-09_r8,5.97875e-09_r8,6.37854e-09_r8,6.80120e-09_r8,7.24782e-09_r8, &
      7.71950e-09_r8,8.21740e-09_r8,8.74271e-09_r8,9.29666e-09_r8,9.88054e-09_r8, &
      1.04956e-08_r8,1.11434e-08_r8,1.18251e-08_r8,1.25422e-08_r8,1.32964e-08_r8, &
      1.40890e-08_r8,1.49217e-08_r8,1.57961e-08_r8,1.67140e-08_r8,1.76771e-08_r8, &
      1.86870e-08_r8,1.97458e-08_r8,2.08553e-08_r8,2.20175e-08_r8,2.32342e-08_r8, &
      2.45077e-08_r8,2.58401e-08_r8,2.72334e-08_r8,2.86900e-08_r8,3.02122e-08_r8, &
      3.18021e-08_r8,3.34624e-08_r8,3.51954e-08_r8,3.70037e-08_r8,3.88899e-08_r8/)
      totplnk(101:150,14) = (/ &
      4.08568e-08_r8,4.29068e-08_r8,4.50429e-08_r8,4.72678e-08_r8,4.95847e-08_r8, &
      5.19963e-08_r8,5.45058e-08_r8,5.71161e-08_r8,5.98309e-08_r8,6.26529e-08_r8, &
      6.55857e-08_r8,6.86327e-08_r8,7.17971e-08_r8,7.50829e-08_r8,7.84933e-08_r8, &
      8.20323e-08_r8,8.57035e-08_r8,8.95105e-08_r8,9.34579e-08_r8,9.75488e-08_r8, &
      1.01788e-07_r8,1.06179e-07_r8,1.10727e-07_r8,1.15434e-07_r8,1.20307e-07_r8, &
      1.25350e-07_r8,1.30566e-07_r8,1.35961e-07_r8,1.41539e-07_r8,1.47304e-07_r8, &
      1.53263e-07_r8,1.59419e-07_r8,1.65778e-07_r8,1.72345e-07_r8,1.79124e-07_r8, &
      1.86122e-07_r8,1.93343e-07_r8,2.00792e-07_r8,2.08476e-07_r8,2.16400e-07_r8, &
      2.24568e-07_r8,2.32988e-07_r8,2.41666e-07_r8,2.50605e-07_r8,2.59813e-07_r8, &
      2.69297e-07_r8,2.79060e-07_r8,2.89111e-07_r8,2.99455e-07_r8,3.10099e-07_r8/)
      totplnk(151:181,14) = (/ &
      3.21049e-07_r8,3.32311e-07_r8,3.43893e-07_r8,3.55801e-07_r8,3.68041e-07_r8, &
      3.80621e-07_r8,3.93547e-07_r8,4.06826e-07_r8,4.20465e-07_r8,4.34473e-07_r8, &
      4.48856e-07_r8,4.63620e-07_r8,4.78774e-07_r8,4.94325e-07_r8,5.10280e-07_r8, &
      5.26648e-07_r8,5.43436e-07_r8,5.60652e-07_r8,5.78302e-07_r8,5.96397e-07_r8, &
      6.14943e-07_r8,6.33949e-07_r8,6.53421e-07_r8,6.73370e-07_r8,6.93803e-07_r8, &
      7.14731e-07_r8,7.36157e-07_r8,7.58095e-07_r8,7.80549e-07_r8,8.03533e-07_r8, &
      8.27050e-07_r8/)
      totplnk(1:50,15) = (/ &
      3.90483e-12_r8,4.47999e-12_r8,5.13122e-12_r8,5.86739e-12_r8,6.69829e-12_r8, &
      7.63467e-12_r8,8.68833e-12_r8,9.87221e-12_r8,1.12005e-11_r8,1.26885e-11_r8, &
      1.43534e-11_r8,1.62134e-11_r8,1.82888e-11_r8,2.06012e-11_r8,2.31745e-11_r8, &
      2.60343e-11_r8,2.92087e-11_r8,3.27277e-11_r8,3.66242e-11_r8,4.09334e-11_r8, &
      4.56935e-11_r8,5.09455e-11_r8,5.67338e-11_r8,6.31057e-11_r8,7.01127e-11_r8, &
      7.78096e-11_r8,8.62554e-11_r8,9.55130e-11_r8,1.05651e-10_r8,1.16740e-10_r8, &
      1.28858e-10_r8,1.42089e-10_r8,1.56519e-10_r8,1.72243e-10_r8,1.89361e-10_r8, &
      2.07978e-10_r8,2.28209e-10_r8,2.50173e-10_r8,2.73999e-10_r8,2.99820e-10_r8, &
      3.27782e-10_r8,3.58034e-10_r8,3.90739e-10_r8,4.26067e-10_r8,4.64196e-10_r8, &
      5.05317e-10_r8,5.49631e-10_r8,5.97347e-10_r8,6.48689e-10_r8,7.03891e-10_r8/)
      totplnk(51:100,15) = (/ &
      7.63201e-10_r8,8.26876e-10_r8,8.95192e-10_r8,9.68430e-10_r8,1.04690e-09_r8, &
      1.13091e-09_r8,1.22079e-09_r8,1.31689e-09_r8,1.41957e-09_r8,1.52922e-09_r8, &
      1.64623e-09_r8,1.77101e-09_r8,1.90401e-09_r8,2.04567e-09_r8,2.19647e-09_r8, &
      2.35690e-09_r8,2.52749e-09_r8,2.70875e-09_r8,2.90127e-09_r8,3.10560e-09_r8, &
      3.32238e-09_r8,3.55222e-09_r8,3.79578e-09_r8,4.05375e-09_r8,4.32682e-09_r8, &
      4.61574e-09_r8,4.92128e-09_r8,5.24420e-09_r8,5.58536e-09_r8,5.94558e-09_r8, &
      6.32575e-09_r8,6.72678e-09_r8,7.14964e-09_r8,7.59526e-09_r8,8.06470e-09_r8, &
      8.55897e-09_r8,9.07916e-09_r8,9.62638e-09_r8,1.02018e-08_r8,1.08066e-08_r8, &
      1.14420e-08_r8,1.21092e-08_r8,1.28097e-08_r8,1.35446e-08_r8,1.43155e-08_r8, &
      1.51237e-08_r8,1.59708e-08_r8,1.68581e-08_r8,1.77873e-08_r8,1.87599e-08_r8/)
      totplnk(101:150,15) = (/ &
      1.97777e-08_r8,2.08423e-08_r8,2.19555e-08_r8,2.31190e-08_r8,2.43348e-08_r8, &
      2.56045e-08_r8,2.69302e-08_r8,2.83140e-08_r8,2.97578e-08_r8,3.12636e-08_r8, &
      3.28337e-08_r8,3.44702e-08_r8,3.61755e-08_r8,3.79516e-08_r8,3.98012e-08_r8, &
      4.17265e-08_r8,4.37300e-08_r8,4.58143e-08_r8,4.79819e-08_r8,5.02355e-08_r8, &
      5.25777e-08_r8,5.50114e-08_r8,5.75393e-08_r8,6.01644e-08_r8,6.28896e-08_r8, &
      6.57177e-08_r8,6.86521e-08_r8,7.16959e-08_r8,7.48520e-08_r8,7.81239e-08_r8, &
      8.15148e-08_r8,8.50282e-08_r8,8.86675e-08_r8,9.24362e-08_r8,9.63380e-08_r8, &
      1.00376e-07_r8,1.04555e-07_r8,1.08878e-07_r8,1.13349e-07_r8,1.17972e-07_r8, &
      1.22751e-07_r8,1.27690e-07_r8,1.32793e-07_r8,1.38064e-07_r8,1.43508e-07_r8, &
      1.49129e-07_r8,1.54931e-07_r8,1.60920e-07_r8,1.67099e-07_r8,1.73473e-07_r8/)
      totplnk(151:181,15) = (/ &
      1.80046e-07_r8,1.86825e-07_r8,1.93812e-07_r8,2.01014e-07_r8,2.08436e-07_r8, &
      2.16082e-07_r8,2.23957e-07_r8,2.32067e-07_r8,2.40418e-07_r8,2.49013e-07_r8, &
      2.57860e-07_r8,2.66963e-07_r8,2.76328e-07_r8,2.85961e-07_r8,2.95868e-07_r8, &
      3.06053e-07_r8,3.16524e-07_r8,3.27286e-07_r8,3.38345e-07_r8,3.49707e-07_r8, &
      3.61379e-07_r8,3.73367e-07_r8,3.85676e-07_r8,3.98315e-07_r8,4.11287e-07_r8, &
      4.24602e-07_r8,4.38265e-07_r8,4.52283e-07_r8,4.66662e-07_r8,4.81410e-07_r8, &
      4.96535e-07_r8/)
      totplnk(1:50,16) = (/ &
      0.28639e-12_r8,0.33349e-12_r8,0.38764e-12_r8,0.44977e-12_r8,0.52093e-12_r8, &
      0.60231e-12_r8,0.69522e-12_r8,0.80111e-12_r8,0.92163e-12_r8,0.10586e-11_r8, &
      0.12139e-11_r8,0.13899e-11_r8,0.15890e-11_r8,0.18138e-11_r8,0.20674e-11_r8, &
      0.23531e-11_r8,0.26744e-11_r8,0.30352e-11_r8,0.34401e-11_r8,0.38936e-11_r8, &
      0.44011e-11_r8,0.49681e-11_r8,0.56010e-11_r8,0.63065e-11_r8,0.70919e-11_r8, &
      0.79654e-11_r8,0.89357e-11_r8,0.10012e-10_r8,0.11205e-10_r8,0.12526e-10_r8, &
      0.13986e-10_r8,0.15600e-10_r8,0.17380e-10_r8,0.19342e-10_r8,0.21503e-10_r8, &
      0.23881e-10_r8,0.26494e-10_r8,0.29362e-10_r8,0.32509e-10_r8,0.35958e-10_r8, &
      0.39733e-10_r8,0.43863e-10_r8,0.48376e-10_r8,0.53303e-10_r8,0.58679e-10_r8, &
      0.64539e-10_r8,0.70920e-10_r8,0.77864e-10_r8,0.85413e-10_r8,0.93615e-10_r8/)
      totplnk(51:100,16) = (/ &
      0.10252e-09_r8,0.11217e-09_r8,0.12264e-09_r8,0.13397e-09_r8,0.14624e-09_r8, &
      0.15950e-09_r8,0.17383e-09_r8,0.18930e-09_r8,0.20599e-09_r8,0.22399e-09_r8, &
      0.24339e-09_r8,0.26427e-09_r8,0.28674e-09_r8,0.31090e-09_r8,0.33686e-09_r8, &
      0.36474e-09_r8,0.39466e-09_r8,0.42676e-09_r8,0.46115e-09_r8,0.49800e-09_r8, &
      0.53744e-09_r8,0.57964e-09_r8,0.62476e-09_r8,0.67298e-09_r8,0.72448e-09_r8, &
      0.77945e-09_r8,0.83809e-09_r8,0.90062e-09_r8,0.96725e-09_r8,0.10382e-08_r8, &
      0.11138e-08_r8,0.11941e-08_r8,0.12796e-08_r8,0.13704e-08_r8,0.14669e-08_r8, &
      0.15694e-08_r8,0.16781e-08_r8,0.17934e-08_r8,0.19157e-08_r8,0.20453e-08_r8, &
      0.21825e-08_r8,0.23278e-08_r8,0.24815e-08_r8,0.26442e-08_r8,0.28161e-08_r8, &
      0.29978e-08_r8,0.31898e-08_r8,0.33925e-08_r8,0.36064e-08_r8,0.38321e-08_r8/)
      totplnk(101:150,16) = (/ &
      0.40700e-08_r8,0.43209e-08_r8,0.45852e-08_r8,0.48636e-08_r8,0.51567e-08_r8, &
      0.54652e-08_r8,0.57897e-08_r8,0.61310e-08_r8,0.64897e-08_r8,0.68667e-08_r8, &
      0.72626e-08_r8,0.76784e-08_r8,0.81148e-08_r8,0.85727e-08_r8,0.90530e-08_r8, &
      0.95566e-08_r8,0.10084e-07_r8,0.10638e-07_r8,0.11217e-07_r8,0.11824e-07_r8, &
      0.12458e-07_r8,0.13123e-07_r8,0.13818e-07_r8,0.14545e-07_r8,0.15305e-07_r8, &
      0.16099e-07_r8,0.16928e-07_r8,0.17795e-07_r8,0.18699e-07_r8,0.19643e-07_r8, &
      0.20629e-07_r8,0.21656e-07_r8,0.22728e-07_r8,0.23845e-07_r8,0.25010e-07_r8, &
      0.26223e-07_r8,0.27487e-07_r8,0.28804e-07_r8,0.30174e-07_r8,0.31600e-07_r8, &
      0.33084e-07_r8,0.34628e-07_r8,0.36233e-07_r8,0.37902e-07_r8,0.39637e-07_r8, &
      0.41440e-07_r8,0.43313e-07_r8,0.45259e-07_r8,0.47279e-07_r8,0.49376e-07_r8/)
      totplnk(151:181,16) = (/ &
      0.51552e-07_r8,0.53810e-07_r8,0.56153e-07_r8,0.58583e-07_r8,0.61102e-07_r8, &
      0.63713e-07_r8,0.66420e-07_r8,0.69224e-07_r8,0.72129e-07_r8,0.75138e-07_r8, &
      0.78254e-07_r8,0.81479e-07_r8,0.84818e-07_r8,0.88272e-07_r8,0.91846e-07_r8, &
      0.95543e-07_r8,0.99366e-07_r8,0.10332e-06_r8,0.10740e-06_r8,0.11163e-06_r8, &
      0.11599e-06_r8,0.12050e-06_r8,0.12515e-06_r8,0.12996e-06_r8,0.13493e-06_r8, &
      0.14005e-06_r8,0.14534e-06_r8,0.15080e-06_r8,0.15643e-06_r8,0.16224e-06_r8, &
      0.16823e-06_r8/)
      totplk16(1:50) = (/ &
      0.28481e-12_r8,0.33159e-12_r8,0.38535e-12_r8,0.44701e-12_r8,0.51763e-12_r8, &
      0.59836e-12_r8,0.69049e-12_r8,0.79549e-12_r8,0.91493e-12_r8,0.10506e-11_r8, &
      0.12045e-11_r8,0.13788e-11_r8,0.15758e-11_r8,0.17984e-11_r8,0.20493e-11_r8, &
      0.23317e-11_r8,0.26494e-11_r8,0.30060e-11_r8,0.34060e-11_r8,0.38539e-11_r8, &
      0.43548e-11_r8,0.49144e-11_r8,0.55387e-11_r8,0.62344e-11_r8,0.70086e-11_r8, &
      0.78692e-11_r8,0.88248e-11_r8,0.98846e-11_r8,0.11059e-10_r8,0.12358e-10_r8, &
      0.13794e-10_r8,0.15379e-10_r8,0.17128e-10_r8,0.19055e-10_r8,0.21176e-10_r8, &
      0.23508e-10_r8,0.26070e-10_r8,0.28881e-10_r8,0.31963e-10_r8,0.35339e-10_r8, &
      0.39034e-10_r8,0.43073e-10_r8,0.47484e-10_r8,0.52299e-10_r8,0.57548e-10_r8, &
      0.63267e-10_r8,0.69491e-10_r8,0.76261e-10_r8,0.83616e-10_r8,0.91603e-10_r8/)
      totplk16(51:100) = (/ &
      0.10027e-09_r8,0.10966e-09_r8,0.11983e-09_r8,0.13084e-09_r8,0.14275e-09_r8, &
      0.15562e-09_r8,0.16951e-09_r8,0.18451e-09_r8,0.20068e-09_r8,0.21810e-09_r8, &
      0.23686e-09_r8,0.25704e-09_r8,0.27875e-09_r8,0.30207e-09_r8,0.32712e-09_r8, &
      0.35400e-09_r8,0.38282e-09_r8,0.41372e-09_r8,0.44681e-09_r8,0.48223e-09_r8, &
      0.52013e-09_r8,0.56064e-09_r8,0.60392e-09_r8,0.65015e-09_r8,0.69948e-09_r8, &
      0.75209e-09_r8,0.80818e-09_r8,0.86794e-09_r8,0.93157e-09_r8,0.99929e-09_r8, &
      0.10713e-08_r8,0.11479e-08_r8,0.12293e-08_r8,0.13157e-08_r8,0.14074e-08_r8, &
      0.15047e-08_r8,0.16079e-08_r8,0.17172e-08_r8,0.18330e-08_r8,0.19557e-08_r8, &
      0.20855e-08_r8,0.22228e-08_r8,0.23680e-08_r8,0.25214e-08_r8,0.26835e-08_r8, &
      0.28546e-08_r8,0.30352e-08_r8,0.32257e-08_r8,0.34266e-08_r8,0.36384e-08_r8/)
      totplk16(101:150) = (/ &
      0.38615e-08_r8,0.40965e-08_r8,0.43438e-08_r8,0.46041e-08_r8,0.48779e-08_r8, &
      0.51658e-08_r8,0.54683e-08_r8,0.57862e-08_r8,0.61200e-08_r8,0.64705e-08_r8, &
      0.68382e-08_r8,0.72240e-08_r8,0.76285e-08_r8,0.80526e-08_r8,0.84969e-08_r8, &
      0.89624e-08_r8,0.94498e-08_r8,0.99599e-08_r8,0.10494e-07_r8,0.11052e-07_r8, &
      0.11636e-07_r8,0.12246e-07_r8,0.12884e-07_r8,0.13551e-07_r8,0.14246e-07_r8, &
      0.14973e-07_r8,0.15731e-07_r8,0.16522e-07_r8,0.17347e-07_r8,0.18207e-07_r8, &
      0.19103e-07_r8,0.20037e-07_r8,0.21011e-07_r8,0.22024e-07_r8,0.23079e-07_r8, &
      0.24177e-07_r8,0.25320e-07_r8,0.26508e-07_r8,0.27744e-07_r8,0.29029e-07_r8, &
      0.30365e-07_r8,0.31753e-07_r8,0.33194e-07_r8,0.34691e-07_r8,0.36246e-07_r8, &
      0.37859e-07_r8,0.39533e-07_r8,0.41270e-07_r8,0.43071e-07_r8,0.44939e-07_r8/)
      totplk16(151:181) = (/ &
      0.46875e-07_r8,0.48882e-07_r8,0.50961e-07_r8,0.53115e-07_r8,0.55345e-07_r8, &
      0.57655e-07_r8,0.60046e-07_r8,0.62520e-07_r8,0.65080e-07_r8,0.67728e-07_r8, &
      0.70466e-07_r8,0.73298e-07_r8,0.76225e-07_r8,0.79251e-07_r8,0.82377e-07_r8, &
      0.85606e-07_r8,0.88942e-07_r8,0.92386e-07_r8,0.95942e-07_r8,0.99612e-07_r8, &
      0.10340e-06_r8,0.10731e-06_r8,0.11134e-06_r8,0.11550e-06_r8,0.11979e-06_r8, &
      0.12421e-06_r8,0.12876e-06_r8,0.13346e-06_r8,0.13830e-06_r8,0.14328e-06_r8, &
      0.14841e-06_r8/)

      end subroutine lwavplank

      end module rrtmg_lw_setcoef

