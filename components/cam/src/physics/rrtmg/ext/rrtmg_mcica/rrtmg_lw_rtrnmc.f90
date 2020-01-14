!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.3 $
!     created:   $Date: 2008/04/24 16:17:28 $
!
      module rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngb, ngs
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrrtc, hnamrtc

      implicit none

      contains

!-----------------------------------------------------------------------------
      subroutine rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                        cldfmc, taucmc, planklay, planklev, plankbnd, &
                        pwvcm, fracs, taut, &
                        totuflux, totdflux, fnet, htr, &
                        totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs ) 
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with the McICA stochastic approach and maximum-random
!  cloud overlap. 
!***************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: nlayers                    ! total number of layers
      integer, intent(in) :: istart                     ! beginning band of calculation
      integer, intent(in) :: iend                       ! ending band of calculation
      integer, intent(in) :: iout                       ! output option flag

! Atmosphere
      real(kind=r8), intent(in) :: pz(0:)               ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: pwvcm                ! precipitable water vapor (cm)
      real(kind=r8), intent(in) :: semiss(:)            ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: planklay(:,:)        ! 
                                                        !    Dimensions: (nlayers,nbndlw)
      real(kind=r8), intent(in) :: planklev(0:,:)       ! 
                                                        !    Dimensions: (0:nlayers,nbndlw)
      real(kind=r8), intent(in) :: plankbnd(:)          ! 
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: fracs(:,:)           ! 
                                                        !    Dimensions: (nlayers,ngptw)
      real(kind=r8), intent(in) :: taut(:,:)            ! gaseous + aerosol optical depths
                                                        !    Dimensions: (nlayers,ngptlw)

! Clouds
      integer, intent(in) :: ncbands                    ! number of cloud spectral bands
      real(kind=r8), intent(in) :: cldfmc(:,:)          ! layer cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: taucmc(:,:)          ! layer cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)

! ----- Output -----
      real(kind=r8), intent(out) :: totuflux(0:)        ! upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdflux(0:)        ! downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnet(0:)            ! net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htr(0:)             ! longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totuclfl(0:)        ! clear sky upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdclfl(0:)        ! clear sky downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnetc(0:)           ! clear sky net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htrc(0:)            ! clear sky longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totufluxs(:,0:)     ! upward longwave flux spectral (w/m2)
                                                        !    Dimensions: (nbndlw, 0:nlayers)
      real(kind=r8), intent(out) :: totdfluxs(:,0:)     ! downward longwave flux spectral (w/m2)
                                                        !    Dimensions: (nbndlw, 0:nlayers)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=r8) :: abscld(nlayers,ngptlw)
      real(kind=r8) :: atot(nlayers)
      real(kind=r8) :: atrans(nlayers)
      real(kind=r8) :: bbugas(nlayers)
      real(kind=r8) :: bbutot(nlayers)
      real(kind=r8) :: clrurad(0:nlayers)
      real(kind=r8) :: clrdrad(0:nlayers)
      real(kind=r8) :: efclfrac(nlayers,ngptlw)
      real(kind=r8) :: uflux(0:nlayers)
      real(kind=r8) :: dflux(0:nlayers)
      real(kind=r8) :: urad(0:nlayers)
      real(kind=r8) :: drad(0:nlayers)
      real(kind=r8) :: uclfl(0:nlayers)
      real(kind=r8) :: dclfl(0:nlayers)
      real(kind=r8) :: odcld(nlayers,ngptlw)


      real(kind=r8) :: secdiff(nbndlw)                   ! secant of diffusivity angle
      real(kind=r8) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)  ! diffusivity angle adjustment coefficients
      real(kind=r8) :: wtdiff, rec_6
      real(kind=r8) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=r8) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=r8) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=r8) :: rad0, reflect, radlu, radclru

      integer :: icldlyr(nlayers)                        ! flag for cloud in layer
      integer :: ibnd, ib, iband, lay, lev, l, ig        ! loop indices
      integer :: igc                                     ! g-point interval counter
      integer :: iclddn                                  ! flag for cloud in down path
      integer :: ittot, itgas, itr                       ! lookup table indices

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)


! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
      data wtdiff /0.5_r8/
      data rec_6 /0.166667_r8/

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
      data a0 / 1.66_r8,  1.55_r8,  1.58_r8,  1.66_r8, &
                1.54_r8, 1.454_r8,  1.89_r8,  1.33_r8, &
               1.668_r8,  1.66_r8,  1.66_r8,  1.66_r8, &
                1.66_r8,  1.66_r8,  1.66_r8,  1.66_r8 /
      data a1 / 0.00_r8,  0.25_r8,  0.22_r8,  0.00_r8, &
                0.13_r8, 0.446_r8, -0.10_r8,  0.40_r8, &
              -0.006_r8,  0.00_r8,  0.00_r8,  0.00_r8, &
                0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /
      data a2 / 0.00_r8, -12.0_r8, -11.7_r8,  0.00_r8, &
               -0.72_r8,-0.243_r8,  0.19_r8,-0.062_r8, &
               0.414_r8,  0.00_r8,  0.00_r8,  0.00_r8, &
                0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /

      hvrrtc = '$Revision: 1.3 $'

      do ibnd = 1,nbndlw
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_r8
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)
         endif
!!
!!  - fix the case where the diffusivity angle was going negative under very 
!!    moist conditions
!!
         if (secdiff(ibnd) .gt. 1.80_r8) secdiff(ibnd) = 1.80_r8
         if (secdiff(ibnd) .lt. 1.50_r8) secdiff(ibnd) = 1.50_r8

      enddo

!!    if (pwvcm.lt.1.0) secdiff(6) = 1.80_r8
!!    if (pwvcm.gt.7.1) secdiff(7) = 1.50_r8

      urad(0) = 0.0_r8
      drad(0) = 0.0_r8
      totuflux(0) = 0.0_r8
      totdflux(0) = 0.0_r8
      clrurad(0) = 0.0_r8
      clrdrad(0) = 0.0_r8
      totuclfl(0) = 0.0_r8
      totdclfl(0) = 0.0_r8

      do lay = 1, nlayers
         urad(lay) = 0.0_r8
         drad(lay) = 0.0_r8
         totuflux(lay) = 0.0_r8
         totdflux(lay) = 0.0_r8
         clrurad(lay) = 0.0_r8
         clrdrad(lay) = 0.0_r8
         totuclfl(lay) = 0.0_r8
         totdclfl(lay) = 0.0_r8
         icldlyr(lay) = 0

! Change to band loop?
         do ig = 1, ngptlw
            if (cldfmc(ig,lay) .eq. 1._r8) then
               ib = ngb(ig)
               odcld(lay,ig) = secdiff(ib) * taucmc(ig,lay)
               transcld = exp(-odcld(lay,ig))
               abscld(lay,ig) = 1._r8 - transcld
               efclfrac(lay,ig) = abscld(lay,ig) * cldfmc(ig,lay)
               icldlyr(lay) = 1
            else
               odcld(lay,ig) = 0.0_r8
               abscld(lay,ig) = 0.0_r8
               efclfrac(lay,ig) = 0.0_r8
            endif
         enddo

      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._r8
         radclrd = 0._r8
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracs(lev,igc)
               blay = planklay(lev,iband)
               dplankup = planklev(lev,iband) - blay
               dplankdn = planklev(lev-1,iband) - blay
               odepth = secdiff(iband) * taut(lev,igc)
               if (odepth .lt. 0.0_r8) odepth = 0.0_r8
!  Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(lev,igc)
                  if (odtot .lt. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_r8*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,igc) * (1. - atrans(lev))) + &
                         gassrc + cldfmc(igc,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld
                  
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)

                  elseif (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(lev,igc)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1. - exp_tbl(ittot)

                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,igc) * (1._r8 - atrans(lev))) + &
                         gassrc + cldfmc(igc,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)

                  else

                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_r8
                     odepth = tau_tbl(itgas)
                     atrans(lev) = 1._r8 - exp_tbl(itgas)
                     tfacgas = tfn_tbl(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(lev,igc)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._r8 - exp_tbl(ittot)

                  radld = radld - radld * (atrans(lev) + &
                    efclfrac(lev,igc) * (1._r8 - atrans(lev))) + &
                    gassrc + cldfmc(igc,lev) * &
                    (bbdtot * atot(lev) - gassrc)
                  drad(lev-1) = drad(lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif
!  Clear layer
               else
                  if (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth-0.5_r8*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_r8
                     transc = exp_tbl(itr)
                     atrans(lev) = 1._r8-transc
                     tausfac = tfn_tbl(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif   
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
               endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.

         rad0 = fracs(1,igc) * plankbnd(iband)
!  Add in specular reflection of surface downward radiance.
         reflect = 1._r8 - semiss(iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru

         do lev = 1, nlayers
!  Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfrac(lev,igc) * (1._r8 - atrans(lev))) + &
                   gassrc + cldfmc(igc,lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               urad(lev) = urad(lev) + radlu
!  Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
            endif
!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (ICLDDN=1), because surface 
!  reflectance is different for each stream.
               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_r8
            drad(lev) = 0.0_r8
            totuflux(lev) = totuflux(lev) + uflux(lev) * delwave(iband)
            totdflux(lev) = totdflux(lev) + dflux(lev) * delwave(iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_r8
            clrdrad(lev) = 0.0_r8
            totuclfl(lev) = totuclfl(lev) + uclfl(lev) * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + dclfl(lev) * delwave(iband)
            totufluxs(iband,lev) = uflux(lev) * delwave(iband)
            totdfluxs(iband,lev) = dflux(lev) * delwave(iband)
         enddo

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) * fluxfac
      totdflux(0) = totdflux(0) * fluxfac
      totufluxs(:,0) = totufluxs(:,0) * fluxfac
      totdfluxs(:,0) = totdfluxs(:,0) * fluxfac
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) * fluxfac
      totdclfl(0) = totdclfl(0) * fluxfac
      fnetc(0) = totuclfl(0) - totdclfl(0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) * fluxfac
         totdflux(lev) = totdflux(lev) * fluxfac
         totufluxs(:,lev) = totufluxs(:,lev) * fluxfac
         totdfluxs(:,lev) = totdfluxs(:,lev) * fluxfac
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) * fluxfac
         totdclfl(lev) = totdclfl(lev) * fluxfac
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_r8
      htrc(nlayers) = 0.0_r8

      end subroutine rtrnmc

!--------------------------------------
!+
! rrtmg_lw_rtr2function.f90
! This program prepares input parameters for two-stream source function
! technique.
!-
! Reference
! Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam (1989),
! Rapid calculation of radiative heating rates and photodissociation
! rates in inhomogeneous multiple scattering atmospheres, J. Geophys.
! Res., 94(D13), 16287–16301, doi:10.1029/JD094iD13p16287.
!+
! History
! Oct. 17, 2016 Develop the code (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

      subroutine rtr2function(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                              cldfrac, taucloud, ssacloud, xmomcloud, &
                              planklay, planklev, plankbnd, &
                              pwvcm, fracs, taut, &
                              totuflux, totdflux, fnet, htr, &
                              totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs)
!*************************************************************************

! ----- Input -----
      integer, intent(in) :: nlayers         ! total number of layers
      integer, intent(in) :: istart          ! beginning band of calculation
      integer, intent(in) :: iend            ! ending band of calculation
      integer, intent(in) :: iout            ! output option flag

! Atmosphere
      real(kind=r8), intent(in) :: pz(0:)             ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=r8), intent(in) :: semiss(:)          ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: fracs(:,:)         ! 
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=r8), intent(in) :: taut(:,:)          ! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)
      real(kind=r8), intent(in) :: planklay(:,:)      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=r8), intent(in) :: planklev(0:,:)     ! 
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=r8), intent(in) :: plankbnd(:)        ! 
                                                      !    Dimensions: (nbndlw)

! Clouds
! yihsuan 2017-03-02 note: make the dimension of cloud variables be consisnt with CESM
      integer, intent(in) :: ncbands         ! number of cloud spectral bands
                                                      ! Planck derivative [0=off, 1=on]
      !original real(kind=r8), intent(in) :: cldfrac(:)         ! layer cloud fraction
      !                                                !    Dimensions: (nlayers)
      real(kind=r8), intent(in) :: cldfrac(:,:)         ! layer cloud fraction
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: taucloud(:,:)      ! layer cloud optical depth
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: ssacloud(:,:)      ! layer cloud single-scattering albedo
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: xmomcloud(0:,:,:)  ! layer cloud expansion coefficients of phase function
                                                      !    Dimensions: (0:16,ngptlw,nlayers)

! ----- Output -----
      real(kind=r8), intent(out) :: totuflux(0:)      ! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdflux(0:)      ! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnet(0:)          ! net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htr(0:)           ! longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totuclfl(0:)      ! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdclfl(0:)      ! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnetc(0:)         ! clear sky net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htrc(0:)          ! clear sky longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      !>>> yihsuan 2017-06-06, add spectral flux >>>
      real(kind=r8), intent(out) :: totufluxs(:,0:)  ! upward longwave spectral flux (w/m2)
                                                          !    Dimensions: (nbndlw,0:nlayers)

      real(kind=r8), intent(out) :: totdfluxs(:,0:)  ! downward longwave spectral flux (W/m2)
                                                          !    Dimensions: (nbndlw,0:nlayers)
      !<<< yihsuan 2017-06-06, add spectral flux <<<
! ----- Local -----
! Declarations for radiative transfer
      integer :: iband, lay, lev, ig   ! loop indices
      integer :: ibcnd

      real(kind=r8) :: wavenumlo, wavenumhi
      real(kind=r8) :: plkavg
      real(kind=r8) :: diffus
      real(kind=r8) :: albedosuf
      real(kind=r8) :: dnftoa
      real(kind=r8) :: scatcld
      real(kind=r8) :: planksuf
      real(kind=r8) :: gama1(nlayers)
      real(kind=r8) :: gama2(nlayers)
      real(kind=r8) :: fluxupcld(0:nlayers) ! Upward flux under cloudy sky
      real(kind=r8) :: fluxupclr(0:nlayers) ! Upward flux under clear sky
      real(kind=r8) :: fluxdncld(0:nlayers) ! Downward flux under cloudy sky
      real(kind=r8) :: fluxdnclr(0:nlayers) ! Downward flux under clear sky
      real(kind=r8) :: taurevcld(nlayers)
      real(kind=r8) :: taurevclr(nlayers)
      real(kind=r8) :: ssarevcld(nlayers)
      real(kind=r8) :: ssarevclr(nlayers)
      real(kind=r8) :: asyrevcld(nlayers)
      real(kind=r8) :: asyrevclr(nlayers)
      real(kind=r8) :: plankrev(0:nlayers)
      real(kind=r8) :: fracsrev(nlayers)

      real(kind=r8), parameter :: pi = 3.1415926535897932_r8

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    ssacloud                     ! layer cloud single-scattering albedo
!    xmomcloud                    ! layer cloud expansion coefficients of phase function
!    semiss                       ! surface emissivities for each band
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    odclr                        ! clear sky (gaseous) optical depth
!    bbdgas                       ! gas-only planck function for downward rt
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)

      !>>> yihsuan 2017-06-06, add spectral flux >>>
      totufluxs = 0.0_r8 
      totdfluxs = 0.0_r8 
      !<<< yihsuan 2017-06-06, add spectral flux <<<

      totuclfl = 0.0_r8
      totdclfl = 0.0_r8
      totuflux = 0.0_r8
      totdflux = 0.0_r8
      fnetc = 0.0_r8
      htrc = 0.0_r8
      fnet = 0.0_r8
      htr = 0.0_r8
      taurevcld = 0.0_r8
      taurevclr = 0.0_r8
      ssarevcld = 0.0_r8
      ssarevclr = 0.0_r8
      asyrevcld = 0.0_r8
      asyrevclr = 0.0_r8
      plankrev = 0.0_r8
      planksuf = 0.0_r8
      albedosuf = 0.0_r8
      fracsrev = 0.0_r8
      scatcld = 0.0_r8
      gama1 = 0.0_r8
      gama2 = 0.0_r8
      dnftoa = 0.0_r8
      
      ig = 1
! *** loop over frequency bands.
      do iband = istart, iend

! ***    planck function for each level
         do lay = 0, nlayers
            plankrev(nlayers-lay) = planklev(lay,iband) * 1.0e4_r8 * delwave(iband)
         enddo

! ***    planck function at the surface
         planksuf = plankbnd(iband) * 1.0e4_r8 * delwave(iband)

! ***    set surface albedo for this band
         albedosuf = 1.0_r8 - semiss(iband)

! ***    loop over g-channels.
         if (iout.gt.0.and.iband.ge.2) ig = ngs(iband-1)+1

 1000    continue

! ***    downward radiative transfer.

         do lay = nlayers, 1, -1

! ***       fraction of planck function
            fracsrev(nlayers-lay+1) = fracs(lay,ig)

! ***       clear sky
            if (taut(lay,ig) .lt. 0.0_r8) then
               taurevclr(nlayers-lay+1) = 0.0_r8
            else
               taurevclr(nlayers-lay+1) = taut(lay,ig)
            endif
            ssarevclr(nlayers-lay+1) = 0.0_r8
            asyrevclr(nlayers-lay+1) = 0.0_r8
            
! ***       mix optical properties of cloud and gas
!yihsuan cancel 2017-02-22            if (taut(lay,ig) .lt. 0.0_r8) then
!yihsuan cancel 2017-02-22               taurevcld(nlayers-lay+1) = taucloud(lay,iband)
!yihsuan cancel 2017-02-22            else
!yihsuan cancel 2017-02-22               taurevcld(nlayers-lay+1) = taut(lay,ig) + taucloud(lay,iband)
!yihsuan cancel 2017-02-22            endif
!yihsuan cancel 2017-02-22            scatcld = ssacloud(lay,iband) * taucloud(lay,iband)
!yihsuan cancel 2017-02-22            asyrevcld(nlayers-lay+1) = xmomcloud(1,lay,iband)

            !>>> yihsuan add 2017-02-22, note that taucloud dimension is
            !different than original TAMU scheme>>>
            if (taut(lay,ig) .lt. 0.0_r8) then
               taurevcld(nlayers-lay+1) = taucloud(ig,lay)
            else
               taurevcld(nlayers-lay+1) = taut(lay,ig) + taucloud(ig,lay)
            endif
            scatcld = ssacloud(ig,lay) * taucloud(ig,lay)
            asyrevcld(nlayers-lay+1) = xmomcloud(1,ig,lay)
            !<<< yihsuan add 2017-02-22 <<<

            if (taurevcld(nlayers-lay+1) .ne. 0.0) &
                 ssarevcld(nlayers-lay+1) = scatcld / taurevcld(nlayers-lay+1)

            !>>> yihsuan 2017-03-03 cancel >>> 
            !if (ssarevcld(nlayers-lay+1) .gt. 1.0) then
            !   write(*,'(A,F10.5,A,I10,A,I10)') '!! warning ssarevcld:
            !   ', &
            !                ssarevcld(nlayers-lay+1),' > 1.0 at layer
            !                ', lay, &
            !                'at g point', ig
            !endif
            !if (ssarevcld(nlayers-lay+1) .lt. 0.0) then
            !   write(*,'(A,F10.5,A,I10,A,I10)') '!! warning ssarevcld:
            !   ', &
            !                ssarevcld(nlayers-lay+1),' < 0.0 at layer
            !                ', lay, &
            !                'at g point', ig
            !endif
            !<<< yihsuan 2017-03-03 cancel <<<

         enddo

! ***    Cloudy sky

! ***    Hemispheric mean
         diffus = 1.66_r8
         gama1 = diffus * 0.5_r8 * (2.0_r8 - ssarevcld*(1.0_r8+asyrevcld))
         gama2 = diffus * 0.5_r8 * ssarevcld * (1.0_r8-asyrevcld)

! ***    boundary conditions
         dnftoa  = 0.0_r8
         
         call TwsFunctionIR(nlayers,diffus,taurevcld,ssarevcld,fracsrev,plankrev,planksuf,&
                            gama1,gama2,albedosuf,dnftoa,fluxupcld,fluxdncld)

         ! Set unphysical values to 0.0
         do lev = 0, nlayers
            if (fluxdncld(lev) < 0.0) fluxdncld(lev) = 0.0_r8
            if (fluxupcld(lev) < 0.0) fluxupcld(lev) = 0.0_r8
         enddo

         do lev = nlayers, 0, -1
            totuflux(lev) = totuflux(lev) + fluxupcld(nlayers-lev)
            totdflux(lev) = totdflux(lev) + fluxdncld(nlayers-lev) 

            !>>> yihsuan 2017-06-06, add spectral flux >>>
            totufluxs(iband,lev)= totufluxs(iband,lev) + fluxupcld(nlayers-lev)  ! upward spectral flux
            totdfluxs(iband,lev)= totdfluxs(iband,lev) + fluxdncld(nlayers-lev)  ! downward spectral flux
            !<<< yihsuan 2017-06-06, add spectral flux <<<
         enddo

         !>>> yihsuan 2017-03-03 cancel >>>
         !if (fluxdncld(0) .gt. 1.e-5) &
         !     write(*,9000) fluxdncld(0), iband, ig
         !do lev = 0, nlayers
         !   if (fluxdncld(lev) < 0.0) &
         !        write(*,9001) fluxdncld(lev), nlayers-lev, iband, ig
         !enddo
         !<<< yihsuan 2017-03-03 cancel <<<

! *** Clear sky

! ***    Hemispheric mean
         diffus = 1.66_r8
         gama1 = diffus * 0.5_r8 * (2.0_r8 - ssarevclr*(1.0_r8+asyrevclr))
         gama2 = diffus * 0.5_r8 * ssarevclr * (1.0_r8-asyrevclr)

! ***    boundary conditions
         dnftoa  = 0.0_r8
         
         call TwsFunctionIR(nlayers,diffus,taurevclr,ssarevclr,fracsrev,plankrev,planksuf,&
                            gama1,gama2,albedosuf,dnftoa,fluxupclr,fluxdnclr)

         ! Set unphysical values to 0.0
         do lev = 0, nlayers
            if (fluxdnclr(lev) < 0.0) fluxdnclr(lev) = 0.0_r8
            if (fluxupclr(lev) < 0.0) fluxupclr(lev) = 0.0_r8
         enddo

         do lev = nlayers, 0, -1
            totuclfl (lev) = totuclfl(lev) + fluxupclr(nlayers-lev)
            totdclfl (lev) = totdclfl(lev) + fluxdnclr(nlayers-lev) 
            !>>> yihsuan 2017-06-06, add spectral flux >>>
            !totufluxsclr(iband,lev)= fluxupclr(nlayers-lev)  ! upward spectral flux
            !totdfluxsclr(iband,lev)= fluxdnclr(nlayers-lev)  ! downward spectral flux
            !<<< yihsuan 2017-06-06, add spectral flux <<<
         enddo

         !>>> yihsuan 2017-03-03 cancel >>>
         !if (fluxdnclr(0) .gt. 1.e-5) &
         !     write(*,9000) fluxdnclr(0), iband, ig
         !do lev = 0, nlayers
         !   if (fluxdnclr(lev) < 0.0) &
         !        write(*,9001) fluxdnclr(lev), nlayers-lev, iband, ig
         !enddo
         !<<< yihsuan 2017-03-03 cancel <<<

         ig = ig + 1
         if (ig .le. ngs(iband)) go to 1000
            
      enddo
      
      fnet(nlayers)  = totuflux(nlayers) - totdflux(nlayers)
      fnetc(nlayers) = totuclfl(nlayers) - totdclfl(nlayers)
      htr(nlayers)  = 0.0_r8
      htrc(nlayers) = 0.0_r8
      do lev = nlayers-1, 0, -1
         fnet(lev)  = totuflux(lev) - totdflux(lev)
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         htr(lev)  = heatfac * (fnet(lev) -fnet(lev+1)) / (pz(lev) - pz(lev+1))
         htrc(lev) = heatfac * (fnetc(lev) -fnetc(lev+1)) / (pz(lev) - pz(lev+1))
      enddo

 9000 format('DOWNWARD FLUX ',ES15.7,' AT TOA GTR THAN 0. IN BAND ',I3, &
      ' AT IG =',I3,'. POSSIBLE', &
      ' INSTABILITY IN TwoStream.')
 9001 format('DOWNWARD FLUX ',ES15.7,' AT LEVEL ',I3,&
      ' GTR THAN 0. IN BAND ',I2, &
      ' AT IG =',I3,'. POSSIBLE',/, &
      ' INSTABILITY IN TwoStream.')

      end subroutine rtr2function
      
!*************************************************************************
      subroutine TwsFunctionIR(nlayers,difactor,taulay,ssalay,planckfracs,planckflev,planckfsuf,&
                               gama1,gama2,sufalb,dnf0,upf,dnf)
!*************************************************************************

!+
! TwsFunction
! This program is two-stream source function technique.
!-
! Reference
! Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam (1989),
! Rapid calculation of radiative heating rates and photodissociation
! rates in inhomogeneous multiple scattering atmospheres, J. Geophys.
! Res., 94(D13), 16287–16301, doi:10.1029/JD094iD13p16287.
!+
! History
! Nov. 11, 2016 Develop the code (Chia-Pang Kuo)
! Mar. 29, 2017 Fix a bug in paramH and paramG (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

!  Structure of atmosphere
!
!  Top of the atmosphere
!  -------------------- level 0
!         layer 1
!  -------------------- level 1
!         layer 2
!  -------------------- level 2
!            .
!            .
!            .
!  -------------------- level N-2
!         layer N-1
!  -------------------- level N-1
!         layer N
!  -------------------- level N
!         Surface

! ----- Input -----
      integer, intent(in) :: nlayers ! number of atmopheric layers

      real(kind=r8), intent(in) :: difactor        ! diffusivity factor
      real(kind=r8), intent(in) :: dnf0            ! downward flux at the top of the atmosphere (boundary condition)
      real(kind=r8), intent(in) :: sufalb          ! surface albedo
      real(kind=r8), intent(in) :: taulay(:)       ! optical thickness in each layer
                                                   ! dimension(layer)
      real(kind=r8), intent(in) :: ssalay(:)       ! single scattering albedo in each layer
                                                   ! dimension(layer)
      real(kind=r8), intent(in) :: gama1(:)        ! gamma 1
                                                   ! dimension(layer)
      real(kind=r8), intent(in) :: gama2(:)        ! gamma 2
                                                   ! dimension(layer)
      real(kind=r8), intent(in) :: planckflev(0:)  ! plank function in each level
                                                   ! dimension(0:layer)
      real(kind=r8), intent(in) :: planckfracs(:)  ! 
                                                   ! dimensions: (layer)
      real(kind=r8), intent(in) :: planckfsuf      ! plank function on the ground
! ----- Output -----
      real(kind=r8), intent(out) :: upf(0:) ! upward flux at the each level
                                            ! dimension(0:layer)
      real(kind=r8), intent(out) :: dnf(0:) ! downward flux at the each level
                                            ! dimension(0:layer)
! ----- Local -----
      integer, parameter :: nangle = 2 ! number of double gauss quadrature
      integer :: iangle
      integer :: lay      ! atmopheric layer
      integer :: lev      ! atmopheric level
      real(kind=r8) :: difactorev  ! inverse diffusivity factor
      real(kind=r8) :: upfn        ! upward flux at the surface (boundary condition)
      real(kind=r8) :: planktop
      real(kind=r8) :: plankbas
      real(kind=r8) :: lamda(nlayers) 
      real(kind=r8) :: gama(nlayers) 
      real(kind=r8) :: exp1(nlayers)
      real(kind=r8) :: exp2(nlayers)
      real(kind=r8) :: e1n(nlayers)
      real(kind=r8) :: e2n(nlayers)
      real(kind=r8) :: e3n(nlayers)
      real(kind=r8) :: e4n(nlayers)
      real(kind=r8) :: mu1
      real(kind=r8) :: b0(nlayers)
      real(kind=r8) :: b1(nlayers)
      real(kind=r8) :: tempc(nlayers)
      real(kind=r8) :: upc0(nlayers)
      real(kind=r8) :: upcn(nlayers)
      real(kind=r8) :: dnc0(nlayers)
      real(kind=r8) :: dncn(nlayers)
      real(kind=r8) :: upc_tmp(nlayers)
      real(kind=r8) :: dnc_tmp(nlayers)
      real(kind=r8) :: aa(2*nlayers)
      real(kind=r8) :: bb(2*nlayers)
      real(kind=r8) :: cc(2*nlayers)
      real(kind=r8) :: dd(2*nlayers)
      real(kind=r8) :: ee(2*nlayers)
      real(kind=r8) :: xx(2*nlayers)
      real(kind=r8) :: k1n(nlayers)
      real(kind=r8) :: k2n(nlayers)
      real(kind=r8) :: dgausa(nangle)
      real(kind=r8) :: dgausw(nangle)
      real(kind=r8) :: upi(0:nlayers)
      real(kind=r8) :: dni(0:nlayers)
      real(kind=r8) :: cosangle
      real(kind=r8) :: cosanglerev
      real(kind=r8) :: param1(nlayers)
      real(kind=r8) :: param2(nlayers)
      real(kind=r8) :: paramG(nlayers)
      real(kind=r8) :: paramK(nlayers)
      real(kind=r8) :: paramJ(nlayers)
      real(kind=r8) :: paramH(nlayers)
      real(kind=r8) :: alpha1(nlayers)
      real(kind=r8) :: alpha2(nlayers)
      real(kind=r8) :: sigma1(nlayers)
      real(kind=r8) :: sigma2(nlayers)
      real(kind=r8) :: exptmp1(nlayers)
      real(kind=r8) :: exptmp2(nlayers)
      real(kind=r8) :: revpar1(nlayers)
      real(kind=r8) :: revpar2(nlayers)
      
      real(kind=r8), parameter :: pi = 3.1415926535897932_r8

      upf     = 0.0_r8
      dnf     = 0.0_r8
      lamda   = 0.0_r8
      gama    = 0.0_r8
      exp1    = 0.0_r8
      exp2    = 0.0_r8
      e1n     = 0.0_r8
      e2n     = 0.0_r8
      e3n     = 0.0_r8
      e4n     = 0.0_r8
      mu1     = 0.0_r8
      b0      = 0.0_r8
      b1      = 0.0_r8
      tempc   = 0.0_r8
      upc0    = 0.0_r8
      upcn    = 0.0_r8
      dnc0    = 0.0_r8
      dncn    = 0.0_r8
      upc_tmp = 0.0_r8
      dnc_tmp = 0.0_r8
      aa      = 0.0_r8
      bb      = 0.0_r8
      cc      = 0.0_r8
      dd      = 0.0_r8
      ee      = 0.0_r8
      xx      = 0.0_r8
      k1n     = 0.0_r8
      k2n     = 0.0_r8
      upi     = 0.0_r8
      dni     = 0.0_r8
      cosangle = 0.0_r8
      cosanglerev = 0.0_r8
      param1 = 0.0_r8
      param2 = 0.0_r8
      paramG = 0.0_r8
      paramK = 0.0_r8
      paramJ = 0.0_r8
      paramH = 0.0_r8
      alpha1 = 0.0_r8
      alpha2 = 0.0_r8
      sigma1 = 0.0_r8
      sigma2 = 0.0_r8
      exptmp1 = 0.0_r8
      exptmp2 = 0.0_r8
      revpar1 = 0.0_r8
      revpar2 = 0.0_r8
      dgausa  = (/0.2113248_r8, 0.7886752_r8/)
      dgausw  = (/0.5_r8, 0.5_r8/)

      ! lamda in eq. 21
      lamda = sqrt(gama1*gama1 - gama2*gama2)

      ! gamma in eq. 22
      gama = gama2 / (gama1 + lamda)
      
      exp1 = exp(lamda*taulay)
      exp2 = exp(-lamda*taulay)
    
      ! exponential functions in eq. 44
      e1n = 1.0_r8 + gama*exp2
      e2n = 1.0_r8 - gama*exp2
      e3n = gama + exp2
      e4n = gama - exp2

      ! constant in eq. 18 for Hemispheric mean 
      ! (mu = diffusivity factor * inv(diffusivity factor) = 1.0)
      mu1 = 1.0_r8

      ! approximate Planck function by first two terms of Taylor expansion in eq. 25
      ! bn(tau) = b0n + b1n*tau
      do lay = 1, nlayers 
         ! weighted planck function for a g interval
         planktop = planckflev(lay-1) * planckfracs(lay)
         plankbas = planckflev(lay) * planckfracs(lay)
         if (taulay(lay) .lt. 1.0e-6_r8) then
            b1(lay) = 0.0_r8
            b0(lay) = planktop
         else
            b1(lay) = (plankbas - planktop) / taulay(lay)
            b0(lay) = planktop
         endif
      enddo

      ! particular solution in eq. 27
      tempc = 1.0_r8 / (gama1 + gama2)
      upc0 = b0 + b1*(0.0_r8+tempc) ! upward direction at the top of the layer
      upcn = b0 + b1*(taulay+tempc) ! upward direction at the bottom of the layer 
      dnc0 = b0 + b1*(0.0_r8-tempc) ! downward direction at the top of the layer
      dncn = b0 + b1*(taulay-tempc) ! downward direction at the bottom of the layer 

      ! upCn+1(0) - upCn(tau) in eq. 41 and 42
      upc_tmp(1:nlayers-1)  = upc0(2:nlayers) - upcn(1:nlayers-1)
      ! dnCn+1(0) - dnCn(tau) in eq. 41 and 42
      dnc_tmp(1:nlayers-1)  = dnc0(2:nlayers) - dncn(1:nlayers-1)

      ! upward flux at the surface (boundary condition)
      upfn = planckfsuf * planckfracs(nlayers)

      ! develop tridiagonal matrix (aa)Yn-1 + (bb)Yn + (dd)Yn+1 = ee in eq. 39
      aa(1) = 0.0_r8
      bb(1) = e1n(1)
      dd(1) = -1.0_r8*e2n(1)
      ee(1) = dnf0/pi - dnc0(1)

      aa(2*nlayers) = e1n(nlayers) - sufalb*e3n(nlayers)
      bb(2*nlayers) = e2n(nlayers) - sufalb*e4n(nlayers)
      dd(2*nlayers) = 0.0_r8
      ee(2*nlayers) = upfn - upcn(nlayers) + sufalb*dncn(nlayers)

      ! even element in a tridiagonal matrix
      aa(2:2*nlayers-2:2) = e1n(1:nlayers-1)*e2n(2:nlayers) - &
                            e3n(1:nlayers-1)*e4n(2:nlayers)
      bb(2:2*nlayers-2:2) = e2n(1:nlayers-1)*e2n(2:nlayers) - &
                            e4n(1:nlayers-1)*e4n(2:nlayers)
      dd(2:2*nlayers-2:2) = e1n(2:nlayers)*e4n(2:nlayers) - &
                            e2n(2:nlayers)*e3n(2:nlayers)
      ee(2:2*nlayers-2:2) = upc_tmp(1:nlayers-1)*e2n(2:nlayers) - &
                            dnc_tmp(1:nlayers-1)*e4n(2:nlayers)
      ! odd element in a tridiagonal matrix
      aa(3:2*nlayers-1:2) = e2n(1:nlayers-1)*e3n(1:nlayers-1) - &
                            e4n(1:nlayers-1)*e1n(1:nlayers-1)
      bb(3:2*nlayers-1:2) = e1n(1:nlayers-1)*e1n(2:nlayers) - &
                            e3n(1:nlayers-1)*e3n(2:nlayers)
      dd(3:2*nlayers-1:2) = e3n(1:nlayers-1)*e4n(2:nlayers) - &
                            e1n(1:nlayers-1)*e2n(2:nlayers)
      ee(3:2*nlayers-1:2) = upc_tmp(1:nlayers-1)*e3n(1:nlayers-1) - &
                            dnc_tmp(1:nlayers-1)*e1n(1:nlayers-1)

      ! use Thomas algorithm to solve tridiagonal matrix
      call TDMA(2*nlayers,aa,bb,dd,ee,xx)

      ! use eq. 29 and 30 to get coefficients
      k1n = xx(1:2*nlayers-1:2) + xx(2:2*nlayers:2)  
      k2n = xx(1:2*nlayers-1:2) - xx(2:2*nlayers:2)

      ! parameter in Table 3
      difactorev = 1.0_r8 / difactor
      param1 = (difactor-lamda) * difactorev
      param2 = gama * (lamda+difactor) * difactorev
      paramH = k2n * param2
      paramG = k1n * param1
      paramK = k2n * param1
      paramJ = k1n * param2
      sigma1 = b0 - b1*(tempc-difactorev)
      sigma2 = b1
      alpha1 = b0 + b1*(tempc-difactorev)
      alpha2 = b1

      ! downward intensity for a given direction
      do iangle = 1, nangle ! loop over double gauss quadrature

         ! double gauss quadrature
         cosangle  = dgausa(iangle)
         cosanglerev = 1.0_r8 / cosangle

         exptmp1 = exp(-taulay*cosanglerev)
         exptmp2 = exp(-taulay*(lamda+cosanglerev))
         revpar1 = 1.0_r8 / (lamda*cosangle + 1.0_r8)
         revpar2 = 1.0_r8 / (lamda*cosangle - 1.0_r8)
        
         ! incident downward intensity 
         dni(0) = dnf0 / pi
         ! isotropic source
         dnf(0) = dnf0

         ! downward radiative transfer in eq. 56
         do lev = 1, nlayers
            dni(lev) = dni(lev-1)*exptmp1(lev) + &
                       paramJ(lev)*revpar1(lev)*(1.0_r8-exptmp2(lev)) + &
                       paramK(lev)*revpar2(lev)*(exptmp1(lev)-exp2(lev)) + &
                       sigma1(lev)*(1.0_r8-exptmp1(lev)) + &
                       sigma2(lev)*(cosangle*exptmp1(lev)+taulay(lev)-cosangle)
            ! downward flux at the bottom of the layer
            dnf(lev) = dnf(lev) + 2.0_r8*pi*dni(lev)*cosangle*dgausw(iangle)
         enddo

         upi(nlayers) = upi(nlayers) + cosangle*dgausw(iangle)*dni(nlayers)

      enddo
      
      ! reflected upward intensity for a given direction at the
      ! lambertian surface
      upi(nlayers) = 2.0_r8*sufalb*upi(nlayers) + planckfsuf*planckfracs(nlayers)

      ! upward intensity for a given direction
      do iangle = 1, nangle ! loop over double gauss quadrature

         ! double gauss quadrature
         cosangle  = dgausa(iangle)
         cosanglerev = 1.0_r8 / cosangle

         exptmp1 = exp(-taulay*cosanglerev)
         exptmp2 = exp(-taulay*(lamda+cosanglerev))
         revpar1 = 1.0_r8 / (lamda*cosangle + 1.0_r8)
         revpar2 = 1.0_r8 / (lamda*cosangle - 1.0_r8)

         ! upward radiative transfer in eq. 55
         do lev = nlayers, 1, -1
            upi(lev-1) = upi(lev)*exptmp1(lev) + &
                         paramG(lev)*revpar2(lev)*(exptmp1(lev)-exp2(lev)) + &
                         paramH(lev)*revpar1(lev)*(1.0_r8-exptmp2(lev)) + &
                         alpha1(lev)*(1.0_r8-exptmp1(lev)) + &
                         alpha2(lev)*(cosangle-(cosangle+taulay(lev))*exptmp1(lev))
            ! upward flux at the bottom of the layer
            upf(lev-1) = upf(lev-1) + 2.0_r8*pi*upi(lev-1)*cosangle*dgausw(iangle)
         enddo

         upf(nlayers) = upf(nlayers) + 2.0_r8*pi*upi(nlayers)*cosangle*dgausw(iangle)

      enddo

      end subroutine TwsFunctionIR

!*************************************************************************
      subroutine TDMA(N,A,B,C,D,X)
!*************************************************************************
!+
! TDMA.f90
! This program is Tridiagonal matrix algorithm (TDMA), also called Thomas 
! algorithm.
! A tridiagonal matrix for N unknowns can be written as
! A(i)X(i-1) + B(i)X(i) + C(i)X(i+1) = D(i), where A(1) = 0, and C(N) =
! 0.
! --                           -- --  --   --  --
! |B1   C1                     0| | X1 |   | D1 |
! |A2   B2   C2                 | | X2 |   | D2 |
! |     A3   B3   C3            | | X3 |   | D3 |
! |      .    .    .            | | :  | = | :  |
! |           .    .    .       | | :  |   | :  |
! |                .    .   CN-1| | :  |   | :  |
! |0                   AN     BN| | XN |   | DN |
! --                           -- --  --   --  --
!-
! Reference
! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
! http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_%28Thomas_algorithm%29
!+
! History
! Oct. 17, 2016 Develop the code (Chia-Pang Kuo)
!-
! Contact Information
!  Chia-Pang Kuo @ Texas A&M University
!  email: intelb52158@tamu.edu

! ----- Input -----
      integer, intent(in) :: N ! N*N of the matrix
      real(kind=r8), intent(in) :: A(N)
      real(kind=r8), intent(in) :: B(N)
      real(kind=r8), intent(in) :: C(N)
      real(kind=r8), intent(in) :: D(N)
! ----- Output -----
      real(kind=r8), intent(out) :: X(N)
! ----- Local -----
      integer :: i 
      real(kind=r8) :: P(0:N)
      real(kind=r8) :: Q(0:N)
      real(kind=r8) :: denominator

      X = 0.0_r8
      P = 0.0_r8
      Q = 0.0_r8
      
      ! forward elimination
      do i = 1, N
         denominator = B(i) + A(i) * P(i-1)
         P(i) = -C(i) / denominator
         Q(i) = (D(i)-A(i)*Q(i-1)) / denominator
      enddo

      ! back substitution
      !yihsuan comment 2017-11-26, do i = N, 1, -1
      !yihsuan comment 2017-11-26   X(i) = P(i) * X(i+1) + Q(i)
      !yihsuan comment 2017-11-26 enddo

      !>>> yihsuan add 2017-11-26, to avoid use X(N+1) which could be
      !NaN and let the CESM crash >>>
      i = N
      X(i) = P(i) * 0.0_r8 + Q(i)

      do i = N-1, 1, -1
        X(i) = P(i) * X(i+1) + Q(i)
      enddo
      !<<< yihsuan add 2017-11-26 <<<

      end subroutine TDMA

      end module rrtmg_lw_rtrnmc

