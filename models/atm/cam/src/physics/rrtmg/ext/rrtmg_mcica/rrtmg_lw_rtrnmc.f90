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
      enddo
      if (pwvcm.lt.1.0) secdiff(6) = 1.80_r8
      if (pwvcm.gt.7.1) secdiff(7) = 1.50_r8

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

      end module rrtmg_lw_rtrnmc

