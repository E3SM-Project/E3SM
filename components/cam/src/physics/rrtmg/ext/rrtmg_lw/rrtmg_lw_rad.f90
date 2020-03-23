!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw.nomcica.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/04/24 16:17:28 $
!

       module rrtmg_lw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                       email:  miacono@aer.com                            *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! -------- Modules --------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb 
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      use rrtmg_lw_cldprop, only: cldprop
! Move call to rrtmg_lw_ini and following use association to 
! GCM initialization area
!      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrn, only: rtrn
      use rrtmg_lw_rtrnmr, only: rtrnmr
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

! public interfaces/functions/subroutines
      public :: rrtmg_lw, inatm

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_lw &
            (lchnk   ,ncol    ,nlay    ,icld    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr, n2ovmr  ,cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cldfr   ,taucld  ,cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc)

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPROP to set cloud optical depth based on input
!       cloud properties
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMR (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation with a maximum-random cloud
!       overlap method, or calls RTRN to use random cloud overlap.
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflglw and liqflglw
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input 
!     flag iaer, see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available. 
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum

      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfr(:,:)           ! Cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cicewp(:,:)          ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cliqwp(:,:)          ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reice(:,:)           ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reliq(:,:)           ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucld(:,:,:)        ! Cloud optical depth
                                                        !    Dimensions: (nbndlw,ncol,nlay)
!      real(kind=r8), intent(in) :: ssacld(:,:,:)       ! Cloud single scattering albedo
                                                        !    Dimensions: (nbndlw,ncol,nlay)
                                                        !   for future expansion 
                                                        !   (lw scattering not yet available)
!      real(kind=r8), intent(in) :: asmcld(:,:,:)       ! Cloud asymmetry parameter
                                                        !    Dimensions: (nbndlw,ncol,nlay)
                                                        !   for future expansion 
                                                        !   (lw scattering not yet available)
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=r8), intent(in) :: ssaaer(:,:,:)       ! aerosol single scattering albedo
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)
!      real(kind=r8), intent(in) :: asmaer(:,:,:)       ! aerosol asymmetry parameter
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)


! ----- Output -----

      real(kind=r8), intent(out) :: uflx(:,:)           ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: dflx(:,:)           ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: hr(:,:)             ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out) :: uflxc(:,:)          ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: dflxc(:,:)          ! Clear sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: hrc(:,:)            ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)

! ----- Local -----

! Control
      integer :: nlayers                        ! total number of layers
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: iplon                          ! column loop index
      integer :: imca                           ! flag for mcica [0=off, 1=on]
      integer :: k                              ! layer loop index
      integer :: ig                             ! g-point loop index

! Atmosphere
      real(kind=r8) :: pavel(nlay+1)            ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay+1)            ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay+1)             ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay+1)             ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: coldry(nlay+1)           ! dry air column density (mol/cm2)
      real(kind=r8) :: wbrodl(nlay+1)           ! broadening gas column density (mol/cm2)
      real(kind=r8) :: wkl(mxmol,nlay+1)        ! molecular amounts (mol/cm-2)
      real(kind=r8) :: wx(maxxsec,nlay+1)       ! cross-section amounts (mol/cm-2)
      real(kind=r8) :: pwvcm                    ! precipitable water vapor (cm)
      real(kind=r8) :: semiss(nbndlw)           ! lw surface emissivity
      real(kind=r8) :: fracs(nlay+1,ngptlw)     ! 
      real(kind=r8) :: taug(nlay+1,ngptlw)      ! gaseous optical depths
      real(kind=r8) :: taut(nlay+1,ngptlw)      ! gaseous + aerosol optical depths

      real(kind=r8) :: taua(nlay+1,nbndlw)      ! aerosol optical depth
!      real(kind=r8) :: ssaa(nlay+1,nbndlw)      ! aerosol single scattering albedo
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)
!      real(kind=r8) :: asma(nlay+1,nbndlw)      ! aerosol asymmetry parameter
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer :: laytrop                          ! tropopause layer index
      integer :: jp(nlay+1)                       ! lookup table index 
      integer :: jt(nlay+1)                       ! lookup table index 
      integer :: jt1(nlay+1)                      ! lookup table index 
      real(kind=r8) :: planklay(nlay+1,nbndlw)    ! 
      real(kind=r8) :: planklev(0:nlay+1,nbndlw)  ! 
      real(kind=r8) :: plankbnd(nbndlw)           ! 

      real(kind=r8) :: colh2o(nlay+1)             ! column amount (h2o)
      real(kind=r8) :: colco2(nlay+1)             ! column amount (co2)
      real(kind=r8) :: colo3(nlay+1)              ! column amount (o3)
      real(kind=r8) :: coln2o(nlay+1)             ! column amount (n2o)
      real(kind=r8) :: colco(nlay+1)              ! column amount (co)
      real(kind=r8) :: colch4(nlay+1)             ! column amount (ch4)
      real(kind=r8) :: colo2(nlay+1)              ! column amount (o2)
      real(kind=r8) :: colbrd(nlay+1)             ! column amount (broadening gases)

      integer :: indself(nlay+1)
      integer :: indfor(nlay+1)
      real(kind=r8) :: selffac(nlay+1)
      real(kind=r8) :: selffrac(nlay+1)
      real(kind=r8) :: forfac(nlay+1)
      real(kind=r8) :: forfrac(nlay+1)

      integer :: indminor(nlay+1)
      real(kind=r8) :: minorfrac(nlay+1)
      real(kind=r8) :: scaleminor(nlay+1)
      real(kind=r8) :: scaleminorn2(nlay+1)

      real(kind=r8) :: &                          !
                         fac00(nlay+1), fac01(nlay+1), &
                         fac10(nlay+1), fac11(nlay+1) 
      real(kind=r8) :: &                          !
                         rat_h2oco2(nlay+1),rat_h2oco2_1(nlay+1), &
                         rat_h2oo3(nlay+1),rat_h2oo3_1(nlay+1), &
                         rat_h2on2o(nlay+1),rat_h2on2o_1(nlay+1), &
                         rat_h2och4(nlay+1),rat_h2och4_1(nlay+1), &
                         rat_n2oco2(nlay+1),rat_n2oco2_1(nlay+1), &
                         rat_o3co2(nlay+1),rat_o3co2_1(nlay+1)

! Atmosphere/clouds - cldprop
      integer :: ncbands                          ! number of cloud spectral bands
      integer :: inflag                           ! flag for cloud property method
      integer :: iceflag                          ! flag for ice cloud properties
      integer :: liqflag                          ! flag for liquid cloud properties

      real(kind=r8) :: cldfrac(nlay+1)            ! layer cloud fraction
      real(kind=r8) :: tauc(nbndlw,nlay+1)        ! cloud optical depth
!      real(kind=r8) :: ssac(nbndlw,nlay+1)        ! cloud single scattering albedo
                                                  !   for future expansion 
                                                  !   (lw scattering not yet available)
!      real(kind=r8) :: asmc(nbndlw,nlay+1)        ! cloud asymmetry parameter
                                                  !   for future expansion 
                                                  !   (lw scattering not yet available)
      real(kind=r8) :: ciwp(nlay+1)               ! cloud ice water path
      real(kind=r8) :: clwp(nlay+1)               ! cloud liquid water path
      real(kind=r8) :: rel(nlay+1)                ! cloud liquid particle effective radius (microns)
      real(kind=r8) :: rei(nlay+1)                ! cloud ice particle effective radius (microns)
      real(kind=r8) :: dge(nlay+1)                ! cloud ice particle generalized effective size (microns)
      real(kind=r8) :: taucloud(nlay+1,nbndlw)    ! layer cloud optical depth

! Output
      real(kind=r8) :: totuflux(0:nlay+1)         ! upward longwave flux (w/m2)
      real(kind=r8) :: totdflux(0:nlay+1)         ! downward longwave flux (w/m2)
      real(kind=r8) :: fnet(0:nlay+1)             ! net longwave flux (w/m2)
      real(kind=r8) :: htr(0:nlay+1)              ! longwave heating rate (k/day)
      real(kind=r8) :: totuclfl(0:nlay+1)         ! clear sky upward longwave flux (w/m2)
      real(kind=r8) :: totdclfl(0:nlay+1)         ! clear sky downward longwave flux (w/m2)
      real(kind=r8) :: fnetc(0:nlay+1)            ! clear sky net longwave flux (w/m2)
      real(kind=r8) :: htrc(0:nlay+1)             ! clear sky longwave heating rate (k/day)

! Local storage arrays for taug and fracs, etc. for use when dotau is false 

      real(kind=r8), allocatable, save :: taugst(:,:,:,:)                  ! Optical depth storage array
      real(kind=r8), allocatable, save :: fracst(:,:,:,:)                  ! Planck fraction storage array
      real(kind=r8), allocatable, save :: planklayst(:,:,:,:)              ! 
      real(kind=r8), allocatable, save :: planklevst(:,:,:,:)              ! 
      real(kind=r8), allocatable, save :: plankbndst(:,:,:)                ! 
      real(kind=r8), allocatable, save :: pwvcmst(:,:)                     ! precipitable water vapor (cm)


! Allocate storage arrays
      if (.not.allocated(taugst)) allocate (taugst(pcols,nlay+1,ngptlw,begchunk:endchunk))
      if (.not.allocated(fracst)) allocate (fracst(pcols,nlay+1,ngptlw,begchunk:endchunk))
      if (.not.allocated(planklayst)) allocate (planklayst(pcols,nlay+1,nbndlw,begchunk:endchunk))
      if (.not.allocated(planklevst)) allocate (planklevst(pcols,0:nlay+1,nbndlw,begchunk:endchunk))
      if (.not.allocated(plankbndst)) allocate (plankbndst(pcols,nbndlw,begchunk:endchunk))
      if (.not.allocated(pwvcmst)) allocate (pwvcmst(pcols,begchunk:endchunk))
       
! Initializations

      oneminus = 1._r8 - 1.e-6_r8
      pi = 2._r8*asin(1._r8)
      fluxfac = pi * 2.e4_r8                    ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0

! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version does not use McICA (imca = 0) ***

! Set default icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 10, input total aerosol optical depth (tauaer) directly 
      iaer = 10

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_lw_ini

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, ncol

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfr, taucld, cicewp, cliqwp, reice, reliq, tauaer, &
              nlayers, pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfrac, tauc, ciwp, clwp, rei, dge, rel, taua)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

         call cldprop(nlayers, inflag, iceflag, liqflag, cldfrac, tauc, &
                      ciwp, clwp, rei, dge, rel, ncbands, taucloud)

! Perform calculation of optical depth and planck fractions only at interval
! specified by dotau and store output arrays taug and fracs in memory.  In
! intervening time steps, obtain taug and fracs from stored arrays. 
! Memory storage should be replaced with I/O to netCDF file.

         if (dotau) then
 
! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

            call setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, &
                      coldry, wkl, wbrodl, &
                      laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                      colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                      colbrd, fac00, fac01, fac10, fac11, &
                      rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                      rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                      selffac, selffrac, indself, forfac, forfrac, indfor, &
                      minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

            call taumol(nlayers, pavel, wx, coldry, &
                     laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                     colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                     colbrd, fac00, fac01, fac10, fac11, &
                     rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                     rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                     rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     minorfrac, scaleminor, scaleminorn2, indminor, &
                     fracs, taug)

! Store taug and fracs for use when dotau is false

            taugst(iplon,:nlayers,:ngptlw,lchnk) = taug(:nlayers,:ngptlw)
            fracst(iplon,:nlayers,:ngptlw,lchnk) = fracs(:nlayers,:ngptlw)
            planklayst(iplon,:nlayers,:nbndlw,lchnk) = planklay(:nlayers,:nbndlw)
            planklevst(iplon,0:nlayers,:nbndlw,lchnk) = planklev(0:nlayers,:nbndlw)
            plankbndst(iplon,:nbndlw,lchnk) = plankbnd(:nbndlw)
            pwvcmst(iplon,lchnk) = pwvcm

         else

! Restore taug, fracs and other arrays from storage when dotau is false

            taug(:nlayers,:ngptlw) = taugst(iplon,:nlayers,:ngptlw,lchnk)
            fracs(:nlayers,:ngptlw) = fracst(iplon,:nlayers,:ngptlw,lchnk)
            planklay(:nlayers,:nbndlw) = planklayst(iplon,:nlayers,:nbndlw,lchnk)
            planklev(0:nlayers,:nbndlw) = planklevst(iplon,0:nlayers,:nbndlw,lchnk)
            plankbnd(:nbndlw) = plankbndst(iplon,:nbndlw,lchnk)
            pwvcm = pwvcmst(iplon,lchnk)

         endif

! Combine gaseous and aerosol optical depths, if aerosol active
         if (iaer .eq. 0) then
            do k = 1, nlayers
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig)
               enddo
            enddo
         elseif (iaer .eq. 10) then
            do k = 1, nlayers
               do ig = 1, ngptlw 
                  taut(k,ig) = taug(k,ig) + taua(k,ngb(ig))
               enddo
            enddo
         endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.

        if (icld .eq. 1) then
           call rtrn(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                  cldfrac, taucloud, planklay, planklev, plankbnd, &
                  pwvcm, fracs, taut, &
                  totuflux, totdflux, fnet, htr, &
                  totuclfl, totdclfl, fnetc, htrc ) 
        else
           call rtrnmr(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                  cldfrac, taucloud, planklay, planklev, plankbnd, &
                  pwvcm, fracs, taut, &
                  totuflux, totdflux, fnet, htr, &
                  totuclfl, totdclfl, fnetc, htrc ) 
        endif

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top

           do k = 0, nlayers
              uflx(iplon,k+1) = totuflux(k)
              dflx(iplon,k+1) = totdflux(k)
              uflxc(iplon,k+1) = totuclfl(k)
              dflxc(iplon,k+1) = totdclfl(k)
           enddo
           do k = 0, nlayers-1
              hr(iplon,k+1) = htr(k)
              hrc(iplon,k+1) = htrc(k)
           enddo

      enddo

      end subroutine rrtmg_lw

!***************************************************************************
      subroutine inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfr, taucld, cicewp, cliqwp, reice, reliq, tauaer, &
              nlayers, pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfrac, tauc, ciwp, clwp, rei, dge, rel, taua)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
!  Set other RRTMG_LW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, nmol, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, ixindx

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfr(:,:)           ! Cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cicewp(:,:)          ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cliqwp(:,:)          ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reice(:,:)           ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reliq(:,:)           ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucld(:,:,:)        ! Cloud optical depth
                                                        !    Dimensions: (nbndlw,ncol,nlay)
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth
                                                        !    Dimensions: (ncol,nlay,nbndlw)

! ----- Output -----
! Atmosphere
      integer, intent(out) :: nlayers                   ! number of layers

      real(kind=r8), intent(out) :: pavel(:)            ! layer pressures (mb) 
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel(:)            ! layer temperatures (K)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:)              ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:)              ! level (interface) temperatures (K)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: coldry(:)           ! dry air column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wbrodl(:)           ! broadening gas column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(:,:)            ! molecular amounts (mol/cm-2)
                                                        !    Dimensions: (mxmol,nlay)
      real(kind=r8), intent(out) :: wx(:,:)             ! cross-section amounts (mol/cm-2)
                                                        !    Dimensions: (maxxsec,nlay)
      real(kind=r8), intent(out) :: pwvcm               ! precipitable water vapor (cm)
      real(kind=r8), intent(out) :: semiss(:)           ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfrac(:)          ! layer cloud fraction
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: ciwp(:)             ! cloud ice water path
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: clwp(:)             ! cloud liquid water path
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: rel(:)              ! cloud liquid particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: rei(:)              ! cloud ice particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: dge(:)              ! cloud ice particle generalized effective size (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: tauc(:,:)           ! cloud optical depth
                                                        !    Dimensions: (nbndlw,nlay)
      real(kind=r8), intent(out) :: taua(:,:)           ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndlw)


! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8      ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8      ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8     ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8    ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8    ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=r8), parameter :: amc11 = 137.3684_r8   ! Molecular weight of CFC11 (g/mol) - CCL3F
!      real(kind=r8), parameter :: amc12 = 120.9138_r8   ! Molecular weight of CFC12 (g/mol) - CCL2F2
!      real(kind=r8), parameter :: amc22 = 86.4688_r8    ! Molecular weight of CFC22 (g/mol) - CHCLF2
!      real(kind=r8), parameter :: amcl4 = 153.823_r8    ! Molecular weight of CCL4 (g/mol) - CCL4

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
      real(kind=r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
      real(kind=r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib                ! Loop indices
      real(kind=r8) :: amm, amttl, wvttl, wvsh, summol  

! Add one to nlayers here to include extra model layer at top of atmosphere
      nlayers = nlay + 1

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.

      wkl(:,:) = 0.0_r8
      wx(:,:) = 0.0_r8
      cldfrac(:) = 0.0_r8
      tauc(:,:) = 0.0_r8
      ciwp(:) = 0.0_r8
      clwp(:) = 0.0_r8
      rei(:) = 0.0_r8
      dge(:) = 0.0_r8
      rel(:) = 0.0_r8
      taua(:,:) = 0.0_r8
      amttl = 0.0_r8
      wvttl = 0.0_r8
 
!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,nlayers)
      tz(0) = tlev(iplon,nlayers)
      do l = 1, nlayers-1
         pavel(l) = play(iplon,nlayers-l)
         tavel(l) = tlay(iplon,nlayers-l)
         pz(l) = plev(iplon,nlayers-l)
         tz(l) = tlev(iplon,nlayers-l)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,nlayers-l)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,nlayers-l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,nlayers-l)/(1._r8 - h2o(iplon,nlayers-l)))*amdw
         wkl(2,l) = co2vmr(iplon,nlayers-l)
         wkl(3,l) = o3vmr(iplon,nlayers-l)
         wkl(4,l) = n2ovmr(iplon,nlayers-l)
         wkl(6,l) = ch4vmr(iplon,nlayers-l)
         wkl(7,l) = o2vmr(iplon,nlayers-l) 
         amm = (1._r8 - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(1,l)))
      enddo

! Set cross section molecule amounts from input; convert to vmr if necessary
      do l=1, nlayers-1
         wx(1,l) = ccl4vmr(iplon,nlayers-l)
         wx(2,l) = cfc11vmr(iplon,nlayers-l)
         wx(3,l) = cfc12vmr(iplon,nlayers-l)
         wx(4,l) = cfc22vmr(iplon,nlayers-l)
      enddo      

! The following section can be used to set values for an additional layer (from
! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (extra) layer values here. 

      pavel(nlayers) = 0.5_r8 * pz(nlayers-1)
      tavel(nlayers) = tavel(nlayers-1)
      pz(nlayers) = 1.e-4_r8
      tz(nlayers-1) = 0.5_r8 * (tavel(nlayers)+tavel(nlayers-1))
      tz(nlayers) = tz(nlayers-1)
      wkl(1,nlayers) = wkl(1,nlayers-1)
      wkl(2,nlayers) = wkl(2,nlayers-1)
      wkl(3,nlayers) = wkl(3,nlayers-1)
      wkl(4,nlayers) = wkl(4,nlayers-1)
      wkl(6,nlayers) = wkl(6,nlayers-1)
      wkl(7,nlayers) = wkl(7,nlayers-1)
      amm = (1._r8 - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlayers-1)))
      wx(1,nlayers) = wx(1,nlayers-1)
      wx(2,nlayers) = wx(2,nlayers-1)
      wx(3,nlayers) = wx(3,nlayers-1)
      wx(4,nlayers) = wx(4,nlayers-1)

! At this point all moleculular amounts in wkl and wx are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.

      do l = 1, nlayers
         summol = 0.0_r8
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         wbrodl(l) = coldry(l) * (1._r8 - summol)
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         do ix = 1,maxxsec
            if (ixindx(ix) .ne. 0) then
               wx(ixindx(ix),l) = coldry(l) * wx(ix,l) * 1.e-20_r8
            endif
         enddo
      enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_r8 * pz(0)) / (1.e2_r8 * grav)

! Set spectral surface emissivity for each longwave band.  

      do n=1,nbndlw
         semiss(n) = emis(iplon,n)
!          semiss(n) = 1.0_r8
      enddo

! Transfer aerosol optical properties to RRTM variable;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlayers-1
            do ib = 1, nbndlw
               taua(l,ib) = tauaer(iplon,nlayers-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlayers-1
            cldfrac(l) = cldfr(iplon,nlayers-l)
            ciwp(l) = cicewp(iplon,nlayers-l)
            clwp(l) = cliqwp(iplon,nlayers-l)
            rei(l) = reice(iplon,nlayers-l)
            if (iceflag .eq. 3) then
               dge(l) = 1.5396_r8 * reice(iplon,nlayers-l)
            endif
            rel(l) = reliq(iplon,nlayers-l)
            do n=1,nbndlw
               tauc(n,l) = taucld(n,iplon,nlayers-l)
!               ssac(n,l) = ssacld(n,iplon,nlayers-l)
!               asmc(n,l) = asmcld(n,iplon,nlayers-l)
            enddo
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

         cldfrac(nlayers) = 0.0_r8
         tauc(:nbndlw,nlayers) = 0.0_r8
         ciwp(nlayers) = 0.0_r8
         clwp(nlayers) = 0.0_r8
         rei(nlayers) = 0.0_r8
         dge(nlayers) = 0.0_r8
         rel(nlayers) = 0.0_r8
         taua(nlayers,:) = 0.0_r8

      endif
      
      end subroutine inatm

      end module rrtmg_lw_rad

