!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw.nomcica.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/01/03 21:35:36 $
!

       module rrtmg_sw_rad

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
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                      email:  miacono@aer.com                             *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
! Move call to rrtmg_sw_ini and following use association to 
! GCM initialization area
!      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw

      implicit none

! public interfaces/functions/subroutines
!      public :: rrtmg_sw, inatm_sw, earth_sun
      public :: rrtmg_sw

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_sw &
            (lchnk   ,ncol    ,nlay    ,icld    ,dotau   , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon, &
             inflgsw ,iceflgsw,liqflgsw, &
             cldfr   ,taucld  ,ssacld  ,asmcld  , &
             cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  ,ssaaer  ,asmaer  , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc)

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP_SW to set cloud optical depth based on input
!       cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVRT to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer;
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky and full or partial cloud conditions.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo 
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag
!     Dec 2008: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2

! ------- Declarations

! ----- Input -----
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns     
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum
      logical, intent(in) :: dotau                      ! True -> do tau calculation

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
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: asdir(:)             ! UV/vis surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldir(:)             ! Near-IR surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: asdif(:)             ! UV/vis surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldif(:)             ! Near-IR surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)

      integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: coszen(:)            ! Cosine of solar zenith angle
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: scon                 ! Solar constant (Wm-2)

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfr(:,:)           ! Cloud fraction
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucld(:,:,:)        ! Cloud optical depth
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacld(:,:,:)        ! Cloud single scattering albedo
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcld(:,:,:)        ! Cloud asymmetry parameter
                                                        !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: cicewp(:,:)          ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cliqwp(:,:)          ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reice(:,:)           ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reliq(:,:)           ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
!      real(kind=r8), intent(in) :: ecaer(:,:,:)         ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                        !    Dimensions: (ncol,nlay,naerec)
                                                        ! (non-delta scaled)      

! ----- Output -----

      real(kind=r8), intent(out) :: swuflx(:,:)         ! Total sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: swdflx(:,:)         ! Total sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: swhr(:,:)           ! Total sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(out) :: swuflxc(:,:)        ! Clear sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: swdflxc(:,:)        ! Clear sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(out) :: swhrc(:,:)          ! Clear sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)

      real(kind=r8), intent(out) :: dirdnuv(:,:)        ! Direct downward shortwave flux, UV/vis
      real(kind=r8), intent(out) :: difdnuv(:,:)        ! Diffuse downward shortwave flux, UV/vis
      real(kind=r8), intent(out) :: dirdnir(:,:)        ! Direct downward shortwave flux, near-IR
      real(kind=r8), intent(out) :: difdnir(:,:)        ! Diffuse downward shortwave flux, near-IR

      real(kind=r8), intent(out) :: ninflx(:,:)         ! Net shortwave flux, near-IR
      real(kind=r8), intent(out) :: ninflxc(:,:)        ! Net clear sky shortwave flux, near-IR

! ----- Local -----

! Control
      integer :: nlayers                        ! total number of layers
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: icpr                           ! cldprop/cldprmc use flag
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: idelm                          ! delta-m scaling flag
                                                ! [0 = direct and diffuse fluxes are unscaled]
                                                ! [1 = direct and diffuse fluxes are scaled]
                                                ! (total downward fluxes are always delta scaled)
      integer :: isccos                         ! instrumental cosine response flag (inactive)
      integer :: iplon                          ! column loop index
      integer :: i                              ! layer loop index                       ! jk
      integer :: ib                             ! band loop index                        ! jsw
      integer :: ia, ig                         ! indices
      integer :: k                              ! layer loop index
      integer :: ims                            ! value for changing mcica permute seed
      integer :: imca                           ! flag for mcica [0=off, 1=on]

      real(kind=r8) :: zepsec, zepzen           ! epsilon
      real(kind=r8) :: zdpgcp                   ! flux to heating conversion ratio

! Atmosphere
      real(kind=r8) :: pavel(nlay+1)            ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay+1)            ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay+1)             ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay+1)             ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: pdp(nlay+1)              ! layer pressure thickness (hPa, mb)
      real(kind=r8) :: coldry(nlay+1)           ! dry air column amount
      real(kind=r8) :: wkl(mxmol,nlay+1)        ! molecular amounts (mol/cm-2)

!      real(kind=r8) :: earth_sun                ! function for Earth/Sun distance factor
      real(kind=r8) :: cossza                   ! Cosine of solar zenith angle
      real(kind=r8) :: adjflux(jpband)          ! adjustment for current Earth/Sun distance
      real(kind=r8) :: solvar(jpband)           ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8) :: albdir(nbndsw)           ! surface albedo, direct          ! zalbp
      real(kind=r8) :: albdif(nbndsw)           ! surface albedo, diffuse         ! zalbd

      real(kind=r8) :: taua(nlay+1,nbndsw)      ! Aerosol optical depth
      real(kind=r8) :: ssaa(nlay+1,nbndsw)      ! Aerosol single scattering albedo
      real(kind=r8) :: asma(nlay+1,nbndsw)      ! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer :: laytrop                        ! tropopause layer index
      integer :: layswtch                       ! 
      integer :: laylow                         ! 
      integer :: jp(nlay+1)                     ! 
      integer :: jt(nlay+1)                     !
      integer :: jt1(nlay+1)                    !

      real(kind=r8) :: colh2o(nlay+1)           ! column amount (h2o)
      real(kind=r8) :: colco2(nlay+1)           ! column amount (co2)
      real(kind=r8) :: colo3(nlay+1)            ! column amount (o3)
      real(kind=r8) :: coln2o(nlay+1)           ! column amount (n2o)
      real(kind=r8) :: colch4(nlay+1)           ! column amount (ch4)
      real(kind=r8) :: colo2(nlay+1)            ! column amount (o2)
      real(kind=r8) :: colmol(nlay+1)           ! column amount
      real(kind=r8) :: co2mult(nlay+1)          ! column amount 

      integer :: indself(nlay+1)
      integer :: indfor(nlay+1)
      real(kind=r8) :: selffac(nlay+1)
      real(kind=r8) :: selffrac(nlay+1)
      real(kind=r8) :: forfac(nlay+1)
      real(kind=r8) :: forfrac(nlay+1)

      real(kind=r8) :: &                        !
                         fac00(nlay+1), fac01(nlay+1), &
                         fac10(nlay+1), fac11(nlay+1) 

! Atmosphere/clouds - cldprop
      integer :: ncbands                        ! number of cloud spectral bands
      integer :: inflag                         ! flag for cloud property method
      integer :: iceflag                        ! flag for ice cloud properties
      integer :: liqflag                        ! flag for liquid cloud properties

      real(kind=r8) :: cldfrac(nlay+1)          ! layer cloud fraction
      real(kind=r8) :: tauc(nbndsw,nlay+1)      ! cloud optical depth (non-delta scaled)
      real(kind=r8) :: ssac(nbndsw,nlay+1)      ! cloud single scattering albedo (non-delta scaled)
      real(kind=r8) :: asmc(nbndsw,nlay+1)      ! cloud asymmetry parameter (non-delta scaled)
      real(kind=r8) :: ciwp(nlay+1)             ! cloud ice water path
      real(kind=r8) :: clwp(nlay+1)             ! cloud liquid water path
      real(kind=r8) :: rel(nlay+1)              ! cloud liquid particle effective radius (microns)
      real(kind=r8) :: rei(nlay+1)              ! cloud ice particle effective radius (microns)
      real(kind=r8) :: dge(nlay+1)              ! cloud ice particle generalized effective size (microns)

      real(kind=r8) :: taucloud(nlay+1,jpband)  ! cloud optical depth
      real(kind=r8) :: taucldorig(nlay+1,jpband)! cloud optical depth (non-delta scaled)
      real(kind=r8) :: ssacloud(nlay+1,jpband)  ! cloud single scattering albedo
      real(kind=r8) :: asmcloud(nlay+1,jpband)  ! cloud asymmetry parameter

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=r8) :: ztauc(nlay+1,nbndsw)     ! cloud optical depth
      real(kind=r8) :: ztaucorig(nlay+1,nbndsw) ! unscaled cloud optical depth
      real(kind=r8) :: zasyc(nlay+1,nbndsw)     ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=r8) :: zomgc(nlay+1,nbndsw)     ! cloud single scattering albedo
      real(kind=r8) :: ztaua(nlay+1,nbndsw)     ! total aerosol optical depth
      real(kind=r8) :: zasya(nlay+1,nbndsw)     ! total aerosol asymmetry parameter 
      real(kind=r8) :: zomga(nlay+1,nbndsw)     ! total aerosol single scattering albedo

      real(kind=r8) :: zbbfu(nlay+2)          ! temporary upward shortwave flux (w/m2)
      real(kind=r8) :: zbbfd(nlay+2)          ! temporary downward shortwave flux (w/m2)
      real(kind=r8) :: zbbcu(nlay+2)          ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=r8) :: zbbcd(nlay+2)          ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=r8) :: zbbfddir(nlay+2)       ! temporary downward direct shortwave flux (w/m2)
      real(kind=r8) :: zbbcddir(nlay+2)       ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=r8) :: zuvfd(nlay+2)          ! temporary UV downward shortwave flux (w/m2)
      real(kind=r8) :: zuvcd(nlay+2)          ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=r8) :: zuvfddir(nlay+2)       ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=r8) :: zuvcddir(nlay+2)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=r8) :: znifd(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znicd(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znifddir(nlay+2)       ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=r8) :: znicddir(nlay+2)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)
! Added for near-IR flux diagnostic
      real(kind=r8) :: znifu(nlay+2)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8) :: znicu(nlay+2)          ! temporary clear sky near-IR downward shortwave flux (w/m2)

! Optional output fields 
      real(kind=r8) :: swnflx(nlay+2)         ! Total sky shortwave net flux (W/m2)
      real(kind=r8) :: swnflxc(nlay+2)        ! Clear sky shortwave net flux (W/m2)
      real(kind=r8) :: dirdflux(nlay+2)       ! Direct downward shortwave surface flux
      real(kind=r8) :: difdflux(nlay+2)       ! Diffuse downward shortwave surface flux
      real(kind=r8) :: uvdflx(nlay+2)         ! Total sky downward shortwave flux, UV/vis   
      real(kind=r8) :: nidflx(nlay+2)         ! Total sky downward shortwave flux, near-IR  

! Output - inactive
!      real(kind=r8) :: zuvfu(nlay+2)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvfd(nlay+2)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcu(nlay+2)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcd(nlay+2)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zvsfu(nlay+2)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvsfd(nlay+2)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscu(nlay+2)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscd(nlay+2)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=r8) :: znifu(nlay+2)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znifd(nlay+2)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicu(nlay+2)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicd(nlay+2)         ! temporary clear sky downward near-IR shortwave flux (w/m2)


! Initializations

      zepsec = 1.e-06_r8
      zepzen = 1.e-10_r8
      oneminus = 1.0_r8 - zepsec
      pi = 2._r8 * asin(1._r8)

      istart = jpb1
      iend = jpb2
      icpr = 0

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version does not use McICA (imca = 0) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! Without McICA, SW calculation is limited to clear or fully overcast conditions. 
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 6, use six ECMWF aerosol types
!           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
      iaer = 10

! Set idelm to select between delta-M scaled or unscaled output direct and diffuse fluxes
! NOTE: total downward fluxes are always delta scaled
! idelm = 0, output direct and diffuse flux components are not delta scaled
!            (direct flux does not include forward scattering peak)
! idelm = 1, output direct and diffuse flux components are delta scaled (default)
!            (direct flux includes part or most of forward scattering peak)
      idelm = 1

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_sw_ini

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, ncol

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, &
              h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, scon, &
              inflgsw, iceflgsw, liqflgsw, &
              cldfr, taucld, ssacld, asmcld, cicewp, cliqwp, &
              reice, reliq, tauaer, ssaaer, asmaer, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
              ssac, asmc, ciwp, clwp, rei, dge, rel, taua, ssaa, asma)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

!  Without McICA, SW calculation is limited to clear or fully overcast conditions. 
!  Stop model if partial cloudiness is present.  

         do i = 1, nlayers
            if (cldfrac(i).gt.zepsec .and. cldfrac(i).lt.oneminus) then
               stop 'PARTIAL CLOUD NOT ALLOWED'
            endif
         enddo
         call cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                         tauc, ssac, asmc, ciwp, clwp, rei, dge, rel, &
                         taucldorig, taucloud, ssacloud, asmcloud)
         icpr = 1

! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

         call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                         laytrop, layswtch, laylow, jp, jt, jt1, &
                         co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                         colo2, colo3, fac00, fac01, fac10, fac11, &
                         selffac, selffrac, indself, forfac, forfrac, indfor)


! Cosine of the solar zenith angle 
!  Prevent using value of zero; ideally, SW model is not called from host model when sun 
!  is below horizon

         cossza = coszen(iplon)
         if (cossza .eq. 0._r8) cossza = zepzen


! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
         do ib=1,9
            albdir(ib) = aldir(iplon)
            albdif(ib) = aldif(iplon)
         enddo
         albdir(nbndsw) = aldir(iplon)
         albdif(nbndsw) = aldif(iplon)
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(ib) = asdir(iplon)
            albdif(ib) = asdif(iplon)
         enddo


! Clouds
         if (icld.eq.0) then

            ztauc(:,:) = 0._r8
            ztaucorig(:,:) = 0._r8
            zasyc(:,:) = 0._r8
            zomgc(:,:) = 1._r8

         elseif (icld.ge.1) then
            do i=1,nlayers
               do ib=1,nbndsw
                  if (cldfrac(i) .ge. zepsec) then
                     ztauc(i,ib) = taucloud(i,jpb1-1+ib)
                     ztaucorig(i,ib) = taucldorig(i,jpb1-1+ib)
                     zasyc(i,ib) = asmcloud(i,jpb1-1+ib)
                     zomgc(i,ib) = ssacloud(i,jpb1-1+ib)
                  endif
               enddo
            enddo

         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:,:) = 0._r8
            zasya(:,:) = 0._r8
            zomga(:,:) = 1._r8

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

!            do i = 1, nlayers
!               do ia = 1, naerec
!                  ecaer(iplon,i,ia) = 1.0e-15_r8
!               enddo
!            enddo

!            do i = 1, nlayers
!               do ib = 1, nbndsw
!                  ztaua(i,ib) = 0._r8
!                  zasya(i,ib) = 0._r8
!                  zomga(i,ib) = 1._r8
!                  do ia = 1, naerec
!                     ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia)
!                     zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
!                                   rsrpiza(ib,ia)
!                     zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
!                                   rsrpiza(ib,ia) * rsrasya(ib,ia)
!                  enddo
!                  if (zomga(i,ib) /= 0._r8) then
!                     zasya(i,ib) = zasya(i,ib) / zomga(i,ib)
!                  endif
!                  if (ztaua(i,ib) /= 0._r8) then
!                     zomga(i,ib) = zomga(i,ib) / ztaua(i,ib)
!                  endif
!               enddo
!            enddo

! IAER=10: Direct specification of aerosol optical properties from GCM
         elseif (iaer.eq.10) then

            do i = 1 ,nlayers
               do ib = 1 ,nbndsw
                  ztaua(i,ib) = taua(i,ib)
                  zasya(i,ib) = asma(i,ib)
                  zomga(i,ib) = ssaa(i,ib)
               enddo
            enddo

         endif


! Call the 2-stream radiation transfer model

         do i=1,nlayers+1
            zbbcu(i) = 0._r8
            zbbcd(i) = 0._r8
            zbbfu(i) = 0._r8
            zbbfd(i) = 0._r8
            zbbcddir(i) = 0._r8
            zbbfddir(i) = 0._r8
            zuvcd(i) = 0._r8
            zuvfd(i) = 0._r8
            zuvcddir(i) = 0._r8
            zuvfddir(i) = 0._r8
            znicd(i) = 0._r8
            znifd(i) = 0._r8
            znicddir(i) = 0._r8
            znifddir(i) = 0._r8
            znicu(i) = 0._r8
            znifu(i) = 0._r8
         enddo

         call spcvrt_sw &
             (nlayers, istart, iend, icpr, idelm, iout, dotau, &
              pavel, tavel, pz, tz, tbound, albdif, albdir, &
              cldfrac, ztauc, zasyc, zomgc, ztaucorig, &
              ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, znifu, znicu, &
              zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top

         do i = 1, nlayers+1
            swuflxc(iplon,i) = zbbcu(i)
            swdflxc(iplon,i) = zbbcd(i)
            swuflx(iplon,i) = zbbfu(i)
            swdflx(iplon,i) = zbbfd(i)
            uvdflx(i) = zuvfd(i)
            nidflx(i) = znifd(i)
!  Direct/diffuse fluxes
            dirdflux(i) = zbbfddir(i)
            difdflux(i) = swdflx(iplon,i) - dirdflux(i)
!  UV/visible direct/diffuse fluxes
            dirdnuv(iplon,i) = zuvfddir(i)
            difdnuv(iplon,i) = zuvfd(i) - dirdnuv(iplon,i)
!  Near-IR direct/diffuse fluxes
            dirdnir(iplon,i) = znifddir(i)
            difdnir(iplon,i) = znifd(i) - dirdnir(iplon,i)
!  Added for net near-IR diagnostic
            ninflx(iplon,i) = znifd(i) - znifu(i)
            ninflxc(iplon,i) = znicd(i) - znicu(i)
         enddo

!  Total and clear sky net fluxes
         do i = 1, nlayers+1
            swnflxc(i) = swdflxc(iplon,i) - swuflxc(iplon,i)
            swnflx(i) = swdflx(iplon,i) - swuflx(iplon,i)
         enddo

!  Total and clear sky heating rates
!  Heating units are in K/d.  Flux units are in W/m2.
         do i = 1, nlayers
            zdpgcp = heatfac / pdp(i)
            swhrc(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
            swhr(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
         enddo
         swhrc(iplon,nlayers) = 0._r8
         swhr(iplon,nlayers) = 0._r8

! End longitude loop
      enddo

      end subroutine rrtmg_sw

!*************************************************************************
      real(kind=r8) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer, intent(in) :: idn

      real(kind=r8) :: gamma

      gamma = 2._r8*pi*(idn-1)/365._r8

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_r8 + .034221_r8 * cos(gamma) + .001289_r8 * sin(gamma) + &
                   .000719_r8 * cos(2._r8*gamma) + .000077_r8 * sin(2._r8*gamma)

      end function earth_sun

!***************************************************************************
      subroutine inatm_sw (iplon, nlay, icld, iaer, &
            play, plev, tlay, tlev, tsfc, &
            h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, scon, &
            inflgsw, iceflgsw, liqflgsw, &
            cldfr, taucld, ssacld, asmcld, cicewp, cliqwp, &
            reice, reliq, tauaer, ssaaer, asmaer, &
            nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
            ssac, asmc, ciwp, clwp, rei, dge, rel, taua, ssaa, asma)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrsw, only : nbndsw, ngptsw, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlay                       ! number of model layers
      integer, intent(in) :: icld                       ! clear/cloud flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        ! Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)

      integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: scon                 ! Solar constant (W/m2)

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfr(:,:)           ! Cloud fraction
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucld(:,:,:)        ! Cloud optical depth (optional)
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacld(:,:,:)        ! Cloud single scattering albedo
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcld(:,:,:)        ! Cloud asymmetry parameter
                                                        ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=r8), intent(in) :: cicewp(:,:)          ! Cloud ice water path (g/m2)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cliqwp(:,:)          ! Cloud liquid water path (g/m2)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reice(:,:)           ! Cloud ice effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: reliq(:,:)           ! Cloud water drop effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)

      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter
                                                        ! Dimensions: (ncol,nlay,nbndsw)

! Atmosphere
      integer, intent(out) :: nlayers                   ! number of layers

      real(kind=r8), intent(out) :: pavel(:)            ! layer pressures (mb) 
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel(:)            ! layer temperatures (K)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:)              ! level (interface) pressures (hPa, mb)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:)              ! level (interface) temperatures (K)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: pdp(:)              ! layer pressure thickness (hPa, mb)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: coldry(:)           ! dry air column density (mol/cm2)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(:,:)            ! molecular amounts (mol/cm-2)
                                                        ! Dimensions: (mxmol,nlay)

      real(kind=r8), intent(out) :: adjflux(:)          ! adjustment for current Earth/Sun distance
                                                        ! Dimensions: (jpband)
      real(kind=r8), intent(out) :: solvar(:)           ! solar constant scaling factor from rrtmg_sw
                                                        ! Dimensions: (jpband)
                                                        !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8), intent(out) :: taua(:,:)           ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: ssaa(:,:)           ! Aerosol single scattering albedo
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: asma(:,:)           ! Aerosol asymmetry parameter
                                                        ! Dimensions: (nlay,nbndsw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfrac(:)          ! layer cloud fraction
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: tauc(:,:)           ! cloud optical depth (non-delta scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=r8), intent(out) :: ssac(:,:)           ! cloud single scattering albedo (non-delta-scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=r8), intent(out) :: asmc(:,:)           ! cloud asymmetry parameter (non-delta scaled)
                                                        ! Dimensions: (nbndsw,nlay)
      real(kind=r8), intent(out) :: ciwp(:)             ! cloud ice water path
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: clwp(:)             ! cloud liquid water path
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: rel(:)              ! cloud liquid particle effective radius (microns)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: rei(:)              ! cloud ice particle effective radius (microns)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: dge(:)              ! cloud ice particle generalized effective size (microns)
                                                        ! Dimensions: (nlay)

! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8     ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8     ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8     ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8     ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8    ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8   ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8   ! Molecular weight of nitrous oxide (g/mol)

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib                ! Loop indices
      real(kind=r8) :: amm, summol                      ! 
      real(kind=r8) :: adjflx                           ! flux adjustment for Earth/Sun distance
!      real(kind=r8) :: earth_sun                        ! function for Earth/Sun distance adjustment

! Add one to nlayers here to include extra model layer at top of atmosphere
      nlayers = nlay + 1

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

      wkl(:,:) = 0.0_r8
      cldfrac(:) = 0.0_r8
      tauc(:,:) = 0.0_r8
      ssac(:,:) = 1.0_r8
      asmc(:,:) = 0.0_r8
      ciwp(:) = 0.0_r8
      clwp(:) = 0.0_r8
      rei(:) = 0.0_r8
      dge(:) = 0.0_r8
      rel(:) = 0.0_r8
      taua(:,:) = 0.0_r8
      ssaa(:,:) = 1.0_r8
      asma(:,:) = 0.0_r8
 
! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
      if (dyofyr .gt. 0) then
         adjflx = earth_sun(dyofyr)
      endif

! Set incoming solar flux adjustment to include adjustment for
! current Earth/Sun distance (ADJFLX) and scaling of default internal
! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 
      do ib = jpb1,jpb2
!         solvar(ib) = 1._r8
         solvar(ib) = scon / rrsw_scon
         adjflux(ib) = adjflx * solvar(ib)
      enddo

!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
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
         pdp(l) = pz(l-1) - pz(l)
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

! The following section can be used to set values for an additional layer (from
! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (nlayers) layer values here. 

      pavel(nlayers) = 0.5_r8 * pz(nlayers-1)
      tavel(nlayers) = tavel(nlayers-1)
      pz(nlayers) = 1.e-4_r8
      tz(nlayers-1) = 0.5_r8 * (tavel(nlayers)+tavel(nlayers-1))
      tz(nlayers) = tz(nlayers-1)
      pdp(nlayers) = pz(nlayers-1) - pz(nlayers)
      wkl(1,nlayers) = wkl(1,nlayers-1)
      wkl(2,nlayers) = wkl(2,nlayers-1)
      wkl(3,nlayers) = wkl(3,nlayers-1)
      wkl(4,nlayers) = wkl(4,nlayers-1)
      wkl(6,nlayers) = wkl(6,nlayers-1)
      wkl(7,nlayers) = wkl(7,nlayers-1)
      amm = (1._r8 - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlayers-1)))

! At this point all molecular amounts in wkl are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  

      do l = 1, nlayers
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
      enddo

! Transfer aerosol optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlayers-1
            do ib = 1, nbndsw
               taua(l,ib) = tauaer(iplon,nlayers-l,ib)
               ssaa(l,ib) = ssaaer(iplon,nlayers-l,ib)
               asma(l,ib) = asmaer(iplon,nlayers-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflgsw
         iceflag = iceflgsw
         liqflag = liqflgsw

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
            do n = 1,nbndsw
               tauc(n,l) = taucld(n,iplon,nlayers-l)
               ssac(n,l) = ssacld(n,iplon,nlayers-l)
               asmc(n,l) = asmcld(n,iplon,nlayers-l)
            enddo
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

         cldfrac(nlayers) = 0.0_r8
         tauc(:,nlayers) = 0.0_r8
         ssac(:,nlayers) = 1.0_r8
         asmc(:,nlayers) = 0.0_r8
         ciwp(nlayers) = 0.0_r8
         clwp(nlayers) = 0.0_r8
         rei(nlayers) = 0.0_r8
         dge(nlayers) = 0.0_r8
         rel(nlayers) = 0.0_r8
      
      endif

      end subroutine inatm_sw

      end module rrtmg_sw_rad


