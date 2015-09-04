!=======================================================================
!BOP
!
! !MODULE: ice_flux - flux variable declarations: coupler, diagnostic and
!          internal
!
! !DESCRIPTION:
!
! Flux variable declarations; these include fields sent from the coupler
! ("in"), sent to the coupler ("out"), written to diagnostic history files
! ("diagnostic"), and used internally ("internal").
!
! !REVISION HISTORY:
!  SVN:$Id: ice_flux.F90 51 2007-01-29 22:25:24Z eclare $
!
! author Elizabeth C. Hunke, LANL
!
! 2004: Block structure added by William Lipscomb
!       Swappped, revised, and added some subroutines
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_flux
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &

       ! in from atmos (if .not.calc_strair)
         strax   , & ! wind stress components (N/m^2)
         stray   , & !

       ! in from ocean
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty , & ! sea surface slope, y-direction

       ! out to atmosphere (if calc_strair)

         strairxT, & ! stress on ice by air, x-direction
         strairyT, & ! stress on ice by air, y-direction

       ! use for super-cycling the dynamics

         strairxT_accum, & ! stress on ice by air, x-direction
         strairyT_accum, & ! stress on ice by air, y-direction

       ! out to ocean          T-cell (kg/m s^2)
       ! Note, CICE_IN_NEMO uses strocnx and strocny for coupling
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

       ! diagnostic

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         sig1    , & ! principal stress component
         sig2    , & ! principal stress component
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strtltx , & ! stress due to sea surface slope, x-direction
         strtlty , & ! stress due to sea surface slope, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         daidtd  , & ! ice area tendency due to transport   (1/s)
         dvidtd  , & ! ice volume tendency due to transport (m/s)
         dardg1dt, & ! rate of area loss by ridging ice (1/s)
         dardg2dt, & ! rate of area gain by new ridges (1/s)
         dvirdgdt, & ! rate of ice volume ridged (m/s)
         opening     ! rate of opening due to divergence/shear (1/s)

       ! restart

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
       ! ice stress tensor in each corner of T cell (kg/s^2)
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,max_blocks) :: &
         iceumask   ! ice extent mask (U-cell)

       ! internal

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         prs_sig  , & ! replacement pressure, for stress calc
         fm           ! Coriolis param. * mass in U-cell (kg/s)

      !-----------------------------------------------------------------
      ! Thermodynamic component
      !-----------------------------------------------------------------

       ! in from atmosphere (if calc_Tsfc)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         zlvl    , & ! atm level height (m)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         wind    , & ! wind speed (m/s)
         potT    , & ! air potential temperature  (K)
         Tair    , & ! air temperature  (K)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         flw         ! incoming longwave radiation (W/m^2)

       ! in from atmosphere (if .not. Tsfc_calc)
       ! required for coupling to HadGEM3
       ! NOTE: when in CICE_IN_NEMO mode, these are gridbox mean fields,
       ! not per ice area. When in standalone mode, these are per ice area.

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         fsurfn_f   , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f, & ! downward cond flux at top surface (W m-2)
         flatn_f        ! latent heat flux (W m-2)

       ! in from atmosphere

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow       ! snowfall rate (kg/m^2 s)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,n_aeromx,max_blocks) :: &
         faero    ! aerosol deposition rate (kg/m^2 s)  MH

       ! in from ocean

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         frzmlt  , & ! freezing/melting potential (W/m^2)
         Tf      , & ! freezing temperature (C)
         qdp     , & ! deep ocean heat flux (W/m^2), negative upward
         hmix        ! mixed layer depth (m)

      character (char_len) :: &
         Tfrzpt      ! ocean freezing temperature formulation
                     ! 'constant' (-1.8C), 'linear_S'

       ! out to atmosphere (if calc_Tsfc)
       ! note Tsfc is in ice_state.F

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         fsens   , & ! sensible heat flux (W/m^2)
         flat    , & ! latent heat flux   (W/m^2)
         fswabs  , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         flwout  , & ! outgoing longwave radiation (W/m^2)
         Uref    , & ! 10m reference wind speed (m/s)
         Tref    , & ! 2m atm reference temperature (K)
         Qref    , & ! 2m atm reference spec humidity (kg/kg)
         evap        ! evaporative water flux (kg/m^2/s)

       ! albedos aggregated over categories (if calc_Tsfc)
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         alvdr   , & ! visible, direct   (fraction)
         alidr   , & ! near-ir, direct   (fraction)
         alvdf   , & ! visible, diffuse  (fraction)
         alidf   , & ! near-ir, diffuse  (fraction)
#ifdef AEROFRC
         dalvdr_noaero   , & ! visible, direct   (fraction) (diag)
         dalidr_noaero   , & ! near-ir, direct   (fraction) (diag)
         dalvdf_noaero   , & ! visible, diffuse  (fraction) (diag)
         dalidf_noaero   , & ! near-ir, diffuse  (fraction) (diag)
#endif
#ifdef CCSM3FRC
         dalvdr_ccsm3   , & ! visible, direct   (fraction) (diag)
         dalidr_ccsm3   , & ! near-ir, direct   (fraction) (diag)
         dalvdf_ccsm3   , & ! visible, diffuse  (fraction) (diag)
         dalidf_ccsm3   , & ! near-ir, diffuse  (fraction) (diag)
#endif
#ifdef PONDFRC
         dalvdr_nopond   , & ! visible, direct   (fraction) (diag)
         dalidr_nopond   , & ! near-ir, direct   (fraction) (diag)
         dalvdf_nopond   , & ! visible, diffuse  (fraction) (diag)
         dalidf_nopond   , & ! near-ir, diffuse  (fraction) (diag)
#endif
         ! grid-box-mean versions
         alvdr_gbm, & ! visible, direct   (fraction)
         alidr_gbm, & ! near-ir, direct   (fraction)
         alvdf_gbm, & ! visible, diffuse  (fraction)
         alidf_gbm, & ! near-ir, diffuse  (fraction)
         ! components for history
#ifdef AEROFRC
         dalbice_noaero   , & ! bare ice albedo (diag)
         dalbsno_noaero   , & ! snow albedo (diag)
         dalbpnd_noaero   , & ! melt pond albedo (diag)
#endif
#ifdef CCSM3FRC
         dalbice_ccsm3   , & ! bare ice albedo (diag)
         dalbsno_ccsm3   , & ! snow albedo (diag)
#endif
#ifdef PONDFRC
         dalbice_nopond   , & ! bare ice albedo (diag)
         dalbsno_nopond   , & ! snow albedo (diag)
         dalbpnd_nopond   , & ! melt pond albedo (diag)
#endif
         albice   , & ! bare ice albedo
         albsno   , & ! snow albedo
         albpnd       ! melt pond albedo

      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,max_blocks,max_nstrm) :: &
         albcnt       ! counter for zenith angle

       ! out to ocean
       ! (Note CICE_IN_NEMO does not use these for coupling.
       !  It uses fresh_gbm,fsalt_gbm,fhocn_gbm and fswthru_gbm)
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         fresh   , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt   , & ! salt flux to ocean (kg/m^2/s)
         fhocn   , & ! net heat flux to ocean (W/m^2)
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
         fswsfc  , & ! shortwave absorbed at the surface (W/m^2)
         fswint  , & ! shortwave absorbed internally (W/m^2)
#endif
#ifdef AEROFRC
         dfswabs_noaero  , & ! shortwave absorbed (diag) (W/m^2)
         dfswsfc_noaero  , & ! shortwave absorbed at the surface (diag) (W/m^2)
         dfswint_noaero  , & ! shortwave absorbed internally (diag) (W/m^2)
         dfswthru_noaero , & ! shortwave penetrating to ocean (diag) (W/m^2)
#endif
#ifdef CCSM3FRC
         dfswabs_ccsm3  , & ! shortwave absorbed (diag) (W/m^2)
         dfswsfc_ccsm3  , & ! shortwave absorbed at the surface (diag) (W/m^2)
         dfswint_ccsm3  , & ! shortwave absorbed internally (diag) (W/m^2)
         dfswthru_ccsm3 , & ! shortwave penetrating to ocean (diag) (W/m^2)
#endif
#ifdef PONDFRC
         dfswabs_nopond  , & ! shortwave absorbed (diag) (W/m^2)
         dfswsfc_nopond  , & ! shortwave absorbed at the surface (diag) (W/m^2)
         dfswint_nopond  , & ! shortwave absorbed internally (diag) (W/m^2)
         dfswthru_nopond , & ! shortwave penetrating to ocean (diag) (W/m^2)
#endif
         fswthru     ! shortwave penetrating to ocean (W/m^2)

      real (kind=dbl_kind), &
        dimension (nx_block,ny_block,n_aeromx,max_blocks) :: &
         fsoot        

       ! internal

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_blocks) :: &
         fswfac  , & ! for history
         scale_factor! scaling factor for shortwave components

      logical (kind=log_kind) :: &
         update_ocn_f ! if true, update fresh water and salt fluxes

      !-----------------------------------------------------------------
      ! quantities passed from ocean mixed layer to atmosphere
      ! (for running with CAM)
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         strairx_ocn , & ! stress on ocean by air, x-direction
         strairy_ocn , & ! stress on ocean by air, y-direction
         fsens_ocn   , & ! sensible heat flux (W/m^2)
         flat_ocn    , & ! latent heat flux   (W/m^2)
         flwout_ocn  , & ! outgoing longwave radiation (W/m^2)
         evap_ocn    , & ! evaporative water flux (kg/m^2/s)
         alvdr_ocn   , & ! visible, direct   (fraction)
         alidr_ocn   , & ! near-ir, direct   (fraction)
         alvdf_ocn   , & ! visible, diffuse  (fraction)
         alidf_ocn   , & ! near-ir, diffuse  (fraction)
         Uref_ocn    , & ! 2m reference wind speed (m/s)
         Tref_ocn    , & ! 2m atm reference temperature (K)
         Qref_ocn        ! 2m atm reference spec humidity (kg/kg)

      !-----------------------------------------------------------------
      ! diagnostic
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         fsurf , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop,&! top surface conductive flux        (W/m^2)
         congel, & ! basal ice growth         (m/step-->cm/day)
         frazil, & ! frazil ice growth        (m/step-->cm/day)
         snoice, & ! snow-ice formation       (m/step-->cm/day)
         meltt , & ! top ice melt             (m/step-->cm/day)
         melts , & ! snow melt                (m/step-->cm/day)
         meltb , & ! basal ice melt           (m/step-->cm/day)
         meltl , & ! lateral ice melt         (m/step-->cm/day)
         daidtt, & ! ice area tendency thermo.   (s^-1)
         dvidtt, & ! ice volume tendency thermo. (m/s)
         mlt_onset, &! day of year that sfc melting begins
         frz_onset   ! day of year that freezing begins (congel or frazil)


      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         fsurfn,   & ! category fsurf
         fcondtopn,& ! category fcondtop
         flatn       ! cagegory latent heat flux

      ! As above but these remain grid box mean values i.e. they are not
      ! divided by aice at end of ice_dynamics.  These are used in
      ! CICE_IN_NEMO for coupling and also for generating
      ! ice diagnostics and history files as these are more accurate.
      ! (The others suffer from problem of incorrect values at grid boxes
      !  that change from an ice free state to an icy state.)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         fresh_gbm, & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_gbm, & ! salt flux to ocean (kg/m^2/s)
         fhocn_gbm, & ! net heat flux to ocean (W/m^2)
         fswthru_gbm  ! shortwave penetrating to ocean (W/m^2)

      !-----------------------------------------------------------------
      ! internal
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         rside   , & ! fraction of ice that melts laterally
         fsw     , & ! incoming shortwave radiation (W/m^2)
         coszen  , & ! cosine solar zenith angle, < 0 for sun below horizon 
         rdg_conv, & ! convergence term for ridging (1/s)
         rdg_shear   ! shear term for ridging (1/s)

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: init_coupler_flux - initialize fluxes exchanged with coupler
!
! !INTERFACE:
!
      subroutine init_coupler_flux
!
! !DESCRIPTION:
!
! Initialize all fluxes exchanged with flux coupler
! and some data-derived fields
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: i, j, n, iblk

      logical (kind=log_kind), parameter ::     &
!         l_winter = .true.   ! winter/summer default switch
         l_winter = .false.   ! winter/summer default switch

      real (kind=dbl_kind) :: fcondtopn_d(6), fsurfn_d(6)

      data fcondtopn_d / -50.0_dbl_kind,-17.0_dbl_kind,-12.0_dbl_kind, &
                          -9.0_dbl_kind, -7.0_dbl_kind, -3.0_dbl_kind /
      data fsurfn_d    /  0.20_dbl_kind, 0.15_dbl_kind, 0.10_dbl_kind, &
                          0.05_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind /

      !-----------------------------------------------------------------
      ! fluxes received from atmosphere
      !-----------------------------------------------------------------
      zlvl  (:,:,:) = c10             ! atm level height (m)
      rhoa  (:,:,:) = 1.3_dbl_kind    ! air density (kg/m^3)
      uatm  (:,:,:) = c5              ! wind velocity    (m/s)
      vatm  (:,:,:) = c5
      strax (:,:,:) = 0.05_dbl_kind
      stray (:,:,:) = 0.05_dbl_kind
      if (l_winter) then
         !typical winter values
         potT  (:,:,:) = 253.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 253.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0006_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         flw   (:,:,:) = c180            ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         fsnow (:,:,:) = 4.0e-6_dbl_kind ! snowfall rate (kg/m2/s)
         do n = 1, ncat                  ! conductive heat flux (W/m^2)
            fcondtopn_f(:,:,n,:) = fcondtopn_d(n)
         enddo
         fsurfn_f = fcondtopn_f          ! surface heat flux (W/m^2)
         flatn_f(:,:,:,:) = c0           ! latent heat flux (kg/m2/s)
      else
         !typical summer values
         potT  (:,:,:) = 273.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 273.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0035_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         flw   (:,:,:) = 280.0_dbl_kind  ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         fsnow (:,:,:) = c0              ! snowfall rate (kg/m2/s)
         do n = 1, ncat                  ! surface heat flux (W/m^2)
            fsurfn_f(:,:,n,:) = fsurfn_d(n)
         enddo
         fcondtopn_f(:,:,:,:) = 0.0_dbl_kind ! conductive heat flux (W/m^2)
         flatn_f(:,:,:,:) = -2.0_dbl_kind    ! latent heat flux (W/m^2)
      endif !     l_winter

      faero (:,:,:,:) = c0            ! aerosol deposition rate (kg/m2/s) MH

      !-----------------------------------------------------------------
      ! fluxes received from ocean
      !-----------------------------------------------------------------

      ss_tltx(:,:,:)= c0              ! sea surface tilt (m/m)
      ss_tlty(:,:,:)= c0
      uocn  (:,:,:) = c0              ! surface ocean currents (m/s)
      vocn  (:,:,:) = c0
      frzmlt(:,:,:) = c0              ! freezing/melting potential (W/m^2)
      sss   (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (o/oo)
      if (trim(Tfrzpt) == 'constant') then
         Tf    (:,:,:) = -1.8_dbl_kind   ! freezing temp (C)
      else ! default:  Tfrzpt = 'linear_S'
         Tf    (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
      endif
#ifndef CICE_IN_NEMO
      sst   (:,:,:) = Tf(:,:,:)       ! sea surface temp (C)
#endif
      qdp   (:,:,:) = c0              ! deep ocean heat flux (W/m^2)
      hmix  (:,:,:) = c20             ! ocean mixed layer depth

      !-----------------------------------------------------------------
      ! fluxes sent to atmosphere
      !-----------------------------------------------------------------

!echmod - for rectangular grid tests without thermo
!      strairxT(:,:,:) = 0.15_dbl_kind
!      strairyT(:,:,:) = 0.15_dbl_kind
      strairxT(:,:,:) = c0            ! wind stress, T grid
      strairyT(:,:,:) = c0
!echmod
      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
      flwout  (:,:,:) = -stefan_boltzmann*Tffresh**4   
                        ! in case atm model diagnoses Tsfc from flwout
      evap    (:,:,:) = c0
      Uref    (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0
      alvdr   (:,:,:) = c0
      alidr   (:,:,:) = c0
      alvdf   (:,:,:) = c0
      alidf   (:,:,:) = c0
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
      fswsfc  (:,:,:) = c0
      fswint  (:,:,:) = c0
#endif
#ifdef AEROFRC
      dfswabs_noaero  (:,:,:) = c0
      dfswsfc_noaero  (:,:,:) = c0
      dfswint_noaero  (:,:,:) = c0
      dfswthru_noaero  (:,:,:) = c0
      dalvdr_noaero   (:,:,:) = c0
      dalidr_noaero   (:,:,:) = c0
      dalvdf_noaero   (:,:,:) = c0
      dalidf_noaero   (:,:,:) = c0
#endif
#ifdef CCSM3FRC
      dfswabs_ccsm3  (:,:,:) = c0
      dfswsfc_ccsm3  (:,:,:) = c0
      dfswint_ccsm3  (:,:,:) = c0
      dfswthru_ccsm3  (:,:,:) = c0
      dalvdr_ccsm3   (:,:,:) = c0
      dalidr_ccsm3   (:,:,:) = c0
      dalvdf_ccsm3   (:,:,:) = c0
      dalidf_ccsm3   (:,:,:) = c0
#endif
#ifdef PONDFRC
      dfswabs_nopond  (:,:,:) = c0
      dfswsfc_nopond  (:,:,:) = c0
      dfswint_nopond  (:,:,:) = c0
      dfswthru_nopond  (:,:,:) = c0
      dalvdr_nopond   (:,:,:) = c0
      dalidr_nopond   (:,:,:) = c0
      dalvdf_nopond   (:,:,:) = c0
      dalidf_nopond   (:,:,:) = c0
#endif

      !-----------------------------------------------------------------
      ! fluxes sent to ocean
      !-----------------------------------------------------------------

      strocnxT(:,:,:) = c0    ! ice-ocean stress, x-direction (T-cell)
      strocnyT(:,:,:) = c0    ! ice-ocean stress, y-direction (T-cell)
      fresh   (:,:,:) = c0
      fsalt   (:,:,:) = c0
      fhocn   (:,:,:) = c0
      fswthru (:,:,:) = c0

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

      fsw     (:,:,:) = c0            ! shortwave radiation (W/m^2)
      scale_factor(:,:,:) = c0        ! shortwave scaling factor
      wind    (:,:,:) = sqrt(uatm(:,:,:)**2 &
                           + vatm(:,:,:)**2)  ! wind speed, (m/s)

      coszen  (:,:,:) = c0

      end subroutine init_coupler_flux

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_atm - initialize all atmospheric fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_atm
!
! !DESCRIPTION:
!
! Initialize some fluxes sent to coupler for use by the atm model
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      !-----------------------------------------------------------------
      ! initialize albedo and fluxes
      !-----------------------------------------------------------------

      strairxT(:,:,:) = c0      ! wind stress, T grid
      strairyT(:,:,:) = c0
      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
      fswsfc  (:,:,:) = c0
      fswint  (:,:,:) = c0
#endif
#ifdef AEROFRC
      dfswabs_noaero  (:,:,:) = c0
      dfswsfc_noaero  (:,:,:) = c0
      dfswint_noaero  (:,:,:) = c0
      dfswthru_noaero  (:,:,:) = c0
#endif
#ifdef CCSM3FRC
      dfswabs_ccsm3  (:,:,:) = c0
      dfswsfc_ccsm3  (:,:,:) = c0
      dfswint_ccsm3  (:,:,:) = c0
      dfswthru_ccsm3  (:,:,:) = c0
#endif
#ifdef PONDFRC
      dfswabs_nopond  (:,:,:) = c0
      dfswsfc_nopond  (:,:,:) = c0
      dfswint_nopond  (:,:,:) = c0
      dfswthru_nopond  (:,:,:) = c0
#endif
      flwout  (:,:,:) = c0
      evap    (:,:,:) = c0
      Uref    (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0

      end subroutine init_flux_atm

!=======================================================================
!BOP
!
! !IROUTINE: init_flux_ocn - initialize ocean fluxes sent to coupler
!
! !INTERFACE:
!
      subroutine init_flux_ocn
!
! !DESCRIPTION:
!
! Initialize some fluxes sent to coupler for use by the ocean model
!
! NOTE: These fluxes should be initialized immediately after the
!       call to the coupler.  The atmospheric fluxes can be initialized
!       at the beginning of the following time step because they are
!       not modified by any subroutines between the call to_coupler
!       and the end of the time step.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------

      fresh  (:,:,:)  = c0
      fsalt  (:,:,:)  = c0
      fhocn  (:,:,:)  = c0
      fsoot  (:,:,:,:)  = c0
      fswthru(:,:,:)  = c0

      end subroutine init_flux_ocn

!=======================================================================
!BOP
!
! !IROUTINE: init_history_therm - initialize thermo history fields
!
! !INTERFACE:
!
      subroutine init_history_therm
!
! !DESCRIPTION:
!
! Initialize thermodynamic fields written to history files.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      use ice_state, only: aice, vice

      fsurf  (:,:,:) = c0
      fcondtop(:,:,:)= c0
      congel (:,:,:) = c0
      frazil (:,:,:) = c0
      snoice (:,:,:) = c0
      meltt  (:,:,:) = c0
      melts  (:,:,:) = c0
      meltb  (:,:,:) = c0
      meltl  (:,:,:) = c0
      daidtt (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtt (:,:,:) = vice(:,:,:) ! temporary initial volume
      fsurfn    (:,:,:,:) = c0
      fcondtopn (:,:,:,:) = c0
      flatn     (:,:,:,:) = c0
      fresh_gbm  (:,:,:) = c0
      fsalt_gbm  (:,:,:) = c0
      fhocn_gbm  (:,:,:) = c0
      fswthru_gbm(:,:,:) = c0
      albice (:,:,:) = c0
      albsno (:,:,:) = c0
      albpnd (:,:,:) = c0
#ifdef AEROFRC
      dalbice_noaero (:,:,:) = c0
      dalbsno_noaero (:,:,:) = c0
      dalbpnd_noaero (:,:,:) = c0
#endif
#ifdef CCSM3FRC
      dalbice_ccsm3 (:,:,:) = c0
      dalbsno_ccsm3 (:,:,:) = c0
#endif
#ifdef PONDFRC
      dalbice_nopond (:,:,:) = c0
      dalbsno_nopond (:,:,:) = c0
      dalbpnd_nopond (:,:,:) = c0
#endif
     
      end subroutine init_history_therm

!=======================================================================
!BOP
!
! !IROUTINE: init_history_dyn - initialize dynamic history fields
!
! !INTERFACE:
!
      subroutine init_history_dyn
!
! !DESCRIPTION:
!
! Initialize dynamic fields written to history files.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_state, only: aice, vice
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      sig1    (:,:,:) = c0
      sig2    (:,:,:) = c0
      strocnx (:,:,:) = c0
      strocny (:,:,:) = c0
      strairx (:,:,:) = c0
      strairy (:,:,:) = c0
      strtltx (:,:,:) = c0
      strtlty (:,:,:) = c0
      strintx (:,:,:) = c0
      strinty (:,:,:) = c0
      dardg1dt(:,:,:) = c0
      dardg2dt(:,:,:) = c0
      dvirdgdt(:,:,:) = c0
      opening (:,:,:) = c0
      daidtd  (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtd  (:,:,:) = vice(:,:,:) ! temporary initial volume
      fm      (:,:,:) = c0
      prs_sig (:,:,:) = c0

      end subroutine init_history_dyn

!=======================================================================
!BOP
!
! !IROUTINE: merge_fluxes - aggregate flux information over ITD
!
! !INTERFACE:
!
      subroutine merge_fluxes (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               aicen,                &    
                               flw,      coszn,      &
                               strairxn, strairyn,   &
                               fsurfn,   fcondtopn,  &
                               fsensn,   flatn,      & 
                               fswabsn,  flwoutn,    &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
                               fswsfcn,  fswintn,    &
#endif
#ifdef AEROFRC
                               dfswabsn_noaero,  dfswsfcn_noaero,    &
                               dfswintn_noaero,  dfswthrun_noaero,    &
#endif
#ifdef CCSM3FRC
                               dfswabsn_ccsm3,  dfswsfcn_ccsm3,    &
                               dfswintn_ccsm3,  dfswthrun_ccsm3,    &
#endif
#ifdef PONDFRC
                               dfswabsn_nopond,  dfswsfcn_nopond,    &
                               dfswintn_nopond,  dfswthrun_nopond,    &
#endif
                               evapn,                &
                               Urefn,                &
                               Trefn,    Qrefn,      &
                               freshn,   fsaltn,     &
                               fhocnn,   fswthrun,   &
                               strairxT, strairyT,   &  
                               fsurf,    fcondtop,   &
                               fsens,    flat,       & 
                               fswabs,   flwout,     &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
                               fswsfc,   fswint,     &
#endif
#ifdef AEROFRC
                               dfswabs_noaero, dfswsfc_noaero, &
                               dfswint_noaero, dfswthru_noaero, &
#endif
#ifdef CCSM3FRC
                               dfswabs_ccsm3, dfswsfc_ccsm3, &
                               dfswint_ccsm3, dfswthru_ccsm3, &
#endif
#ifdef PONDFRC
                               dfswabs_nopond, dfswsfc_nopond, &
                               dfswint_nopond, dfswthru_nopond, &
#endif
                               evap,                 & 
                               Uref,                 &
                               Tref,     Qref,       &
                               fresh,    fsalt,      &
                               fhocn,    fswthru,    &
                               melttn, meltsn, meltbn, congeln, snoicen, &
                               meltt,  melts,  &
                               meltb,                       &
                               congel,  snoice)

!
! !DESCRIPTION:
!
! Aggregate flux information from all ice thickness categories
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj    ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &
          aicen   , & ! concentration of ice
          flw     , & ! downward longwave flux          (W/m**2)
          coszn   , & ! cosine of solar zenith angle 
          strairxn, & ! air/ice zonal  strss,           (N/m**2)
          strairyn, & ! air/ice merdnl strss,           (N/m**2)
          fsurfn  , & ! net heat flux to top surface    (W/m**2)
          fcondtopn,& ! downward cond flux at top sfc   (W/m**2)
          fsensn  , & ! sensible heat flx               (W/m**2)
          flatn   , & ! latent   heat flx               (W/m**2)
          fswabsn , & ! shortwave absorbed heat flx     (W/m**2)
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
          fswsfcn , & ! shortwave surface absorbed heat flx   (W/m**2)
          fswintn , & ! shortwave internal absorbed heat flx  (W/m**2)
#endif
#ifdef AEROFRC
          dfswabsn_noaero, & ! shortwave absorbed heat flx     (W/m**2)
          dfswsfcn_noaero, & ! shortwave absorbed heat flx     (W/m**2)
          dfswintn_noaero, & ! shortwave absorbed heat flx     (W/m**2)
          dfswthrun_noaero, & ! shortwave absorbed heat flx     (W/m**2)
#endif
#ifdef CCSM3FRC
          dfswabsn_ccsm3, & ! shortwave absorbed heat flx     (W/m**2)
          dfswsfcn_ccsm3, & ! shortwave absorbed heat flx     (W/m**2)
          dfswintn_ccsm3, & ! shortwave absorbed heat flx     (W/m**2)
          dfswthrun_ccsm3, & ! shortwave absorbed heat flx     (W/m**2)
#endif
#ifdef PONDFRC
          dfswabsn_nopond, & ! shortwave absorbed heat flx     (W/m**2)
          dfswsfcn_nopond, & ! shortwave absorbed heat flx     (W/m**2)
          dfswintn_nopond, & ! shortwave absorbed heat flx     (W/m**2)
          dfswthrun_nopond, & ! shortwave absorbed heat flx     (W/m**2)
#endif
          flwoutn , & ! upwd lw emitted heat flx        (W/m**2)
          evapn   , & ! evaporation                     (kg/m2/s)
          Urefn   , & ! wind speed reference level  (m/s)
          Trefn   , & ! air tmp reference level         (K)
          Qrefn   , & ! air sp hum reference level      (kg/kg)
          freshn  , & ! fresh water flux to ocean       (kg/m2/s)
          fsaltn  , & ! salt flux to ocean              (kg/m2/s)
          fhocnn  , & ! actual ocn/ice heat flx         (W/m**2)
          fswthrun, & ! sw radiation through ice bot    (W/m**2)
          melttn  , & ! top ice melt                    (m)
          meltbn  , & ! bottom ice melt                 (m)
          meltsn  , & ! snow melt                       (m)
          congeln , & ! congelation ice growth          (m)
          snoicen     ! snow-ice growth                 (m)

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          strairxT, & ! air/ice zonal  strss,           (N/m**2)
          strairyT, & ! air/ice merdnl strss,           (N/m**2)
          fsurf   , & ! net heat flux to top surface    (W/m**2)
          fcondtop, & ! downward cond flux at top sfc   (W/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
          fswsfc  , & ! shortwave surface absorbed heat flx         (W/m**2)
          fswint  , & ! shortwave internal absorbed heat flx        (W/m**2)
#endif
#ifdef AEROFRC
          dfswabs_noaero  , & ! shortwave absorbed heat flx          (W/m**2)
          dfswsfc_noaero  , & ! shortwave surface absorbed heat flx  (W/m**2)
          dfswint_noaero  , & ! shortwave internal absorbed heat flx (W/m**2)
          dfswthru_noaero  , & ! shortwave penetrating heat flux     (W/m**2)
#endif
#ifdef CCSM3FRC
          dfswabs_ccsm3  , & ! shortwave absorbed heat flx          (W/m**2)
          dfswsfc_ccsm3  , & ! shortwave surface absorbed heat flx  (W/m**2)
          dfswint_ccsm3  , & ! shortwave internal absorbed heat flx (W/m**2)
          dfswthru_ccsm3  , & ! shortwave penetrating heat flux     (W/m**2)
#endif
#ifdef PONDFRC
          dfswabs_nopond  , & ! shortwave absorbed heat flx          (W/m**2)
          dfswsfc_nopond  , & ! shortwave surface absorbed heat flx  (W/m**2)
          dfswint_nopond  , & ! shortwave internal absorbed heat flx (W/m**2)
          dfswthru_nopond  , & ! shortwave penetrating heat flux     (W/m**2)
#endif
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Uref    , & ! wind speed reference level      (m/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          meltt   , & ! top ice melt                    (m)
          meltb   , & ! bottom ice melt                 (m)
          melts   , & ! snow melt                       (m)
          congel  , & ! congelation ice growth          (m)
          snoice      ! snow-ice growth                 (m)

!
!EOP
!
      integer (kind=int_kind) :: &
          ij, i, j    ! horizontal indices

      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! atmo fluxes

         strairxT (i,j)  = strairxT(i,j) + strairxn(i,j)*aicen(i,j)
         strairyT (i,j)  = strairyT(i,j) + strairyn(i,j)*aicen(i,j)
         fsurf    (i,j)  = fsurf   (i,j) + fsurfn  (i,j)*aicen(i,j)
         fcondtop (i,j)  = fcondtop(i,j) + fcondtopn(i,j)*aicen(i,j)
         fsens    (i,j)  = fsens   (i,j) + fsensn  (i,j)*aicen(i,j)
         flat     (i,j)  = flat    (i,j) + flatn   (i,j)*aicen(i,j)
         fswabs   (i,j)  = fswabs  (i,j) + fswabsn (i,j)*aicen(i,j)
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
         fswint(i,j)=fswint(i,j) + fswintn(i,j)*aicen(i,j)
         fswsfc(i,j)=fswsfc(i,j) + fswsfcn(i,j)*aicen(i,j)
#endif
#ifdef AEROFRC
         dfswsfc_noaero(i,j)=dfswsfc_noaero(i,j) + dfswsfcn_noaero(i,j)*aicen(i,j)
         dfswabs_noaero(i,j)=dfswabs_noaero(i,j) + dfswabsn_noaero(i,j)*aicen(i,j)
         dfswint_noaero(i,j)=dfswint_noaero(i,j) + dfswintn_noaero(i,j)*aicen(i,j)
         dfswthru_noaero(i,j)=dfswthru_noaero(i,j)+dfswthrun_noaero(i,j)*aicen(i,j)
#endif
#ifdef CCSM3FRC
         dfswsfc_ccsm3(i,j)=dfswsfc_ccsm3(i,j) + dfswsfcn_ccsm3(i,j)*aicen(i,j)
         dfswabs_ccsm3(i,j)=dfswabs_ccsm3(i,j) + dfswabsn_ccsm3(i,j)*aicen(i,j)
         dfswint_ccsm3(i,j)=dfswint_ccsm3(i,j) + dfswintn_ccsm3(i,j)*aicen(i,j)
         dfswthru_ccsm3(i,j)=dfswthru_ccsm3(i,j)+dfswthrun_ccsm3(i,j)*aicen(i,j)
#endif
#ifdef PONDFRC
         dfswsfc_nopond(i,j)=dfswsfc_nopond(i,j) + dfswsfcn_nopond(i,j)*aicen(i,j)
         dfswabs_nopond(i,j)=dfswabs_nopond(i,j) + dfswabsn_nopond(i,j)*aicen(i,j)
         dfswint_nopond(i,j)=dfswint_nopond(i,j) + dfswintn_nopond(i,j)*aicen(i,j)
         dfswthru_nopond(i,j)=dfswthru_nopond(i,j)+dfswthrun_nopond(i,j)*aicen(i,j)
#endif
         flwout   (i,j)  = flwout  (i,j) &
             + (flwoutn(i,j) - (c1-emissivity)*flw(i,j))*aicen(i,j)
         evap     (i,j)  = evap    (i,j) + evapn   (i,j)*aicen(i,j)
         Uref     (i,j)  = Uref    (i,j) + Urefn   (i,j)*aicen(i,j)
         Tref     (i,j)  = Tref    (i,j) + Trefn   (i,j)*aicen(i,j)
         Qref     (i,j)  = Qref    (i,j) + Qrefn   (i,j)*aicen(i,j)

         ! ocean fluxes

         fresh    (i,j) = fresh     (i,j) + freshn  (i,j)*aicen(i,j)
         fsalt    (i,j) = fsalt     (i,j) + fsaltn  (i,j)*aicen(i,j)
         fhocn    (i,j) = fhocn     (i,j) + fhocnn  (i,j)*aicen(i,j)
         fswthru  (i,j) = fswthru   (i,j) + fswthrun(i,j)*aicen(i,j)

         ! ice/snow thickness

         meltt    (i,j) = meltt    (i,j) + melttn  (i,j)*aicen(i,j)
         meltb    (i,j) = meltb    (i,j) + meltbn  (i,j)*aicen(i,j)
         melts    (i,j) = melts    (i,j) + meltsn  (i,j)*aicen(i,j)
         congel   (i,j) = congel   (i,j) + congeln (i,j)*aicen(i,j)
         snoice   (i,j) = snoice   (i,j) + snoicen (i,j)*aicen(i,j)

      enddo                     ! ij
      
      end subroutine merge_fluxes

!=======================================================================
!BOP
!
! !IROUTINE: scale_fluxes
!
! !DESCRIPTION:
!
!  Divide ice fluxes by ice area before sending them to the
!  coupler, since the coupler multiplies by ice area. This
!  is the ice area at the beginning of the timestep, i.e.
!  coupler, since the coupler multiplies by ice area.
!
! !INTERFACE:
!
      subroutine scale_fluxes (nx_block, ny_block, &
                               tmask,              &
                               aice,     Tf,       &
                               Tair,     Qa,       &
                               wind,               &
                               strairxT, strairyT, &
                               fsens,    flat,     &
                               fswabs,   flwout,   &
                               evap,               &
                               Uref,               &
                               Tref,     Qref,     &
                               fresh,    fsalt,    &
                               fhocn,    fswthru,  &
                               fsoot,              &
                               alvdr,    alidr,    &
                               alvdf,    alidf)
!
! !REVISION HISTORY:
!
! authors: C.M.Bitz, William H. Lipscomb
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block   ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask     ! land/boundary mask, thickness (T-cell)


      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice    , & ! fractional ice area
          Tf      , & ! freezing temperature            (C)
          wind    , & ! wind speed                      (m/s)
          Tair    , & ! surface air temperature         (K)
          Qa          ! sfc air specific humidity       (kg/kg)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          strairxT, & ! air/ice zonal  stress           (N/m**2)
          strairyT, & ! air/ice merdnl stress           (N/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Uref    , & ! wind speed reference level      (m/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          alvdr   , & ! visible, direct   (fraction)
          alidr   , & ! near-ir, direct   (fraction)
          alvdf   , & ! visible, diffuse  (fraction)
          alidf       ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx), &
          intent(inout):: &
          fsoot       ! 
!
!EOP
!
      real (kind=dbl_kind) :: ar   ! 1/aice

      integer (kind=int_kind) :: &
          i, j    ! horizontal indices


!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j) .and. aice(i,j) > c0) then
            ar = c1 / aice(i,j)
            strairxT(i,j) = strairxT(i,j) * ar
            strairyT(i,j) = strairyT(i,j) * ar
            fsens   (i,j) = fsens   (i,j) * ar
            flat    (i,j) = flat    (i,j) * ar
            fswabs  (i,j) = fswabs  (i,j) * ar
            flwout  (i,j) = flwout  (i,j) * ar
            evap    (i,j) = evap    (i,j) * ar
            Uref    (i,j) = Uref    (i,j) * ar
            Tref    (i,j) = Tref    (i,j) * ar
            Qref    (i,j) = Qref    (i,j) * ar
            fresh   (i,j) = fresh   (i,j) * ar
            fsalt   (i,j) = fsalt   (i,j) * ar
            fhocn   (i,j) = fhocn   (i,j) * ar
            fswthru (i,j) = fswthru (i,j) * ar
            alvdr   (i,j) = alvdr   (i,j) * ar
            alidr   (i,j) = alidr   (i,j) * ar
            alvdf   (i,j) = alvdf   (i,j) * ar
            alidf   (i,j) = alidf   (i,j) * ar
            fsoot   (i,j,:) = fsoot (i,j,:) * ar
         else                   ! zero out fluxes
            strairxT(i,j) = c0
            strairyT(i,j) = c0
            fsens   (i,j) = c0
            flat    (i,j) = c0
            fswabs  (i,j) = c0
            flwout  (i,j) = -stefan_boltzmann *(Tf(i,j) + Tffresh)**4
               ! to make upward longwave over ocean reasonable for history file
            evap    (i,j) = c0
            Uref    (i,j) = wind(i,j)
            Tref    (i,j) = Tair(i,j)
            Qref    (i,j) = Qa  (i,j)
            fresh   (i,j) = c0
            fsalt   (i,j) = c0
            fhocn   (i,j) = c0
            fsoot   (i,j,:) = c0
            fswthru (i,j) = c0
            alvdr   (i,j) = c0  ! zero out albedo where ice is absent
            alidr   (i,j) = c0
            alvdf   (i,j) = c0 
            alidf   (i,j) = c0
         endif                  ! tmask and aice > 0
      enddo                     ! i
      enddo                     ! j
      
      end subroutine scale_fluxes

!=======================================================================

      end module ice_flux

!=======================================================================
