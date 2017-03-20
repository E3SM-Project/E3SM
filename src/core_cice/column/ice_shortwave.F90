!  SVN:$Id: ice_shortwave.F90 1182 2017-03-16 19:29:26Z njeffery $
!=======================================================================
!
! The albedo and absorbed/transmitted flux parameterizations for
! snow over ice, bare ice and ponded ice.
!
! Presently, two methods are included:
!   (1) CCSM3 
!   (2) Delta-Eddington 
! as two distinct routines.
! Either can be called from the ice driver.
!
! The Delta-Eddington method is described here:
!
! Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple 
!    Scattering Parameterization for Solar Radiation in the Sea Ice 
!    Component of the Community Climate System Model, NCAR Technical 
!    Note  NCAR/TN-472+STR  February 2007
!
! name: originally ice_albedo
!
! authors:  Bruce P. Briegleb, NCAR
!           Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2005, WHL: Moved absorbed_solar from ice_therm_vertical to this 
!            module and changed name from ice_albedo
! 2006, WHL: Added Delta Eddington routines from Bruce Briegleb
! 2006, ECH: Changed data statements in Delta Eddington routines (no 
!            longer hardwired)
!            Converted to free source form (F90)
! 2007, BPB: Completely updated Delta-Eddington code, so that:
!            (1) multiple snow layers enabled (i.e. nslyr > 1)
!            (2) included SSL for snow surface absorption
!            (3) added Sswabs for internal snow layer absorption
!            (4) variable sea ice layers allowed (i.e. not hardwired)
!            (5) updated all inherent optical properties
!            (6) included algae absorption for sea ice lowest layer
!            (7) very complete internal documentation included
! 2007, ECH: Improved efficiency
! 2008, BPB: Added aerosols to Delta Eddington code
! 2013, ECH: merged with NCAR version, cleaned up

      module ice_shortwave

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c1p5, c2, c3, c4, c10, &
          p01, p1, p15, p25, p5, p75, puny, &
          albocn, Timelt, snowpatch, awtvdr, awtidr, awtvdf, awtidf, &
          kappav, hs_min, rhofresh, rhos, nspint
      use ice_colpkg_shared, only: hi_ssl, hs_ssl, modal_aero
      use ice_warnings, only: add_warning

      implicit none

      private
      public :: run_dEdd, shortwave_ccsm3, compute_shortwave_trcr

      real (kind=dbl_kind), parameter :: &
         hpmin  = 0.005_dbl_kind, & ! minimum allowed melt pond depth (m)
         hp0    = 0.200_dbl_kind    ! pond depth below which transition to bare ice

      real (kind=dbl_kind) :: &
         exp_min                    ! minimum exponential value

!=======================================================================

      contains

!=======================================================================
!
! Driver for basic solar radiation from CCSM3.  Albedos and absorbed solar.

      subroutine shortwave_ccsm3 (aicen,    vicen,    &
                                  vsnon,    Tsfcn,    &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  heat_capacity,      &
                                  albedo_type,        &
                                  albicev,  albicei,  &
                                  albsnowv, albsnowi, &
                                  ahmax,              &
                                  alvdrn,   alidrn,   &
                                  alvdfn,   alidfn,   &
                                  fswsfc,   fswint,   &
                                  fswthru,  fswpenl,  &
                                  Iswabs,   SSwabs,   &
                                  albin,    albsn,    &
                                  coszen,   ncat)

      integer (kind=int_kind), intent(in) :: &
         ncat         ! number of ice thickness categories

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen    , & ! concentration of ice per category
         vicen    , & ! volume of ice per category
         vsnon    , & ! volume of ice per category
         Tsfcn        ! surface temperature

      real (kind=dbl_kind), intent(in) :: &
         swvdr    , & ! sw down, visible, direct  (W/m^2)
         swvdf    , & ! sw down, visible, diffuse (W/m^2)
         swidr    , & ! sw down, near IR, direct  (W/m^2)
         swidf        ! sw down, near IR, diffuse (W/m^2)

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), intent(in) :: &
         albicev , & ! visible ice albedo for h > ahmax
         albicei , & ! near-ir ice albedo for h > ahmax
         albsnowv, & ! cold snow albedo, visible
         albsnowi, & ! cold snow albedo, near IR
         ahmax       ! thickness above which ice albedo is constant (m)

      logical(kind=log_kind), intent(in) :: &
         heat_capacity! if true, ice has nonzero heat capacity

      character (len=char_len), intent(in) :: &
         albedo_type  ! albedo parameterization, 'default' ('ccsm3') or 'constant'

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         fswsfc   , & ! SW absorbed at ice/snow surface (W m-2)
         fswint   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthru  , & ! SW through ice to ocean (W m-2)
         albin    , & ! bare ice albedo
         albsn        ! snow albedo

      real (kind=dbl_kind), intent(inout) :: &
         coszen       ! cosine(zenith angle)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         fswpenl  , & ! SW entering ice layers (W m-2)
         Iswabs   , & ! SW absorbed in particular layer (W m-2)
         Sswabs       ! SW absorbed in particular layer (W m-2)

      ! local variables

      integer (kind=int_kind) :: &
         n                  ! thickness category index

      ! ice and snow albedo for each category

      real (kind=dbl_kind) :: &
         alvdrni, & ! visible, direct, ice    (fraction)
         alidrni, & ! near-ir, direct, ice    (fraction)
         alvdfni, & ! visible, diffuse, ice   (fraction)
         alidfni, & ! near-ir, diffuse, ice   (fraction)
         alvdrns, & ! visible, direct, snow   (fraction)
         alidrns, & ! near-ir, direct, snow   (fraction)
         alvdfns, & ! visible, diffuse, snow  (fraction)
         alidfns    ! near-ir, diffuse, snow  (fraction)

      !-----------------------------------------------------------------
      ! Solar radiation: albedo and absorbed shortwave
      !-----------------------------------------------------------------

      ! For basic shortwave, set coszen to a constant between 0 and 1.
      coszen = p5 ! sun above the horizon

      do n = 1, ncat

         Sswabs(:,n) = c0

      alvdrni = albocn
      alidrni = albocn
      alvdfni = albocn
      alidfni = albocn
      
      alvdrns = albocn
      alidrns = albocn
      alvdfns = albocn
      alidfns = albocn
      
      alvdrn(n) = albocn
      alidrn(n) = albocn
      alvdfn(n) = albocn
      alidfn(n) = albocn
      
      albin(n) = c0
      albsn(n) = c0    

      fswsfc(n)  = c0
      fswint(n)  = c0
      fswthru(n) = c0
      fswpenl(:,n)  = c0
      Iswabs (:,n) = c0

      if (aicen(n) > puny) then

      !-----------------------------------------------------------------
      ! Compute albedos for ice and snow.
      !-----------------------------------------------------------------

         if (trim(albedo_type) == 'constant') then

            call constant_albedos (aicen(n),             &
                                   vsnon(n),             &
                                   Tsfcn(n),             &
                                   alvdrni,    alidrni,  &
                                   alvdfni,    alidfni,  &
                                   alvdrns,    alidrns,  &
                                   alvdfns,    alidfns,  &
                                   alvdrn(n),            &
                                   alidrn(n),            &
                                   alvdfn(n),            &
                                   alidfn(n),            &
                                   albin(n),             &
                                   albsn(n))
         else ! default

            call compute_albedos (aicen(n),             &
                                  vicen(n),             &
                                  vsnon(n),             &
                                  Tsfcn(n),             &
                                  albicev,    albicei,  &
                                  albsnowv,   albsnowi, &
                                  ahmax,                &
                                  alvdrni,    alidrni,  &
                                  alvdfni,    alidfni,  &
                                  alvdrns,    alidrns,  &
                                  alvdfns,    alidfns,  &
                                  alvdrn(n),            &
                                  alidrn(n),            &
                                  alvdfn(n),            &
                                  alidfn(n),            &
                                  albin(n),             &
                                  albsn(n))
         endif

      !-----------------------------------------------------------------
      ! Compute solar radiation absorbed in ice and penetrating to ocean.
      !-----------------------------------------------------------------

         call absorbed_solar  (heat_capacity,        &
                               ncat,                 &
                               aicen(n),             &
                               vicen(n),             &
                               vsnon(n),             &
                               swvdr,      swvdf,    &
                               swidr,      swidf,    &
                               alvdrni,    alvdfni,  &
                               alidrni,    alidfni,  &
                               alvdrns,    alvdfns,  &
                               alidrns,    alidfns,  &
                               fswsfc(n),            &
                               fswint(n),            &
                               fswthru(n),           &
                               fswpenl(:,n),         &
                               Iswabs(:,n))

      endif ! aicen > puny

      enddo                  ! ncat

      end subroutine shortwave_ccsm3

!=======================================================================
!
! Compute albedos for each thickness category

      subroutine compute_albedos (aicen,    vicen,    &
                                  vsnon,    Tsfcn,    &
                                  albicev,  albicei,  &
                                  albsnowv, albsnowi, &
                                  ahmax,              &
                                  alvdrni,  alidrni,  &
                                  alvdfni,  alidfni,  &
                                  alvdrns,  alidrns,  &
                                  alvdfns,  alidfns,  &
                                  alvdrn,   alidrn,   &
                                  alvdfn,   alidfn,   &
                                  albin,    albsn)

      real (kind=dbl_kind), intent(in) :: &
         aicen   , & ! concentration of ice per category
         vicen   , & ! volume of ice per category
         vsnon   , & ! volume of ice per category
         Tsfcn       ! surface temperature

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), intent(in) :: &
         albicev , & ! visible ice albedo for h > ahmax
         albicei , & ! near-ir ice albedo for h > ahmax
         albsnowv, & ! cold snow albedo, visible
         albsnowi, & ! cold snow albedo, near IR
         ahmax       ! thickness above which ice albedo is constant (m)

      real (kind=dbl_kind), intent(out) :: &
         alvdrni  , & ! visible, direct, ice   (fraction)
         alidrni  , & ! near-ir, direct, ice   (fraction)
         alvdfni  , & ! visible, diffuse, ice  (fraction)
         alidfni  , & ! near-ir, diffuse, ice  (fraction)
         alvdrns  , & ! visible, direct, snow  (fraction)
         alidrns  , & ! near-ir, direct, snow  (fraction)
         alvdfns  , & ! visible, diffuse, snow (fraction)
         alidfns  , & ! near-ir, diffuse, snow (fraction)
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         albin    , & ! bare ice 
         albsn        ! snow 

      ! local variables

      real (kind=dbl_kind), parameter :: &
         dT_melt    = c1          , & ! change in temp to give dalb_mlt 
                                     ! albedo change
         dalb_mlt  = -0.075_dbl_kind, & ! albedo change per dT_melt change
                                     ! in temp for ice
         dalb_mltv = -p1         , & ! albedo vis change per dT_melt change
                                     ! in temp for snow
         dalb_mlti = -p15            ! albedo nir change per dT_melt change
                                     ! in temp for snow

      real (kind=dbl_kind) :: &
         hi  , & ! ice thickness  (m)
         hs  , & ! snow thickness  (m)
         albo, & ! effective ocean albedo, function of ice thickness
         fh  , & ! piecewise linear function of thickness
         fT  , & ! piecewise linear function of surface temperature
         dTs , & ! difference of Tsfc and Timelt
         fhtan,& ! factor used in albedo dependence on ice thickness
         asnow   ! fractional area of snow cover

      fhtan = atan(ahmax*c4)

      !-----------------------------------------------------------------
      ! Compute albedo for each thickness category.
      !-----------------------------------------------------------------

         hi = vicen / aicen
         hs = vsnon / aicen            

         ! bare ice, thickness dependence
         fh = min(atan(hi*c4)/fhtan,c1)
         albo = albocn*(c1-fh)
         alvdfni = albicev*fh + albo
         alidfni = albicei*fh + albo

         ! bare ice, temperature dependence
         dTs = Timelt - Tsfcn
         fT = min(dTs/dT_melt-c1,c0)
         alvdfni = alvdfni - dalb_mlt*fT
         alidfni = alidfni - dalb_mlt*fT

         ! avoid negative albedos for thin, bare, melting ice
         alvdfni = max (alvdfni, albocn)
         alidfni = max (alidfni, albocn)

         if (hs > puny) then

            alvdfns = albsnowv
            alidfns = albsnowi

            ! snow on ice, temperature dependence
            alvdfns = alvdfns - dalb_mltv*fT
            alidfns = alidfns - dalb_mlti*fT

         endif                  ! hs > puny

         ! direct albedos (same as diffuse for now)
         alvdrni = alvdfni
         alidrni = alidfni
         alvdrns = alvdfns
         alidrns = alidfns

         ! fractional area of snow cover
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

         ! combine ice and snow albedos (for coupler)
         alvdfn = alvdfni*(c1-asnow) + &
                  alvdfns*asnow
         alidfn = alidfni*(c1-asnow) + &
                  alidfns*asnow
         alvdrn = alvdrni*(c1-asnow) + &
                  alvdrns*asnow
         alidrn = alidrni*(c1-asnow) + &
                  alidrns*asnow

         ! save ice and snow albedos (for history)
         albin = awtvdr*alvdrni + awtidr*alidrni &
               + awtvdf*alvdfni + awtidf*alidfni 
         albsn = awtvdr*alvdrns + awtidr*alidrns &
               + awtvdf*alvdfns + awtidf*alidfns 

      end subroutine compute_albedos

!=======================================================================
!
! Compute albedos for each thickness category

      subroutine constant_albedos (aicen,              &
                                   vsnon,    Tsfcn,    &
                                   alvdrni,  alidrni,  &
                                   alvdfni,  alidfni,  &
                                   alvdrns,  alidrns,  &
                                   alvdfns,  alidfns,  &
                                   alvdrn,   alidrn,   &
                                   alvdfn,   alidfn,   &
                                   albin,    albsn)

      real (kind=dbl_kind), intent(in) :: &
         aicen   , & ! concentration of ice per category
         vsnon   , & ! volume of ice per category
         Tsfcn       ! surface temperature

      real (kind=dbl_kind), intent(out) :: &
         alvdrni  , & ! visible, direct, ice   (fraction)
         alidrni  , & ! near-ir, direct, ice   (fraction)
         alvdfni  , & ! visible, diffuse, ice  (fraction)
         alidfni  , & ! near-ir, diffuse, ice  (fraction)
         alvdrns  , & ! visible, direct, snow  (fraction)
         alidrns  , & ! near-ir, direct, snow  (fraction)
         alvdfns  , & ! visible, diffuse, snow (fraction)
         alidfns  , & ! near-ir, diffuse, snow (fraction)
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         albin    , & ! bare ice 
         albsn        ! snow 

      ! local variables

      real (kind=dbl_kind), parameter :: &
         warmice  = 0.68_dbl_kind, &
         coldice  = 0.70_dbl_kind, &
         warmsnow = 0.77_dbl_kind, &
         coldsnow = 0.81_dbl_kind

      real (kind=dbl_kind) :: &
         hs      ! snow thickness  (m)

      !-----------------------------------------------------------------
      ! Compute albedo for each thickness category.
      !-----------------------------------------------------------------

         hs = vsnon / aicen

         if (hs > puny) then
            ! snow, temperature dependence
            if (Tsfcn >= -c2*puny) then
               alvdfn = warmsnow
               alidfn = warmsnow
            else
               alvdfn = coldsnow
               alidfn = coldsnow
            endif
         else      ! hs < puny
            ! bare ice, temperature dependence
            if (Tsfcn >= -c2*puny) then
               alvdfn = warmice
               alidfn = warmice
            else
               alvdfn = coldice
               alidfn = coldice
            endif
         endif                  ! hs > puny

         ! direct albedos (same as diffuse for now)
         alvdrn  = alvdfn
         alidrn  = alidfn

         alvdrni = alvdrn
         alidrni = alidrn
         alvdrns = alvdrn
         alidrns = alidrn
         alvdfni = alvdfn
         alidfni = alidfn
         alvdfns = alvdfn
         alidfns = alidfn

         ! save ice and snow albedos (for history)
         albin = awtvdr*alvdrni + awtidr*alidrni &
               + awtvdf*alvdfni + awtidf*alidfni 
         albsn = awtvdr*alvdrns + awtidr*alidrns &
               + awtvdf*alvdfns + awtidf*alidfns 

      end subroutine constant_albedos

!=======================================================================
!
! Compute solar radiation absorbed in ice and penetrating to ocean
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine absorbed_solar (heat_capacity,      &
                                 nilyr,    aicen,    &
                                 vicen,    vsnon,    &
                                 swvdr,    swvdf,    &
                                 swidr,    swidf,    &
                                 alvdrni,  alvdfni,  &
                                 alidrni,  alidfni,  &
                                 alvdrns,  alvdfns,  &
                                 alidrns,  alidfns,  &
                                 fswsfc,   fswint,   &
                                 fswthru,  fswpenl,  &
                                 Iswabs)

      logical(kind=log_kind), intent(in) :: &
         heat_capacity   ! if true, ice has nonzero heat capacity

      integer (kind=int_kind), intent(in) :: & 
         nilyr           ! number of ice layers

      real (kind=dbl_kind), intent(in) :: &
         aicen       , & ! fractional ice area
         vicen       , & ! ice volume
         vsnon       , & ! snow volume
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         alvdrni     , & ! visible, direct albedo,ice
         alidrni     , & ! near-ir, direct albedo,ice
         alvdfni     , & ! visible, diffuse albedo,ice
         alidfni     , & ! near-ir, diffuse albedo,ice
         alvdrns     , & ! visible, direct albedo, snow
         alidrns     , & ! near-ir, direct albedo, snow
         alvdfns     , & ! visible, diffuse albedo, snow
         alidfns         ! near-ir, diffuse albedo, snow

      real (kind=dbl_kind), intent(out):: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthru         ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         Iswabs      , & ! SW absorbed in particular layer (W m-2)
         fswpenl         ! visible SW entering ice layers (W m-2)

      ! local variables

      real (kind=dbl_kind), parameter :: &
         i0vis = 0.70_dbl_kind  ! fraction of penetrating solar rad (visible)

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         fswpen      , & ! SW penetrating beneath surface (W m-2)
         trantop     , & ! transmitted frac of penetrating SW at layer top
         tranbot         ! transmitted frac of penetrating SW at layer bot

      real (kind=dbl_kind) :: &
         swabs       , & ! net SW down at surface (W m-2)
         swabsv      , & ! swabs in vis (wvlngth < 700nm)  (W/m^2)
         swabsi      , & ! swabs in nir (wvlngth > 700nm)  (W/m^2)
         fswpenvdr   , & ! penetrating SW, vis direct
         fswpenvdf   , & ! penetrating SW, vis diffuse
         hi          , & ! ice thickness (m)
         hs          , & ! snow thickness (m)
         hilyr       , & ! ice layer thickness
         asnow           ! fractional area of snow cover

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      trantop = c0
      tranbot = c0

         hs  = vsnon / aicen

      !-----------------------------------------------------------------
      ! Fractional snow cover
      !-----------------------------------------------------------------
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

      !-----------------------------------------------------------------
      ! Shortwave flux absorbed at surface, absorbed internally,
      !  and penetrating to mixed layer.
      ! This parameterization assumes that all IR is absorbed at the
      !  surface; only visible is absorbed in the ice interior or
      !  transmitted to the ocean.
      !-----------------------------------------------------------------

         swabsv  = swvdr * ( (c1-alvdrni)*(c1-asnow) &
                           + (c1-alvdrns)*asnow ) &
                 + swvdf * ( (c1-alvdfni)*(c1-asnow) &
                           + (c1-alvdfns)*asnow )

         swabsi  = swidr * ( (c1-alidrni)*(c1-asnow) &
                           + (c1-alidrns)*asnow ) &
                 + swidf * ( (c1-alidfni)*(c1-asnow) &
                           + (c1-alidfns)*asnow )

         swabs   = swabsv + swabsi

         fswpenvdr = swvdr * (c1-alvdrni) * (c1-asnow) * i0vis
         fswpenvdf = swvdf * (c1-alvdfni) * (c1-asnow) * i0vis

          ! no penetrating radiation in near IR
!         fswpenidr = swidr * (c1-alidrni) * (c1-asnow) * i0nir
!         fswpenidf = swidf * (c1-alidfni) * (c1-asnow) * i0nir  

         fswpen = fswpenvdr + fswpenvdf
                      
         fswsfc = swabs - fswpen

         trantop = c1  ! transmittance at top of ice

      !-----------------------------------------------------------------
      ! penetrating SW absorbed in each ice layer
      !-----------------------------------------------------------------

         do k = 1, nilyr

            hi  = vicen / aicen
            hilyr = hi / real(nilyr,kind=dbl_kind)

            tranbot = exp (-kappav * hilyr * real(k,kind=dbl_kind))
            Iswabs(k) = fswpen * (trantop-tranbot)

            ! bottom of layer k = top of layer k+1
            trantop = tranbot

            ! bgc layer model
            if (k == 1) then   ! surface flux
               fswpenl(k)   = fswpen
               fswpenl(k+1) = fswpen * tranbot
            else
               fswpenl(k+1) = fswpen * tranbot
            endif
         enddo                     ! nilyr

         ! SW penetrating thru ice into ocean
         fswthru = fswpen * tranbot

         ! SW absorbed in ice interior
         fswint  = fswpen - fswthru

      !----------------------------------------------------------------
      ! if zero-layer model (no heat capacity), no SW is absorbed in ice
      ! interior, so add to surface absorption
      !----------------------------------------------------------------
         
         if (.not. heat_capacity) then

            ! SW absorbed at snow/ice surface
            fswsfc = fswsfc + fswint

            ! SW absorbed in ice interior (nilyr = 1)
            fswint    = c0
            Iswabs(1) = c0

         endif                       ! heat_capacity

      end subroutine absorbed_solar

! End ccsm3 shortwave method
!=======================================================================
! Begin Delta-Eddington shortwave method

! Compute initial data for Delta-Eddington method, specifically, 
! the approximate exponential look-up table.
!
! author:  Bruce P. Briegleb, NCAR
! 2011 ECH modified for melt pond tracers
! 2013 ECH merged with NCAR version

      subroutine run_dEdd(dt,       tr_aero,   &
                          tr_pond_cesm,        &
                          tr_pond_lvl,         &
                          tr_pond_topo,        &
                          ncat,     n_aero,    &
                          n_zaero,  dEdd_algae,&
                          nlt_chl_sw,          &
                          nlt_zaero_sw,        &
                          tr_bgc_N, tr_zaero,  &
                          nilyr,    nslyr,     &
                          aicen,    vicen,     &
                          vsnon,    Tsfcn,     &
                          alvln,    apndn,     &
                          hpndn,    ipndn,     &
                          aeron,    kalg,      &
                          zbion,               &
                          heat_capacity,       &
                          tlat,     tlon,      & 
                          calendar_type,       &
                          days_per_year,       &
                          nextsw_cday,   yday, &
                          sec,      R_ice,     &
                          R_pnd,    R_snw,     &
                          dT_mlt,   rsnw_mlt,  &
                          hs0,      hs1,  hp1, &
                          pndaspect,           &
                          kaer_tab, waer_tab,  &
                          gaer_tab,            &
                          kaer_bc_tab,         &
                          waer_bc_tab,         &
                          gaer_bc_tab,         &
                          bcenh,               &
                          modal_aero,          &
                          swvdr,    swvdf,     &
                          swidr,    swidf,     &
                          coszen,   fsnow,     &
                          alvdrn,   alvdfn,    &
                          alidrn,   alidfn,    &
                          fswsfcn,  fswintn,   &
                          fswthrun, fswpenln,  &
                          Sswabsn,  Iswabsn,   &
                          albicen,  albsnon,   &
                          albpndn,  apeffn,    &
                          snowfracn,           &
                          dhsn,     ffracn,    &
                          l_print_point,       &
                          initonly)

      use ice_orbital, only: compute_coszen

      integer (kind=int_kind), intent(in) :: &
         ncat   , & ! number of ice thickness categories
         nilyr  , & ! number of ice layers
         nslyr  , & ! number of snow layers
         n_aero , & ! number of aerosol tracers
         n_zaero, & ! number of zaerosol tracers 
         nlt_chl_sw ! index for chla

      integer (kind=int_kind), dimension(:), intent(in) :: &
        nlt_zaero_sw   ! index for zaerosols

      logical(kind=log_kind), intent(in) :: &
         heat_capacity,& ! if true, ice has nonzero heat capacity
         tr_aero     , & ! if .true., use aerosol tracers
         tr_pond_cesm, & ! if .true., use explicit topography-based ponds
         tr_pond_lvl , & ! if .true., use explicit topography-based ponds
         tr_pond_topo, & ! if .true., use explicit topography-based ponds
         dEdd_algae,   & ! .true. use prognostic chla in dEdd
         tr_bgc_N,     & ! .true. active bgc (skl or z)
         tr_zaero,     & ! .true. use zaerosols
         modal_aero      ! .true. use modal aerosol treatment

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), intent(in) :: &
         R_ice , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt, & ! change in temp for non-melt to melt snow grain radius change (C)
         rsnw_mlt, & ! maximum melting snow grain radius (10^-6 m)
         hs0      , & ! snow depth for transition to bare sea ice (m)
         pndaspect, & ! ratio of pond depth to pond fraction
         hs1      , & ! tapering parameter for snow on pond ice
         hp1      , & ! critical parameter for pond ice thickness
         kalg         ! algae absorption coefficient

      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), dimension(:,:), intent(in) :: & ! Modal aerosol treatment
         kaer_bc_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), dimension(:,:,:), intent(in) :: & ! Modal aerosol treatment
         bcenh          ! BC absorption enhancement factor

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=dbl_kind), intent(in) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real(kind=dbl_kind), intent(in) :: &
           dt,    & ! time step (s)
           tlat,  & ! latitude of temp pts (radians)
           tlon,  & ! longitude of temp pts (radians)
           swvdr, & ! sw down, visible, direct  (W/m^2)
           swvdf, & ! sw down, visible, diffuse (W/m^2)
           swidr, & ! sw down, near IR, direct  (W/m^2)
           swidf, & ! sw down, near IR, diffuse (W/m^2)
           fsnow    ! snowfall rate (kg/m^2 s)

      real(kind=dbl_kind), dimension(:), intent(in) :: &
           aicen, & ! concentration of ice
           vicen, & ! volume per unit area of ice (m)
           vsnon, & ! volume per unit area of snow (m)
           ffracn,& ! fraction of fsurfn used to melt ipond
           Tsfcn, & ! surface temperature (deg C)
           alvln, & ! level-ice area fraction
           apndn, & ! pond area fraction
           hpndn, & ! pond depth (m)
           ipndn    ! pond refrozen lid thickness (m)

      real(kind=dbl_kind), dimension(:,:), intent(in) :: &
           aeron, & ! aerosols (kg/m^3)
           zbion    ! zaerosols (kg/m^3) + chlorophyll on shorthwave grid

      real(kind=dbl_kind), dimension(:), intent(inout) :: &
           dhsn     ! depth difference for snow on sea ice and pond ice

      real(kind=dbl_kind), intent(inout) :: &
           coszen   ! cosine solar zenith angle, < 0 for sun below horizon 

      real(kind=dbl_kind), dimension(:), intent(inout) :: &
           alvdrn,   & ! visible direct albedo (fraction)
           alvdfn,   & ! near-ir direct albedo (fraction)
           alidrn,   & ! visible diffuse albedo (fraction)
           alidfn,   & ! near-ir diffuse albedo (fraction)
           fswsfcn,  & ! SW absorbed at ice/snow surface (W m-2)
           fswintn,  & ! SW absorbed in ice interior, below surface (W m-2)
           fswthrun, & ! SW through ice to ocean (W/m^2) 
           albicen,  & ! albedo bare ice 
           albsnon,  & ! albedo snow 
           albpndn,  & ! albedo pond 
           apeffn,   & ! effective pond area used for radiation calculation
           snowfracn   ! snow fraction on each category used for radiation

      real(kind=dbl_kind), dimension(:,:), intent(inout) :: &
           Sswabsn , & ! SW radiation absorbed in snow layers (W m-2)
           Iswabsn , & ! SW radiation absorbed in ice layers (W m-2) 
           fswpenln    ! visible SW entering ice layers (W m-2)

      logical (kind=log_kind), intent(in) :: &
           l_print_point

      logical (kind=log_kind), optional :: &
           initonly    ! flag to indicate init only, default is false

      ! local temporary variables

      ! other local variables
      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind) :: &
         fsn         , & ! snow horizontal fraction
         hsn             ! snow depth (m)

      real (kind=dbl_kind), dimension (nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micrometers)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind) :: &
         fpn         , & ! pond fraction of ice cover
         hpn             ! actual pond depth (m)

      integer (kind=int_kind) :: &
         n           , & ! thickness category index
         na              ! aerosol index               

      real (kind=dbl_kind) :: &
         ipn         , & ! refrozen pond ice thickness (m), mean over ice fraction
         hp          , & ! pond depth
         hs          , & ! snow depth
         asnow       , & ! fractional area of snow cover
         rp          , & ! volume fraction of retained melt water to total liquid content
         hmx         , & ! maximum available snow infiltration equivalent depth
         dhs         , & ! local difference in snow depth on sea ice and pond ice
         spn         , & ! snow depth on refrozen pond (m)
         tmp             ! 0 or 1

      logical (kind=log_kind) :: &
         linitonly       ! local initonly value

      real (kind=dbl_kind), parameter :: & 
         argmax = c10    ! maximum argument of exponential

      linitonly = .false.
      if (present(initonly)) then
         linitonly = initonly
      endif

      exp_min = exp(-argmax)

      ! cosine of the zenith angle
         call compute_coszen (tlat,          tlon, &
                              calendar_type, days_per_year, &
                              nextsw_cday,   yday,  sec, &
                              coszen,        dt)

      do n = 1, ncat

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

         ! set snow properties
         fsn        = c0
         hsn        = c0
         rhosnwn(:) = c0
         rsnwn(:)   = c0
         apeffn(n)    = c0 ! for history
         snowfracn(n) = c0 ! for history

         if (aicen(n) > puny) then

            call shortwave_dEdd_set_snow(nslyr,      R_snw,    &
                                         dT_mlt,     rsnw_mlt, &
                                         aicen(n),   vsnon(n), &
                                         Tsfcn(n),   fsn,      &
                                         hs0,        hsn,      &
                                         rhosnwn,    rsnwn)    

            ! set pond properties
            if (tr_pond_cesm) then
               ! fraction of ice area
               fpn = apndn(n)
               ! pond depth over fraction fpn 
               hpn = hpndn(n)
               ! snow infiltration
               if (hsn >= hs_min .and. hs0 > puny) then
                  asnow = min(hsn/hs0, c1) ! delta-Eddington formulation
                  fpn = (c1 - asnow) * fpn
                  hpn = pndaspect * fpn
               endif
               ! Zero out fraction of thin ponds for radiation only
               if (hpn < hpmin) fpn = c0
               fsn = min(fsn, c1-fpn)
               apeffn(n) = fpn ! for history
            elseif (tr_pond_lvl) then
               fpn = c0  ! fraction of ice covered in pond
               hpn = c0  ! pond depth over fpn
               ! refrozen pond lid thickness avg over ice
               ! allow snow to cover pond ice
               ipn = alvln(n) * apndn(n) * ipndn(n)
               dhs = dhsn(n) ! snow depth difference, sea ice - pond
               if (.not. linitonly .and. ipn > puny .and. &
                    dhs < puny .and. fsnow*dt > hs_min) &
                    dhs = hsn - fsnow*dt ! initialize dhs>0
               spn = hsn - dhs   ! snow depth on pond ice
               if (.not. linitonly .and. ipn*spn < puny) dhs = c0
               dhsn(n) = dhs ! save: constant until reset to 0
                  
               ! not using ipn assumes that lid ice is perfectly clear
               ! if (ipn <= 0.3_dbl_kind) then
               
               ! fraction of ice area
               fpn = apndn(n) * alvln(n) 
               ! pond depth over fraction fpn
               hpn = hpndn(n)
               
               ! reduce effective pond area absorbing surface heat flux
               ! due to flux already having been used to melt pond ice
               fpn = (c1 - ffracn(n)) * fpn
               
               ! taper pond area with snow on pond ice
               if (dhs > puny .and. spn >= puny .and. hs1 > puny) then
                  asnow = min(spn/hs1, c1)
                  fpn = (c1 - asnow) * fpn
               endif
               
               ! infiltrate snow
               hp = hpn
               if (hp > puny) then
                  hs = hsn
                  rp = rhofresh*hp/(rhofresh*hp + rhos*hs)
                  if (rp < p15) then
                     fpn = c0
                     hpn = c0
                  else
                     hmx = hs*(rhofresh - rhos)/rhofresh
                     tmp = max(c0, sign(c1, hp-hmx)) ! 1 if hp>=hmx, else 0
                     hp = (rhofresh*hp + rhos*hs*tmp) &
                          / (rhofresh    - rhos*(c1-tmp))
                     hsn = hs - hp*fpn*(c1-tmp)
                     hpn = hp * tmp
                     fpn = fpn * tmp
                  endif
               endif ! hp > puny

               ! Zero out fraction of thin ponds for radiation only
               if (hpn < hpmin) fpn = c0
               fsn = min(fsn, c1-fpn)

               ! endif    ! masking by lid ice
               apeffn(n) = fpn ! for history

            elseif (tr_pond_topo) then
               ! Lid effective if thicker than hp1
               if (apndn(n)*aicen(n) > puny .and. ipndn(n) < hp1) then
                  fpn = apndn(n)
               else
                  fpn = c0
               endif
               if (apndn(n) > puny) then
                  hpn = hpndn(n) 
               else
                  fpn = c0
                  hpn = c0
               endif

               ! Zero out fraction of thin ponds for radiation only
               if (hpn < hpmin) fpn = c0

               ! If ponds are present snow fraction reduced to
               ! non-ponded part dEdd scheme 
               fsn = min(fsn, c1-fpn)

               apeffn(n) = fpn
            else
               fpn = c0
               hpn = c0
               call shortwave_dEdd_set_pond(Tsfcn(n), &
                                            fsn, fpn,   &
                                            hpn)
               
               apeffn(n) = fpn ! for history
               fpn = c0
               hpn = c0
            endif ! pond type
            
         snowfracn(n) = fsn ! for history

         call shortwave_dEdd(n_aero,        n_zaero,        &
                             dEdd_algae,    nlt_chl_sw,     &
                             nlt_zaero_sw(:),               &
                             tr_bgc_N,      tr_zaero,       &
                             nslyr,         nilyr,          &
                             coszen,        heat_capacity,  &
                             aicen(n),      vicen(n),       &
                             hsn,           fsn,            &
                             rhosnwn,       rsnwn,          &
                             fpn,           hpn,            &
                             aeron(:,n),    tr_aero,        &
                             R_ice,         R_pnd,          &
                             kaer_tab,      waer_tab,       &
                             gaer_tab,                      &
                             kaer_bc_tab,                   &
                             waer_bc_tab,                   &
                             gaer_bc_tab,                   &
                             bcenh,         modal_aero,     &   
                             kalg,                          &
                             swvdr,         swvdf,          &
                             swidr,         swidf,          &
                             alvdrn(n),     alvdfn(n),      &
                             alidrn(n),     alidfn(n),      &
                             fswsfcn(n),    fswintn(n),     &
                             fswthrun(n),                   &
                             Sswabsn(:,n),                  &
                             Iswabsn(:,n),                  &
                             albicen(n),                    &
                             albsnon(n),    albpndn(n),     &
                             fswpenln(:,n), zbion(:,n),     &
                             l_print_point)

         endif ! aicen > puny

      enddo  ! ncat

      end subroutine run_dEdd
 
!=======================================================================
!
!   Compute snow/bare ice/ponded ice shortwave albedos, absorbed and transmitted 
!   flux using the Delta-Eddington solar radiation method as described in:
!
!   "A Delta-Eddington Multiple Scattering Parameterization for Solar Radiation
!        in the Sea Ice Component of the Community Climate System Model"
!            B.P.Briegleb and B.Light   NCAR/TN-472+STR  February 2007
!
!   Compute shortwave albedos and fluxes for three surface types: 
!   snow over ice, bare ice and ponded ice. 
!   
!   Albedos and fluxes are output for later use by thermodynamic routines. 
!   Invokes three calls to compute_dEdd, which sets inherent optical properties 
!   appropriate for the surface type. Within compute_dEdd, a call to solution_dEdd 
!   evaluates the Delta-Eddington solution. The final albedos and fluxes are then
!   evaluated in compute_dEdd. Albedos and fluxes are transferred to output in 
!   this routine.
!
!   NOTE regarding albedo diagnostics:  This method yields zero albedo values
!   if there is no incoming solar and thus the albedo diagnostics are masked
!   out when the sun is below the horizon.  To estimate albedo from the history 
!   output (post-processing), compute ice albedo using
!   (1 - albedo)*swdn = swabs. -ECH
!
! author:  Bruce P. Briegleb, NCAR 
!   2013:  E Hunke merged with NCAR version
!
      subroutine shortwave_dEdd  (n_aero,   n_zaero,     &
                                  dEdd_algae,            &
                                  nlt_chl_sw,            &
                                  nlt_zaero_sw,          &
                                  tr_bgc_N, tr_zaero,    &
                                  nslyr,    nilyr,       &
                                  coszen,   heat_capacity,&
                                  aice,     vice,        &
                                  hs,       fs,          & 
                                  rhosnw,   rsnw,        &
                                  fp,       hp,          &
                                  aero,     tr_aero,     &
                                  R_ice,    R_pnd,       &
                                  kaer_tab, waer_tab,    &
                                  gaer_tab,              &
                                  kaer_bc_tab,           &
                                  waer_bc_tab,           &
                                  gaer_bc_tab,           &
                                  bcenh,  modal_aero,    &   
                                  kalg,                  &
                                  swvdr,    swvdf,       &
                                  swidr,    swidf,       &
                                  alvdr,    alvdf,       &
                                  alidr,    alidf,       &
                                  fswsfc,   fswint,      &
                                  fswthru,  Sswabs,      &
                                  Iswabs,   albice,      &
                                  albsno,   albpnd,      &
                                  fswpenl,  zbio,        &
                                  l_print_point)

      integer (kind=int_kind), intent(in) :: &
         nilyr   , & ! number of ice layers
         nslyr   , & ! number of snow layers
         n_aero  , & ! number of aerosol tracers in use
         n_zaero , & ! number of zaerosol tracers in use
         nlt_chl_sw  ! index for chla

      integer (kind=int_kind), dimension(:), intent(in) :: &
        nlt_zaero_sw   ! index for zaerosols

      logical (kind=log_kind), intent(in) :: &
         heat_capacity, & ! if true, ice has nonzero heat capacity
         tr_aero,       & ! if .true., use aerosol tracers
         dEdd_algae,    & ! .true. use prognostic chla in dEdd
         tr_bgc_N,      & ! .true. active bgc (skl or z)
         tr_zaero,      & ! .true. use zaerosols
         modal_aero       ! .true. use modal aerosol treatment
 
      real (kind=dbl_kind), dimension(:,:), intent(in) :: & ! Modal aerosol treatment
         kaer_bc_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), dimension(:,:,:), intent(in) :: & ! Modal aerosol treatment
         bcenh          ! BC absorption enhancement factor

      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), intent(in) :: &
         kalg    , & ! algae absorption coefficient
         R_ice , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         aice    , & ! concentration of ice 
         vice    , & ! volume of ice 
         hs      , & ! snow depth
         fs          ! horizontal coverage of snow

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         rhosnw  , & ! density in snow layer (kg/m3)
         rsnw    , & ! grain radius in snow layer (m)
         aero    , & ! aerosol tracers
         zbio        ! shortwave tracers (zaero+chla)

      real (kind=dbl_kind), intent(in) :: &
         fp      , & ! pond fractional coverage (0 to 1) 
         hp      , & ! pond depth (m) 
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf       ! sw down, near IR, diffuse (W/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         coszen  , & ! cosine of solar zenith angle 
         alvdr   , & ! visible, direct, albedo (fraction) 
         alvdf   , & ! visible, diffuse, albedo (fraction) 
         alidr   , & ! near-ir, direct, albedo (fraction) 
         alidf   , & ! near-ir, diffuse, albedo (fraction) 
         fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
         fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
         fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)
 
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         fswpenl , & ! visible SW entering ice layers (W m-2)
         Sswabs  , & ! SW absorbed in snow layer (W m-2)
         Iswabs      ! SW absorbed in ice layer (W m-2)

      real (kind=dbl_kind), intent(out) :: &
         albice  , & ! bare ice albedo, for history  
         albsno  , & ! snow albedo, for history  
         albpnd      ! pond albedo, for history  

      logical (kind=log_kind) , intent(in) :: &
         l_print_point

      ! local variables

      real (kind=dbl_kind) :: &
         netsw    , & ! net shortwave
         fnidr    , & ! fraction of direct to total down surface flux in nir
         hstmp    , & ! snow thickness (set to 0 for bare ice case)
         hi       , & ! ice thickness (all sea ice layers, m)
         fi           ! snow/bare ice fractional coverage (0 to 1)

      real (kind=dbl_kind), dimension (4*n_aero) :: &
         aero_mp      ! aerosol mass path in kg/m2

      integer (kind=int_kind) :: &
         srftyp       ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      integer (kind=int_kind) :: &
         k        , & ! level index
         na       , & ! aerosol index
         klev     , & ! number of radiation layers - 1
         klevp        ! number of radiation interfaces - 1
                      ! (0 layer is included also)

      real (kind=dbl_kind) :: &
         vsno         ! volume of snow 

      ! for printing points
      integer (kind=int_kind) :: &
         n            ! point number for prints

      real (kind=dbl_kind) :: &
         swdn  , & ! swvdr(i,j)+swvdf(i,j)+swidr(i,j)+swidf(i,j)
         swab  , & ! fswsfc(i,j)+fswint(i,j)+fswthru(i,j)
         swalb     ! (1.-swab/(swdn+.0001))

      ! for history
      real (kind=dbl_kind) :: &
         avdrl   , & ! visible, direct, albedo (fraction) 
         avdfl   , & ! visible, diffuse, albedo (fraction) 
         aidrl   , & ! near-ir, direct, albedo (fraction) 
         aidfl       ! near-ir, diffuse, albedo (fraction) 

      character(len=char_len_long) :: &
         warning ! warning message
      
!-----------------------------------------------------------------------

         klev    = nslyr + nilyr + 1   ! number of radiation layers - 1
         klevp   = klev  + 1           ! number of radiation interfaces - 1
                                       ! (0 layer is included also)

      ! zero storage albedos and fluxes for accumulation over surface types:
      hstmp    = c0
      hi       = c0
      fi       = c0
      alvdr    = c0
      alvdf    = c0
      alidr    = c0
      alidf    = c0
      avdrl    = c0
      avdfl    = c0
      aidrl    = c0
      aidfl    = c0
      fswsfc   = c0
      fswint   = c0
      fswthru  = c0
      ! compute fraction of nir down direct to total over all points:
      fnidr = c0
      if( swidr + swidf > puny ) then
         fnidr = swidr/(swidr+swidf)
      endif
      albice    = c0
      albsno    = c0
      albpnd    = c0
      fswpenl(:) = c0
      Sswabs(:)  = c0
      Iswabs(:)  = c0

      ! compute aerosol mass path
      
         aero_mp(:) = c0
         if( tr_aero ) then
            ! assume 4 layers for each aerosol, a snow SSL, snow below SSL,
            ! sea ice SSL, and sea ice below SSL, in that order.
            do na = 1, 4*n_aero, 4
               vsno = hs * aice
               netsw = swvdr + swidr + swvdf + swidf
               if (netsw > puny) then ! sun above horizon
                  aero_mp(na  ) = aero(na  )*vsno
                  aero_mp(na+1) = aero(na+1)*vsno
                  aero_mp(na+2) = aero(na+2)*vice
                  aero_mp(na+3) = aero(na+3)*vice
               endif                  ! aice > 0 and netsw > 0
            enddo      ! na
         endif      ! if aerosols

         ! compute shortwave radiation accounting for snow/ice (both snow over 
         ! ice and bare ice) and ponded ice (if any):
         
         ! sea ice points with sun above horizon
         netsw = swvdr + swidr + swvdf + swidf
         if (netsw > puny) then ! sun above horizon
            coszen = max(puny,coszen)
            ! evaluate sea ice thickness and fraction
            hi  = vice / aice
            fi  = c1 - fs - fp
            ! bare sea ice points
            if(fi > c0) then
               ! calculate bare sea ice

               srftyp = 0
               call compute_dEdd(nilyr,       nslyr,   klev,   klevp,   & 
                      n_zaero,   zbio,        dEdd_algae,               &
                      nlt_chl_sw,nlt_zaero_sw,         tr_bgc_N,        &
                      tr_zaero,                                         &
                      heat_capacity,          fnidr,   coszen,          &
                      n_aero,    tr_aero,     R_ice,   R_pnd,           &
                      kaer_tab,  waer_tab,    gaer_tab,                 &
                      kaer_bc_tab, waer_bc_tab, gaer_bc_tab,            &
                      bcenh,     modal_aero,  kalg,                     &
                      swvdr,     swvdf,       swidr,   swidf,  srftyp,  &
                      hstmp,     rhosnw,      rsnw,    hi,     hp,      &
                      fi,        aero_mp,     avdrl,   avdfl,           &
                      aidrl,     aidfl,                                 &
                      fswsfc,    fswint,                                &
                      fswthru,   Sswabs,                                &
                      Iswabs,    fswpenl)
               
               alvdr   = alvdr   + avdrl *fi
               alvdf   = alvdf   + avdfl *fi
               alidr   = alidr   + aidrl *fi
               alidf   = alidf   + aidfl *fi
               ! for history
               albice = albice &
                      + awtvdr*avdrl + awtidr*aidrl &
                      + awtvdf*avdfl + awtidf*aidfl 
            endif
         endif
         
         ! sea ice points with sun above horizon
         netsw = swvdr + swidr + swvdf + swidf
         if (netsw > puny) then ! sun above horizon
            coszen = max(puny,coszen)
            ! snow-covered sea ice points
            if(fs > c0) then
               ! calculate snow covered sea ice

               srftyp = 1
               call compute_dEdd(nilyr,       nslyr,   klev,   klevp,   & 
                      n_zaero,   zbio,        dEdd_algae,               &
                      nlt_chl_sw,nlt_zaero_sw,         tr_bgc_N,        &
                      tr_zaero,                                         &
                      heat_capacity,          fnidr,   coszen,          &
                      n_aero,    tr_aero,     R_ice,   R_pnd,           &
                      kaer_tab,  waer_tab,    gaer_tab,                 &
                      kaer_bc_tab, waer_bc_tab, gaer_bc_tab,            &
                      bcenh,     modal_aero,  kalg,                     &
                      swvdr,     swvdf,       swidr,   swidf,  srftyp,  &
                      hs,        rhosnw,      rsnw,    hi,     hp,      &
                      fs,        aero_mp,     avdrl,   avdfl,           &
                      aidrl,     aidfl,                                 &
                      fswsfc,    fswint,                                &
                      fswthru,   Sswabs,                                &
                      Iswabs,    fswpenl)
               
               alvdr   = alvdr   + avdrl *fs
               alvdf   = alvdf   + avdfl *fs
               alidr   = alidr   + aidrl *fs
               alidf   = alidf   + aidfl *fs
               ! for history
               albsno = albsno &
                      + awtvdr*avdrl + awtidr*aidrl &
                      + awtvdf*avdfl + awtidf*aidfl 
            endif
         endif
         
         hi = c0

         ! sea ice points with sun above horizon
         netsw = swvdr + swidr + swvdf + swidf
         if (netsw > puny) then ! sun above horizon
            coszen = max(puny,coszen)
            hi  = vice / aice
            ! if nonzero pond fraction and sufficient pond depth
            ! if( fp > puny .and. hp > hpmin ) then
            if (fp > puny) then
               
               ! calculate ponded ice

               srftyp = 2
               call compute_dEdd(nilyr,       nslyr,   klev,   klevp,   & 
                      n_zaero,   zbio,        dEdd_algae,               &
                      nlt_chl_sw,nlt_zaero_sw,         tr_bgc_N,        &
                      tr_zaero,                                         &
                      heat_capacity,          fnidr,   coszen,          &
                      n_aero,    tr_aero,     R_ice,   R_pnd,           &
                      kaer_tab,  waer_tab,    gaer_tab,                 &
                      kaer_bc_tab, waer_bc_tab, gaer_bc_tab,            &
                      bcenh,     modal_aero,  kalg,                     &
                      swvdr,     swvdf,       swidr,   swidf,  srftyp,  & 
                      hs,        rhosnw,      rsnw,    hi,     hp,      &
                      fp,        aero_mp,     avdrl,   avdfl,           &
                      aidrl,     aidfl,                                 &
                      fswsfc,    fswint,                                &
                      fswthru,   Sswabs,                                &
                      Iswabs,    fswpenl)
               
               alvdr   = alvdr   + avdrl *fp
               alvdf   = alvdf   + avdfl *fp
               alidr   = alidr   + aidrl *fp
               alidf   = alidf   + aidfl *fp
               ! for history
               albpnd = albpnd &
                      + awtvdr*avdrl + awtidr*aidrl &
                      + awtvdf*avdfl + awtidf*aidfl 
            endif
         endif

         ! if no incoming shortwave, set albedos to 1
         netsw = swvdr + swidr + swvdf + swidf
         if (netsw <= puny) then ! sun above horizon
            alvdr = c1
            alvdf = c1
            alidr = c1
            alidf = c1
         endif

      if (l_print_point .and. netsw > puny) then

         write(warning,*) ' printing point = ',n
         call add_warning(warning)
         write(warning,*) ' coszen = ', &
                            coszen
         call add_warning(warning)
         write(warning,*) ' swvdr  swvdf = ', &
                            swvdr,swvdf
         call add_warning(warning)
         write(warning,*) ' swidr  swidf = ', &
                            swidr,swidf
         call add_warning(warning)
         write(warning,*) ' aice = ', &
                            aice
         call add_warning(warning)
         write(warning,*) ' hs = ', &
                            hs
         call add_warning(warning)
         write(warning,*) ' hp = ', &
                            hp
         call add_warning(warning)
         write(warning,*) ' fs = ', &
                            fs
         call add_warning(warning)
         write(warning,*) ' fi = ', &
                            fi
         call add_warning(warning)
         write(warning,*) ' fp = ', &
                            fp
         call add_warning(warning)
         write(warning,*) ' hi = ', &
                            hi
         call add_warning(warning)
         write(warning,*) ' alvdr  alvdf = ', &
                            alvdr,alvdf
         call add_warning(warning)
         write(warning,*) ' alidr  alidf = ', &
                            alidr,alidf
         call add_warning(warning)
         write(warning,*) ' fswsfc fswint fswthru = ', &
                            fswsfc,fswint,fswthru
         call add_warning(warning)
         swdn  = swvdr+swvdf+swidr+swidf
         swab  = fswsfc+fswint+fswthru
         swalb = (1.-swab/(swdn+.0001))
         write(warning,*) ' swdn swab swalb = ',swdn,swab,swalb
         do k = 1, nslyr               
            write(warning,*) ' snow layer k    = ', k, &
                             ' rhosnw = ', &
                               rhosnw(k), &
                             ' rsnw = ', &
                               rsnw(k)
            call add_warning(warning)
         enddo
         do k = 1, nslyr               
            write(warning,*) ' snow layer k    = ', k, &
                             ' Sswabs(k)       = ', Sswabs(k)
            call add_warning(warning)
         enddo
         do k = 1, nilyr               
            write(warning,*) ' sea ice layer k = ', k, &
                             ' Iswabs(k)       = ', Iswabs(k)
            call add_warning(warning)
         enddo

      endif  ! l_print_point .and. coszen > .01

      end subroutine shortwave_dEdd

!=======================================================================
!
! Evaluate snow/ice/ponded ice inherent optical properties (IOPs), and 
! then calculate the multiple scattering solution by calling solution_dEdd.
!
! author:  Bruce P. Briegleb, NCAR 
!   2013:  E Hunke merged with NCAR version

      subroutine compute_dEdd (nilyr,    nslyr,    klev,  klevp,  &
                    n_zaero,   zbio,     dEdd_algae,              &
                    nlt_chl_sw,nlt_zaero_sw,       tr_bgc_N,      &
                    tr_zaero,                                     &
                    heat_capacity,       fnidr,    coszen,        &
                    n_aero,    tr_aero,  R_ice,    R_pnd,         &
                    kaer_tab,  waer_tab,           gaer_tab,      &
                    kaer_bc_tab, waer_bc_tab,      gaer_bc_tab,   &
                    bcenh,     modal_aero,         kalg,          &
                    swvdr,     swvdf,    swidr,    swidf, srftyp, &
                    hs,        rhosnw,   rsnw,     hi,    hp,     &
                    fi,        aero_mp,  alvdr,    alvdf,         &
                               alidr,    alidf,         &
                               fswsfc,   fswint,        &
                               fswthru,  Sswabs,        &
                               Iswabs,   fswpenl)

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         n_zaero , & ! number of zaerosol tracers in use
         nlt_chl_sw, &! index for chla
         klev  , & ! number of radiation layers - 1
         klevp     ! number of radiation interfaces - 1
                   ! (0 layer is included also)
 
      integer (kind=int_kind), dimension(:), intent(in) :: &
        nlt_zaero_sw   ! index for zaerosols

      logical (kind=log_kind), intent(in) :: &
         heat_capacity,& ! if true, ice has nonzero heat capacity
         tr_aero,      & ! if .true., use aerosol tracers
         dEdd_algae,   & ! .true. use prognostic chla in dEdd
         tr_bgc_N,     & ! .true. active bgc (skl or z)
         tr_zaero,     & ! .true. use zaerosols
         modal_aero      ! .true. use modal aerosol treatment
 
      real (kind=dbl_kind), dimension(:,:), intent(in) :: & ! Modal aerosol treatment
         kaer_bc_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), dimension(:,:,:), intent(in) :: & ! Modal aerosol treatment
         bcenh          ! BC absorption enhancement factor

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), intent(in) :: &
         R_ice , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd     ! ponded ice tuning parameter; +1 > 1sig increase in albedo

      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), intent(in) :: &
         kalg    , & ! algae absorption coefficient
         fnidr   , & ! fraction of direct to total down flux in nir
         coszen  , & ! cosine solar zenith angle
         swvdr   , & ! shortwave down at surface, visible, direct  (W/m^2)
         swvdf   , & ! shortwave down at surface, visible, diffuse (W/m^2)
         swidr   , & ! shortwave down at surface, near IR, direct  (W/m^2)
         swidf       ! shortwave down at surface, near IR, diffuse (W/m^2)
 
      integer (kind=int_kind), intent(in) :: &
         srftyp      ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      real (kind=dbl_kind), intent(in) :: &
         hs          ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         rhosnw  , & ! snow density in snow layer (kg/m3)
         rsnw    , & ! snow grain radius in snow layer (m)
         zbio    , & ! zaerosol + chla shortwave tracers kg/m^3
         aero_mp     ! aerosol mass path in kg/m2

      real (kind=dbl_kind), intent(in) :: &
         hi      , & ! ice thickness (m)
         hp      , & ! pond depth (m)
         fi          ! snow/bare ice fractional coverage (0 to 1)
 
      real (kind=dbl_kind), intent(inout) :: &
         alvdr   , & ! visible, direct, albedo (fraction) 
         alvdf   , & ! visible, diffuse, albedo (fraction) 
         alidr   , & ! near-ir, direct, albedo (fraction) 
         alidf   , & ! near-ir, diffuse, albedo (fraction) 
         fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
         fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
         fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)
 
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         fswpenl , & ! visible SW entering ice layers (W m-2)
         Sswabs  , & ! SW absorbed in snow layer (W m-2)
         Iswabs      ! SW absorbed in ice layer (W m-2)

!-----------------------------------------------------------------------
!
! Set up optical property profiles, based on snow, sea ice and ponded 
! ice IOPs from:
!
! Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple 
!    Scattering Parameterization for Solar Radiation in the Sea Ice 
!    Component of the Community Climate System Model, NCAR Technical 
!    Note  NCAR/TN-472+STR  February 2007
!
! Computes column Delta-Eddington radiation solution for specific
! surface type: either snow over sea ice, bare sea ice, or ponded sea ice.
!
! Divides solar spectrum into 3 intervals: 0.2-0.7, 0.7-1.19, and
! 1.19-5.0 micro-meters. The latter two are added (using an assumed
! partition of incident shortwave in the 0.7-5.0 micro-meter band between
! the 0.7-1.19 and 1.19-5.0 micro-meter band) to give the final output 
! of 0.2-0.7 visible and 0.7-5.0 near-infrared albedos and fluxes.
!
! Specifies vertical layer optical properties based on input snow depth,
! density and grain radius, along with ice and pond depths, then computes
! layer by layer Delta-Eddington reflectivity, transmissivity and combines
! layers (done by calling routine solution_dEdd). Finally, surface albedos
! and internal fluxes/flux divergences are evaluated.
!
!  Description of the level and layer index conventions. This is
!  for the standard case of one snow layer and four sea ice layers.
!
!  Please read the following; otherwise, there is 99.9% chance you
!  will be confused about indices at some point in time........ :)
!
!  CICE4.0 snow treatment has one snow layer above the sea ice. This 
!  snow layer has finite heat capacity, so that surface absorption must
!  be distinguished from internal. The Delta-Eddington solar radiation
!  thus adds extra surface scattering layers to both snow and sea ice.
!  Note that in the following, we assume a fixed vertical layer structure
!  for the radiation calculation. In other words, we always have the 
!  structure shown below for one snow and four sea ice layers, but for 
!  ponded ice the pond fills "snow" layer 1 over the sea ice, and for 
!  bare sea ice the top layers over sea ice are treated as transparent air.
!
!  SSL = surface scattering layer for either snow or sea ice
!  DL  = drained layer for sea ice immediately under sea ice SSL
!  INT = interior layers for sea ice below the drained layer.
!
!  Notice that the radiation level starts with 0 at the top. Thus,
!  the total number radiation layers is klev+1, where klev is the
!  sum of nslyr, the number of CCSM snow layers, and nilyr, the
!  number of CCSM sea ice layers, plus the sea ice SSL:
!  klev = 1 + nslyr + nilyr
!
!  For the standard case illustrated below, nslyr=1, nilyr=4,
!  and klev=6, with the number of layer interfaces klevp=klev+1.
!  Layer interfaces are the surfaces on which reflectivities,
!  transmissivities and fluxes are evaluated.
!
!  CCSM3 Sea Ice Model            Delta-Eddington Solar Radiation
!                                     Layers and Interfaces
!                             Layer Index             Interface Index
!    ---------------------            ---------------------  0
!                                  0  \\\   snow SSL    \\\
!       snow layer 1                  ---------------------  1
!                                  1    rest of snow layer
!    +++++++++++++++++++++            +++++++++++++++++++++  2
!                                  2  \\\ sea ice SSL   \\\
!      sea ice layer 1                ---------------------  3
!                                  3      sea ice  DL
!    ---------------------            ---------------------  4
!
!      sea ice layer 2             4      sea ice INT
!
!    ---------------------            ---------------------  5
!
!      sea ice layer 3             5      sea ice INT
!
!    ---------------------            ---------------------  6
!
!      sea ice layer 4             6      sea ice INT
!
!    ---------------------            ---------------------  7
!
! When snow lies over sea ice, the radiation absorbed in the
! snow SSL is used for surface heating, and that in the rest
! of the snow layer for its internal heating. For sea ice in
! this case, all of the radiant heat absorbed in both the
! sea ice SSL and the DL are used for sea ice layer 1 heating.
!
! When pond lies over sea ice, and for bare sea ice, all of the
! radiant heat absorbed within and above the sea ice SSL is used
! for surface heating, and that absorbed in the sea ice DL is
! used for sea ice layer 1 heating.
!
! Basically, vertical profiles of the layer extinction optical depth (tau), 
! single scattering albedo (w0) and asymmetry parameter (g) are required over
! the klev+1 layers, where klev+1 = 2 + nslyr + nilyr. All of the surface type
! information and snow/ice iop properties are evaulated in this routine, so
! the tau,w0,g profiles can be passed to solution_dEdd for multiple scattering
! evaluation. Snow, bare ice and ponded ice iops are contained in data arrays
! in this routine.
!
!-----------------------------------------------------------------------

      ! local variables

      integer (kind=int_kind) :: &
         k       , & ! level index
         ns      , & ! spectral index
         nr      , & ! index for grain radius tables
         ki      , & ! index for internal absorption
         km      , & ! k starting index for snow, sea ice internal absorption
         kp      , & ! k+1 or k+2 index for snow, sea ice internal absorption
         ksrf    , & ! level index for surface absorption
         ksnow   , & ! level index for snow density and grain size
         kii         ! level starting index for sea ice (nslyr+1)

      integer (kind=int_kind), parameter :: & 
         nmbrad  = 32        ! number of snow grain radii in tables
 
      real (kind=dbl_kind) :: & 
         avdr    , & ! visible albedo, direct   (fraction)
         avdf    , & ! visible albedo, diffuse  (fraction)
         aidr    , & ! near-ir albedo, direct   (fraction)
         aidf        ! near-ir albedo, diffuse  (fraction)
 
      real (kind=dbl_kind) :: & 
         fsfc    , & ! shortwave absorbed at snow/bare ice/ponded ice surface (W m-2)
         fint    , & ! shortwave absorbed in interior (W m-2)
         fthru       ! shortwave through snow/bare ice/ponded ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(nslyr) :: & 
         Sabs        ! shortwave absorbed in snow layer (W m-2)

      real (kind=dbl_kind), dimension(nilyr) :: & 
         Iabs        ! shortwave absorbed in ice layer (W m-2)
 
      real (kind=dbl_kind), dimension(nilyr+1) :: & 
         fthrul      ! shortwave through to ice layers (W m-2)

      real (kind=dbl_kind), dimension (nspint) :: &
         wghtns              ! spectral weights
 
      real (kind=dbl_kind), parameter :: & 
         cp67    = 0.67_dbl_kind   , & ! nir band weight parameter
         cp33    = 0.33_dbl_kind   , & ! nir band weight parameter
         cp78    = 0.78_dbl_kind   , & ! nir band weight parameter
         cp22    = 0.22_dbl_kind   , & ! nir band weight parameter
         cp01    = 0.01_dbl_kind       ! for ocean visible albedo
 
      real (kind=dbl_kind), dimension (0:klev) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (0:klevp) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above
 
      real (kind=dbl_kind), dimension (0:klevp) :: &
         dfdir   , & ! down-up flux at interface due to direct beam at top surface
         dfdif       ! down-up flux at interface due to diffuse beam at top surface

      real (kind=dbl_kind) :: &
         refk    , & ! interface k multiple scattering term
         delr    , & ! snow grain radius interpolation parameter
      ! inherent optical properties (iop) for snow
         Qs      , & ! Snow extinction efficiency
         ks      , & ! Snow extinction coefficient (/m)
         ws      , & ! Snow single scattering albedo
         gs          ! Snow asymmetry parameter

      real (kind=dbl_kind), dimension(nslyr) :: & 
         frsnw       ! snow grain radius in snow layer * adjustment factor (m)

      ! actual used ice and ponded ice IOPs, allowing for tuning 
      ! modifications of the above "_mn" value
      real (kind=dbl_kind), dimension (nspint) :: &
         ki_ssl       , & ! Surface-scattering-layer ice extinction coefficient (/m)
         wi_ssl       , & ! Surface-scattering-layer ice single scattering albedo
         gi_ssl       , & ! Surface-scattering-layer ice asymmetry parameter
         ki_dl        , & ! Drained-layer ice extinction coefficient (/m)
         wi_dl        , & ! Drained-layer ice single scattering albedo
         gi_dl        , & ! Drained-layer ice asymmetry parameter
         ki_int       , & ! Interior-layer ice extinction coefficient (/m)
         wi_int       , & ! Interior-layer ice single scattering albedo
         gi_int       , & ! Interior-layer ice asymmetry parameter
         ki_p_ssl     , & ! Ice under pond srf scat layer extinction coefficient (/m)
         wi_p_ssl     , & ! Ice under pond srf scat layer single scattering albedo
         gi_p_ssl     , & ! Ice under pond srf scat layer asymmetry parameter
         ki_p_int     , & ! Ice under pond extinction coefficient (/m)
         wi_p_int     , & ! Ice under pond single scattering albedo
         gi_p_int         ! Ice under pond asymmetry parameter

      real (kind=dbl_kind), dimension(0:klev) :: &
         dzk              ! layer thickness

      real (kind=dbl_kind) :: &
         dz           , & ! snow, sea ice or pond water layer thickness
         dz_ssl       , & ! snow or sea ice surface scattering layer thickness
         fs               ! scaling factor to reduce (nilyr<4) or increase (nilyr>4) DL
                          ! extinction coefficient to maintain DL optical depth constant
                          ! with changing number of sea ice layers, to approximately 
                          ! conserve computed albedo for constant physical depth of sea
                          ! ice when the number of sea ice layers vary
      real (kind=dbl_kind) :: &
         sig          , & ! scattering coefficient for tuning
         kabs         , & ! absorption coefficient for tuning
         sigp             ! modified scattering coefficient for tuning

      real (kind=dbl_kind), dimension(nspint, 0:klev) :: &
         kabs_chl    , & ! absorption coefficient for chlorophyll (/m)
         tzaer       , & ! total aerosol extinction optical depth
         wzaer       , & ! total aerosol single scatter albedo
         gzaer           ! total aerosol asymmetry parameter

      real (kind=dbl_kind) :: &
         albodr       , & ! spectral ocean albedo to direct rad
         albodf           ! spectral ocean albedo to diffuse rad
      
      ! for melt pond transition to bare sea ice for small pond depths 
      real (kind=dbl_kind) :: &
         sig_i        , & ! ice scattering coefficient (/m)
         sig_p        , & ! pond scattering coefficient (/m)
         kext             ! weighted extinction coefficient (/m)

      ! aerosol optical properties from Mark Flanner, 26 June 2008
      ! order assumed: hydrophobic black carbon, hydrophilic black carbon,
      ! four dust aerosols by particle size range:
      ! dust1(.05-0.5 micron), dust2(0.5-1.25 micron),
      ! dust3(1.25-2.5 micron), dust4(2.5-5.0 micron)
      ! spectral bands same as snow/sea ice: (0.3-0.7 micron, 0.7-1.19 micron
      ! and 1.19-5.0 micron in wavelength)

      integer (kind=int_kind) :: &
         na , n                        ! aerosol index

      real (kind=dbl_kind) :: &
         taer                     , & ! total aerosol extinction optical depth
         waer                     , & ! total aerosol single scatter albedo
         gaer                     , & ! total aerosol asymmetry parameter
         swdr                     , & ! shortwave down at surface, direct  (W/m^2)
         swdf                     , & ! shortwave down at surface, diffuse (W/m^2)
         rnilyr                   , & ! real(nilyr)
         rnslyr                   , & ! real(nslyr)
         rns                      , & ! real(ns)
         tmp_0, tmp_ks, tmp_kl        ! temp variables

      integer(kind=int_kind), dimension(0:klev) :: &
         k_bcini        , & !
         k_bcins        , &
         k_bcexs
      real(kind=dbl_kind)::  &
          tmp_gs, tmp1               ! temp variables

      ! snow grain radii (micro-meters) for table
      real (kind=dbl_kind), dimension(nmbrad), parameter :: &
         rsnw_tab = (/ &   ! snow grain radius for each table entry (micro-meters)
          5._dbl_kind,    7._dbl_kind,   10._dbl_kind,   15._dbl_kind, &
         20._dbl_kind,   30._dbl_kind,   40._dbl_kind,   50._dbl_kind, &
         65._dbl_kind,   80._dbl_kind,  100._dbl_kind,  120._dbl_kind, &
        140._dbl_kind,  170._dbl_kind,  200._dbl_kind,  240._dbl_kind, &
        290._dbl_kind,  350._dbl_kind,  420._dbl_kind,  500._dbl_kind, &
        570._dbl_kind,  660._dbl_kind,  760._dbl_kind,  870._dbl_kind, &
       1000._dbl_kind, 1100._dbl_kind, 1250._dbl_kind, 1400._dbl_kind, &
       1600._dbl_kind, 1800._dbl_kind, 2000._dbl_kind, 2500._dbl_kind/)

      ! snow extinction efficiency (unitless)
      real (kind=dbl_kind), dimension (nspint,nmbrad), parameter :: &
         Qs_tab = reshape((/ &
          2.131798_dbl_kind,  2.187756_dbl_kind,  2.267358_dbl_kind, &
          2.104499_dbl_kind,  2.148345_dbl_kind,  2.236078_dbl_kind, &
          2.081580_dbl_kind,  2.116885_dbl_kind,  2.175067_dbl_kind, &
          2.062595_dbl_kind,  2.088937_dbl_kind,  2.130242_dbl_kind, &
          2.051403_dbl_kind,  2.072422_dbl_kind,  2.106610_dbl_kind, &
          2.039223_dbl_kind,  2.055389_dbl_kind,  2.080586_dbl_kind, &
          2.032383_dbl_kind,  2.045751_dbl_kind,  2.066394_dbl_kind, &
          2.027920_dbl_kind,  2.039388_dbl_kind,  2.057224_dbl_kind, &
          2.023444_dbl_kind,  2.033137_dbl_kind,  2.048055_dbl_kind, &
          2.020412_dbl_kind,  2.028840_dbl_kind,  2.041874_dbl_kind, &
          2.017608_dbl_kind,  2.024863_dbl_kind,  2.036046_dbl_kind, &
          2.015592_dbl_kind,  2.022021_dbl_kind,  2.031954_dbl_kind, &
          2.014083_dbl_kind,  2.019887_dbl_kind,  2.028853_dbl_kind, &
          2.012368_dbl_kind,  2.017471_dbl_kind,  2.025353_dbl_kind, &
          2.011092_dbl_kind,  2.015675_dbl_kind,  2.022759_dbl_kind, &
          2.009837_dbl_kind,  2.013897_dbl_kind,  2.020168_dbl_kind, &
          2.008668_dbl_kind,  2.012252_dbl_kind,  2.017781_dbl_kind, &
          2.007627_dbl_kind,  2.010813_dbl_kind,  2.015678_dbl_kind, &
          2.006764_dbl_kind,  2.009577_dbl_kind,  2.013880_dbl_kind, &
          2.006037_dbl_kind,  2.008520_dbl_kind,  2.012382_dbl_kind, &
          2.005528_dbl_kind,  2.007807_dbl_kind,  2.011307_dbl_kind, &
          2.005025_dbl_kind,  2.007079_dbl_kind,  2.010280_dbl_kind, &
          2.004562_dbl_kind,  2.006440_dbl_kind,  2.009333_dbl_kind, &
          2.004155_dbl_kind,  2.005898_dbl_kind,  2.008523_dbl_kind, &
          2.003794_dbl_kind,  2.005379_dbl_kind,  2.007795_dbl_kind, &
          2.003555_dbl_kind,  2.005041_dbl_kind,  2.007329_dbl_kind, &
          2.003264_dbl_kind,  2.004624_dbl_kind,  2.006729_dbl_kind, &
          2.003037_dbl_kind,  2.004291_dbl_kind,  2.006230_dbl_kind, &
          2.002776_dbl_kind,  2.003929_dbl_kind,  2.005700_dbl_kind, &
          2.002590_dbl_kind,  2.003627_dbl_kind,  2.005276_dbl_kind, &
          2.002395_dbl_kind,  2.003391_dbl_kind,  2.004904_dbl_kind, &
          2.002071_dbl_kind,  2.002922_dbl_kind,  2.004241_dbl_kind/), &
          (/nspint,nmbrad/))

      ! snow single scattering albedo (unitless)
      real (kind=dbl_kind), dimension (nspint,nmbrad), parameter :: &
        ws_tab = reshape((/ &
         0.9999994_dbl_kind,  0.9999673_dbl_kind,  0.9954589_dbl_kind, &
         0.9999992_dbl_kind,  0.9999547_dbl_kind,  0.9938576_dbl_kind, &
         0.9999990_dbl_kind,  0.9999382_dbl_kind,  0.9917989_dbl_kind, &
         0.9999985_dbl_kind,  0.9999123_dbl_kind,  0.9889724_dbl_kind, &
         0.9999979_dbl_kind,  0.9998844_dbl_kind,  0.9866190_dbl_kind, &
         0.9999970_dbl_kind,  0.9998317_dbl_kind,  0.9823021_dbl_kind, &
         0.9999960_dbl_kind,  0.9997800_dbl_kind,  0.9785269_dbl_kind, &
         0.9999951_dbl_kind,  0.9997288_dbl_kind,  0.9751601_dbl_kind, &
         0.9999936_dbl_kind,  0.9996531_dbl_kind,  0.9706974_dbl_kind, &
         0.9999922_dbl_kind,  0.9995783_dbl_kind,  0.9667577_dbl_kind, &
         0.9999903_dbl_kind,  0.9994798_dbl_kind,  0.9621007_dbl_kind, &
         0.9999885_dbl_kind,  0.9993825_dbl_kind,  0.9579541_dbl_kind, &
         0.9999866_dbl_kind,  0.9992862_dbl_kind,  0.9541924_dbl_kind, &
         0.9999838_dbl_kind,  0.9991434_dbl_kind,  0.9490959_dbl_kind, &
         0.9999810_dbl_kind,  0.9990025_dbl_kind,  0.9444940_dbl_kind, &
         0.9999772_dbl_kind,  0.9988171_dbl_kind,  0.9389141_dbl_kind, &
         0.9999726_dbl_kind,  0.9985890_dbl_kind,  0.9325819_dbl_kind, &
         0.9999670_dbl_kind,  0.9983199_dbl_kind,  0.9256405_dbl_kind, &
         0.9999605_dbl_kind,  0.9980117_dbl_kind,  0.9181533_dbl_kind, &
         0.9999530_dbl_kind,  0.9976663_dbl_kind,  0.9101540_dbl_kind, &
         0.9999465_dbl_kind,  0.9973693_dbl_kind,  0.9035031_dbl_kind, &
         0.9999382_dbl_kind,  0.9969939_dbl_kind,  0.8953134_dbl_kind, &
         0.9999289_dbl_kind,  0.9965848_dbl_kind,  0.8865789_dbl_kind, &
         0.9999188_dbl_kind,  0.9961434_dbl_kind,  0.8773350_dbl_kind, &
         0.9999068_dbl_kind,  0.9956323_dbl_kind,  0.8668233_dbl_kind, &
         0.9998975_dbl_kind,  0.9952464_dbl_kind,  0.8589990_dbl_kind, &
         0.9998837_dbl_kind,  0.9946782_dbl_kind,  0.8476493_dbl_kind, &
         0.9998699_dbl_kind,  0.9941218_dbl_kind,  0.8367318_dbl_kind, &
         0.9998515_dbl_kind,  0.9933966_dbl_kind,  0.8227881_dbl_kind, &
         0.9998332_dbl_kind,  0.9926888_dbl_kind,  0.8095131_dbl_kind, &
         0.9998148_dbl_kind,  0.9919968_dbl_kind,  0.7968620_dbl_kind, &
         0.9997691_dbl_kind,  0.9903277_dbl_kind,  0.7677887_dbl_kind/), &
         (/nspint,nmbrad/))

      ! snow asymmetry parameter (unitless)
      real (kind=dbl_kind), dimension (nspint,nmbrad), parameter :: &
         gs_tab = reshape((/ &
          0.859913_dbl_kind,  0.848003_dbl_kind,  0.824415_dbl_kind, &
          0.867130_dbl_kind,  0.858150_dbl_kind,  0.848445_dbl_kind, &
          0.873381_dbl_kind,  0.867221_dbl_kind,  0.861714_dbl_kind, &
          0.878368_dbl_kind,  0.874879_dbl_kind,  0.874036_dbl_kind, &
          0.881462_dbl_kind,  0.879661_dbl_kind,  0.881299_dbl_kind, &
          0.884361_dbl_kind,  0.883903_dbl_kind,  0.890184_dbl_kind, &
          0.885937_dbl_kind,  0.886256_dbl_kind,  0.895393_dbl_kind, &
          0.886931_dbl_kind,  0.887769_dbl_kind,  0.899072_dbl_kind, &
          0.887894_dbl_kind,  0.889255_dbl_kind,  0.903285_dbl_kind, &
          0.888515_dbl_kind,  0.890236_dbl_kind,  0.906588_dbl_kind, &
          0.889073_dbl_kind,  0.891127_dbl_kind,  0.910152_dbl_kind, &
          0.889452_dbl_kind,  0.891750_dbl_kind,  0.913100_dbl_kind, &
          0.889730_dbl_kind,  0.892213_dbl_kind,  0.915621_dbl_kind, &
          0.890026_dbl_kind,  0.892723_dbl_kind,  0.918831_dbl_kind, &
          0.890238_dbl_kind,  0.893099_dbl_kind,  0.921540_dbl_kind, &
          0.890441_dbl_kind,  0.893474_dbl_kind,  0.924581_dbl_kind, &
          0.890618_dbl_kind,  0.893816_dbl_kind,  0.927701_dbl_kind, &
          0.890762_dbl_kind,  0.894123_dbl_kind,  0.930737_dbl_kind, &
          0.890881_dbl_kind,  0.894397_dbl_kind,  0.933568_dbl_kind, &
          0.890975_dbl_kind,  0.894645_dbl_kind,  0.936148_dbl_kind, &
          0.891035_dbl_kind,  0.894822_dbl_kind,  0.937989_dbl_kind, &
          0.891097_dbl_kind,  0.895020_dbl_kind,  0.939949_dbl_kind, &
          0.891147_dbl_kind,  0.895212_dbl_kind,  0.941727_dbl_kind, &
          0.891189_dbl_kind,  0.895399_dbl_kind,  0.943339_dbl_kind, &
          0.891225_dbl_kind,  0.895601_dbl_kind,  0.944915_dbl_kind, &
          0.891248_dbl_kind,  0.895745_dbl_kind,  0.945950_dbl_kind, &
          0.891277_dbl_kind,  0.895951_dbl_kind,  0.947288_dbl_kind, &
          0.891299_dbl_kind,  0.896142_dbl_kind,  0.948438_dbl_kind, &
          0.891323_dbl_kind,  0.896388_dbl_kind,  0.949762_dbl_kind, &
          0.891340_dbl_kind,  0.896623_dbl_kind,  0.950916_dbl_kind, &
          0.891356_dbl_kind,  0.896851_dbl_kind,  0.951945_dbl_kind, &
          0.891386_dbl_kind,  0.897399_dbl_kind,  0.954156_dbl_kind/), &
          (/nspint,nmbrad/))

      ! inherent optical property (iop) arrays for ice and ponded ice
      ! mn = specified mean (or base) value
      ! ki = extinction coefficient (/m)
      ! wi = single scattering albedo
      ! gi = asymmetry parameter

      ! ice surface scattering layer (ssl) iops
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         ki_ssl_mn = (/ 1000.1_dbl_kind, 1003.7_dbl_kind, 7042._dbl_kind/), &
         wi_ssl_mn = (/ .9999_dbl_kind,  .9963_dbl_kind,  .9088_dbl_kind/), &
         gi_ssl_mn = (/  .94_dbl_kind,     .94_dbl_kind,    .94_dbl_kind/)

      ! ice drained layer (dl) iops
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         ki_dl_mn = (/ 100.2_dbl_kind, 107.7_dbl_kind,  1309._dbl_kind /), &
         wi_dl_mn = (/ .9980_dbl_kind,  .9287_dbl_kind, .0305_dbl_kind /), &
         gi_dl_mn = (/ .94_dbl_kind,     .94_dbl_kind,    .94_dbl_kind /)

      ! ice interior layer (int) iops
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         ki_int_mn = (/  20.2_dbl_kind,  27.7_dbl_kind,  1445._dbl_kind /), &
         wi_int_mn = (/ .9901_dbl_kind, .7223_dbl_kind,  .0277_dbl_kind /), &
         gi_int_mn = (/ .94_dbl_kind,    .94_dbl_kind,     .94_dbl_kind /)

      ! ponded ice surface scattering layer (ssl) iops
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         ki_p_ssl_mn = (/ 70.2_dbl_kind,  77.7_dbl_kind,  1309._dbl_kind/), &
         wi_p_ssl_mn = (/ .9972_dbl_kind, .9009_dbl_kind, .0305_dbl_kind/), &
         gi_p_ssl_mn = (/ .94_dbl_kind,   .94_dbl_kind,   .94_dbl_kind  /)

      ! ponded ice interior layer (int) iops
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         ki_p_int_mn = (/  20.2_dbl_kind,  27.7_dbl_kind, 1445._dbl_kind/), &
         wi_p_int_mn = (/ .9901_dbl_kind, .7223_dbl_kind, .0277_dbl_kind/), &
         gi_p_int_mn = (/ .94_dbl_kind,   .94_dbl_kind,   .94_dbl_kind  /)

      ! inherent optical property (iop) arrays for pond water and underlying ocean
      ! kw = Pond water extinction coefficient (/m)
      ! ww = Pond water single scattering albedo
      ! gw = Pond water asymmetry parameter
      real (kind=dbl_kind), dimension (nspint), parameter :: &
         kw = (/ 0.20_dbl_kind,   12.0_dbl_kind,   729._dbl_kind /), &
         ww = (/ 0.00_dbl_kind,   0.00_dbl_kind,   0.00_dbl_kind /), &
         gw = (/ 0.00_dbl_kind,   0.00_dbl_kind,   0.00_dbl_kind /)

      real (kind=dbl_kind), parameter :: &
         rhoi   = 917.0_dbl_kind,& ! pure ice mass density (kg/m3)
         fr_max = 1.00_dbl_kind, & ! snow grain adjustment factor max
         fr_min = 0.80_dbl_kind, & ! snow grain adjustment factor min
      ! tuning parameters
      ! ice and pond scat coeff fractional change for +- one-sigma in albedo
         fp_ice = 0.15_dbl_kind, & ! ice fraction of scat coeff for + stn dev in alb
         fm_ice = 0.15_dbl_kind, & ! ice fraction of scat coeff for - stn dev in alb
         fp_pnd = 2.00_dbl_kind, & ! ponded ice fraction of scat coeff for + stn dev in alb
         fm_pnd = 0.50_dbl_kind    ! ponded ice fraction of scat coeff for - stn dev in alb

      real (kind=dbl_kind),  parameter :: &   !chla-specific absorption coefficient
         kchl_tab = 0.01 !0.0023-0.0029 Perovich 1993, also 0.0067 m^2 (mg Chl)^-1
                         ! found values of 0.006 to 0.023 m^2/ mg  (676 nm)  Neukermans 2014
                         ! and averages over the 300-700nm of 0.0075 m^2/mg in ice Fritsen (2011)
                         ! at 440nm values as high as 0.2 m^2/mg in under ice bloom (Balch 2014)
                         ! Grenfell 1991 uses 0.004 (m^2/mg) which is (0.0078 * spectral weighting)
                         !chlorophyll mass extinction cross section (m^2/mg chla)

      character(len=char_len_long) :: &
         warning ! warning message

!-----------------------------------------------------------------------
! Initialize and tune bare ice/ponded ice iops

      k_bcini(:) = c0
      k_bcins(:) = c0
      k_bcexs(:) = c0

      rnilyr = c1/real(nilyr,kind=dbl_kind)
      rnslyr = c1/real(nslyr,kind=dbl_kind)
      kii = nslyr + 1

      ! initialize albedos and fluxes to 0
      fthrul = c0                
      Iabs = c0
      kabs_chl(:,:) = c0
      tzaer(:,:) = c0
      wzaer(:,:) = c0
      gzaer(:,:) = c0

      avdr   = c0
      avdf   = c0
      aidr   = c0
      aidf   = c0
      fsfc   = c0
      fint   = c0
      fthru  = c0
 
      ! spectral weights
      ! weights 2 (0.7-1.19 micro-meters) and 3 (1.19-5.0 micro-meters) 
      ! are chosen based on 1D calculations using ratio of direct to total 
      ! near-infrared solar (0.7-5.0 micro-meter) which indicates clear/cloudy 
      ! conditions: more cloud, the less 1.19-5.0 relative to the 
      ! 0.7-1.19 micro-meter due to cloud absorption.
      wghtns(1) = c1
      wghtns(2) = cp67 + (cp78-cp67)*(c1-fnidr)
!      wghtns(3) = cp33 + (cp22-cp33)*(c1-fnidr)
      wghtns(3) = c1 - wghtns(2)

      ! find snow grain adjustment factor, dependent upon clear/overcast sky
      ! estimate. comparisons with SNICAR show better agreement with DE when
      ! this factor is included (clear sky near 1 and overcast near 0.8 give
      ! best agreement).  Multiply by rnsw here for efficiency.
      do k = 1, nslyr
         frsnw(k) = (fr_max*fnidr + fr_min*(c1-fnidr))*rsnw(k)
         Sabs(k) = c0
      enddo

      ! layer thicknesses
      ! snow
      dz = hs*rnslyr
      ! for small enough snow thickness, ssl thickness half of top snow layer
!ech: note this is highly resolution dependent!
      dzk(0) = min(hs_ssl, dz/c2)
      dzk(1) = dz - dzk(0)
      if (nslyr > 1) then
         do k = 2, nslyr
            dzk(k) = dz
         enddo
      endif
      
      ! ice
      dz = hi*rnilyr
      ! empirical reduction in sea ice ssl thickness for ice thinner than 1.5m;
      ! factor of 30 gives best albedo comparison with limited observations
      dz_ssl = hi_ssl
!ech: note hardwired parameters
!         if( hi < 1.5_dbl_kind ) dz_ssl = hi/30._dbl_kind
      dz_ssl = min(hi_ssl, hi/30._dbl_kind)
      ! set sea ice ssl thickness to half top layer if sea ice thin enough
!ech: note this is highly resolution dependent!
      dz_ssl = min(dz_ssl, dz/c2)
      
      dzk(kii)   = dz_ssl
      dzk(kii+1) = dz - dz_ssl
      if (kii+2 <= klev) then
         do k = kii+2, klev
            dzk(k) = dz
         enddo
      endif

      ! adjust sea ice iops with tuning parameters; tune only the
      ! scattering coefficient by factors of R_ice, R_pnd, where
      ! R values of +1 correspond approximately to +1 sigma changes in albedo, and
      ! R values of -1 correspond approximately to -1 sigma changes in albedo
      ! Note: the albedo change becomes non-linear for R values > +1 or < -1
      if( R_ice >= c0 ) then
        do ns = 1, nspint
          sigp       = ki_ssl_mn(ns)*wi_ssl_mn(ns)*(c1+fp_ice*R_ice)
          ki_ssl(ns) = sigp+ki_ssl_mn(ns)*(c1-wi_ssl_mn(ns))
          wi_ssl(ns) = sigp/ki_ssl(ns)
          gi_ssl(ns) = gi_ssl_mn(ns)

          sigp       = ki_dl_mn(ns)*wi_dl_mn(ns)*(c1+fp_ice*R_ice)
          ki_dl(ns)  = sigp+ki_dl_mn(ns)*(c1-wi_dl_mn(ns))
          wi_dl(ns)  = sigp/ki_dl(ns)
          gi_dl(ns)  = gi_dl_mn(ns)

          sigp       = ki_int_mn(ns)*wi_int_mn(ns)*(c1+fp_ice*R_ice)
          ki_int(ns) = sigp+ki_int_mn(ns)*(c1-wi_int_mn(ns))
          wi_int(ns) = sigp/ki_int(ns)
          gi_int(ns) = gi_int_mn(ns)
        enddo
      else !if( R_ice < c0 ) then
        do ns = 1, nspint
          sigp       = ki_ssl_mn(ns)*wi_ssl_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_ssl(ns) = sigp+ki_ssl_mn(ns)*(c1-wi_ssl_mn(ns))
          wi_ssl(ns) = sigp/ki_ssl(ns)
          gi_ssl(ns) = gi_ssl_mn(ns)

          sigp       = ki_dl_mn(ns)*wi_dl_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_dl(ns)  = sigp+ki_dl_mn(ns)*(c1-wi_dl_mn(ns))
          wi_dl(ns)  = sigp/ki_dl(ns)
          gi_dl(ns)  = gi_dl_mn(ns)

          sigp       = ki_int_mn(ns)*wi_int_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_int(ns) = sigp+ki_int_mn(ns)*(c1-wi_int_mn(ns))
          wi_int(ns) = sigp/ki_int(ns)
          gi_int(ns) = gi_int_mn(ns)
        enddo
      endif          ! adjust ice iops

      ! adjust ponded ice iops with tuning parameters
      if( R_pnd >= c0 ) then
        do ns = 1, nspint
          sigp         = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)*(c1+fp_pnd*R_pnd)
          ki_p_ssl(ns) = sigp+ki_p_ssl_mn(ns)*(c1-wi_p_ssl_mn(ns))
          wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
          gi_p_ssl(ns) = gi_p_ssl_mn(ns)

          sigp         = ki_p_int_mn(ns)*wi_p_int_mn(ns)*(c1+fp_pnd*R_pnd)
          ki_p_int(ns) = sigp+ki_p_int_mn(ns)*(c1-wi_p_int_mn(ns))
          wi_p_int(ns) = sigp/ki_p_int(ns)
          gi_p_int(ns) = gi_p_int_mn(ns)
        enddo
      else !if( R_pnd < c0 ) then
        do ns = 1, nspint
          sigp         = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)*(c1+fm_pnd*R_pnd)
          sigp         = max(sigp, c0)
          ki_p_ssl(ns) = sigp+ki_p_ssl_mn(ns)*(c1-wi_p_ssl_mn(ns))
          wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
          gi_p_ssl(ns) = gi_p_ssl_mn(ns)

          sigp         = ki_p_int_mn(ns)*wi_p_int_mn(ns)*(c1+fm_pnd*R_pnd)
          sigp         = max(sigp, c0)
          ki_p_int(ns) = sigp+ki_p_int_mn(ns)*(c1-wi_p_int_mn(ns))
          wi_p_int(ns) = sigp/ki_p_int(ns)
          gi_p_int(ns) = gi_p_int_mn(ns)
        enddo
      endif            ! adjust ponded ice iops

      ! use srftyp to determine interface index of surface absorption
      if (srftyp == 1) then
         ! snow covered sea ice
         ksrf = 1
      else
         ! bare sea ice or ponded ice
         ksrf = nslyr + 2 
      endif

      if (tr_bgc_N .and. dEdd_algae) then ! compute kabs_chl for chlorophyll
          do k = 0, klev
             kabs_chl(1,k) = kchl_tab*zbio(nlt_chl_sw+k)
          enddo
      else
            k = klev
            kabs_chl(1,k) = kalg*(0.50_dbl_kind/dzk(k)) 
      endif        ! kabs_chl

!mgf++
      if (modal_aero) then
           do k=0,klev   
              if (k < nslyr+1) then ! define indices for snow layer
                 ! use top rsnw, rhosnw for snow ssl and rest of top layer
                 ksnow = k - min(k-1,0)
                 tmp_gs = frsnw(ksnow)
                
                 ! get grain size index:
                 ! works for 25 < snw_rds < 1625 um:
                 if (tmp_gs < 125) then
                   tmp1 = tmp_gs/50
                   k_bcini(k) = nint(tmp1)
                 elseif (tmp_gs < 175) then
                   k_bcini(k) = 2
                 else
                   tmp1 = (tmp_gs/250)+2
                   k_bcini(k) = nint(tmp1)
                 endif
              else                  ! use the largest snow grain size for ice
                 k_bcini(k) = 8
              endif
              ! Set index corresponding to BC effective radius.  Here,
              ! asssume constant BC effective radius of 100nm
              ! (corresponding to index 2)
              k_bcins(k) = 2
              k_bcexs(k) = 2

              ! check bounds:
              if (k_bcini(k) < 1)  k_bcini(k) = 1
              if (k_bcini(k) > 8)  k_bcini(k) = 8
              if (k_bcins(k) < 1)  k_bcins(k) = 1
              if (k_bcins(k) > 10) k_bcins(k) = 10
              if (k_bcexs(k) < 1)  k_bcexs(k) = 1
              if (k_bcexs(k) > 10) k_bcexs(k) = 10

              ! print ice radius index:
              ! write(warning,*) "MGFICE2:k, ice index= ",k,  k_bcini(k)
              ! call add_warning(warning)
            enddo   ! k

        if (tr_zaero .and. dEdd_algae) then ! compute kzaero for chlorophyll  
        do n = 1,n_zaero        
           if (n == 1) then ! interstitial BC
            do k = 0, klev
               do ns = 1,nspint   ! not weighted by aice
                  tzaer(ns,k) = tzaer(ns,k)+kaer_bc_tab(ns,k_bcexs(k))* &
                                zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  wzaer(ns,k) = wzaer(ns,k)+kaer_bc_tab(ns,k_bcexs(k))* &
                                waer_bc_tab(ns,k_bcexs(k))* &
                                zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  gzaer(ns,k) = gzaer(ns,k)+kaer_bc_tab(ns,k_bcexs(k))* &
                                waer_bc_tab(ns,k_bcexs(k))* &
                                gaer_bc_tab(ns,k_bcexs(k))*zbio(nlt_zaero_sw(n)+k)*dzk(k)
               enddo  ! nspint
            enddo
           elseif (n==2) then ! within-ice BC
            do k = 0, klev
               do ns = 1,nspint  
                  tzaer(ns,k) = tzaer(ns,k)+kaer_bc_tab(ns,k_bcins(k))  * &
                                bcenh(ns,k_bcins(k),k_bcini(k))* &
                                zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  wzaer(ns,k) = wzaer(ns,k)+kaer_bc_tab(ns,k_bcins(k))* &
                                waer_bc_tab(ns,k_bcins(k))* &
                                zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  gzaer(ns,k) = gzaer(ns,k)+kaer_bc_tab(ns,k_bcins(k))* &
                                waer_bc_tab(ns,k_bcins(k))* &
                                gaer_bc_tab(ns,k_bcins(k))*zbio(nlt_zaero_sw(n)+k)*dzk(k)
               enddo  ! nspint
            enddo
           else                ! dust
            do k = 0, klev
               do ns = 1,nspint   ! not weighted by aice
                  tzaer(ns,k) = tzaer(ns,k)+kaer_tab(ns,n)* &
                                   zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  wzaer(ns,k) = wzaer(ns,k)+kaer_tab(ns,n)*waer_tab(ns,n)* &
                                   zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  gzaer(ns,k) = gzaer(ns,k)+kaer_tab(ns,n)*waer_tab(ns,n)* &
                                   gaer_tab(ns,n)*zbio(nlt_zaero_sw(n)+k)*dzk(k)
               enddo  ! nspint
            enddo
           endif      !(n=1)
        enddo         ! n_zaero
        endif         ! tr_zaero and dEdd_algae

     else  ! Bulk aerosol treatment
        if (tr_zaero .and. dEdd_algae) then ! compute kzaero for chlorophyll
        do n = 1,n_zaero          ! multiply by aice?
            do k = 0, klev
               do ns = 1,nspint   ! not weighted by aice
                  tzaer(ns,k) = tzaer(ns,k)+kaer_tab(ns,n)* &
                                   zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  wzaer(ns,k) = wzaer(ns,k)+kaer_tab(ns,n)*waer_tab(ns,n)* &
                                   zbio(nlt_zaero_sw(n)+k)*dzk(k)
                  gzaer(ns,k) = gzaer(ns,k)+kaer_tab(ns,n)*waer_tab(ns,n)* &
                                   gaer_tab(ns,n)*zbio(nlt_zaero_sw(n)+k)*dzk(k)
               enddo  ! nspint
            enddo
        enddo
        endif  !tr_zaero  

     endif  ! modal_aero

!-----------------------------------------------------------------------
 
      ! begin spectral loop
      do ns = 1, nspint
         
         ! set optical properties of air/snow/pond overlying sea ice
         ! air
         if( srftyp == 0 ) then
            do k=0,nslyr 
               tau(k) = c0
               w0(k)  = c0
               g(k)   = c0
            enddo
            ! snow
         else if( srftyp == 1 ) then
            ! interpolate snow iops using input snow grain radius,
            ! snow density and tabular data
            do k=0,nslyr
               ! use top rsnw, rhosnw for snow ssl and rest of top layer
               ksnow = k - min(k-1,0)
               ! find snow iops using input snow density and snow grain radius:
               if( frsnw(ksnow) < rsnw_tab(1) ) then
                  Qs     = Qs_tab(ns,1)
                  ws     = ws_tab(ns,1)
                  gs     = gs_tab(ns,1)
               else if( frsnw(ksnow) >= rsnw_tab(nmbrad) ) then
                  Qs     = Qs_tab(ns,nmbrad)
                  ws     = ws_tab(ns,nmbrad)
                  gs     = gs_tab(ns,nmbrad)
               else
                  ! linear interpolation in rsnw
                  do nr=2,nmbrad
                     if( rsnw_tab(nr-1) <= frsnw(ksnow) .and. &
                         frsnw(ksnow) < rsnw_tab(nr)) then
                        delr = (frsnw(ksnow) - rsnw_tab(nr-1)) / &
                               (rsnw_tab(nr) - rsnw_tab(nr-1))
                        Qs   = Qs_tab(ns,nr-1)*(c1-delr) + &
                               Qs_tab(ns,nr)*delr
                        ws   = ws_tab(ns,nr-1)*(c1-delr) + &
                               ws_tab(ns,nr)*delr
                        gs   = gs_tab(ns,nr-1)*(c1-delr) + &
                               gs_tab(ns,nr)*delr
                     endif
                  enddo       ! nr
               endif
               ks = Qs*((rhosnw(ksnow)/rhoi)*3._dbl_kind / &
                       (4._dbl_kind*frsnw(ksnow)*1.0e-6_dbl_kind))

               tau(k) = (ks + kabs_chl(ns,k))*dzk(k)
               w0(k)  = ks/(ks + kabs_chl(ns,k)) *ws
               g(k)   = gs
            enddo       ! k


            ! aerosol in snow
             if (tr_zaero .and. dEdd_algae) then 
               do k = 0,nslyr
                  gzaer(ns,k) = gzaer(ns,k)/(wzaer(ns,k)+puny)
                  wzaer(ns,k) = wzaer(ns,k)/(tzaer(ns,k)+puny)
                  g(k)   = (g(k)*w0(k)*tau(k) + gzaer(ns,k)*wzaer(ns,k)*tzaer(ns,k)) / &
                                  (w0(k)*tau(k) + wzaer(ns,k)*tzaer(ns,k))
                  w0(k)  = (w0(k)*tau(k) + wzaer(ns,k)*tzaer(ns,k)) / &
                                   (tau(k) + tzaer(ns,k))
                  tau(k) = tau(k) + tzaer(ns,k)
               enddo
             elseif (tr_aero) then
                  k = 0  ! snow SSL
               taer = c0
               waer = c0
               gaer = c0

               do na=1,4*n_aero,4
! mgf++
               if (modal_aero) then
                  if (na == 1) then  
                  !interstitial BC
                     taer = taer + &
                          aero_mp(na)*kaer_bc_tab(ns,k_bcexs(k))
                     waer = waer + &
                          aero_mp(na)*kaer_bc_tab(ns,k_bcexs(k))* &
                          waer_bc_tab(ns,k_bcexs(k))
                     gaer = gaer + &
                          aero_mp(na)*kaer_bc_tab(ns,k_bcexs(k))* &
                           waer_bc_tab(ns,k_bcexs(k))*gaer_bc_tab(ns,k_bcexs(k))
                  elseif (na == 5)then  
                  !within-ice BC
                      taer = taer + &
                           aero_mp(na)*kaer_bc_tab(ns,k_bcins(k))* &
                           bcenh(ns,k_bcins(k),k_bcini(k))
                      waer = waer + &
                           aero_mp(na)*kaer_bc_tab(ns,k_bcins(k))* &
                           waer_bc_tab(ns,k_bcins(k))
                      gaer = gaer + &
                           aero_mp(na)*kaer_bc_tab(ns,k_bcins(k))* &
                           waer_bc_tab(ns,k_bcins(k))*gaer_bc_tab(ns,k_bcins(k))
                   else
                      ! other species (dust)
                      taer = taer + &
                           aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))
                      waer = waer + &
                           aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                      gaer = gaer + &
                           aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                   endif
               else
                  taer = taer + &
                       aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))
                  waer = waer + &
                       aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))* &
                       waer_tab(ns,(1+(na-1)/4))
                  gaer = gaer + &
                       aero_mp(na)*kaer_tab(ns,(1+(na-1)/4))* &
                       waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
               endif  !modal_aero
!mgf--
               enddo       ! na
               gaer = gaer/(waer+puny)
               waer = waer/(taer+puny)

               do k=1,nslyr
                  taer = c0
                  waer = c0
                  gaer = c0
                  do na=1,4*n_aero,4
                  if (modal_aero) then
!mgf++
                    if (na==1) then
                      ! interstitial BC
                      taer = taer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcexs(k))
                      waer = waer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcexs(k))* &
                           waer_bc_tab(ns,k_bcexs(k))
                      gaer = gaer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcexs(k))* &
                           waer_bc_tab(ns,k_bcexs(k))*gaer_bc_tab(ns,k_bcexs(k))
                    elseif (na==5) then
                      ! within-ice BC
                      taer = taer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcins(k))*&
                           bcenh(ns,k_bcins(k),k_bcini(k))
                      waer = waer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcins(k))* &
                           waer_bc_tab(ns,k_bcins(k))
                      gaer = gaer + &
                           (aero_mp(na+1)/rnslyr)*kaer_bc_tab(ns,k_bcins(k))* &
                           waer_bc_tab(ns,k_bcins(k))*gaer_bc_tab(ns,k_bcins(k))
                      
                    else
                      ! other species (dust)
                      taer = taer + &
                           (aero_mp(na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))
                      waer = waer + &
                           (aero_mp(na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                      gaer = gaer + &
                           (aero_mp(na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                    endif   !(na==1)

                  else
                     taer = taer + &
                            (aero_mp(na+1)*rnslyr)*kaer_tab(ns,(1+(na-1)/4))
                     waer = waer + &
                            (aero_mp(na+1)*rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                            waer_tab(ns,(1+(na-1)/4))
                     gaer = gaer + &
                            (aero_mp(na+1)*rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                            waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                  endif       ! modal_aero
!mgf--
                  enddo       ! na
                  gaer = gaer/(waer+puny)
                  waer = waer/(taer+puny)
                  g(k)   = (g(k)*w0(k)*tau(k) + gaer*waer*taer) / &
                           (w0(k)*tau(k) + waer*taer)
                  w0(k)  = (w0(k)*tau(k) + waer*taer) / &
                           (tau(k) + taer)
                  tau(k) = tau(k) + taer
               enddo       ! k
            endif     ! tr_aero

            ! pond
         else !if( srftyp == 2 ) then
            ! pond water layers evenly spaced
            dz = hp/(c1/rnslyr+c1)
            do k=0,nslyr
               tau(k) = kw(ns)*dz
               w0(k)  = ww(ns)
               g(k)   = gw(ns)
               ! no aerosol in pond
            enddo       ! k
         endif        ! srftyp
         
         ! set optical properties of sea ice
         
         ! bare or snow-covered sea ice layers
         if( srftyp <= 1 ) then
            ! ssl
            k = kii
            tau(k) = (ki_ssl(ns)+kabs_chl(ns,k))*dzk(k)
            w0(k)  = ki_ssl(ns)/(ki_ssl(ns) + kabs_chl(ns,k))*wi_ssl(ns)
            g(k)   = gi_ssl(ns)
            ! dl
            k = kii + 1
            ! scale dz for dl relative to 4 even-layer-thickness 1.5m case
            fs = p25/rnilyr
            tau(k) = (ki_dl(ns) + kabs_chl(ns,k)) *dzk(k)*fs
            w0(k)  = ki_dl(ns)/(ki_dl(ns) + kabs_chl(ns,k)) *wi_dl(ns)
            g(k)   = gi_dl(ns)
            ! int above lowest layer
            if (kii+2 <= klev-1) then
               do k = kii+2, klev-1
                  tau(k) = (ki_int(ns) + kabs_chl(ns,k))*dzk(k)
                  w0(k)  = ki_int(ns)/(ki_int(ns) + kabs_chl(ns,k)) *wi_int(ns)
                  g(k)   = gi_int(ns)
               enddo
            endif
            ! lowest layer
            k = klev
            ! add algae to lowest sea ice layer, visible only:
            kabs = ki_int(ns)*(c1-wi_int(ns))
            if( ns == 1 ) then
               ! total layer absorption optical depth fixed at value
               ! of kalg*0.50m, independent of actual layer thickness
               kabs = kabs + kabs_chl(ns,k) 
            endif
            sig        = ki_int(ns)*wi_int(ns)
            tau(k) = (kabs+sig)*dzk(k)
            w0(k)  = (sig/(sig+kabs))
            g(k)   = gi_int(ns)
            ! aerosol in sea ice
            if (tr_zaero .and. dEdd_algae) then
               do k = kii, klev                  
                  gzaer(ns,k) = gzaer(ns,k)/(wzaer(ns,k)+puny)
                  wzaer(ns,k) = wzaer(ns,k)/(tzaer(ns,k)+puny)
                  g(k)   = (g(k)*w0(k)*tau(k) + gzaer(ns,k)*wzaer(ns,k)*tzaer(ns,k)) / &
                                  (w0(k)*tau(k) + wzaer(ns,k)*tzaer(ns,k))
                  w0(k)  = (w0(k)*tau(k) + wzaer(ns,k)*tzaer(ns,k)) / &
                                   (tau(k) + tzaer(ns,k))
                  tau(k) = tau(k) + tzaer(ns,k)
               enddo
            elseif (tr_aero) then 
               k = kii   ! sea ice SSL
               taer = c0
               waer = c0
               gaer = c0
               do na=1,4*n_aero,4
!mgf++
               if (modal_aero) then
                  if (na==1) then
                  ! interstitial BC
                     taer = taer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcexs(k))
                     waer = waer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcexs(k))* &
                          waer_bc_tab(ns,k_bcexs(k))
                     gaer = gaer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcexs(k))* &
                          waer_bc_tab(ns,k_bcexs(k))*gaer_bc_tab(ns,k_bcexs(k))
                  elseif (na==5) then
                  ! within-ice BC
                     taer = taer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcins(k))* &
                          bcenh(ns,k_bcins(k),k_bcini(k))
                     waer = waer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcins(k))* &
                          waer_bc_tab(ns,k_bcins(k))
                     gaer = gaer + &
                          aero_mp(na+2)*kaer_bc_tab(ns,k_bcins(k))* &
                          waer_bc_tab(ns,k_bcins(k))*gaer_bc_tab(ns,k_bcins(k))    
                  else
                  ! other species (dust)
                     taer = taer + &
                          aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))
                     waer = waer + &
                          aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                          waer_tab(ns,(1+(na-1)/4))
                     gaer = gaer + &
                          aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                          waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                  endif
               else      !bulk
                  taer = taer + &
                         aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))
                  waer = waer + &
                         aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))
                  gaer = gaer + &
                         aero_mp(na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                endif     ! modal_aero
!mgf--
               enddo      ! na

               gaer = gaer/(waer+puny)
               waer = waer/(taer+puny)
               g(k)   = (g(k)*w0(k)*tau(k) + gaer*waer*taer) / &
                    (w0(k)*tau(k) + waer*taer)
               w0(k)  = (w0(k)*tau(k) + waer*taer) / &
                    (tau(k) + taer)
               tau(k) = tau(k) + taer
               do k = kii+1, klev
                  taer = c0
                  waer = c0
                  gaer = c0
                  do na=1,4*n_aero,4
!mgf++
                  if (modal_aero) then
                     if (na==1) then
                        ! interstitial BC
                        taer = taer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcexs(k))
                        waer = waer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcexs(k))* &
                             waer_bc_tab(ns,k_bcexs(k))
                        gaer = gaer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcexs(k))* &
                             waer_bc_tab(ns,k_bcexs(k))*gaer_bc_tab(ns,k_bcexs(k))
                     elseif (na==5) then
                        ! within-ice BC
                        taer = taer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcins(k))* &
                             bcenh(ns,k_bcins(k),k_bcini(k))
                        waer = waer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcins(k))* &
                             waer_bc_tab(ns,k_bcins(k))
                        gaer = gaer + &
                             (aero_mp(na+3)/rnilyr)*kaer_bc_tab(ns,k_bcins(k))* &
                             waer_bc_tab(ns,k_bcins(k))*gaer_bc_tab(ns,k_bcins(k))
                        
                     else
                        ! other species (dust)
                        taer = taer + &
                             (aero_mp(na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))
                        waer = waer + &
                             (aero_mp(na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))
                        gaer = gaer + &
                             (aero_mp(na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                     endif
                  else       !bulk

                     taer = taer + &
                          (aero_mp(na+3)*rnilyr)*kaer_tab(ns,(1+(na-1)/4))
                     waer = waer + &
                          (aero_mp(na+3)*rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                          waer_tab(ns,(1+(na-1)/4))
                     gaer = gaer + &
                          (aero_mp(na+3)*rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                          waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                  endif       ! modal_aero
!mgf--
                  enddo       ! na
                  gaer = gaer/(waer+puny)
                  waer = waer/(taer+puny)
                  g(k)   = (g(k)*w0(k)*tau(k) + gaer*waer*taer) / &
                           (w0(k)*tau(k) + waer*taer)
                  w0(k)  = (w0(k)*tau(k) + waer*taer) / &
                           (tau(k) + taer)
                  tau(k) = tau(k) + taer
               enddo ! k
            endif ! tr_aero

            ! sea ice layers under ponds
         else !if( srftyp == 2 ) then
            k = kii
            tau(k) = ki_p_ssl(ns)*dzk(k)
            w0(k)  = wi_p_ssl(ns)
            g(k)   = gi_p_ssl(ns)
            k = kii + 1
            tau(k) = ki_p_int(ns)*dzk(k)
            w0(k)  = wi_p_int(ns)
            g(k)   = gi_p_int(ns)
            if (kii+2 <= klev) then
               do k = kii+2, klev
                  tau(k) = ki_p_int(ns)*dzk(k)
                  w0(k)  = wi_p_int(ns)
                  g(k)   = gi_p_int(ns)
               enddo       ! k
            endif
            ! adjust pond iops if pond depth within specified range
            if( hpmin <= hp .and. hp <= hp0 ) then
               k = kii
               sig_i      = ki_ssl(ns)*wi_ssl(ns)
               sig_p      = ki_p_ssl(ns)*wi_p_ssl(ns)
               sig        = sig_i + (sig_p-sig_i)*(hp/hp0)
               kext       = sig + ki_p_ssl(ns)*(c1-wi_p_ssl(ns))
               tau(k) = kext*dzk(k)
               w0(k) = sig/kext
               g(k)  = gi_p_int(ns)
               k = kii + 1
               ! scale dz for dl relative to 4 even-layer-thickness 1.5m case
               fs = p25/rnilyr
               sig_i      = ki_dl(ns)*wi_dl(ns)*fs
               sig_p      = ki_p_int(ns)*wi_p_int(ns)
               sig        = sig_i + (sig_p-sig_i)*(hp/hp0)
               kext       = sig + ki_p_int(ns)*(c1-wi_p_int(ns))
               tau(k) = kext*dzk(k)
               w0(k) = sig/kext
               g(k)  = gi_p_int(ns)
               if (kii+2 <= klev) then
                  do k = kii+2, klev
                     sig_i      = ki_int(ns)*wi_int(ns)
                     sig_p      = ki_p_int(ns)*wi_p_int(ns)
                     sig        = sig_i + (sig_p-sig_i)*(hp/hp0)
                     kext       = sig + ki_p_int(ns)*(c1-wi_p_int(ns))
                     tau(k) = kext*dzk(k)
                     w0(k) = sig/kext
                     g(k)  = gi_p_int(ns)
                  enddo       ! k
               endif
            endif        ! small pond depth transition to bare sea ice
         endif         ! srftyp  
         
         ! set reflectivities for ocean underlying sea ice
         rns = real(ns-1, kind=dbl_kind)
         albodr = cp01 * (c1 - min(rns, c1))
         albodf = cp01 * (c1 - min(rns, c1))
         
         ! layer input properties now completely specified: tau, w0, g,
         ! albodr, albodf; now compute the Delta-Eddington solution 
         ! reflectivities and transmissivities for each layer; then,
         ! combine the layers going downwards accounting for multiple
         ! scattering between layers, and finally start from the 
         ! underlying ocean and combine successive layers upwards to
         ! the surface; see comments in solution_dEdd for more details.
         
         call solution_dEdd &
               (coszen,     srftyp,     klev,       klevp,      nslyr,     &
                tau,        w0,         g,          albodr,     albodf,    &
                trndir,     trntdr,     trndif,     rupdir,     rupdif,    &
                rdndif)   

         ! the interface reflectivities and transmissivities required
         ! to evaluate interface fluxes are returned from solution_dEdd;
         ! now compute up and down fluxes for each interface, using the 
         ! combined layer properties at each interface:
         !
         !              layers       interface
         !
         !       ---------------------  k
         !                 k
         !       --------------------- 
         
         do k = 0, klevp 
            ! interface scattering
            refk          = c1/(c1 - rdndif(k)*rupdif(k))
            ! dir tran ref from below times interface scattering, plus diff
            ! tran and ref from below times interface scattering
            ! fdirup(k) = (trndir(k)*rupdir(k) + &
            !                 (trntdr(k)-trndir(k))  &
            !                 *rupdif(k))*refk
            ! dir tran plus total diff trans times interface scattering plus
            ! dir tran with up dir ref and down dif ref times interface scattering 
            ! fdirdn(k) = trndir(k) + (trntdr(k) &
            !               - trndir(k) + trndir(k)  &
            !               *rupdir(k)*rdndif(k))*refk
            ! diffuse tran ref from below times interface scattering
            ! fdifup(k) = trndif(k)*rupdif(k)*refk
            ! diffuse tran times interface scattering
            ! fdifdn(k) = trndif(k)*refk

            ! dfdir = fdirdn - fdirup
            dfdir(k) = trndir(k) &
                        + (trntdr(k)-trndir(k)) * (c1 - rupdif(k)) * refk &
                        -  trndir(k)*rupdir(k)  * (c1 - rdndif(k)) * refk
            if (dfdir(k) < puny) dfdir(k) = c0 !echmod necessary?
            ! dfdif = fdifdn - fdifup
            dfdif(k) = trndif(k) * (c1 - rupdif(k)) * refk
            if (dfdif(k) < puny) dfdif(k) = c0 !echmod necessary?
         enddo       ! k 
         
         ! calculate final surface albedos and fluxes-
         ! all absorbed flux above ksrf is included in surface absorption
         
         if( ns == 1) then      ! visible
            
            swdr = swvdr
            swdf = swvdf
            avdr  = rupdir(0)
            avdf  = rupdif(0)
            
            tmp_0  = dfdir(0    )*swdr + dfdif(0    )*swdf
            tmp_ks = dfdir(ksrf )*swdr + dfdif(ksrf )*swdf
            tmp_kl = dfdir(klevp)*swdr + dfdif(klevp)*swdf

            ! for layer biology: save visible only
            do k = nslyr+2, klevp ! Start at DL layer of ice after SSL scattering
               fthrul(k-nslyr-1) = dfdir(k)*swdr + dfdif(k)*swdf
            enddo
            
            fsfc  = fsfc  + tmp_0  - tmp_ks
            fint  = fint  + tmp_ks - tmp_kl
            fthru = fthru + tmp_kl
            
            ! if snow covered ice, set snow internal absorption; else, Sabs=0
            if( srftyp == 1 ) then
               ki = 0
               do k=1,nslyr
                  ! skip snow SSL, since SSL absorption included in the surface
                  ! absorption fsfc above
                  km  = k
                  kp  = km + 1
                  ki  = ki + 1
                  Sabs(ki) = Sabs(ki) &
                           +  dfdir(km)*swdr + dfdif(km)*swdf &
                           - (dfdir(kp)*swdr + dfdif(kp)*swdf)
               enddo       ! k
            endif
            
            ! complex indexing to insure proper absorptions for sea ice
            ki = 0
            do k=nslyr+2,nslyr+1+nilyr
               ! for bare ice, DL absorption for sea ice layer 1
               km = k  
               kp = km + 1
               ! modify for top sea ice layer for snow over sea ice
               if( srftyp == 1 ) then
                  ! must add SSL and DL absorption for sea ice layer 1
                  if( k == nslyr+2 ) then
                     km = k  - 1
                     kp = km + 2
                  endif
               endif
               ki = ki + 1
               Iabs(ki) = Iabs(ki) &
                        +  dfdir(km)*swdr + dfdif(km)*swdf &
                        - (dfdir(kp)*swdr + dfdif(kp)*swdf)
            enddo       ! k
            
         else !if(ns > 1) then  ! near IR
            
            swdr = swidr
            swdf = swidf

            ! let fr1 = alb_1*swd*wght1 and fr2 = alb_2*swd*wght2 be the ns=2,3
            ! reflected fluxes respectively, where alb_1, alb_2 are the band
            ! albedos, swd = nir incident shortwave flux, and wght1, wght2 are
            ! the 2,3 band weights. thus, the total reflected flux is:
            ! fr = fr1 + fr2 = alb_1*swd*wght1 + alb_2*swd*wght2  hence, the
            ! 2,3 nir band albedo is alb = fr/swd = alb_1*wght1 + alb_2*wght2

            aidr   = aidr + rupdir(0)*wghtns(ns)
            aidf   = aidf + rupdif(0)*wghtns(ns)

            tmp_0  = dfdir(0    )*swdr + dfdif(0    )*swdf
            tmp_ks = dfdir(ksrf )*swdr + dfdif(ksrf )*swdf
            tmp_kl = dfdir(klevp)*swdr + dfdif(klevp)*swdf

            tmp_0  = tmp_0  * wghtns(ns)
            tmp_ks = tmp_ks * wghtns(ns)
            tmp_kl = tmp_kl * wghtns(ns)

            fsfc  = fsfc  + tmp_0  - tmp_ks
            fint  = fint  + tmp_ks - tmp_kl
            fthru = fthru + tmp_kl

            ! if snow covered ice, set snow internal absorption; else, Sabs=0
            if( srftyp == 1 ) then
               ki = 0
               do k=1,nslyr
                  ! skip snow SSL, since SSL absorption included in the surface
                  ! absorption fsfc above
                  km  = k
                  kp  = km + 1
                  ki  = ki + 1
                  Sabs(ki) = Sabs(ki) &
                           + (dfdir(km)*swdr + dfdif(km)*swdf   &
                           - (dfdir(kp)*swdr + dfdif(kp)*swdf)) & 
                           * wghtns(ns)
               enddo       ! k
            endif
            
            ! complex indexing to insure proper absorptions for sea ice
            ki = 0
            do k=nslyr+2,nslyr+1+nilyr
               ! for bare ice, DL absorption for sea ice layer 1
               km = k  
               kp = km + 1
               ! modify for top sea ice layer for snow over sea ice
               if( srftyp == 1 ) then
                  ! must add SSL and DL absorption for sea ice layer 1
                  if( k == nslyr+2 ) then
                     km = k  - 1
                     kp = km + 2
                  endif
               endif
               ki = ki + 1
               Iabs(ki) = Iabs(ki) &
                        + (dfdir(km)*swdr + dfdif(km)*swdf &
                        - (dfdir(kp)*swdr + dfdif(kp)*swdf)) &
                        * wghtns(ns)
            enddo       ! k
            
         endif        ! ns = 1, ns > 1
         
      enddo         ! end spectral loop  ns

      ! accumulate fluxes over bare sea ice
      alvdr   = avdr
      alvdf   = avdf
      alidr   = aidr
      alidf   = aidf
      fswsfc  = fswsfc  + fsfc *fi
      fswint  = fswint  + fint *fi
      fswthru = fswthru + fthru*fi

      do k = 1, nslyr
         Sswabs(k) = Sswabs(k) + Sabs(k)*fi
      enddo                     ! k
      
      do k = 1, nilyr
         Iswabs(k) = Iswabs(k) + Iabs(k)*fi
         
         ! bgc layer 
         fswpenl(k) = fswpenl(k) + fthrul(k)* fi
         if (k == nilyr) then
            fswpenl(k+1) = fswpenl(k+1) + fthrul(k+1)*fi
         endif
      enddo                     ! k
      
      !----------------------------------------------------------------
      ! if ice has zero heat capacity, no SW can be absorbed 
      ! in the ice/snow interior, so add to surface absorption.
      ! Note: nilyr = nslyr = 1 for this case
      !----------------------------------------------------------------

      if (.not. heat_capacity) then
         
         ! SW absorbed at snow/ice surface
         fswsfc = fswsfc + Iswabs(1) + Sswabs(1)
         
         ! SW absorbed in ice interior
         fswint   = c0
         Iswabs(1) = c0
         Sswabs(1) = c0

      endif                       ! heat_capacity

      end subroutine compute_dEdd

!=======================================================================
!
! Given input vertical profiles of optical properties, evaluate the 
! monochromatic Delta-Eddington solution.
!
! author:  Bruce P. Briegleb, NCAR
!   2013:  E Hunke merged with NCAR version
      subroutine solution_dEdd                                 &
            (coszen,     srftyp,    klev,      klevp,  nslyr,  &
             tau,        w0,        g,         albodr, albodf, &
             trndir,     trntdr,    trndif,    rupdir, rupdif, &
             rdndif)

      real (kind=dbl_kind), intent(in) :: &
         coszen      ! cosine solar zenith angle

      integer (kind=int_kind), intent(in) :: &
         srftyp   , & ! surface type over ice: (0=air, 1=snow, 2=pond)
         klev     , & ! number of radiation layers - 1
         klevp    , & ! number of radiation interfaces - 1
                      ! (0 layer is included also)
         nslyr        ! number of snow layers
 
      real (kind=dbl_kind), dimension(0:klev), intent(in) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter
 
      real (kind=dbl_kind), intent(in) :: &
         albodr  , & ! ocean albedo to direct rad
         albodf      ! ocean albedo to diffuse rad
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (0:klevp), intent(out) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above

!-----------------------------------------------------------------------
!
! Delta-Eddington solution for snow/air/pond over sea ice
!
! Generic solution for a snow/air/pond input column of klev+1 layers,
! with srftyp determining at what interface fresnel refraction occurs.
!
! Computes layer reflectivities and transmissivities, from the top down
! to the lowest interface using the Delta-Eddington solutions for each
! layer; combines layers from top down to lowest interface, and from the
! lowest interface (underlying ocean) up to the top of the column.
!
! Note that layer diffuse reflectivity and transmissivity are computed
! by integrating the direct over several gaussian angles. This is
! because the diffuse reflectivity expression sometimes is negative,
! but the direct reflectivity is always well-behaved. We assume isotropic
! radiation in the upward and downward hemispheres for this integration.
!
! Assumes monochromatic (spectrally uniform) properties across a band
! for the input optical parameters.
!
! If total transmission of the direct beam to the interface above a particular 
! layer is less than trmin, then no further Delta-Eddington solutions are
! evaluated for layers below.
!
! The following describes how refraction is handled in the calculation.
!
! First, we assume that radiation is refracted when entering either
! sea ice at the base of the surface scattering layer, or water (i.e. melt
! pond); we assume that radiation does not refract when entering snow, nor 
! upon entering sea ice from a melt pond, nor upon entering the underlying 
! ocean from sea ice.
!
! To handle refraction, we define a "fresnel" layer, which physically
! is of neglible thickness and is non-absorbing, which can be combined to 
! any sea ice layer or top of melt pond. The fresnel layer accounts for 
! refraction of direct beam and associated reflection and transmission for
! solar radiation. A fresnel layer is combined with the top of a melt pond 
! or to the surface scattering layer of sea ice if no melt pond lies over it. 
!
! Some caution must be exercised for the fresnel layer, because any layer
! to which it is combined is no longer a homogeneous layer, as are all other
! individual layers. For all other layers for example, the direct and diffuse
! reflectivities/transmissivities (R/T) are the same for radiation above or
! below the layer. This is the meaning of homogeneous! But for the fresnel
! layer this is not so. Thus, the R/T for this layer must be distinguished
! for radiation above from that from radiation below. For generality, we
! treat all layers to be combined as inhomogeneous.
!
!-----------------------------------------------------------------------

      ! local variables

      integer (kind=int_kind) :: &
         kfrsnl      ! radiation interface index for fresnel layer
 
      ! following variables are defined for each layer; 0 refers to the top
      ! layer. In general we must distinguish directions above and below in 
      ! the diffuse reflectivity and transmissivity, as layers are not assumed
      ! to be homogeneous (apart from the single layer Delta-Edd solutions); 
      ! the direct is always from above.
      real (kind=dbl_kind), dimension (0:klev) :: &
         rdir    , & ! layer reflectivity to direct radiation
         rdif_a  , & ! layer reflectivity to diffuse radiation from above
         rdif_b  , & ! layer reflectivity to diffuse radiation from below
         tdir    , & ! layer transmission to direct radiation (solar beam + diffuse)
         tdif_a  , & ! layer transmission to diffuse radiation from above
         tdif_b  , & ! layer transmission to diffuse radiation from below
         trnlay      ! solar beam transm for layer (direct beam only)

      integer (kind=int_kind) :: & 
         k           ! level index
 
      real (kind=dbl_kind), parameter :: &
         trmin = 0.001_dbl_kind   ! minimum total transmission allowed
      ! total transmission is that due to the direct beam; i.e. it includes
      ! both the directly transmitted solar beam and the diffuse downwards
      ! transmitted radiation resulting from scattering out of the direct beam 
      real (kind=dbl_kind) :: &
         tautot   , & ! layer optical depth
         wtot     , & ! layer single scattering albedo
         gtot     , & ! layer asymmetry parameter
         ftot     , & ! layer forward scattering fraction
         ts       , & ! layer scaled extinction optical depth
         ws       , & ! layer scaled single scattering albedo
         gs       , & ! layer scaled asymmetry parameter
         rintfc   , & ! reflection (multiple) at an interface
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref 
         tdndif       ! total down diffuse = tot tran - direct tran
 
      ! perpendicular and parallel relative to plane of incidence and scattering
      real (kind=dbl_kind) :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below
 
      ! refractive index for sea ice, water; pre-computed, band-independent,
      ! diffuse fresnel reflectivities
      real (kind=dbl_kind), parameter :: & 
         refindx = 1.310_dbl_kind  , & ! refractive index of sea ice (water also)
         cp063   = 0.063_dbl_kind  , & ! diffuse fresnel reflectivity from above
         cp455   = 0.455_dbl_kind      ! diffuse fresnel reflectivity from below
 
      real (kind=dbl_kind) :: &
         mu0      , & ! cosine solar zenith angle incident
         mu0nij       ! cosine solar zenith angle in medium below fresnel level
 
      real (kind=dbl_kind) :: &
         mu0n         ! cosine solar zenith angle in medium
 
      real (kind=dbl_kind) :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         agamm    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,agamm,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         mu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n
 
      real (kind=dbl_kind) :: &
         w        , & ! dummy argument for statement function
         uu       , & ! dummy argument for statement function
         gg       , & ! dummy argument for statement function
         e        , & ! dummy argument for statement function
         f        , & ! dummy argument for statement function
         t        , & ! dummy argument for statement function
         et           ! dummy argument for statement function
 
      real (kind=dbl_kind) :: &
         alp      , & ! temporary for alpha
         gam      , & ! temporary for agamm
         ue       , & ! temporary for u
         extins   , & ! extinction
         amg      , & ! alp - gam
         apg          ! alp + gam
 
      integer (kind=int_kind), parameter :: &
         ngmax = 8    ! number of gaussian angles in hemisphere
 
      real (kind=dbl_kind), dimension (ngmax), parameter :: &
         gauspt     & ! gaussian angles (radians)
            = (/ .9894009_dbl_kind,  .9445750_dbl_kind, &
                 .8656312_dbl_kind,  .7554044_dbl_kind, &
                 .6178762_dbl_kind,  .4580168_dbl_kind, &
                 .2816036_dbl_kind,  .0950125_dbl_kind/), &
         gauswt     & ! gaussian weights
            = (/ .0271525_dbl_kind,  .0622535_dbl_kind, &
                 .0951585_dbl_kind,  .1246290_dbl_kind, &
                 .1495960_dbl_kind,  .1691565_dbl_kind, &
                 .1826034_dbl_kind,  .1894506_dbl_kind/)
  
      integer (kind=int_kind) :: &
         ng           ! gaussian integration index
 
      real (kind=dbl_kind) :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt          ! accumulator for tdif gaussian integration
 
      ! Delta-Eddington solution expressions
      alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
      agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
      n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
      u(w,gg,e)        = c1p5*(c1 - w*gg)/e
      el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
      taus(w,f,t)      = (c1 - w*f)*t
      omgs(w,f)        = (c1 - f)*w/(c1 - w*f)
      asys(gg,f)       = (gg - f)/(c1 - f)
 
!-----------------------------------------------------------------------

      do k = 0, klevp 
         trndir(k) = c0
         trntdr(k) = c0
         trndif(k) = c0
         rupdir(k) = c0
         rupdif(k) = c0
         rdndif(k) = c0
      enddo
 
      ! initialize top interface of top layer 
      trndir(0) =   c1
      trntdr(0) =   c1
      trndif(0) =   c1
      rdndif(0) =   c0

      ! mu0 is cosine solar zenith angle above the fresnel level; make 
      ! sure mu0 is large enough for stable and meaningful radiation
      ! solution: .01 is like sun just touching horizon with its lower edge
      mu0  = max(coszen,p01)

      ! mu0n is cosine solar zenith angle used to compute the layer
      ! Delta-Eddington solution; it is initially computed to be the
      ! value below the fresnel level, i.e. the cosine solar zenith 
      ! angle below the fresnel level for the refracted solar beam:
      mu0nij = sqrt(c1-((c1-mu0**2)/(refindx*refindx)))
 
      ! compute level of fresnel refraction
      ! if ponded sea ice, fresnel level is the top of the pond.
      kfrsnl = 0
      ! if snow over sea ice or bare sea ice, fresnel level is
      ! at base of sea ice SSL (and top of the sea ice DL); the
      ! snow SSL counts for one, then the number of snow layers,
      ! then the sea ice SSL which also counts for one:
      if( srftyp < 2 ) kfrsnl = nslyr + 2 

      ! proceed down one layer at a time; if the total transmission to
      ! the interface just above a given layer is less than trmin, then no
      ! Delta-Eddington computation for that layer is done.

      ! begin main level loop
      do k = 0, klev
         
         ! initialize all layer apparent optical properties to 0
         rdir  (k) = c0
         rdif_a(k) = c0
         rdif_b(k) = c0
         tdir  (k) = c0
         tdif_a(k) = c0
         tdif_b(k) = c0
         trnlay(k) = c0

         ! compute next layer Delta-eddington solution only if total transmission
         ! of radiation to the interface just above the layer exceeds trmin.
         
         if (trntdr(k) > trmin ) then

            ! calculation over layers with penetrating radiation

            tautot  = tau(k)
            wtot    = w0(k)
            gtot    = g(k)
            ftot    = gtot*gtot

            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            ue   = u(ws,gs,lm)

            mu0n = mu0nij
            ! if level k is above fresnel level and the cell is non-pond, use the
            ! non-refracted beam instead
            if( srftyp < 2 .and. k < kfrsnl ) mu0n = mu0

            extins = max(exp_min, exp(-lm*ts))
            ne = n(ue,extins)

            ! first calculation of rdif, tdif using Delta-Eddington formulas
!            rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
            rdif_a(k) = (ue**2-c1)*(c1/extins - extins)/ne
            tdif_a(k) = c4*ue/ne

            ! evaluate rdir,tdir for direct beam
            trnlay(k) = max(exp_min, exp(-ts/mu0n))
            alp = alpha(ws,mu0n,gs,lm)
            gam = agamm(ws,mu0n,gs,lm)
            apg = alp + gam
            amg = alp - gam
            rdir(k) = apg*rdif_a(k) +  amg*(tdif_a(k)*trnlay(k) - c1)
            tdir(k) = apg*tdif_a(k) + (amg* rdif_a(k)-apg+c1)*trnlay(k)
            
            ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
            ! since Delta-Eddington rdif formula is not well-behaved (it is usually
            ! biased low and can even be negative); use ngmax angles and gaussian
            ! integration for most accuracy:
            R1 = rdif_a(k) ! use R1 as temporary
            T1 = tdif_a(k) ! use T1 as temporary
            swt = c0
            smr = c0
            smt = c0
            do ng=1,ngmax
               mu  = gauspt(ng)
               gwt = gauswt(ng)
               swt = swt + mu*gwt
               trn = max(exp_min, exp(-ts/mu))
               alp = alpha(ws,mu,gs,lm)
               gam = agamm(ws,mu,gs,lm)
               apg = alp + gam
               amg = alp - gam
               rdr = apg*R1 + amg*T1*trn - amg
               tdr = apg*T1 + amg*R1*trn - apg*trn + trn
               smr = smr + mu*rdr*gwt
               smt = smt + mu*tdr*gwt
            enddo      ! ng
            rdif_a(k) = smr/swt
            tdif_a(k) = smt/swt
            
            ! homogeneous layer
            rdif_b(k) = rdif_a(k)
            tdif_b(k) = tdif_a(k)

            ! add fresnel layer to top of desired layer if either 
            ! air or snow overlies ice; we ignore refraction in ice 
            ! if a melt pond overlies it:

            if( k == kfrsnl ) then
               ! compute fresnel reflection and transmission amplitudes
               ! for two polarizations: 1=perpendicular and 2=parallel to
               ! the plane containing incident, reflected and refracted rays.
               R1 = (mu0 - refindx*mu0n) / & 
                    (mu0 + refindx*mu0n)
               R2 = (refindx*mu0 - mu0n) / &
                    (refindx*mu0 + mu0n)
               T1 = c2*mu0 / &
                    (mu0 + refindx*mu0n)
               T2 = c2*mu0 / &
                    (refindx*mu0 + mu0n)

               ! unpolarized light for direct beam
               Rf_dir_a = p5 * (R1*R1 + R2*R2)
               Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0

               ! precalculated diffuse reflectivities and transmissivities
               ! for incident radiation above and below fresnel layer, using
               ! the direct albedos and accounting for complete internal
               ! reflection from below; precalculated because high order
               ! number of gaussian points (~256) is required for convergence:

               ! above
               Rf_dif_a = cp063
               Tf_dif_a = c1 - Rf_dif_a
               ! below
               Rf_dif_b = cp455
               Tf_dif_b = c1 - Rf_dif_b

               ! the k = kfrsnl layer properties are updated to combined 
               ! the fresnel (refractive) layer, always taken to be above
               ! the present layer k (i.e. be the top interface):

               rintfc   = c1 / (c1-Rf_dif_b*rdif_a(k))
               tdir(k)   = Tf_dir_a*tdir(k) + &
                    Tf_dir_a*rdir(k) * &
                    Rf_dif_b*rintfc*tdif_a(k)
               rdir(k)   = Rf_dir_a + &
                    Tf_dir_a*rdir(k) * &
                    rintfc*Tf_dif_b
               rdif_a(k) = Rf_dif_a + &
                    Tf_dif_a*rdif_a(k) * &
                    rintfc*Tf_dif_b
               rdif_b(k) = rdif_b(k) + &
                    tdif_b(k)*Rf_dif_b * &
                    rintfc*tdif_a(k)
               tdif_a(k) = tdif_a(k)*rintfc*Tf_dif_a
               tdif_b(k) = tdif_b(k)*rintfc*Tf_dif_b

               ! update trnlay to include fresnel transmission
               trnlay(k) = Tf_dir_a*trnlay(k)

            endif      ! k = kfrsnl

         endif ! trntdr(k) > trmin
         
         ! initialize current layer properties to zero; only if total
         ! transmission to the top interface of the current layer exceeds the
         ! minimum, will these values be computed below:
         ! Calculate the solar beam transmission, total transmission, and
         ! reflectivity for diffuse radiation from below at interface k, 
         ! the top of the current layer k:
         !
         !              layers       interface
         !         
         !       ---------------------  k-1 
         !                k-1
         !       ---------------------  k
         !                 k
         !       ---------------------  
         !       For k = klevp
         ! note that we ignore refraction between sea ice and underlying ocean:
         !
         !              layers       interface
         !
         !       ---------------------  k-1 
         !                k-1
         !       ---------------------  k
         !       \\\\\\\ ocean \\\\\\\
         
         trndir(k+1) = trndir(k)*trnlay(k)
         refkm1         = c1/(c1 - rdndif(k)*rdif_a(k))
         tdrrdir        = trndir(k)*rdir(k)
         tdndif         = trntdr(k) - trndir(k)
         trntdr(k+1) = trndir(k)*tdir(k) + &
              (tdndif + tdrrdir*rdndif(k))*refkm1*tdif_a(k)
         rdndif(k+1) = rdif_b(k) + &
              (tdif_b(k)*rdndif(k)*refkm1*tdif_a(k))
         trndif(k+1) = trndif(k)*refkm1*tdif_a(k)

      enddo       ! k   end main level loop

      ! compute reflectivity to direct and diffuse radiation for layers 
      ! below by adding succesive layers starting from the underlying 
      ! ocean and working upwards:
      !
      !              layers       interface
      !
      !       ---------------------  k
      !                 k
      !       ---------------------  k+1
      !                k+1
      !       ---------------------

      rupdir(klevp) = albodr
      rupdif(klevp) = albodf 

      do k=klev,0,-1
         ! interface scattering
         refkp1        = c1/( c1 - rdif_b(k)*rupdif(k+1))
         ! dir from top layer plus exp tran ref from lower layer, interface
         ! scattered and tran thru top layer from below, plus diff tran ref
         ! from lower layer with interface scattering tran thru top from below
         rupdir(k) = rdir(k) &
              + (        trnlay(k)  *rupdir(k+1) &
              +  (tdir(k)-trnlay(k))*rupdif(k+1))*refkp1*tdif_b(k)
         ! dif from top layer from above, plus dif tran upwards reflected and
         ! interface scattered which tran top from below
         rupdif(k) = rdif_a(k) + tdif_a(k)*rupdif(k+1)*refkp1*tdif_b(k)
      enddo       ! k

      end subroutine solution_dEdd

!=======================================================================
!
!   Set snow horizontal coverage, density and grain radius diagnostically 
!   for the Delta-Eddington solar radiation method.
!
! author:  Bruce P. Briegleb, NCAR 
!   2013:  E Hunke merged with NCAR version

      subroutine shortwave_dEdd_set_snow(nslyr,    R_snw,    &
                                         dT_mlt,   rsnw_mlt, &
                                         aice,     vsno,     &
                                         Tsfc,     fs,       &
                                         hs0,      hs,       &
                                         rhosnw,   rsnw)

      integer (kind=int_kind), intent(in) :: & 
         nslyr      ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         R_snw , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt, & ! change in temp for non-melt to melt snow grain radius change (C)
         rsnw_mlt  ! maximum melting snow grain radius (10^-6 m)

      real (kind=dbl_kind), intent(in) :: &
         aice   , & ! concentration of ice
         vsno   , & ! volume of snow
         Tsfc   , & ! surface temperature 
         hs0        ! snow depth for transition to bare sea ice (m)

      real (kind=dbl_kind), intent(out) :: &
         fs     , & ! horizontal coverage of snow
         hs         ! snow depth

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         rhosnw , & ! density in snow layer (kg/m3)
         rsnw       ! grain radius in snow layer (micro-meters)

      ! local variables

      integer (kind=int_kind) :: &
         ks           ! snow vertical index

      real (kind=dbl_kind) :: &
         fT  , & ! piecewise linear function of surface temperature
         dTs , & ! difference of Tsfc and Timelt
         rsnw_nm ! actual used nonmelt snow grain radius (micro-meters)

      real (kind=dbl_kind), parameter :: &
         ! units for the following are 1.e-6 m (micro-meters)
         rsnw_fresh    =  100._dbl_kind, & ! freshly-fallen snow grain radius 
         rsnw_nonmelt  =  500._dbl_kind, & ! nonmelt snow grain radius
         rsnw_sig      =  250._dbl_kind    ! assumed sigma for snow grain radius

!-----------------------------------------------------------------------

      ! set snow horizontal fraction
      hs = vsno / aice
      
      if (hs >= hs_min) then
         fs = c1
         if (hs0 > puny) fs = min(hs/hs0, c1)
      endif
      
      ! bare ice, temperature dependence
      dTs = Timelt - Tsfc
      fT  = -min(dTs/dT_mlt-c1,c0)
      ! tune nonmelt snow grain radius if desired: note that
      ! the sign is negative so that if R_snw is 1, then the
      ! snow grain radius is reduced and thus albedo increased.
      rsnw_nm = rsnw_nonmelt - R_snw*rsnw_sig
      rsnw_nm = max(rsnw_nm, rsnw_fresh)
      rsnw_nm = min(rsnw_nm, rsnw_mlt) 
      do ks = 1, nslyr
         ! snow density ccsm3 constant value
         rhosnw(ks) = rhos
         ! snow grain radius between rsnw_nonmelt and rsnw_mlt
         rsnw(ks) = rsnw_nm + (rsnw_mlt-rsnw_nm)*fT
         rsnw(ks) = max(rsnw(ks), rsnw_fresh)
         rsnw(ks) = min(rsnw(ks), rsnw_mlt)
      enddo        ! ks

      end subroutine shortwave_dEdd_set_snow

!=======================================================================
!
!   Set pond fraction and depth diagnostically for
!   the Delta-Eddington solar radiation method.
!
! author:  Bruce P. Briegleb, NCAR 
!   2013:  E Hunke merged with NCAR version

      subroutine shortwave_dEdd_set_pond(Tsfc,               &
                                         fs,       fp,       &
                                         hp)

      real (kind=dbl_kind), intent(in) :: &
         Tsfc   , & ! surface temperature
         fs         ! horizontal coverage of snow

      real (kind=dbl_kind), intent(out) :: &
         fp     , & ! pond fractional coverage (0 to 1)
         hp         ! pond depth (m)

      ! local variables

      real (kind=dbl_kind) :: &
         fT  , & ! piecewise linear function of surface temperature
         dTs     ! difference of Tsfc and Timelt

      real (kind=dbl_kind), parameter :: &
         dT_pnd = c1   ! change in temp for pond fraction and depth

!-----------------------------------------------------------------------

      ! bare ice, temperature dependence
      dTs = Timelt - Tsfc
      fT  = -min(dTs/dT_pnd-c1,c0)
      ! pond
      fp = 0.3_dbl_kind*fT*(c1-fs)
      hp = 0.3_dbl_kind*fT*(c1-fs)
      
      end subroutine shortwave_dEdd_set_pond

! End Delta-Eddington shortwave method

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine compute_shortwave_trcr(n_algae,  nslyr,  &
                                    trcrn,        trcrn_sw,  &
                                    sw_grid,      hin,       &
                                    hbri,         ntrcr,     &
                                    nilyr,        nblyr,     &
                                    i_grid,                  &
                                    nbtrcr_sw,    n_zaero,   &
                                    skl_bgc,      z_tracers, &
                                    l_stop,       stop_label)
      
      use ice_constants_colpkg, only: c0, c1, c2, p5
      use ice_colpkg_tracers, only: nt_bgc_N, nt_zaero, tr_bgc_N, &
          tr_zaero, nlt_chl_sw, nlt_zaero_sw
      use ice_colpkg_shared, only: dEdd_algae, bgc_flux_type, sk_l, &
           R_chl2N, min_bgc, F_abs_chl,  hi_ssl
      use ice_zbgc_shared, only: remap_zbgc

      integer (kind=int_kind), intent(in) :: &
         nslyr, & ! number of snow layers
         n_zaero    , & ! number of cells with aicen > puny 
         nbtrcr_sw, n_algae, & ! nilyr+nslyr+2 for chlorophyll
         ntrcr

      integer (kind=int_kind), intent(in) :: &
         nblyr      , & ! number of bio layers
         nilyr          ! number of ice layers 

      real (kind=dbl_kind), dimension (ntrcr), intent(in) ::       &
         trcrn          ! aerosol or chlorophyll

      real (kind=dbl_kind), dimension (nbtrcr_sw), &
         intent(out) ::    &
         trcrn_sw       ! ice on shortwave grid tracers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         sw_grid     , & ! 
         i_grid          ! CICE bio grid 
         
      real(kind=dbl_kind), intent(in) :: &
         hin          , & ! CICE ice thickness
         hbri             ! brine height 

      logical (kind=log_kind), intent(in) :: &
         skl_bgc, & ! skeletal layer bgc  
         z_tracers  ! zbgc   

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      !  local variables

      integer (kind=int_kind) :: k, n, nn

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0, &      ! temporary, remapped tracers
         trtmp

      real (kind=dbl_kind), dimension (nilyr+1):: &
         icegrid        ! correct for large ice surface layers

      real (kind=dbl_kind):: &
         top_conc       ! 1% (min_bgc) of surface concentration 
                        ! when hin > hbri:  just used in sw calculation

      !-----------------------------------------------------------------
      ! Compute aerosols and algal chlorophyll on shortwave grid
      !-----------------------------------------------------------------

      trtmp0(:) = c0
      trtmp(:) = c0
      trcrn_sw(:) = c0

      do k = 1,nilyr+1
         icegrid(k) = sw_grid(k)
      enddo    
      if (sw_grid(1)*hin*c2 > hi_ssl) then
         icegrid(1) = hi_ssl/c2/hin
      endif

      if (z_tracers) then
      if (tr_bgc_N)  then
         do k = 1, nblyr+1
            do n = 1, n_algae
               trtmp0(nt_bgc_N(1) + k-1) = trtmp0(nt_bgc_N(1) + k-1) + &
                                R_chl2N(n)*F_abs_chl(n)*trcrn(nt_bgc_N(n)+k-1)
            enddo ! n
         enddo    ! k
 
         top_conc = trtmp0(nt_bgc_N(1))*min_bgc
         call remap_zbgc (ntrcr,             nilyr+1, &
                          nt_bgc_N(1),                &
                          trtmp0(1:ntrcr  ),          &
                          trtmp (1:ntrcr+2),          &
                          1,                 nblyr+1, &
                          hin,               hbri,    &
                          icegrid(1:nilyr+1),         &
                          i_grid(1:nblyr+1), top_conc, & 
                          l_stop,            stop_label) 

         if (l_stop) return

         do k = 1, nilyr+1
            trcrn_sw(nlt_chl_sw+nslyr+k) = trtmp(nt_bgc_N(1) + k-1)
         enddo       ! k

         do n = 1, n_algae   ! snow contribution
            trcrn_sw(nlt_chl_sw)= trcrn_sw(nlt_chl_sw) &
                     + R_chl2N(n)*F_abs_chl(n)*trcrn(nt_bgc_N(n)+nblyr+1) 
                              ! snow surface layer
            trcrn_sw(nlt_chl_sw+1:nlt_chl_sw+nslyr) = &
                     trcrn_sw(nlt_chl_sw+1:nlt_chl_sw+nslyr) &
                     + R_chl2N(n)*F_abs_chl(n)*trcrn(nt_bgc_N(n)+nblyr+2) 
                              ! only 1 snow layer in zaero
         enddo ! n
      endif    ! tr_bgc_N

      if (tr_zaero) then
         do n = 1, n_zaero

            trtmp0(:) = c0
            trtmp(:) = c0

            do k = 1, nblyr+1
               trtmp0(nt_zaero(n) + k-1) = trcrn(nt_zaero(n)+k-1)
            enddo

            top_conc = trtmp0(nt_zaero(n))*min_bgc
            call remap_zbgc (ntrcr,             nilyr+1, &
                             nt_zaero(n),                &
                             trtmp0(1:ntrcr  ),          &
                             trtmp (1:ntrcr+2),          &
                             1,                 nblyr+1, &
                             hin,               hbri,    &
                             icegrid(1:nilyr+1),         &
                             i_grid(1:nblyr+1), top_conc, &
                             l_stop,            stop_label) 

            if (l_stop) return

            do k = 1,nilyr+1
               trcrn_sw(nlt_zaero_sw(n)+nslyr+k) = trtmp(nt_zaero(n) + k-1)
            enddo
            trcrn_sw(nlt_zaero_sw(n))= trcrn(nt_zaero(n)+nblyr+1) !snow ssl
            trcrn_sw(nlt_zaero_sw(n)+1:nlt_zaero_sw(n)+nslyr)= trcrn(nt_zaero(n)+nblyr+2)
         enddo ! n
      endif    ! tr_zaero
      elseif (skl_bgc) then

         do nn = 1,n_algae
            trcrn_sw(nbtrcr_sw) = trcrn_sw(nbtrcr_sw) &
                                + F_abs_chl(nn)*R_chl2N(nn) &
                                * trcrn(nt_bgc_N(nn))*sk_l/hin &
                                * real(nilyr,kind=dbl_kind)
         enddo 

      endif
      end subroutine compute_shortwave_trcr

!=======================================================================

      end module ice_shortwave

!=======================================================================
