!  SVN:$Id: ice_therm_vertical.F90 1175 2017-03-02 19:53:26Z akt $
!=========================================================================
!
! Update ice and snow internal temperatures and compute
! thermodynamic growth rates and atmospheric fluxes.
!
! NOTE: The thermodynamic calculation is split in two for load balancing.
!       First ice_therm_vertical computes vertical growth rates and coupler
!       fluxes.  Then ice_therm_itd does thermodynamic calculations not
!       needed for coupling.
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW
!          Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_therm_vertical

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c3, p001, p5, puny, &
          pi, depressT, Lvap, hs_min, cp_ice, &
          cp_ocn, rhow, rhoi, rhos, Lfresh, rhofresh, ice_ref_salinity
      use ice_colpkg_shared, only: ktherm, heat_capacity, calc_Tsfc, min_salin
      use ice_therm_shared, only: ferrmax, l_brine, &
                                  calculate_tin_from_qin, Tmin
      use ice_therm_bl99, only: temperature_changes
      use ice_therm_0layer, only: zerolayer_temperature
      use ice_warnings, only: add_warning

      implicit none
      save

      private
      public :: frzmlt_bottom_lateral, thermo_vertical

!=======================================================================

      contains

!=======================================================================
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and atmospheric fluxes.
!
! authors: William H. Lipscomb, LANL
!          C. M. Bitz, UW

      subroutine thermo_vertical (nilyr,       nslyr,     &
                                  dt,          aicen,     &
                                  vicen,       vsnon,     &
                                  Tsf,         zSin,      &
                                  zqin,        zqsn,      &
                                  apond,       hpond,     &
                                  iage,        tr_pond_topo,&
                                  flw,         potT,      &
                                  Qa,          rhoa,      &
                                  fsnow,       fpond,     &
                                  fbot,        Tbot,      &
                                  sss,                    &
                                  lhcoef,      shcoef,    &
                                  fswsfc,      fswint,    &
                                  Sswabs,      Iswabs,    &
                                  fsurfn,      fcondtopn, &
                                  fsensn,      flatn,     &
                                  flwoutn,     evapn,     &
                                  freshn,      fsaltn,    &
                                  fhocnn,      meltt,     &
                                  melts,       meltb,     &
                                  congel,      snoice,    &
                                  mlt_onset,   frz_onset, &
                                  yday,        dsnow,     &
                                  l_stop,      stop_label,&
                                  prescribed_ice)

      use ice_therm_mushy, only: temperature_changes_salinity

      integer (kind=int_kind), intent(in) :: &
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step

      ! ice state variables
      real (kind=dbl_kind), intent(inout) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      ! tracers
      real (kind=dbl_kind), intent(inout) :: &
         Tsf     , & ! ice/snow top surface temp, same as Tsfcn (deg C)
         apond   , & ! melt pond area fraction
         hpond   , & ! melt pond depth (m)
         iage        ! ice age (s)

      logical (kind=log_kind), intent(in) :: &
         tr_pond_topo    ! if .true., use melt pond tracer

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqsn    , & ! snow layer enthalpy, zqsn < 0 (J m-3)
         zqin    , & ! ice layer enthalpy, zqin < 0 (J m-3)
         zSin        ! internal ice layer salinities

      ! input from atmosphere
      real (kind=dbl_kind), &
         intent(in) :: &
         flw     , & ! incoming longwave radiation (W/m^2)
         potT    , & ! air potential temperature  (K) 
         Qa      , & ! specific humidity (kg/kg) 
         rhoa    , & ! air density (kg/m^3) 
         fsnow   , & ! snowfall rate (kg m-2 s-1)
         shcoef  , & ! transfer coefficient for sensible heat
         lhcoef      ! transfer coefficient for latent heat

      real (kind=dbl_kind), &
         intent(inout) :: &
         fswsfc  , & ! SW absorbed at ice/snow surface (W m-2)
         fswint  , & ! SW absorbed in ice interior, below surface (W m-2)
         fpond       ! fresh water flux to ponds (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         Sswabs  , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabs      ! SW radiation absorbed in ice layers (W m-2)

      ! input from ocean
      real (kind=dbl_kind), intent(in) :: &
         fbot    , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot    , & ! ice bottom surface temperature (deg C)
         sss         ! ocean salinity

      ! coupler fluxes to atmosphere
      real (kind=dbl_kind), intent(out):: &
         flwoutn , & ! outgoing longwave radiation (W/m^2) 
         evapn       ! evaporative water flux (kg/m^2/s) 

      ! Note: these are intent out if calc_Tsfc = T, otherwise intent in
      real (kind=dbl_kind), intent(inout):: &
         fsensn   , & ! sensible heat flux (W/m^2) 
         flatn    , & ! latent heat flux   (W/m^2) 
         fsurfn   , & ! net flux to top surface, excluding fcondtopn
         fcondtopn    ! downward cond flux at top surface (W m-2)

      ! coupler fluxes to ocean
      real (kind=dbl_kind), intent(out):: &
         freshn  , & ! fresh water flux to ocean (kg/m^2/s)
         fsaltn  , & ! salt flux to ocean (kg/m^2/s)
         fhocnn      ! net heat flux to ocean (W/m^2) 

      ! diagnostic fields
      real (kind=dbl_kind), &
         intent(inout):: &
         meltt    , & ! top ice melt             (m/step-->cm/day) 
         melts    , & ! snow melt                (m/step-->cm/day) 
         meltb    , & ! basal ice melt           (m/step-->cm/day) 
         congel   , & ! basal ice growth         (m/step-->cm/day) 
         snoice   , & ! snow-ice formation       (m/step-->cm/day) 
         dsnow    , & ! change in snow thickness (m/step-->cm/day) 
         mlt_onset, & ! day of year that sfc melting begins 
         frz_onset    ! day of year that freezing begins (congel or frazil) 

      real (kind=dbl_kind), intent(in) :: &
         yday         ! day of year

      logical (kind=log_kind), intent(out) :: &
         l_stop       ! if true, print diagnostics and abort on return

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local variables

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         dhi         , & ! change in ice thickness
         dhs             ! change in snow thickness

! 2D state variables (thickness, temperature)

      real (kind=dbl_kind) :: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hsn_new     , & ! thickness of new snow (m)
         worki       , & ! local work array
         works           ! local work array

      real (kind=dbl_kind), dimension (nilyr) :: &
         zTin            ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (nslyr) :: &
         zTsn            ! internal snow layer temperatures

! other 2D flux and energy variables

      real (kind=dbl_kind) :: &
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         einit       , & ! initial energy of melting (J m-2)
         efinal      , & ! final energy of melting (J m-2)
         einter          ! intermediate energy

      real (kind=dbl_kind) :: &
         fadvocn ! advective heat flux to ocean

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      flwoutn = c0
      evapn   = c0
      freshn  = c0
      fsaltn  = c0
      fhocnn  = c0
      fadvocn = c0

      meltt   = c0
      meltb   = c0
      melts   = c0
      congel  = c0
      snoice  = c0
      dsnow   = c0

      if (calc_Tsfc) then
         fsensn  = c0
         flatn     = c0
         fsurfn    = c0
         fcondtopn = c0
      endif

      !-----------------------------------------------------------------
      ! Compute variables needed for vertical thermo calculation
      !-----------------------------------------------------------------

      call init_vertical_profile (nilyr,    nslyr,   &
                                  aicen,             &
                                  vicen,    vsnon,   &
                                  hin,      hilyr,   &
                                  hsn,      hslyr,   &
                                  zqin,     zTin,    &
                                  zqsn,     zTsn,    &
                                  zSin,              &
                                  einit,    Tbot,    &
                                  l_stop,   stop_label)

      if (l_stop) return

      ! Save initial ice and snow thickness (for fresh and fsalt)
      worki = hin
      works = hsn

      !-----------------------------------------------------------------
      ! Compute new surface temperature and internal ice and snow
      !  temperatures.
      !-----------------------------------------------------------------

      if (heat_capacity) then   ! usual case

         if (ktherm == 2) then

            call temperature_changes_salinity(dt,                   & 
                                              nilyr,     nslyr,     &
                                              rhoa,      flw,       &
                                              potT,      Qa,        &
                                              shcoef,    lhcoef,    &
                                              fswsfc,    fswint,    &
                                              Sswabs,    Iswabs,    &
                                              hilyr,     hslyr,     &
                                              apond,     hpond,     &
                                              zqin,      zTin,      &
                                              zqsn,      zTsn,      &
                                              zSin,                 &
                                              Tsf,       Tbot,      &
                                              sss,                  &
                                              fsensn,    flatn,     &
                                              flwoutn,   fsurfn,    &
                                              fcondtopn, fcondbot,  &
                                              fadvocn,   snoice,    &
                                              einit,                &
                                              l_stop,    stop_label)
               
            if (l_stop) return

         else ! ktherm

            call temperature_changes(dt,                   &  
                                     nilyr,     nslyr,     &
                                     rhoa,      flw,       &
                                     potT,      Qa,        &
                                     shcoef,    lhcoef,    &
                                     fswsfc,    fswint,    &
                                     Sswabs,    Iswabs,    &
                                     hilyr,     hslyr,     &
                                     zqin,      zTin,      &
                                     zqsn,      zTsn,      &
                                     zSin,                 &
                                     Tsf,       Tbot,      &
                                     fsensn,    flatn,     &
                                     flwoutn,   fsurfn,    &
                                     fcondtopn, fcondbot,  &
                                     einit,     l_stop,    &
                                     stop_label)

            if (l_stop) return

         endif ! ktherm
            
      else

         if (calc_Tsfc) then       

            call zerolayer_temperature(dt,                  & 
                                       nilyr,     nslyr,    &
                                       rhoa,      flw,      &
                                       potT,      Qa,       &
                                       shcoef,    lhcoef,   &
                                       fswsfc,              &
                                       hilyr,     hslyr,    &
                                       Tsf,       Tbot,     &
                                       fsensn,    flatn,    &
                                       flwoutn,   fsurfn,   &
                                       fcondtopn, fcondbot, &
                                       l_stop,    stop_label)

            if (l_stop) return

         else

            !------------------------------------------------------------
            ! Set fcondbot = fcondtop for zero layer thermodynamics
            ! fcondtop is set in call to set_sfcflux in step_therm1
            !------------------------------------------------------------

            fcondbot  = fcondtopn   ! zero layer         
      
         endif      ! calc_Tsfc

      endif         ! heat_capacity

      ! intermediate energy for error check
      
      einter = c0
      do k = 1, nslyr
         einter = einter + hslyr * zqsn(k)
      enddo ! k
      do k = 1, nilyr
         einter = einter + hilyr * zqin(k)
      enddo ! k

      if (l_stop) return

      !-----------------------------------------------------------------
      ! Compute growth and/or melting at the top and bottom surfaces.
      ! Add new snowfall.
      ! Repartition ice into equal-thickness layers, conserving energy.
      !----------------------------------------------------------------- 

      call thickness_changes(nilyr,       nslyr,     &
                             dt,          yday,      &
                             efinal,                 &
                             hin,         hilyr,     &
                             hsn,         hslyr,     &
                             zqin,        zqsn,      &
                             fbot,        Tbot,      &
                             flatn,       fsurfn,    &
                             fcondtopn,   fcondbot,  &
                             fsnow,       hsn_new,   &
                             fhocnn,      evapn,     &
                             meltt,       melts,     &
                             meltb,       iage,      &
                             congel,      snoice,    &
                             mlt_onset,   frz_onset, &
                             zSin,        sss,       &
                             dsnow)

      !-----------------------------------------------------------------
      ! Check for energy conservation by comparing the change in energy
      ! to the net energy input
      !-----------------------------------------------------------------

      call conservation_check_vthermo(dt,                  &
                                      fsurfn,    flatn,    &
                                      fhocnn,    fswint,   &
                                      fsnow,     einit,    &
                                      einter,    efinal,   &
                                      fcondtopn, fcondbot, &
                                      fadvocn,   fbot,     &
                                      l_stop,    stop_label)
      
      if (l_stop) return

      !-----------------------------------------------------------------
      ! If prescribed ice, set hi back to old values
      !-----------------------------------------------------------------

#ifdef CCSMCOUPLED
      if (present(prescribed_ice)) then
          if (prescribed_ice) then
            hin    = worki
            fhocnn = c0             ! for diagnostics
          endif
      endif
#endif

      !-----------------------------------------------------------------
      ! Compute fluxes of water and salt from ice to ocean.
      ! evapn < 0 => sublimation, evapn > 0 => condensation
      ! aerosol flux is accounted for in ice_aerosol.F90
      !-----------------------------------------------------------------
            
      dhi = hin - worki
      dhs = hsn - works - hsn_new
      
      freshn = freshn + evapn - (rhoi*dhi + rhos*dhs) / dt
      fsaltn = fsaltn - rhoi*dhi*ice_ref_salinity*p001/dt
      fhocnn = fhocnn + fadvocn ! for ktherm=2 

      if (hin == c0) then
         if (tr_pond_topo) fpond = fpond - aicen * apond * hpond
      endif

      !-----------------------------------------------------------------
      !  Given the vertical thermo state variables, compute the new ice 
      !   state variables.
      !-----------------------------------------------------------------

      call update_state_vthermo(nilyr,   nslyr,   &
                                Tbot,    Tsf,     &     
                                hin,     hsn,     &
                                zqin,    zSin,    &
                                zqsn,             &
                                aicen,            &
                                vicen,   vsnon)

      !-----------------------------------------------------------------
      ! Reload passive tracer array
      !-----------------------------------------------------------------

    end subroutine thermo_vertical

!=======================================================================
!
! Compute heat flux to bottom surface.
! Compute fraction of ice that melts laterally.
!
! authors C. M. Bitz, UW
!         William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL

      subroutine frzmlt_bottom_lateral (dt,       ncat,     &
                                        nilyr,    nslyr,    &
                                        aice,     frzmlt,   &
                                        vicen,    vsnon,    &
                                        qicen,    qsnon,    &
                                        sst,      Tf,       & 
                                        ustar_min,          &
                                        fbot_xfer_type,     &
                                        strocnxT, strocnyT, &
                                        Tbot,     fbot,     &
                                        rside,    Cdn_ocn)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt                  ! time step

      real (kind=dbl_kind), intent(in) :: &
         aice    , & ! ice concentration
         frzmlt  , & ! freezing/melting potential (W/m^2)
         sst     , & ! sea surface temperature (C)
         Tf      , & ! freezing temperature (C)
         ustar_min,& ! minimum friction velocity for ice-ocean heat flux
         Cdn_ocn , & ! ocean-ice neutral drag coefficient
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

      character (char_len), intent(in) :: &
         fbot_xfer_type  ! transfer coefficient type for ice-ocean heat flux

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen   , & ! ice volume (m)
         vsnon       ! snow volume (m)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         qicen   , & ! ice layer enthalpy (J m-3)
         qsnon       ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), intent(out) :: &
         Tbot    , & ! ice bottom surface temperature (deg C)
         fbot    , & ! heat flux to ice bottom  (W/m^2)
         rside       ! fraction of ice that melts laterally

      ! local variables

      integer (kind=int_kind) :: &
         n              , & ! thickness category index
         k                  ! layer index

      real (kind=dbl_kind) :: &
         etot    , & ! total energy in column
         fside       ! lateral heat flux (W/m^2)

      real (kind=dbl_kind) :: &
         deltaT    , & ! SST - Tbot >= 0
         ustar     , & ! skin friction velocity for fbot (m/s)
         wlat      , & ! lateral melt rate (m/s)
         xtmp          ! temporary variable

      ! Parameters for bottom melting

      real (kind=dbl_kind) :: &
         cpchr         ! -cp_ocn*rhow*exchange coefficient

      ! Parameters for lateral melting

      real (kind=dbl_kind), parameter :: &
         floediam = 300.0_dbl_kind, & ! effective floe diameter (m)
         floeshape = 0.66_dbl_kind , & ! constant from Steele (unitless)
         m1 = 1.6e-6_dbl_kind     , & ! constant from Maykut & Perovich
                                      ! (m/s/deg^(-m2))
         m2 = 1.36_dbl_kind           ! constant from Maykut & Perovich
                                      ! (unitless)

      !-----------------------------------------------------------------
      ! Identify grid cells where ice can melt.
      !-----------------------------------------------------------------

      rside = c0
      Tbot  = Tf
      fbot  = c0
      
      if (aice > puny .and. frzmlt < c0) then ! ice can melt
         
         fside = c0

      !-----------------------------------------------------------------
      ! Use boundary layer theory for fbot.
      ! See Maykut and McPhee (1995): JGR, 100, 24,691-24,703.
      !-----------------------------------------------------------------

         deltaT = max((sst-Tbot),c0)
         
         ! strocnx has units N/m^2 so strocnx/rho has units m^2/s^2
         ustar = sqrt (sqrt(strocnxT**2+strocnyT**2)/rhow)
         ustar = max (ustar,ustar_min)

         if (trim(fbot_xfer_type) == 'Cdn_ocn') then
            ! Note: Cdn_ocn has already been used for calculating ustar 
            ! (formdrag only) --- David Schroeder (CPOM)
            cpchr = -cp_ocn*rhow*Cdn_ocn
         else ! fbot_xfer_type == 'constant'
            ! 0.006 = unitless param for basal heat flx ala McPhee and Maykut
            cpchr = -cp_ocn*rhow*0.006_dbl_kind
         endif

         fbot = cpchr * deltaT * ustar ! < 0
         fbot = max (fbot, frzmlt) ! frzmlt < fbot < 0
            
!!! uncomment to use all frzmlt for standalone runs
   !     fbot = min (c0, frzmlt)

      !-----------------------------------------------------------------
      ! Compute rside.  See these references:
      !    Maykut and Perovich (1987): JGR, 92, 7032-7044
      !    Steele (1992): JGR, 97, 17,729-17,738
      !-----------------------------------------------------------------

         wlat = m1 * deltaT**m2 ! Maykut & Perovich
         rside = wlat*dt*pi/(floeshape*floediam) ! Steele
         rside = max(c0,min(rside,c1))

      !-----------------------------------------------------------------
      ! Compute heat flux associated with this value of rside.
      !-----------------------------------------------------------------

         do n = 1, ncat
            
            etot = c0
            
            ! melting energy/unit area in each column, etot < 0
            
            do k = 1, nslyr
               etot = etot + qsnon(k,n) * vsnon(n)/real(nslyr,kind=dbl_kind)
            enddo
            
            do k = 1, nilyr
               etot = etot + qicen(k,n) * vicen(n)/real(nilyr,kind=dbl_kind)
            enddo                  ! nilyr
            
            ! lateral heat flux
            fside = fside + rside*etot/dt ! fside < 0
            
         enddo                     ! n
         
      !-----------------------------------------------------------------
      ! Limit bottom and lateral heat fluxes if necessary.
      !-----------------------------------------------------------------
            
         xtmp = frzmlt/(fbot + fside + puny) 
         xtmp = min(xtmp, c1)
         fbot  = fbot  * xtmp
         rside = rside * xtmp
         
      endif

      end subroutine frzmlt_bottom_lateral

!=======================================================================
!
! Given the state variables (vicen, vsnon, zqin, etc.),
! compute variables needed for the vertical thermodynamics
! (hin, hsn, zTin, etc.)
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine init_vertical_profile(nilyr,    nslyr,    &
                                       aicen,    vicen,    &
                                       vsnon,              &
                                       hin,      hilyr,    &
                                       hsn,      hslyr,    &
                                       zqin,     zTin,     &
                                       zqsn,     zTsn,     &
                                       zSin,               &
                                       einit,    Tbot,     &
                                       l_stop,   stop_label)

      use ice_mushy_physics, only: temperature_mush, &
                                   liquidus_temperature_mush, &
                                   enthalpy_of_melting

      use ice_constants_colpkg, only: p1 !!!AKT Column!!!

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)
 
      real (kind=dbl_kind), intent(out):: &
         hilyr       , & ! ice layer thickness
         hslyr       , & ! snow layer thickness
         einit           ! initial energy of melting (J m-2)
 
      real (kind=dbl_kind), intent(in):: &
         Tbot            ! bottom ice temp  (C)

      real (kind=dbl_kind), intent(out):: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zTin            ! internal ice layer temperatures

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zSin            ! internal ice layer salinities
        
      real (kind=dbl_kind), dimension (:), &
         intent(out) :: &
         zqsn        , & ! snow enthalpy
         zTsn            ! snow temperature

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label      ! abort error message

      ! local variables
      real (kind=dbl_kind), dimension(nilyr) :: &
         Tmlts           ! melting temperature

      integer (kind=int_kind) :: &
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         rnslyr,       & ! real(nslyr)
         Tmax            ! maximum allowed snow/ice temperature (deg C)

      logical (kind=log_kind) :: &   ! for vector-friendly error checks
         tsno_high   , & ! flag for zTsn > Tmax
         tice_high   , & ! flag for zTin > Tmlt
         tsno_low    , & ! flag for zTsn < Tmin
         tice_low        ! flag for zTin < Tmin

      character(len=char_len_long) :: &
         warning ! warning message

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rnslyr = real(nslyr,kind=dbl_kind)

      tsno_high = .false.
      tice_high = .false.
      tsno_low  = .false.
      tice_low  = .false.

      einit = c0
 
      !-----------------------------------------------------------------
      ! Surface temperature, ice and snow thickness
      ! Initialize internal energy
      !-----------------------------------------------------------------

      hin    = vicen / aicen
      hsn    = vsnon / aicen
      hilyr    = hin / real(nilyr,kind=dbl_kind)
      hslyr    = hsn / rnslyr

      !-----------------------------------------------------------------
      ! Snow enthalpy and maximum allowed snow temperature
      ! If heat_capacity = F, zqsn and zTsn are used only for checking
      ! conservation.
      !-----------------------------------------------------------------

      do k = 1, nslyr

      !-----------------------------------------------------------------
      ! Tmax based on the idea that dT ~ dq / (rhos*cp_ice)
      !                             dq ~ q dv / v
      !                             dv ~ puny = eps11
      ! where 'd' denotes an error due to roundoff.
      !-----------------------------------------------------------------

         if (hslyr > hs_min/rnslyr .and. heat_capacity) then
            ! zqsn < 0              
            Tmax = -zqsn(k)*puny*rnslyr / &
                 (rhos*cp_ice*vsnon)
         else
            zqsn  (k) = -rhos * Lfresh
            Tmax = puny
         endif

      !-----------------------------------------------------------------
      ! Compute snow temperatures from enthalpies.
      ! Note: zqsn <= -rhos*Lfresh, so zTsn <= 0.
      !-----------------------------------------------------------------
         zTsn(k) = (Lfresh + zqsn(k)/rhos)/cp_ice

      !-----------------------------------------------------------------
      ! Check for zTsn > Tmax (allowing for roundoff error) and zTsn < Tmin.
      !-----------------------------------------------------------------
         if (zTsn(k) > Tmax) then
            tsno_high = .true.
         elseif (zTsn(k) < Tmin) then
            tsno_low  = .true.
         endif

      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! If zTsn is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

      if (tsno_high .and. heat_capacity) then
         do k = 1, nslyr

            if (hslyr > hs_min/rnslyr) then
               Tmax = -zqsn(k)*puny*rnslyr / &
                    (rhos*cp_ice*vsnon)
            else
               Tmax = puny
            endif

            if (zTsn(k) > Tmax) then
               write(warning,*) ' '
               call add_warning(warning)
               write(warning,*) 'Starting thermo, zTsn > Tmax'
               call add_warning(warning)
               write(warning,*) 'zTsn=',zTsn(k)
               call add_warning(warning)
               write(warning,*) 'Tmax=',Tmax
               call add_warning(warning)
               write(warning,*) 'zqsn',zqsn(k),-Lfresh*rhos,zqsn(k)+Lfresh*rhos
               call add_warning(warning)
               l_stop = .true.
               stop_label = "init_vertical_profile: Starting thermo, zTsn > Tmax"
               return
            endif

         enddo                  ! nslyr
      endif                     ! tsno_high

      if (tsno_low .and. heat_capacity) then
         do k = 1, nslyr

            if (zTsn(k) < Tmin) then ! allowing for roundoff error
               write(warning,*) ' '
               call add_warning(warning)
               write(warning,*) 'Starting thermo, zTsn < Tmin'
               call add_warning(warning)
               write(warning,*) 'zTsn=', zTsn(k)
               call add_warning(warning)
               write(warning,*) 'Tmin=', Tmin
               call add_warning(warning)
               write(warning,*) 'zqsn', zqsn(k)
               call add_warning(warning)
               write(warning,*) hin
               call add_warning(warning)
               write(warning,*) hsn
               call add_warning(warning)
               l_stop = .true.
               stop_label = "init_vertical_profile: Starting thermo, zTsn < Tmin"
               return
            endif

         enddo                  ! nslyr
      endif                     ! tsno_low

      do k = 1, nslyr

         if (zTsn(k) > c0) then   ! correct roundoff error
            zTsn(k) = c0
            zqsn(k) = -rhos*Lfresh
         endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------
         einit = einit + hslyr*zqsn(k)

      enddo                     ! nslyr

      do k = 1, nilyr

      !---------------------------------------------------------------------
      !  Use initial salinity profile for thin ice
      !---------------------------------------------------------------------

         if (ktherm == 1 .and. zSin(k) < min_salin-puny) then
            write(warning,*) ' '
            call add_warning(warning)
            write(warning,*) 'Starting zSin < min_salin, layer', k
            call add_warning(warning)
            write(warning,*) 'zSin =', zSin(k)
            call add_warning(warning)
            write(warning,*) 'min_salin =', min_salin
            call add_warning(warning)
            l_stop = .true.
            stop_label = "init_vertical_profile: Starting zSin < min_salin, layer"
            return
         endif
         
         if (ktherm == 2) then
            Tmlts(k) = liquidus_temperature_mush(zSin(k))
         else
            Tmlts(k) = -zSin(k) * depressT
         endif

      !-----------------------------------------------------------------
      ! Compute ice enthalpy
      ! If heat_capacity = F, zqin and zTin are used only for checking
      ! conservation.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------
         
         if (ktherm == 2) then
            zTin(k) = temperature_mush(zqin(k),zSin(k))
         else
            zTin(k) = calculate_Tin_from_qin(zqin(k),Tmlts(k))
         endif

         if (l_brine) then
            Tmax = Tmlts(k)
         else                ! fresh ice
            Tmax = -zqin(k)*puny/(rhos*cp_ice*vicen)
         endif

      !-----------------------------------------------------------------
      ! Check for zTin > Tmax and zTin < Tmin
      !-----------------------------------------------------------------
         if (zTin(k) > Tmax) then
            tice_high = .true.
         elseif (zTin(k) < Tmin) then
            tice_low  = .true.
         endif

      !-----------------------------------------------------------------
      ! If zTin is out of bounds, print diagnostics and exit.
      !-----------------------------------------------------------------

         if (tice_high .and. heat_capacity) then

            if (l_brine) then
               Tmax = Tmlts(k)
            else             ! fresh ice
               Tmax = -zqin(k)*puny/(rhos*cp_ice*vicen)
            endif

            if (zTin(k) > Tmax) then
               write(warning,*) ' '
               call add_warning(warning)
               write(warning,*) 'Starting thermo, T > Tmax, layer', k
               call add_warning(warning)
               write(warning,*) 'k:', k
               call add_warning(warning)
               write(warning,*) 'zTin =',zTin(k),', Tmax=',Tmax
               call add_warning(warning)
               write(warning,*) 'zSin =',zSin(k)
               call add_warning(warning)
               write(warning,*) 'hin =',hin
               call add_warning(warning)
               write(warning,*) 'zqin =',zqin(k)
               call add_warning(warning)
               write(warning,*) 'qmlt=',enthalpy_of_melting(zSin(k))
               call add_warning(warning)
               write(warning,*) 'Tmlt=',Tmlts(k)
               call add_warning(warning)
               
               if (ktherm == 2) then
                  zqin(k) = enthalpy_of_melting(zSin(k)) - c1
                  zTin(k) = temperature_mush(zqin(k),zSin(k))
                  write(warning,*) 'Corrected quantities'
                  call add_warning(warning)
                  write(warning,*) 'zqin=',zqin(k)
                  call add_warning(warning)
                  write(warning,*) 'zTin=',zTin(k)
                  call add_warning(warning)
               else
                  l_stop = .true.
                  stop_label = "init_vertical_profile: Starting thermo, T > Tmax, layer"
                  return
               endif
            endif
         endif                  ! tice_high

         if (tice_low .and. heat_capacity) then

            if (zTin(k) < Tmin) then
               write(warning,*) ' '
               call add_warning(warning)
               write(warning,*) 'Starting thermo T < Tmin, layer', k
               call add_warning(warning)
               write(warning,*) 'zTin =', zTin(k)
               call add_warning(warning)
               write(warning,*) 'Tmin =', Tmin
               call add_warning(warning)
               l_stop = .true.
               stop_label = "init_vertical_profile: Starting thermo, T < Tmin, layer"
               return
            endif
         endif                  ! tice_low

      !-----------------------------------------------------------------
      ! correct roundoff error
      !-----------------------------------------------------------------

         if (ktherm /= 2) then

            if (zTin(k) > c0) then
               zTin(k) = c0
               zqin(k) = -rhoi*Lfresh
            endif

         endif

! echmod: is this necessary?
!         if (ktherm == 1) then
!               if (zTin(k)>= -zSin(k)*depressT) then
!                   zTin(k) = -zSin(k)*depressT - puny
!                   zqin(k) = -rhoi*cp_ocn*zSin(k)*depressT
!               endif
!         endif

      !-----------------------------------------------------------------
      ! initial energy per unit area of ice/snow, relative to 0 C
      !-----------------------------------------------------------------

         einit = einit + hilyr*zqin(k) 

      enddo                     ! nilyr

      end subroutine init_vertical_profile

!=======================================================================
!
! Compute growth and/or melting at the top and bottom surfaces.
! Convert snow to ice if necessary.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine thickness_changes (nilyr,     nslyr,    &
                                    dt,        yday,     &
                                    efinal,              & 
                                    hin,       hilyr,    &
                                    hsn,       hslyr,    &
                                    zqin,      zqsn,     &
                                    fbot,      Tbot,     &
                                    flatn,     fsurfn,   &
                                    fcondtopn, fcondbot, &
                                    fsnow,     hsn_new,  &
                                    fhocnn,    evapn,    &
                                    meltt,     melts,    &
                                    meltb,     iage,     &
                                    congel,    snoice,   &  
                                    mlt_onset, frz_onset,&
                                    zSin,      sss,      &
                                    dsnow)

      use ice_colpkg_shared, only: phi_i_mushy
      use ice_mushy_physics, only: enthalpy_mush, enthalpy_of_melting, &
                           temperature_mush, liquidus_temperature_mush, &
                           liquid_fraction

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         yday            ! day of the year

      real (kind=dbl_kind), intent(in) :: &
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         Tbot        , & ! ice bottom surface temperature (deg C)
         fsnow       , & ! snowfall rate (kg m-2 s-1)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), intent(inout) :: &
         fcondbot        ! downward cond flux at bottom surface (W m-2)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zqsn            ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr           ! snow layer thickness (m)

      real (kind=dbl_kind), intent(inout) :: &
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         dsnow       , & ! snow  formation          (m/step-->cm/day)
         iage        , & ! ice age (s)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(inout) :: &
         hin         , & ! total ice thickness (m)
         hsn             ! total snow thickness (m)

      real (kind=dbl_kind), intent(out):: &
         efinal          ! final energy of melting (J m-2)

      real (kind=dbl_kind), intent(out):: &
         fhocnn      , & ! fbot, corrected for any surplus energy (W m-2)
         evapn           ! ice/snow mass sublimated/condensed (kg m-2 s-1)

      real (kind=dbl_kind), intent(out):: &
         hsn_new         ! thickness of new snow (m)

      ! changes to zSin in this subroutine are not reloaded into the
      ! trcrn array for ktherm /= 2, so we could remove ktherm=2 conditionals
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zSin            ! ice layer salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         sss             ! ocean salinity (PSU) 

      ! local variables

      real (kind=dbl_kind), parameter :: &
         qbotmax = -p5*rhoi*Lfresh  ! max enthalpy of ice growing at bottom

      integer (kind=int_kind) :: &
         k               ! vertical index

      real (kind=dbl_kind) :: &
         esub        , & ! energy for sublimation, > 0    (J m-2)
         econ        , & ! energy for condensation, < 0   (J m-2)
         etop_mlt    , & ! energy for top melting, > 0    (J m-2)
         ebot_mlt    , & ! energy for bottom melting, > 0 (J m-2)
         ebot_gro    , & ! energy for bottom growth, < 0  (J m-2)
         emlt_atm    , & ! total energy of brine, from atmosphere (J m-2)
         emlt_ocn        ! total energy of brine, to ocean        (J m-2)

      real (kind=dbl_kind) :: &
         dhi         , & ! change in ice thickness
         dhs         , & ! change in snow thickness
         Ti          , & ! ice temperature
         Ts          , & ! snow temperature
         qbot        , & ! enthalpy of ice growing at bottom surface (J m-3)
         qsub        , & ! energy/unit volume to sublimate ice/snow (J m-3)
         hqtot       , & ! sum of h*q for two layers
         wk1         , & ! temporary variable
         zqsnew      , & ! enthalpy of new snow (J m-3)
         hstot       , & ! snow thickness including new snow (m)
         Tmlts           ! melting temperature

      real (kind=dbl_kind), dimension (nilyr+1) :: &
         zi1         , & ! depth of ice layer boundaries (m)
         zi2             ! adjusted depths, with equal hilyr (m)

      real (kind=dbl_kind), dimension (nslyr+1) :: &
         zs1         , & ! depth of snow layer boundaries (m)
         zs2             ! adjusted depths, with equal hslyr (m)

      real (kind=dbl_kind), dimension (nilyr) :: &
         dzi             ! ice layer thickness after growth/melting

      real (kind=dbl_kind), dimension (nslyr) :: &
         dzs             ! snow layer thickness after growth/melting

      real (kind=dbl_kind), dimension (nilyr) :: &
         qm          , & ! energy of melting (J m-3) = zqin in BL99 formulation
         qmlt            ! enthalpy of melted ice (J m-3) = zero in BL99 formulation

      real (kind=dbl_kind) :: &
         qbotm       , &
         qbotp       , &
         qbot0

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      hsn_new  = c0

      do k = 1, nilyr
         dzi(k) = hilyr
      enddo

      do k = 1, nslyr
         dzs(k) = hslyr
      enddo

      do k = 1, nilyr
         if (ktherm == 2) then
            qmlt(k) = enthalpy_of_melting(zSin(k))
         else
            qmlt(k) = c0
         endif
         qm(k) = zqin(k) - qmlt(k)
         emlt_atm = c0
         emlt_ocn = c0
      enddo

      !-----------------------------------------------------------------
      ! For l_brine = false (fresh ice), check for temperatures > 0.
      !  Melt ice or snow as needed to bring temperatures back to 0.
      ! For l_brine = true, this should not be necessary.
      !-----------------------------------------------------------------

      if (.not. l_brine) then 

         do k = 1, nslyr
            Ts = (Lfresh + zqsn(k)/rhos) / cp_ice
            if (Ts > c0) then
               dhs = cp_ice*Ts*dzs(k) / Lfresh
               dzs(k) = dzs(k) - dhs
               zqsn(k) = -rhos*Lfresh
            endif
         enddo

         do k = 1, nilyr
            Ti = (Lfresh + zqin(k)/rhoi) / cp_ice
            if (Ti > c0) then
               dhi = cp_ice*Ti*dzi(k) / Lfresh
               dzi(k) = dzi(k) - dhi
               zqin(k) = -rhoi*Lfresh
            endif
         enddo                  ! k

      endif                     ! .not. l_brine

      !-----------------------------------------------------------------
      ! Compute energy available for sublimation/condensation, top melt,
      ! and bottom growth/melt.
      !-----------------------------------------------------------------

      wk1 = -flatn * dt
      esub = max(wk1, c0)     ! energy for sublimation, > 0
      econ = min(wk1, c0)     ! energy for condensation, < 0

      wk1 = (fsurfn - fcondtopn) * dt
      etop_mlt = max(wk1, c0)           ! etop_mlt > 0

      wk1 = (fcondbot - fbot) * dt
      ebot_mlt = max(wk1, c0)           ! ebot_mlt > 0
      ebot_gro = min(wk1, c0)           ! ebot_gro < 0

      !--------------------------------------------------------------
      ! Condensation (evapn > 0)
      ! Note: evapn here has unit of kg/m^2.  Divide by dt later.
      ! This is the only case in which energy from the atmosphere
      ! is used for changes in the brine energy (emlt_atm).
      !--------------------------------------------------------------

      evapn = c0          ! initialize

      if (hsn > puny) then    ! add snow with enthalpy zqsn(1)
         dhs = econ / (zqsn(1) - rhos*Lvap) ! econ < 0, dhs > 0
         dzs(1) = dzs(1) + dhs
         evapn = evapn + dhs*rhos
      else                        ! add ice with enthalpy zqin(1)
         dhi = econ / (qm(1) - rhoi*Lvap) ! econ < 0, dhi > 0
         dzi(1) = dzi(1) + dhi
         evapn = evapn + dhi*rhoi
         ! enthalpy of melt water
         emlt_atm = emlt_atm - qmlt(1) * dhi 
      endif

      !--------------------------------------------------------------
      ! Grow ice (bottom)
      !--------------------------------------------------------------

      if (ktherm == 2) then

         qbotm = enthalpy_mush(Tbot, sss)
         qbotp = -Lfresh * rhoi * (c1 - phi_i_mushy)
         qbot0 = qbotm - qbotp

         dhi = ebot_gro / qbotp     ! dhi > 0

         hqtot = dzi(nilyr)*zqin(nilyr) + dhi*qbotm
         hstot = dzi(nilyr)*zSin(nilyr) + dhi*sss
         emlt_ocn = emlt_ocn - qbot0 * dhi

      else

         Tmlts = -zSin(nilyr) * depressT 

         ! enthalpy of new ice growing at bottom surface
         if (heat_capacity) then
            if (l_brine) then
               qbot = -rhoi * (cp_ice * (Tmlts-Tbot) &
                    + Lfresh * (c1-Tmlts/Tbot) &
                    - cp_ocn * Tmlts)
               qbot = min (qbot, qbotmax)      ! in case Tbot is close to Tmlt
            else
               qbot = -rhoi * (-cp_ice * Tbot + Lfresh)
            endif
         else   ! zero layer
            qbot = -rhoi * Lfresh
         endif

         dhi = ebot_gro / qbot     ! dhi > 0

         hqtot = dzi(nilyr)*zqin(nilyr) + dhi*qbot
         hstot = c0
      endif ! ktherm

      dzi(nilyr) = dzi(nilyr) + dhi
      if (dzi(nilyr) > puny) then
         zqin(nilyr) = hqtot / dzi(nilyr)
         if (ktherm == 2) then
            zSin(nilyr) = hstot / dzi(nilyr)
            qmlt(nilyr) = enthalpy_of_melting(zSin(nilyr))
         else
            qmlt(nilyr) = c0
         endif
         qm(nilyr) = zqin(nilyr) - qmlt(nilyr)
      endif

      ! update ice age due to freezing (new ice age = dt)
      !         if (tr_iage) &
      !            iage = (iage*hin + dt*dhi) / (hin + dhi)

      ! history diagnostics
      congel = congel + dhi
      if (dhi > puny .and. frz_onset < puny) &
           frz_onset = yday

      do k = 1, nslyr

         !--------------------------------------------------------------
         ! Remove internal snow melt 
         !--------------------------------------------------------------
         
         if (ktherm == 2 .and. zqsn(k) > -rhos * Lfresh) then

            dhs = max(-dzs(k), &
                -((zqsn(k) + rhos*Lfresh) / (rhos*Lfresh)) * dzs(k))
            dzs(k) = dzs(k) + dhs
            zqsn(k) = -rhos * Lfresh
            melts = melts - dhs
            ! delta E = zqsn(k) + rhos * Lfresh

         endif

         !--------------------------------------------------------------
         ! Sublimation of snow (evapn < 0)
         !--------------------------------------------------------------

         qsub = zqsn(k) - rhos*Lvap ! qsub < 0
         dhs  = max (-dzs(k), esub/qsub)  ! esub > 0, dhs < 0
         dzs(k) = dzs(k) + dhs
         esub = esub - dhs*qsub
         esub = max(esub, c0)   ! in case of roundoff error
         evapn = evapn + dhs*rhos

         !--------------------------------------------------------------
         ! Melt snow (top)
         !--------------------------------------------------------------

         dhs = max(-dzs(k), etop_mlt/zqsn(k))
         dzs(k) = dzs(k) + dhs         ! zqsn < 0, dhs < 0
         etop_mlt = etop_mlt - dhs*zqsn(k)
         etop_mlt = max(etop_mlt, c0) ! in case of roundoff error

         ! history diagnostics
         if (dhs < -puny .and. mlt_onset < puny) &
              mlt_onset = yday
         melts = melts - dhs

      enddo                     ! nslyr

      do k = 1, nilyr

         !--------------------------------------------------------------
         ! Sublimation of ice (evapn < 0)
         !--------------------------------------------------------------

         qsub = qm(k) - rhoi*Lvap              ! qsub < 0
         dhi  = max (-dzi(k), esub/qsub) ! esub < 0, dhi < 0
         dzi(k) = dzi(k) + dhi
         esub = esub - dhi*qsub
         esub = max(esub, c0)
         evapn = evapn + dhi*rhoi
         emlt_ocn = emlt_ocn - qmlt(k) * dhi 

         !--------------------------------------------------------------
         ! Melt ice (top)
         !--------------------------------------------------------------
   
         if (qm(k) < c0) then
            dhi = max(-dzi(k), etop_mlt/qm(k))
         else
            qm(k) = c0
            dhi = -dzi(k)
         endif
         emlt_ocn = emlt_ocn - max(zqin(k),qmlt(k)) * dhi

         dzi(k) = dzi(k) + dhi         ! zqin < 0, dhi < 0
         etop_mlt = max(etop_mlt - dhi*qm(k), c0)

         ! history diagnostics
         if (dhi < -puny .and. mlt_onset < puny) &
              mlt_onset = yday
         meltt = meltt - dhi

      enddo                     ! nilyr

      do k = nilyr, 1, -1

         !--------------------------------------------------------------
         ! Melt ice (bottom)
         !--------------------------------------------------------------

         if (qm(k) < c0) then
            dhi = max(-dzi(k), ebot_mlt/qm(k))
         else
            qm(k) = c0
            dhi = -dzi(k)
         endif
         emlt_ocn = emlt_ocn - max(zqin(k),qmlt(k)) * dhi

         dzi(k) = dzi(k) + dhi         ! zqin < 0, dhi < 0
         ebot_mlt = max(ebot_mlt - dhi*qm(k), c0)

         ! history diagnostics 
         meltb = meltb -dhi

      enddo                     ! nilyr

      do k = nslyr, 1, -1

         !--------------------------------------------------------------
         ! Melt snow (only if all the ice has melted)
         !--------------------------------------------------------------
         
         dhs = max(-dzs(k), ebot_mlt/zqsn(k))
         dzs(k) = dzs(k) + dhs         ! zqsn < 0, dhs < 0
         ebot_mlt = ebot_mlt - dhs*zqsn(k)
         ebot_mlt = max(ebot_mlt, c0)

      enddo                     ! nslyr

      !-----------------------------------------------------------------
      ! Compute heat flux used by the ice (<=0).
      ! fhocn is the available ocean heat that is left after use by ice
      !-----------------------------------------------------------------

      fhocnn = fbot &
             + (esub + etop_mlt + ebot_mlt)/dt

!---!-----------------------------------------------------------------
!---! Add new snowfall at top surface.
!---!-----------------------------------------------------------------

      !----------------------------------------------------------------
      ! NOTE: If heat flux diagnostics are to work, new snow should
      !       have T = 0 (i.e. q = -rhos*Lfresh) and should not be
      !       converted to rain.
      !----------------------------------------------------------------

      if (fsnow > c0) then

         hsn_new = fsnow/rhos * dt
         zqsnew = -rhos*Lfresh
         hstot = dzs(1) + hsn_new

         if (hstot > c0) then
            zqsn(1) =  (dzs(1) * zqsn(1) &
                    + hsn_new * zqsnew) / hstot
            ! avoid roundoff errors
            zqsn(1) = min(zqsn(1), -rhos*Lfresh)

            dzs(1) = hstot
         endif
      endif

    !-----------------------------------------------------------------
    ! Find the new ice and snow thicknesses.
    !-----------------------------------------------------------------

      hin = c0
      hsn = c0

      do k = 1, nilyr
         hin = hin + dzi(k)
      enddo                     ! k

      do k = 1, nslyr
         hsn = hsn + dzs(k)
         dsnow = dsnow + dzs(k) - hslyr  
      enddo                     ! k

    !-------------------------------------------------------------------
    ! Convert snow to ice if snow lies below freeboard.
    !-------------------------------------------------------------------

      if (ktherm /= 2) &
           call freeboard (nslyr,    dt,       &
                           snoice,   iage,     &
                           hin,      hsn,      &
                           zqin,     zqsn,     &
                           dzi,      dzs,      &
                           dsnow)

!---!-------------------------------------------------------------------
!---! Repartition the ice and snow into equal-thickness layers,
!---! conserving energy.
!---!-------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Compute desired layer thicknesses.
      !-----------------------------------------------------------------
 
      if (hin > c0) then
         hilyr = hin / real(nilyr,kind=dbl_kind)
      else
         hin = c0
         hilyr = c0
      endif
      if (hsn > c0) then
         hslyr = hsn / real(nslyr,kind=dbl_kind)
      else
         hsn = c0
         hslyr = c0
      endif

      !-----------------------------------------------------------------
      ! Compute depths zi1 of old layers (unequal thickness).
      ! Compute depths zi2 of new layers (equal thickness).
      !-----------------------------------------------------------------

      zi1(1) = c0
      zi1(1+nilyr) = hin

      zi2(1) = c0
      zi2(1+nilyr) = hin

      if (heat_capacity) then

         do k = 1, nilyr-1
            zi1(k+1) = zi1(k) + dzi(k)
            zi2(k+1) = zi2(k) + hilyr
         enddo

        !-----------------------------------------------------------------
        ! Conserving energy, compute the enthalpy of the new equal layers.
        !-----------------------------------------------------------------
            
         call adjust_enthalpy (nilyr,              &
                               zi1,      zi2,      &
                               hilyr,    hin,      &
                               zqin)   

         if (ktherm == 2) &
              call adjust_enthalpy (nilyr,              &
                                    zi1,      zi2,      &
                                    hilyr,    hin,      &
                                    zSin)   

      else ! zero layer (nilyr=1)

         zqin(1) = -rhoi * Lfresh
         zqsn(1) = -rhos * Lfresh
       
      endif

      if (nslyr > 1) then

      !-----------------------------------------------------------------
      ! Compute depths zs1 of old layers (unequal thickness).
      ! Compute depths zs2 of new layers (equal thickness).
      !-----------------------------------------------------------------

         zs1(1) = c0
         zs1(1+nslyr) = hsn
         
         zs2(1) = c0
         zs2(1+nslyr) = hsn
         
         do k = 1, nslyr-1
            zs1(k+1) = zs1(k) + dzs(k)
            zs2(k+1) = zs2(k) + hslyr
         enddo

      !-----------------------------------------------------------------
      ! Conserving energy, compute the enthalpy of the new equal layers.
      !-----------------------------------------------------------------

         call adjust_enthalpy (nslyr,              &
                               zs1,      zs2,      &
                               hslyr,    hsn,      &
                               zqsn)   

      endif   ! nslyr > 1

      !-----------------------------------------------------------------
      ! Remove very thin snow layers (ktherm = 2)
      !-----------------------------------------------------------------

      if (ktherm == 2) then
         do k = 1, nslyr
            if (hsn <= puny) then
               fhocnn = fhocnn &
                      + zqsn(k)*hsn/(real(nslyr,kind=dbl_kind)*dt)
               zqsn(k) = -rhos*Lfresh
               hslyr = c0
            endif
         enddo
      endif

      !-----------------------------------------------------------------
      ! Compute final ice-snow energy, including the energy of
      !  sublimated/condensed ice.
      !-----------------------------------------------------------------

      efinal = -evapn*Lvap
      evapn =  evapn/dt

      do k = 1, nslyr
         efinal = efinal + hslyr*zqsn(k)
      enddo

      do k = 1, nilyr
         efinal = efinal + hilyr*zqin(k)
      enddo                     ! k

      if (ktherm < 2) then
         emlt_atm = c0
         emlt_ocn = c0
      endif

      ! melt water is no longer zero enthalpy with ktherm=2
      fhocnn = fhocnn + emlt_ocn/dt
      efinal = efinal + emlt_atm ! for conservation check

      end subroutine thickness_changes

!=======================================================================
!
! If there is enough snow to lower the ice/snow interface below
! sea level, convert enough snow to ice to bring the interface back
! to sea level.
!
! authors William H. Lipscomb, LANL
!         Elizabeth C. Hunke, LANL

      subroutine freeboard (nslyr,    dt,       &
                            snoice,             &
                            iage,               &
                            hin,      hsn,      &
                            zqin,     zqsn,     &
                            dzi,      dzs,      &
                            dsnow)

      integer (kind=int_kind), intent(in) :: &
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      real (kind=dbl_kind), &
         intent(inout) :: &
         snoice  , & ! snow-ice formation       (m/step-->cm/day)
         dsnow   , & ! change in snow thickness after snow-ice formation (m)
         iage        ! ice age (s)

      real (kind=dbl_kind), &
         intent(inout) :: &
         hin     , & ! ice thickness (m)
         hsn         ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         zqsn        ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin    , & ! ice layer enthalpy (J m-3)
         dzi     , & ! ice layer thicknesses (m)
         dzs         ! snow layer thicknesses (m)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! vertical index

      real (kind=dbl_kind) :: &
         dhin        , & ! change in ice thickness (m)
         dhsn        , & ! change in snow thickness (m)
         hqs             ! sum of h*q for snow (J m-2)

      real (kind=dbl_kind) :: &
         wk1         , & ! temporary variable
         dhs             ! snow to remove from layer (m)

      !-----------------------------------------------------------------
      ! Determine whether snow lies below freeboard.
      !-----------------------------------------------------------------
      
      dhin = c0
      dhsn = c0
      hqs  = c0

      wk1 = hsn - hin*(rhow-rhoi)/rhos

      if (wk1 > puny .and. hsn > puny) then  ! snow below freeboard
         dhsn = min(wk1*rhoi/rhow, hsn) ! snow to remove
         dhin = dhsn * rhos/rhoi        ! ice to add
      endif

      !-----------------------------------------------------------------
      ! Adjust snow layer thickness.
      ! Compute energy to transfer from snow to ice.
      !-----------------------------------------------------------------

      do k = nslyr, 1, -1
         if (dhin > puny) then
            dhs = min(dhsn, dzs(k)) ! snow to remove from layer
            hsn = hsn - dhs
            dsnow = dsnow -dhs   !new snow addition term 
            dzs(k) = dzs(k) - dhs
            dhsn = dhsn - dhs
            dhsn = max(dhsn,c0)
            hqs = hqs + dhs * zqsn(k)
         endif               ! dhin > puny
      enddo

      !-----------------------------------------------------------------
      ! Transfer volume and energy from snow to top ice layer.
      !-----------------------------------------------------------------

      if (dhin > puny) then
         ! update ice age due to freezing (new ice age = dt)
         !            if (tr_iage) &
         !               iage = (iage*hin+dt*dhin)/(hin+dhin)
         
         wk1 = dzi(1) + dhin
         hin = hin + dhin
         zqin(1) = (dzi(1)*zqin(1) + hqs) / wk1
         dzi(1) = wk1

         ! history diagnostic
         snoice = snoice + dhin
      endif               ! dhin > puny

      end subroutine freeboard

!=======================================================================
!
! Conserving energy, compute the new enthalpy of equal-thickness ice
! or snow layers.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW

      subroutine adjust_enthalpy (nlyr,               &
                                  z1,       z2,       &
                                  hlyr,     hn,       &
                                  qn)

      integer (kind=int_kind), intent(in) :: &
         nlyr            ! number of layers (nilyr or nslyr)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         z1          , & ! interface depth for old, unequal layers (m)
         z2              ! interface depth for new, equal layers (m)

      real (kind=dbl_kind), intent(in) :: &
         hlyr            ! new layer thickness (m)

      real (kind=dbl_kind), intent(in) :: &
         hn              ! total thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         qn              ! layer quantity (enthalpy, salinity...)

      ! local variables

      integer (kind=int_kind) :: &
         k, k1, k2       ! vertical indices

      real (kind=dbl_kind) :: &
         hovlp           ! overlap between old and new layers (m)

      real (kind=dbl_kind) :: &
         rhlyr           ! 1./hlyr

      real (kind=dbl_kind), dimension (nlyr) :: &
         hq              ! h * q for a layer

      !-----------------------------------------------------------------
      ! Compute reciprocal layer thickness.
      !-----------------------------------------------------------------

      rhlyr = c0
      if (hn > puny) rhlyr = c1 / hlyr

      !-----------------------------------------------------------------
      ! Compute h*q for new layers (k2) given overlap with old layers (k1)
      !-----------------------------------------------------------------

      do k2 = 1, nlyr
         hq(k2) = c0
      enddo                     ! k
      k1 = 1
      k2 = 1
      do while (k1 <= nlyr .and. k2 <= nlyr)
         hovlp = min (z1(k1+1), z2(k2+1)) &
               - max (z1(k1),   z2(k2))
         hovlp = max (hovlp, c0)
         hq(k2) = hq(k2) + hovlp*qn(k1)
         if (z1(k1+1) > z2(k2+1)) then
            k2 = k2 + 1
         else
            k1 = k1 + 1
         endif
      enddo                  ! while

      !-----------------------------------------------------------------
      ! Compute new enthalpies.
      !-----------------------------------------------------------------

      do k = 1, nlyr
         qn(k) = hq(k) * rhlyr
      enddo                     ! k

      end subroutine adjust_enthalpy

!=======================================================================
!
! Check for energy conservation by comparing the change in energy
! to the net energy input.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Adrian K. Turner, LANL

      subroutine conservation_check_vthermo(dt,                 &
                                            fsurfn,   flatn,    &
                                            fhocnn,   fswint,   &
                                            fsnow,              &
                                            einit,    einter,   &
                                            efinal,             &
                                            fcondtopn,fcondbot, &
                                            fadvocn,  fbot,     &
                                            l_stop,   stop_label)

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step

      real (kind=dbl_kind), intent(in) :: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         flatn       , & ! surface downward latent heat (W m-2)
         fhocnn      , & ! fbot, corrected for any surplus energy
         fswint      , & ! SW absorbed in ice interior, below surface (W m-2)
         fsnow       , & ! snowfall rate (kg m-2 s-1)
         fcondtopn   , &
         fadvocn     , &
         fbot           

      real (kind=dbl_kind), intent(in) :: &
         einit       , & ! initial energy of melting (J m-2)
         einter      , & ! intermediate energy of melting (J m-2)
         efinal      , & ! final energy of melting (J m-2)
         fcondbot

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local variables

      real (kind=dbl_kind) :: &
         einp        , & ! energy input during timestep (J m-2)
         ferr            ! energy conservation error (W m-2)

      character(len=char_len_long) :: &
         warning ! warning message

      !----------------------------------------------------------------
      ! If energy is not conserved, print diagnostics and exit.
      !----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Note that fsurf - flat = fsw + flw + fsens; i.e., the latent
      ! heat is not included in the energy input, since (efinal - einit)
      ! is the energy change in the system ice + vapor, and the latent
      ! heat lost by the ice is equal to that gained by the vapor.
      !-----------------------------------------------------------------
       
      einp = (fsurfn - flatn + fswint - fhocnn &
           - fsnow*Lfresh - fadvocn) * dt
      ferr = abs(efinal-einit-einp) / dt

      if (ferr > ferrmax) then
         l_stop = .true.
         stop_label = "conservation_check_vthermo: Thermo energy conservation error"

         write(warning,*) 'Thermo energy conservation error'
         call add_warning(warning)
         write(warning,*) 'Flux error (W/m^2) =', ferr
         call add_warning(warning)
         write(warning,*) 'Energy error (J) =', ferr*dt
         call add_warning(warning)
         write(warning,*) 'Initial energy =', einit
         call add_warning(warning)
         write(warning,*) 'Final energy   =', efinal
         call add_warning(warning)
         write(warning,*) 'efinal - einit  =', efinal-einit
         call add_warning(warning)
         write(warning,*) 'fsurfn,flatn,fswint,fhocn, fsnow*Lfresh:'
         call add_warning(warning)
         write(warning,*) fsurfn,flatn,fswint,fhocnn, fsnow*Lfresh
         call add_warning(warning)
         write(warning,*) 'Input energy =', einp
         call add_warning(warning)
         write(warning,*) 'fbot,fcondbot:'
         call add_warning(warning)
         write(warning,*) fbot,fcondbot
         call add_warning(warning)

         !         if (ktherm == 2) then
         write(warning,*) 'Intermediate energy =', einter
         call add_warning(warning)
         write(warning,*) 'efinal - einter =', &
              efinal-einter
         call add_warning(warning)
         write(warning,*) 'einter - einit  =', &
              einter-einit
         call add_warning(warning)
         write(warning,*) 'Conduction Error =', (einter-einit) &
              - (fcondtopn*dt - fcondbot*dt + fswint*dt)
         call add_warning(warning)
         write(warning,*) 'Melt/Growth Error =', (einter-einit) &
              + ferr*dt - (fcondtopn*dt - fcondbot*dt + fswint*dt)
         call add_warning(warning)
         write(warning,*) 'Advection Error =', fadvocn*dt
         call add_warning(warning)
         !         endif

         !         write(warning,*) fsurfn,flatn,fswint,fhocnn
         !         call add_warning(warning)
         
         write(warning,*) 'dt*(fsurfn, flatn, fswint, fhocn, fsnow*Lfresh, fadvocn):'
         call add_warning(warning)
         write(warning,*) fsurfn*dt, flatn*dt, &
              fswint*dt, fhocnn*dt, &
              fsnow*Lfresh*dt, fadvocn*dt
         call add_warning(warning)
         return
      endif

      end subroutine conservation_check_vthermo

!=======================================================================
!
! Given the vertical thermo state variables (hin, hsn),
! compute the new ice state variables (vicen, vsnon).
! Zero out state variables if ice has melted entirely.
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!         Elizabeth C. Hunke, LANL

      subroutine update_state_vthermo(nilyr,    nslyr,    &
                                      Tf,       Tsf,      &
                                      hin,      hsn,      &
                                      zqin,     zSin,     &
                                      zqsn,               &
                                      aicen,    vicen,    &
                                      vsnon)

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         Tf              ! freezing temperature (C)

      real (kind=dbl_kind), intent(inout) :: &
         Tsf             ! ice/snow surface temperature, Tsfcn

      real (kind=dbl_kind), intent(in) :: &
         hin         , & ! ice thickness (m)
         hsn             ! snow thickness (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! ice salinity    (ppt)
         zqsn            ! snow layer enthalpy (J m-3)

      real (kind=dbl_kind), intent(inout) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon           ! volume per unit area of snow         (m)

      ! local variables

      integer (kind=int_kind) :: &
         k               ! ice layer index

      if (hin <= c0) then
         aicen = c0
         vicen = c0
         vsnon = c0
         Tsf = Tf
         do k = 1, nilyr
            zqin(k) = c0
         enddo
         if (ktherm == 2) then
            do k = 1, nilyr
               zSin(k) = c0
            enddo
         endif
         do k = 1, nslyr
            zqsn(k) = c0
         enddo
      else
         ! aicen is already up to date
         vicen = aicen * hin
         vsnon = aicen * hsn
      endif
      
      end subroutine update_state_vthermo

!=======================================================================

      end module ice_therm_vertical

!=======================================================================
