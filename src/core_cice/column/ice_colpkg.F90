!  SVN:$Id: ice_colpkg.F90 1012 2015-06-26 12:34:09Z eclare $
!=========================================================================
!
! flags and interface routines for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module ice_colpkg

      use ice_kinds_mod
      use ice_colpkg_shared ! namelist and other parameters

      implicit none

      private

      ! initialization
      public :: &
           colpkg_init_itd, &
           colpkg_init_itd_hist, &
           colpkg_init_thermo, &
           colpkg_init_orbit, &
           colpkg_init_trcr, &
           colpkg_init_parameters, &
           colpkg_init_tracer_flags, &
           colpkg_init_tracer_indices, &
           colpkg_init_tracer_numbers

      ! time stepping
      public :: &
           colpkg_step_therm1, &
           colpkg_step_therm2, &
           colpkg_prep_radiation, &
           colpkg_step_radiation, &
           colpkg_step_ridge

      ! other column routines
      public :: &
           colpkg_aggregate, &
           colpkg_ice_strength, &
           colpkg_atm_boundary, &
           colpkg_ocn_mixed_layer

      ! temperature inquiry functions
      public :: &
           colpkg_ice_temperature, &
           colpkg_snow_temperature, &
           colpkg_liquidus_temperature, &
           colpkg_sea_freezing_temperature, &
           colpkg_enthalpy_snow

!=======================================================================

      contains

!=======================================================================
!     Initialization routines
!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd(ncat, hin_max, nu_diag)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_therm_shared, only: hi_min
      use ice_constants_colpkg, only: p01, p1, c0, c1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat, & ! number of thickness categories
           nu_diag ! diagnostic output file unit number

      real (kind=dbl_kind), intent(out) :: &
           hin_max(0:ncat)  ! category limits (m)

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2           , & !
           b1           , & ! parameters for kcatbound = 3
           b2           , & !
           b3

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat
      b1 = p1         ! asymptotic category width (m)
      b2 = c3         ! thickness for which participation function is small (m)
      b3 = max(rncat*(rncat-1), c2*b2/b1)

      hi_min = p01    ! minimum ice thickness allowed (m) for thermo
                      ! note hi_min is reset to 0.1 for kitd=0, below

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of four options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5 
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !
      ! The fourth formula asymptotes to a particular category width as
      ! the number of categories increases, given by the parameter b1.
      ! The parameter b3 is computed so that the category boundaries
      ! are even numbers.
      !
      !    H(n) = b1 * [n + b3*n*(n+1)/(2*N*(N-1))] for N=ncat
      !
      ! kcatbound=-1 is available only for 1-category runs, with
      ! boundaries 0 and 100 m.
      !-----------------------------------------------------------------

      if (kcatbound == -1) then ! single category
         hin_max(0) = c0
         hin_max(1) = c100

      elseif (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
            hi_min = p1    ! minimum ice thickness allowed (m) for thermo
            cc1 = max(1.1_dbl_kind/rncat,hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         ! thinnest 3 categories combined
         data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         ! thinnest 2 categories combined
         data wmo6 / 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind /
!echmod wmo6a
!         data wmo6 /0.30_dbl_kind, 0.70_dbl_kind,  &
!                    1.20_dbl_kind, 2.00_dbl_kind,  &
!                    4.56729_dbl_kind, &
!                    999._dbl_kind /

         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         ! all thickness categories 
         data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         write (nu_diag,*) 'kcatbound=2 (WMO) must have ncat=5, 6 or 7'
         stop
       endif

      elseif (kcatbound == 3) then  ! asymptotic scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = b1 * (rn + b3*rn*(rn+c1)/(c2*rncat*(rncat-c1)))
         enddo

      endif ! kcatbound

      end subroutine colpkg_init_itd

!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd_hist (ncat, hin_max, nu_diag, c_hi_range)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_constants_colpkg, only: p01, p1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat, & ! number of thickness categories
           nu_diag ! diagnostic output file unit number

      real (kind=dbl_kind), intent(in) :: &
           hin_max(0:ncat)  ! category limits (m)

      character (len=35), intent(out) :: &
           c_hi_range(ncat) ! string for history output

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

         write (nu_diag,*) ' '
         write (nu_diag,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         do n = 1, ncat
            write (nu_diag,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            ! Write integer n to character string
            write (c_nc, '(i2)') n    

            ! Write hin_max to character string
            write (c_hinmax1, '(f6.3)') hin_max(n-1)
            write (c_hinmax2, '(f6.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '//c_hinmax2//'m'
         enddo
         write (nu_diag,*) ' '

      end subroutine colpkg_init_itd_hist

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine colpkg_init_thermo(nilyr, sprofile)

      use ice_colpkg_shared, only: saltmax, ktherm, heat_capacity
      use ice_constants_colpkg, only: p1, p5, c0, c1, c2, pi!!!AKT Column!!!
      use ice_therm_shared, only: l_brine
      !!!AKT Column!!!use ice_zbgc_shared, only: min_salin
      real(kind=dbl_kind), parameter :: min_salin = p1!!!AKT Column!!!

      integer (kind=int_kind), intent(in) :: &
         nilyr                            ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         sprofile                         ! vertical salinity profile

      real (kind=dbl_kind), parameter :: &
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      ! Set l_brine to false for zero layer thermodynamics
      !-----------------------------------------------------------------

      heat_capacity = .true.      
      if (ktherm == 0) heat_capacity = .false. ! 0-layer thermodynamics

      l_brine = .false.
      if (saltmax > min_salin .and. heat_capacity) l_brine = .true.

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            sprofile(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
            sprofile(k) = max(sprofile(k), min_salin)
         enddo ! k
         sprofile(nilyr+1) = saltmax

      else ! .not. l_brine
         do k = 1, nilyr+1
            sprofile(k) = c0
         enddo
      endif ! l_brine

      end subroutine colpkg_init_thermo

!=======================================================================

! Compute orbital parameters for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 

      subroutine colpkg_init_orbit(nu_diag, l_stop, stop_label)

      use ice_constants_colpkg, only: iyear_AD, eccen, obliqr, lambm0, &
         mvelpp, obliq, mvelp, decln, eccf, log_print
#ifdef CCSMCOUPLED
      use shr_orb_mod, only: shr_orb_params
#else
      use ice_orbital, only: shr_orb_params
#endif

      integer (kind=int_kind), intent(in) :: &
         nu_diag         ! diagnostic file unit number

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort the model

      character (char_len), intent(out) :: stop_label

      iyear_AD  = 1950
      log_print = .false.   ! if true, write out orbital parameters

      call shr_orb_params( iyear_AD, eccen , obliq , mvelp    , &
                           obliqr  , lambm0, mvelpp, log_print, &
                           nu_diag , l_stop, stop_label)
 
      end subroutine colpkg_init_orbit
 
!=======================================================================

      subroutine colpkg_init_trcr(Tair,     Tf,       &
                                  Sprofile, Tprofile, &
                                  Tsfc,               &
                                  nilyr,    nslyr,    &
                                  qin,      qsn)

      use ice_colpkg_shared, only: calc_Tsfc
      use ice_constants_colpkg, only: Tsmelt, Tffresh, p5, cp_ice, cp_ocn, &
          Lfresh, rhoi, rhos, c0, c1
      use ice_mushy_physics, only: enthalpy_mush

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         Tair, &     ! air temperature (C)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         Sprofile, & ! vertical salinity profile (ppt)
         Tprofile    ! vertical temperature profile (C)

      real (kind=dbl_kind), intent(out) :: &
         Tsfc        ! surface temperature (C)

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         qin, &      ! ice enthalpy profile (J/m3)
         qsn         ! snow enthalpy profile (J/m3)

      ! local variables

      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: &
         slope, Ti

            ! surface temperature
            Tsfc = Tf ! default
            if (calc_Tsfc) Tsfc = min(Tsmelt, Tair - Tffresh) ! deg C

            if (heat_capacity) then

               ! ice enthalpy
               do k = 1, nilyr
                  ! assume linear temp profile and compute enthalpy
                  slope = Tf - Tsfc
                  Ti = Tsfc + slope*(real(k,kind=dbl_kind)-p5) &
                                    /real(nilyr,kind=dbl_kind)
                  if (ktherm == 2) then
                     qin(k) = enthalpy_mush(Ti, Sprofile(k))
                  else
                     qin(k) = -(rhoi * (cp_ice*(Tprofile(k)-Ti) &
                         + Lfresh*(c1-Tprofile(k)/Ti) - cp_ocn*Tprofile(k)))
                  endif
               enddo               ! nilyr

               ! snow enthalpy
               do k = 1, nslyr
                  Ti = min(c0, Tsfc)
                  qsn(k) = -rhos*(Lfresh - cp_ice*Ti)
               enddo               ! nslyr

            else  ! one layer with zero heat capacity

               ! ice energy
               qin(1) = -rhoi * Lfresh 

               ! snow energy
               qsn(1) = -rhos * Lfresh 

            endif               ! heat_capacity

      end subroutine colpkg_init_trcr

!=======================================================================
!     Temperature functions
!=======================================================================

      function colpkg_liquidus_temperature(Sin) result(Tmlt)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: liquidus_temperature_mush

        real(dbl_kind), intent(in) :: Sin
        real(dbl_kind) :: Tmlt

        if (ktherm == 2) then

           Tmlt = liquidus_temperature_mush(Sin)

        else

           Tmlt = -depressT * Sin

        endif

      end function colpkg_liquidus_temperature

!=======================================================================

      function colpkg_sea_freezing_temperature(sss) result(Tf)

        use ice_colpkg_shared, only: tfrz_option
        use ice_constants_colpkg, only: depressT

        real(dbl_kind), intent(in) :: sss
        real(dbl_kind) :: Tf

        if (trim(tfrz_option) == 'mushy') then

           Tf = colpkg_liquidus_temperature(sss) ! deg C
           
        elseif (trim(tfrz_option) == 'linear_salt') then

           Tf = -depressT * sss ! deg C

        else

           Tf = -1.8_dbl_kind

        endif

      end function colpkg_sea_freezing_temperature

!=======================================================================

      function colpkg_ice_temperature(qin, Sin) result(Tin)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: temperature_mush
        use ice_therm_shared, only: calculate_Tin_from_qin

        real(kind=dbl_kind), intent(in) :: qin, Sin
        real(kind=dbl_kind) :: Tin

        real(kind=dbl_kind) :: Tmlts

        if (ktherm == 2) then

           Tin = temperature_mush(qin, Sin)

        else

           Tmlts = -depressT * Sin
           Tin = calculate_Tin_from_qin(qin,Tmlts)

        endif

      end function colpkg_ice_temperature

!=======================================================================

      function colpkg_snow_temperature(qin) result(Tsn)

        use ice_colpkg_shared, only: ktherm
        use ice_mushy_physics, only: temperature_snow
        use ice_constants_colpkg, only: Lfresh, rhos, cp_ice

        real(kind=dbl_kind), intent(in) :: qin
        real(kind=dbl_kind) :: Tsn

        if (ktherm == 2) then

           Tsn = temperature_snow(qin)

        else

           Tsn = (Lfresh + qin/rhos)/cp_ice

        endif

      end function colpkg_snow_temperature

!=======================================================================

      function colpkg_enthalpy_snow(zTsn) result(qsn)

        use ice_mushy_physics, only: enthalpy_snow

        real(kind=dbl_kind), intent(in) :: zTsn
        real(kind=dbl_kind) :: qsn

        qsn = enthalpy_snow(zTsn)

      end function colpkg_enthalpy_snow

!=======================================================================
!     Time-stepping routines
!=======================================================================

! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm1(dt, ncat, nilyr, nslyr,     &
                                    aicen_init  ,               &
                                    vicen_init  , vsnon_init  , &
                                    aice        , aicen       , &
                                    vice        , vicen       , &
                                    vsno        , vsnon       , &
                                    uvel        , vvel        , &
                                    Tsfc        , zqsn        , &
                                    zqin        , zSin        , &
                                    alvl        , vlvl        , &
                                    apnd        , hpnd        , &
                                    ipnd        ,               &
                                    iage        , FY          , &
                                    aerosno     , aeroice     , &
                                    uatm        , vatm        , &
                                    wind        , zlvl        , &
                                    Qa          , rhoa        , &
                                    Tair        , Tref        , &
                                    Qref        , Uref        , &
                                    Cdn_atm_ratio,              &
                                    Cdn_ocn     , Cdn_ocn_skin, &
                                    Cdn_ocn_floe, Cdn_ocn_keel, &
                                    Cdn_atm     , Cdn_atm_skin, &
                                    Cdn_atm_floe, Cdn_atm_pond, &
                                    Cdn_atm_rdg , hfreebd     , &
                                    hdraft      , hridge      , &
                                    distrdg     , hkeel       , &
                                    dkeel       , lfloe       , &
                                    dfloe       ,               &
                                    strax       , stray       , &
                                    strairxT    , strairyT    , &
                                    potT        , sst         , &
                                    sss         , Tf          , &
                                    strocnxT    , strocnyT    , &
                                    frzmlt      , rside       , &
                                    fsnow       , frain       , &
                                    fpond       ,               &
                                    fsurf       , fsurfn      , &
                                    fcondtop    , fcondtopn   , &
                                    fswsfcn     , fswintn     , &
                                    fswthrun    , fswabs      , &
                                    flwout      ,               &
                                    Sswabsn     , Iswabsn     , &
                                    flw         , coszen      , & 
                                    fsens       , fsensn      , &
                                    flat        , flatn       , &
                                    evap        ,               &
                                    fresh       , fsalt       , &
                                    fhocn       , fswthru     , &
                                    flatn_f     , fsensn_f    , &
                                    fsurfn_f    , fcondtopn_f , &
                                    faero_atm   , faero_ocn   , &
                                    dhsn        , ffracn      , &
                                    meltt       , melttn      , &
                                    meltb       , meltbn      , &
                                    meltl       ,               &
                                    melts       , meltsn      , &
                                    congel      , congeln     , &
                                    snoice      , snoicen     , &
                                    dsnown      , frazil      , &
                                    lmask_n     , lmask_s     , &
                                    mlt_onset   , frz_onset   , &
                                    yday        , l_stop      , &
                                    stop_label  , nu_diag)

      !!!AKT Column!!!use ice_aerosol, only: update_aerosol
      use ice_atmo, only: neutral_drag_coeffs
      use ice_age, only: increment_age
      use ice_constants_colpkg, only: rhofresh, rhoi, rhos, c0, c1, puny
      use ice_firstyear, only: update_FYarea
      use ice_flux_colpkg, only: set_sfcflux, merge_fluxes
      use ice_meltpond_cesm, only: compute_ponds_cesm
      use ice_meltpond_lvl, only: compute_ponds_lvl
      use ice_meltpond_topo, only: compute_ponds_topo
      use ice_therm_shared, only: hi_min
      use ice_therm_vertical, only: frzmlt_bottom_lateral, thermo_vertical
      use ice_colpkg_tracers, only: tr_iage, tr_FY, tr_aero, tr_pond, tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         uvel        , & ! x-component of velocity (m/s)
         vvel        , & ! y-component of velocity (m/s)
         strax       , & ! wind stress components (N/m^2)
         stray       , & ! 
         yday            ! day of year

      logical (kind=log_kind), intent(in) :: &
         lmask_n     , & ! northern hemisphere mask
         lmask_s         ! southern hemisphere mask

      real (kind=dbl_kind), intent(inout) :: &
         aice        , & ! sea ice concentration
         vice        , & ! volume per unit area of ice          (m)
         vsno        , & ! volume per unit area of snow         (m)
         zlvl        , & ! atm level height (m)
         uatm        , & ! wind velocity components (m/s)
         vatm        , &
         wind        , & ! wind speed (m/s)
         potT        , & ! air potential temperature  (K)
         Tair        , & ! air temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         rhoa        , & ! air density (kg/m^3)
         frain       , & ! rainfall rate (kg/m^2 s)
         fsnow       , & ! snowfall rate (kg/m^2 s)
         fpond       , & ! fresh water flux to ponds (kg/m^2/s)
         fresh       , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt       , & ! salt flux to ocean (kg/m^2/s)
         fhocn       , & ! net heat flux to ocean (W/m^2)
         fswthru     , & ! shortwave penetrating to ocean (W/m^2)
         fsurf       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop    , & ! top surface conductive flux        (W/m^2)
         fsens       , & ! sensible heat flux (W/m^2)
         flat        , & ! latent heat flux   (W/m^2)
         fswabs      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         coszen      , & ! cosine solar zenith angle, < 0 for sun below horizon 
         flw         , & ! incoming longwave radiation (W/m^2)
         flwout      , & ! outgoing longwave radiation (W/m^2)
         evap        , & ! evaporative water flux (kg/m^2/s)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         frazil      , & ! frazil ice growth        (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         Tref        , & ! 2m atm reference temperature (K)
         Qref        , & ! 2m atm reference spec humidity (kg/kg)
         Uref        , & ! 10m atm reference wind speed (m/s)
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
         hfreebd     , & ! freeboard (m)
         hdraft      , & ! draft of ice + snow column (Stoessel1993)
         hridge      , & ! ridge height
         distrdg     , & ! distance between ridges
         hkeel       , & ! keel depth
         dkeel       , & ! distance between keels
         lfloe       , & ! floe length
         dfloe       , & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg , & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio,& ! ratio drag atm / neutral drag atm
         strairxT    , & ! stress on ice by air, x-direction
         strairyT    , & ! stress on ice by air, y-direction
         strocnxT    , & ! ice-ocean stress, x-direction
         strocnyT    , & ! ice-ocean stress, y-direction
         frzmlt      , & ! freezing/melting potential (W/m^2)
         rside       , & ! fraction of ice that melts laterally
         sst         , & ! sea surface temperature (C)
         Tf          , & ! freezing temperature (C)
         sss         , & ! sea surface salinity (ppt)
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         meltl       , & ! lateral ice melt         (m/step-->cm/day)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init  , & ! fractional area of ice
         vicen_init  , & ! volume per unit area of ice (m)
         vsnon_init  , & ! volume per unit area of snow (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         alvl        , & ! level ice area fraction
         vlvl        , & ! level ice volume fraction
         apnd        , & ! melt pond area fraction
         hpnd        , & ! melt pond depth (m)
         ipnd        , & ! melt pond refrozen lid thickness (m)
         iage        , & ! volume-weighted ice age
         FY          , & ! area-weighted first-year ice area
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         flatn       , & ! latent heat flux (W m-2)
         fsensn      , & ! sensible heat flux (W m-2)
         fsurfn_f    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f , & ! downward cond flux at top surface (W m-2)
         flatn_f     , & ! latent heat flux (W m-2)
         fsensn_f    , & ! sensible heat flux (W m-2)
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm   , & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn   , & ! aerosol flux to ocean  (kg/m^2/s)
         dhsn        , & ! depth difference for snow on sea ice and pond ice
         ffracn      , & ! fraction of fsurfn used to melt ipond
         meltsn      , & ! snow melt                       (m)
         melttn      , & ! top ice melt                    (m)
         meltbn      , & ! bottom ice melt                 (m)
         congeln     , & ! congelation ice growth          (m)
         snoicen     , & ! snow-ice growth                 (m)
         dsnown          ! change in snow thickness (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! internal ice layer salinities
         Sswabsn     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: &
         aerosno    , &  ! snow aerosol tracer (kg/m^2)
         aeroice         ! ice aerosol tracer (kg/m^2)

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=char_len), intent(out) :: &
         stop_label

    integer (kind=int_kind), intent(in) :: &
         nu_diag         ! file unit number (diagnostic only)

      ! local variables

      integer (kind=int_kind) :: &
         n               ! category index

      real (kind=dbl_kind) :: &
         worka, workb    ! temporary variables

      ! 2D coupler variables (computed for each category, then aggregated)
      real (kind=dbl_kind) :: &
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Cdn_atm_ratio_n, & ! drag coefficient ratio
         Trefn       , & ! air tmp reference level                (K)
         Urefn       , & ! air speed reference level            (m/s)
         Qrefn       , & ! air sp hum reference level         (kg/kg)
         Tbot        , & ! ice bottom surface temperature (deg C)
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         rfrac           ! water fraction retained for melt ponds

      real (kind=dbl_kind) :: &
         raice       , & ! 1/aice
         pond            ! water retained in ponds (m)

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

      call frzmlt_bottom_lateral (dt,        ncat,      &
                                  nilyr,     nslyr,     &
                                  aice,      frzmlt,    &
                                  vicen (:), vsnon (:), &
                                  zqin(:,:), zqsn(:,:), &
                                  sst,       Tf,        &
                                  ustar_min,            &
                                  fbot_xfer_type,       &
                                  strocnxT,  strocnyT,  &
                                  Tbot,      fbot,      &
                                  rside,     Cdn_ocn)
      
      !-----------------------------------------------------------------
      ! Update the neutral drag coefficients to account for form drag
      ! Oceanic and atmospheric drag coefficients
      !-----------------------------------------------------------------

      if (formdrag) then
         call neutral_drag_coeffs (apnd      (:), &
                                   hpnd     (:), ipnd      (:), &
                                   alvl     (:), vlvl      (:), &
                                   aice        , vice,          &
                                   vsno        , aicen     (:), &
                                   vicen    (:), vsnon     (:), &
                                   Cdn_ocn     , Cdn_ocn_skin, &
                                   Cdn_ocn_floe, Cdn_ocn_keel, &
                                   Cdn_atm     , Cdn_atm_skin, &
                                   Cdn_atm_floe, Cdn_atm_pond, &
                                   Cdn_atm_rdg , hfreebd     , &
                                   hdraft      , hridge      , &
                                   distrdg     , hkeel       , &
                                   dkeel       , lfloe       , &
                                   dfloe       , ncat)
      endif

      do n = 1, ncat

         meltsn (n) = c0
         melttn (n) = c0
         meltbn (n) = c0
         congeln(n) = c0
         snoicen(n) = c0
         dsnown (n) = c0

         Trefn  = c0
         Qrefn  = c0
         Urefn  = c0
         lhcoef = c0
         shcoef = c0
         worka  = c0
         workb  = c0

         if (aicen_init(n) > puny) then

            if (calc_Tsfc .or. calc_strair) then 

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               call colpkg_atm_boundary( 'ice',                  &
                                        Tsfc(n),  potT,          &
                                        uatm,     vatm,          &
                                        wind,     zlvl,          &
                                        Qa,       rhoa,          &
                                        strairxn, strairyn,      &
                                        Trefn,    Qrefn,         &
                                        worka,    workb,         &
                                        lhcoef,   shcoef,        &
                                        Cdn_atm,                 &
                                        Cdn_atm_ratio_n,         &
                                        uvel,     vvel,          &
                                        Uref=Urefn)

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then
#ifndef CICE_IN_NEMO
               ! Set to data values (on T points)
               strairxn = strax
               strairyn = stray
#else
               ! NEMO wind stress is supplied on u grid, multipied 
               ! by ice concentration and set directly in evp, so
               ! strairxT/yT = 0. Zero u-components here for safety.
               strairxn = c0
               strairyn = c0
#endif
            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) call increment_age (dt, iage(n))
            if (tr_FY)   call update_FYarea (dt,               &
                                             lmask_n, lmask_s, &
                                             yday,    FY(n))

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set 
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
               ! in thickness_changes
 
               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(aicen      (n),                 &
                                flatn_f    (n), fsensn_f   (n), &
                                fcondtopn_f(n),                 &
                                fsurfn_f   (n),                 &
                                flatn      (n), fsensn     (n), &
                                fsurfn     (n),                 &
                                fcondtopn  (n), nu_diag)
            endif

            call thermo_vertical(nilyr,        nslyr,        &
                                 dt,           aicen    (n), &
                                 vicen    (n), vsnon    (n), &
                                 Tsfc     (n), zSin   (:,n), &
                                 zqin   (:,n), zqsn   (:,n), &
                                 apnd     (n), hpnd     (n), &
                                 iage     (n), tr_pond_topo, &
                                 flw,          potT,         &
                                 Qa,           rhoa,         &
                                 fsnow,        fpond,        &
                                 fbot,         Tbot,         &
                                 sss,                        &
                                 lhcoef,       shcoef,       &
                                 fswsfcn  (n), fswintn  (n), &
                                 Sswabsn(:,n), Iswabsn(:,n), &
                                 fsurfn   (n), fcondtopn(n), &
                                 fsensn   (n), flatn    (n), &
                                 flwoutn,      evapn,        &
                                 freshn,       fsaltn,       &
                                 fhocnn,                     &
                                 melttn   (n), meltsn   (n), &
                                 meltbn   (n),               &
                                 congeln  (n), snoicen  (n), &
                                 mlt_onset,    frz_onset,    &
                                 yday,         dsnown   (n), &
                                 l_stop,       nu_diag)
               
            if (l_stop) then
               stop_label = 'ice: Vertical thermo error'
               return
            endif
               
      !-----------------------------------------------------------------
      ! Total absorbed shortwave radiation
      !-----------------------------------------------------------------

            fswabsn = fswsfcn(n) + fswintn(n) + fswthrun(n)

      !-----------------------------------------------------------------
      ! Aerosol update
      !-----------------------------------------------------------------

            !!!AKT Column!!!if (tr_aero) then
            !!!AKT Column!!!   call update_aerosol (dt,                             &
            !!!AKT Column!!!                        melttn     (n), meltsn     (n), &
            !!!AKT Column!!!                        meltbn     (n), congeln    (n), &
            !!!AKT Column!!!                        snoicen    (n), fsnow,          &
            !!!AKT Column!!!                        aerosno(:,:,n), aeroice(:,:,n), &
            !!!AKT Column!!!                        aicen_init (n), vicen_init (n), &
            !!!AKT Column!!!                        vsnon_init (n),                 &
            !!!AKT Column!!!                        vicen      (n), vsnon      (n), &
            !!!AKT Column!!!                        aicen      (n),                 &
            !!!AKT Column!!!                        faero_atm  (:),  faero_ocn(:))
            !!!AKT Column!!!endif

         endif   ! aicen_init

      !-----------------------------------------------------------------
      ! Melt ponds
      ! If using tr_pond_cesm, the full calculation is performed here.
      ! If using tr_pond_topo, the rest of the calculation is done after
      ! the surface fluxes are merged, below.
      !-----------------------------------------------------------------

         !call ice_timer_start(timer_ponds)
         if (tr_pond) then
               
            if (tr_pond_cesm) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n) 
               call compute_ponds_cesm(dt,        hi_min,    &
                                       pndaspect, rfrac,     &
                                       melttn(n), meltsn(n), &
                                       frain,                &
                                       aicen (n), vicen (n), &
                                       vsnon (n), Tsfc  (n), &
                                       apnd  (n), hpnd  (n))
                  
            elseif (tr_pond_lvl) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
               call compute_ponds_lvl(dt,        nilyr,     &
                                      ktherm,               &
                                      hi_min,               &
                                      dpscale,   frzpnd,    &
                                      pndaspect, rfrac,     &
                                      melttn(n), meltsn(n), &
                                      frain,     Tair,      &
                                      fsurfn(n),            &
                                      dhsn  (n), ffracn(n), &
                                      aicen (n), vicen (n), &
                                      vsnon (n),            &
                                      zqin(:,n), zSin(:,n), &
                                      Tsfc  (n), alvl  (n), &
                                      apnd  (n), hpnd  (n), &
                                      ipnd  (n))
                  
            elseif (tr_pond_topo) then
               if (aicen_init(n) > puny) then
                     
                  ! collect liquid water in ponds
                  ! assume salt still runs off
                  rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
                  pond = rfrac/rhofresh * (melttn(n)*rhoi &
                       +                   meltsn(n)*rhos &
                       +                   frain *dt)

                  ! if pond does not exist, create new pond over full ice area
                  ! otherwise increase pond depth without changing pond area
                  if (apnd(n) < puny) then
                     hpnd(n) = c0
                     apnd(n) = c1
                  endif
                  hpnd(n) = (pond + hpnd(n)*apnd(n)) / apnd(n)
                  fpond = fpond + pond * aicen(n) ! m
               endif ! aicen_init
            endif

         endif ! tr_pond
         !call ice_timer_stop(timer_ponds)

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         if (aicen_init(n) > puny) &
            call merge_fluxes (aicen_init(n),            &
                               flw,        coszen,       & 
                               strairxn,   strairyn,     &
                               Cdn_atm_ratio_n,          &
                               fsurfn(n),  fcondtopn(n), &
                               fsensn(n),  flatn(n),     &
                               fswabsn,    flwoutn,      &
                               evapn,                    &
                               Trefn,      Qrefn,        &
                               freshn,     fsaltn,       &
                               fhocnn,     fswthrun(n),  &
                               strairxT,   strairyT,     &
                               Cdn_atm_ratio,            &
                               fsurf,      fcondtop,     &
                               fsens,      flat,         &
                               fswabs,     flwout,       &
                               evap,                     &
                               Tref,       Qref,         &
                               fresh,      fsalt,        &
                               fhocn,      fswthru,      &
                               melttn (n), meltsn(n),    &
                               meltbn (n), congeln(n),   &
                               snoicen(n),               &
                               meltt,      melts,        &
                               meltb,      congel,       &
                               snoice,                   &
                               Uref,       Urefn)

      enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
      !call ice_timer_start(timer_ponds)
      if (tr_pond_topo) then
         call compute_ponds_topo(dt,       ncat,      nilyr,     &
                                 ktherm,   heat_capacity,        &
                                 aice,     aicen (:),            &
                                 vice,     vicen (:),            &
                                 vsno,     vsnon (:),            &
                                 potT,     meltt,                &
                                 fsurf,    fpond,                &
                                 Tsfc (:), Tf,                   &
                                 zqin(:,:),zSin(:,:),            &
                                 apnd (:), hpnd  (:), ipnd(:),   &
                                 l_stop,   stop_label)
      endif
      !call ice_timer_stop(timer_ponds)

      end subroutine colpkg_step_therm1

!=======================================================================
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm2 (dt, ncat, n_aero,            &
                                     nilyr,        nslyr,         &
                                     hin_max,                     &
                                     aicen,                       &
                                     vicen,        vsnon,         &
                                     aicen_init,   vicen_init,    &
                                     trcrn,                       &
                                     aice0,        aice,          &
                                     trcr_depend,                 &
                                     trcr_base,    n_trcr_strata, &
                                     nt_strata,                   &
                                     Tf,           sss,           &
                                     salinz,                      &
                                     rside,        meltl,         &
                                     frzmlt,       frazil,        &
                                     frain,        fpond,         &
                                     fresh,        fsalt,         &
                                     fhocn,        update_ocn_f,  &
                                     faero_ocn,                   &
                                     first_ice,                   &
                                     flux_bio,     ocean_bio,     &
                                     l_stop,       stop_label,    &
                                     nu_diag,                     &
                                     frz_onset,    yday)

      use ice_constants_colpkg, only: puny
      use ice_itd, only: aggregate_area, reduce_area, cleanup_itd
      use ice_therm_itd, only: linear_itd, add_new_ice, lateral_melt
      use ice_colpkg_tracers, only: ntrcr, nbtrcr, tr_aero, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nslyr    , & ! number of snow layers
         n_aero   , & ! number of aerosol tracers
         nu_diag      ! diagnostic file unit number

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f     ! if true, update fresh water and salt fluxes

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind), intent(in) :: &
         dt       , & ! time step
         Tf       , & ! freezing temperature (C)
         sss      , & ! sea surface salinity (ppt)
         rside    , & ! fraction of ice that melts laterally
         frzmlt       ! freezing/melting potential (W/m^2)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         salinz   , & ! initial salinity profile
         ocean_bio    ! ocean concentration of biological tracer

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         frain    , & ! rainfall rate (kg/m^2 s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         meltl    , & ! lateral ice melt         (m/step-->cm/day)
         frazil       ! frazil ice growth        (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init,& ! initial concentration of ice
         vicen_init,& ! initial volume per unit area of ice          (m)
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        ! tracers
 
      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice      ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop         ! if true, abort model

      character (char_len), intent(out) :: stop_label

      real (kind=dbl_kind), intent(inout), optional :: &
         frz_onset    ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday         ! day of year

      l_stop = .false.
      
      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

      fresh  = fresh + frain * aice

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen(:), aice, aice0)

      if (kitd == 1) then

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

         if (aice > puny) then

            call linear_itd (ncat,     hin_max,        &
                             nilyr,    nslyr,          &
                             ntrcr,    trcr_depend(:), &
                             trcr_base (:,:),  & 
                             n_trcr_strata(:),&
                             nt_strata (:,:),          &
                             aicen_init(:),         &
                             vicen_init(:),         &
                             aicen     (:),         &
                             trcrn     (:,:), & 
                             vicen     (:),         &
                             vsnon     (:),         &
                             aice      ,         &
                             aice0     ,         &
                             fpond,       l_stop,      &
                             stop_label,  nu_diag)

            if (l_stop) return

         endif ! aice > puny

      endif  ! kitd = 1

!      call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

      ! identify ice-ocean cells

         call add_new_ice (ncat,          nilyr,        &
                           n_aero,        dt,           &
                           hin_max   (:), ktherm,       &
                           aicen     (:), trcrn (:,:),  &
                           vicen     (:),               &
                           aice0,         aice,         &
                           frzmlt,        frazil,       &
                           frz_onset,     yday,         &
                           update_ocn_f,                &
                           fresh,         fsalt,        &
                           Tf,            sss,          &
                           salinz    (:), phi_init,     &
                           dSin0_frazil,                &
                           nbtrcr,        flux_bio (:), &
                           ocean_bio (:), nu_diag,      &
                           l_stop,        stop_label)

         if (l_stop) return

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------

      call lateral_melt (dt,        ncat,          &
                         nilyr,     nslyr,         &
                         n_aero,    fpond,         &
                         fresh,     fsalt,         &
                         fhocn,     faero_ocn (:), &
                         rside,     meltl,         &
                         aicen (:), vicen     (:), &
                         vsnon (:), trcrn   (:,:))

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!echmod: test this
      if (ncat==1) &
          call reduce_area (hin_max   (0),                &
                            aicen     (1), vicen     (1), &
                            aicen_init(1), vicen_init(1))

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      call cleanup_itd (dt,                   ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max(:),       &
                        aicen(:),             trcrn(1:ntrcr,:), &
                        vicen(:),             vsnon(:),         &
                        aice0,                aice,             &
                        trcr_depend(:),       trcr_base(:,:),   &
                        n_trcr_strata(:),     nt_strata(:,:),   &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn(:),         tr_aero,          &
                        tr_pond_topo,         heat_capacity,    &
                        first_ice(:),                           &
                        flux_bio(1:nbtrcr),   n_aero,           &
                        l_stop,               stop_label,       &
                        nu_diag)

      end subroutine colpkg_step_therm2

!=======================================================================
!
! Scales radiation fields computed on the previous time step.
!
! authors: Elizabeth Hunke, LANL

      subroutine colpkg_prep_radiation (ncat, nilyr, nslyr,    &
                                        aice,        aicen,    &
                                        swvdr,       swvdf,    &
                                        swidr,       swidf,    &
                                        alvdr_ai,    alvdf_ai, &
                                        alidr_ai,    alidf_ai, &
                                        scale_factor,          &
                                        fswsfcn,     fswintn,  &
                                        fswthrun,    fswpenln, &
                                        Sswabsn,     Iswabsn)

      use ice_constants_colpkg, only: c0, c1, puny

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of ice thickness categories
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aice        , & ! ice area fraction
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         ! grid-box-mean albedos aggregated over categories (if calc_Tsfc)
         alvdr_ai    , & ! visible, direct   (fraction)
         alidr_ai    , & ! near-ir, direct   (fraction)
         alvdf_ai    , & ! visible, diffuse  (fraction)
         alidf_ai        ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen           ! ice area fraction in each category

      real (kind=dbl_kind), intent(inout) :: &
         scale_factor    ! shortwave scaling factor, ratio new:old

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun        ! SW through ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         fswpenln    , & ! visible SW entering ice layers (W m-2)
         Iswabsn     , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)

      ! local variables

      integer (kind=int_kind) :: &
         k           , & ! vertical index       
         n               ! thickness category index

      real (kind=dbl_kind) :: netsw 

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

         if (aice > c0 .and. scale_factor > puny) then
            netsw = swvdr*(c1 - alvdr_ai) &
                  + swvdf*(c1 - alvdf_ai) &
                  + swidr*(c1 - alidr_ai) &
                  + swidf*(c1 - alidf_ai)
            scale_factor = netsw / scale_factor
         else
            scale_factor = c1
         endif

         do n = 1, ncat

            if (aicen(n) > puny) then

      !-----------------------------------------------------------------
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

               fswsfcn(n)  = scale_factor*fswsfcn (n)
               fswintn(n)  = scale_factor*fswintn (n)
               fswthrun(n) = scale_factor*fswthrun(n)
               do k = 1,nilyr+1
                  fswpenln(k,n) = scale_factor*fswpenln(k,n)
               enddo       !k
               do k=1,nslyr
                  Sswabsn(k,n) = scale_factor*Sswabsn(k,n)
               enddo
               do k=1,nilyr
                  Iswabsn(k,n) = scale_factor*Iswabsn(k,n)
               enddo

            endif
         enddo                  ! ncat

      end subroutine colpkg_prep_radiation

!=======================================================================
!
! Computes radiation fields
!
! authors: William H. Lipscomb, LANL
!          David Bailey, NCAR
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_radiation (dt,       ncat,      & 
                                        nilyr,    nslyr,     &
                                        n_aero,              &
                                        aicen,    vicen,     &
                                        vsnon,    Tsfcn,     &
                                        alvln,    apndn,     &
                                        hpndn,    ipndn,     &
                                        aeron,               &
                                        TLAT,     TLON,      &
                                        calendar_type,       &
                                        days_per_year,       &
                                        nextsw_cday,         &
                                        yday,     sec,       &
                                        kaer_tab, waer_tab,  &
                                        gaer_tab,            &
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
                                        dhsn,     ffracn,    &
                                        nu_diag,  l_print_point, &
                                        initonly)

      use ice_constants_colpkg, only: c0, puny
      use ice_shortwave, only: run_dEdd, shortwave_ccsm3
      use ice_colpkg_tracers, only: tr_aero, tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat      , & ! number of ice thickness categories
         nilyr     , & ! number of ice layers
         nslyr     , & ! number of snow layers
         n_aero    , & ! number of aerosols
         nu_diag       ! diagnostic file unit number

      real (kind=dbl_kind), intent(in) :: &
         dt        , & ! time step (s)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         fsnow     , & ! snowfall rate (kg/m^2 s)
         TLAT, TLON    ! latitude and longitude (radian)

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=dbl_kind), intent(in) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real (kind=dbl_kind), intent(out) :: &
         coszen        ! cosine solar zenith angle, < 0 for sun below horizon 

      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))
   
      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen     , & ! ice area fraction in each category
         vicen     , & ! ice volume in each category (m)
         vsnon     , & ! snow volume in each category (m)
         Tsfcn     , & ! surface temperature (deg C)
         alvln     , & ! level-ice area fraction
         apndn     , & ! pond area fraction
         hpndn     , & ! pond depth (m)
         ipndn         ! pond refrozen lid thickness (m)

      real(kind=dbl_kind), dimension(:,:), intent(in) :: &
         aeron         ! aerosols (kg/m^3)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         alvdrn    , & ! visible, direct  albedo (fraction)
         alidrn    , & ! near-ir, direct   (fraction)
         alvdfn    , & ! visible, diffuse  (fraction)
         alidfn    , & ! near-ir, diffuse  (fraction)
         fswsfcn   , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun  , & ! SW through ice to ocean (W/m^2)
         dhsn      , & ! depth difference for snow on sea ice and pond ice
         ffracn    , & ! fraction of fsurfn used to melt ipond
                       ! albedo components for history
         albicen   , & ! bare ice 
         albsnon   , & ! snow 
         albpndn   , & ! pond 
         apeffn        ! effective pond area used for radiation calculation

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         fswpenln  , & ! visible SW entering ice layers (W m-2)
         Iswabsn   , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn       ! SW radiation absorbed in snow layers (W m-2)

      logical (kind=log_kind), intent(in) :: &
         l_print_point ! flag for printing diagnostics

      logical (kind=log_kind), optional :: &
           initonly    ! flag to indicate init only, default is false

      ! local variables

      integer (kind=int_kind) :: &
         n                  ! thickness category index

         ! Initialize
         do n = 1, ncat
            alvdrn  (n) = c0
            alidrn  (n) = c0
            alvdfn  (n) = c0
            alidfn  (n) = c0
            fswsfcn (n) = c0
            fswintn (n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         fswpenln (:,:) = c0
         Iswabsn  (:,:) = c0
         Sswabsn  (:,:) = c0
         coszen         = c0

         if (calc_Tsfc) then
         if (trim(shortwave) == 'dEdd') then ! delta Eddington
            
            call run_dEdd(dt,           tr_aero,        &
                          tr_pond_cesm,                 &
                          tr_pond_lvl,                  &
                          tr_pond_topo,                 &
                          ncat,         n_aero,         &
                          nilyr,        nslyr,          &
                          aicen(:),     vicen(:),       &
                          vsnon(:),     Tsfcn(:),       &
                          alvln(:),     apndn(:),       &
                          hpndn(:),     ipndn(:),       &
                          aeron(:,:),   kalg,           &
                          heat_capacity,                &
                          TLAT,         TLON,           &
                          calendar_type,days_per_year,  &
                          nextsw_cday,  yday,           &
                          sec,          R_ice,          &
                          R_pnd,        R_snw,          &
                          dT_mlt,       rsnw_mlt,       &
                          hs0,          hs1,            &
                          hp1,          pndaspect,      &
                          kaer_tab,     waer_tab,       &
                          gaer_tab,                     &
                          swvdr,        swvdf,          &
                          swidr,        swidf,          &
                          coszen,       fsnow,          &
                          alvdrn(:),    alvdfn(:),      &
                          alidrn(:),    alidfn(:),      &
                          fswsfcn(:),   fswintn(:),     &
                          fswthrun(:),  fswpenln(:,:),  &
                          Sswabsn(:,:), Iswabsn(:,:),   &
                          albicen(:),   albsnon(:),     &
                          albpndn(:),   apeffn(:),      &
                          dhsn(:),      ffracn(:),      &
                          nu_diag,      l_print_point,  &
                          initonly)
 
         else  ! .not. dEdd

            call shortwave_ccsm3(aicen(:),   vicen(:),   &
                                 vsnon(:),               &
                                 Tsfcn(:),               &
                                 swvdr,      swvdf,      &
                                 swidr,      swidf,      &
                                 heat_capacity,          &
                                 albedo_type,            &
                                 albicev,    albicei,    &
                                 albsnowv,   albsnowi,   &
                                 ahmax,                  &
                                 alvdrn(:),  alidrn(:),  &
                                 alvdfn(:),  alidfn(:),  &
                                 fswsfcn(:), fswintn(:), &
                                 fswthrun(:),            &
                                 fswpenln(:,:),          &
                                 Iswabsn(:,:),           &
                                 Sswabsn(:,:),           &
                                 albicen(:), albsnon(:), &
                                 coszen,     ncat)

         endif   ! shortwave

      else    ! .not. calc_Tsfc

      ! Calculate effective pond area for HadGEM

         if (tr_pond_topo) then
            do n = 1, ncat
               apeffn(n) = c0 
               if (aicen(n) > puny) then
               ! Lid effective if thicker than hp1
                 if (apndn(n)*aicen(n) > puny .and. ipndn(n) < hp1) then
                    apeffn(n) = apndn(n)
                 else
                    apeffn(n) = c0
                 endif
                 if (apndn(n) < puny) apeffn(n) = c0
               endif
            enddo  ! ncat
 
         endif ! tr_pond_topo

         ! Initialize for safety
         do n = 1, ncat
            alvdrn(n) = c0
            alidrn(n) = c0
            alvdfn(n) = c0
            alidfn(n) = c0
            fswsfcn(n) = c0
            fswintn(n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         Iswabsn(:,:) = c0
         Sswabsn(:,:) = c0

      endif    ! calc_Tsfc

      end subroutine colpkg_step_radiation

!=======================================================================
!
! Computes sea ice mechanical deformation
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_ridge (dt,           ndtd,          &
                                    nilyr,        nslyr,         &
                                    ncat,         hin_max,       &
                                    rdg_conv,     rdg_shear,     &
                                    aicen,                       &
                                    trcrn,                       &
                                    vicen,        vsnon,         &
                                    aice0,        trcr_depend,   &
                                    trcr_base,    n_trcr_strata, &
                                    nt_strata,                   &
                                    dardg1dt,     dardg2dt,      &
                                    dvirdgdt,     opening,       &
                                    fpond,                       &
                                    fresh,        fhocn,         &
                                    n_aero,                      &
                                    faero_ocn,                   &
                                    aparticn,     krdgn,         &
                                    aredistn,     vredistn,      &
                                    dardg1ndt,    dardg2ndt,     &
                                    dvirdgndt,                   &
                                    araftn,       vraftn,        &
                                    aice,         fsalt,         &
                                    first_ice,                   &
                                    flux_bio,     nu_diag,       &
                                    l_stop,       stop_label)

      use ice_mechred, only: ridge_ice
      use ice_itd, only: cleanup_itd
      use ice_colpkg_tracers, only: tr_pond_topo, tr_aero, tr_brine, ntrcr, nbtrcr

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ndtd  , & ! number of dynamics supercycles
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         nu_diag   ! diagnostic file unit number

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear, & ! shear term for ridging (1/s)
         dardg1dt , & ! rate of area loss by ridging ice (1/s)
         dardg2dt , & ! rate of area gain by new ridges (1/s)
         dvirdgdt , & ! rate of ice volume ridged (m/s)
         opening  , & ! rate of opening due to divergence/shear (1/s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn        ! net heat flux to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         dardg1ndt, & ! rate of area loss by ridging ice (1/s)
         dardg2ndt, & ! rate of area gain by new ridges (1/s)
         dvirdgndt, & ! rate of ice volume ridged (m/s)
         aparticn , & ! participation function
         krdgn    , & ! mean ridge thickness/thickness of ridging ice
         araftn   , & ! rafting ice area
         vraftn   , & ! rafting ice volume 
         aredistn , & ! redistribution function: fraction of new ridge area
         vredistn , & ! redistribution function: fraction of new ridge volume
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        ! tracers

      !logical (kind=log_kind), intent(in) :: &
         !tr_pond_topo,& ! if .true., use explicit topography-based ponds
         !tr_aero     ,& ! if .true., use aerosol tracers
         !tr_brine    !,& ! if .true., brine height differs from ice thickness
         !heat_capacity  ! if true, ice has nonzero heat capacity

      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice    ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop       ! if true, abort the model

      character (char_len), intent(out) :: &
         stop_label   ! diagnostic information for abort

      ! local variables

      real (kind=dbl_kind) :: &
         dtt      , & ! thermo time step
         atmp     , & ! temporary ice area
         atmp0        ! temporary open water area

      l_stop = .false.

      !-----------------------------------------------------------------
      ! Identify ice-ocean cells.
      ! Note:  We can not limit the loop here using aice>puny because
      !        aice has not yet been updated since the transport (and
      !        it may be out of whack, which the ridging helps fix).-ECH
      !-----------------------------------------------------------------
           
         call ridge_ice (dt,           ndtd,           &
                         ncat,         n_aero,         &
                         nilyr,        nslyr,          &
                         ntrcr,        hin_max(:),     &
                         rdg_conv,     rdg_shear,      &
                         aicen    (:),                 &
                         trcrn    (:,:),         &
                         vicen    (:), vsnon     (:),  &
                         aice0,                        &
                         trcr_depend  (:),       &
                         trcr_base    (:,:),       &
                         n_trcr_strata(:),       &
                         nt_strata    (:,:),       &
                         l_stop,                       &
                         stop_label,   nu_diag,        &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         dardg1dt,     dardg2dt,       &
                         dvirdgdt,     opening,        &
                         fpond,                        &
                         fresh,        fhocn,          &
                         tr_brine,     faero_ocn(:),   &
                         aparticn (:), krdgn     (:),  &
                         aredistn (:), vredistn  (:),  &
                         dardg1ndt(:), dardg2ndt (:),  &
                         dvirdgndt(:),                 &
                         araftn   (:), vraftn    (:))

         if (l_stop) return

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      dtt = dt * ndtd  ! for proper averaging over thermo timestep
      call cleanup_itd (dtt,                  ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max(:),       &
                        aicen(:),             trcrn(:,:),       &
                        vicen(:),             vsnon(:),         &
                        aice0,                aice,             &
                        trcr_depend(:),       trcr_base(:,:),   &
                        n_trcr_strata(:),     nt_strata(:,:),   &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn(:),         tr_aero,          &
                        tr_pond_topo,         heat_capacity,    &
                        first_ice(:),                           &
                        flux_bio(1:nbtrcr),   n_aero,           &
                        l_stop,               stop_label,       &
                        nu_diag)

      if (l_stop) then
         stop_label = 'ice: ITD cleanup error in colpkg_step_ridge'
      endif

      end subroutine colpkg_step_ridge

!=======================================================================

! Aggregate ice state variables over thickness categories.
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL

      subroutine colpkg_aggregate (ncat,               &
                                   aicen,    trcrn,    &
                                   vicen,    vsnon,    &
                                   aice,     trcr,     &
                                   vice,     vsno,     &
                                   aice0,              &
                                   ntrcr,              &
                                   trcr_depend,        &
                                   trcr_base,          & 
                                   n_trcr_strata,      &
                                   nt_strata)

      use ice_constants_colpkg, only: c0, c1
      use ice_colpkg_tracers, only: colpkg_compute_tracers

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn     ! ice tracers

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent(out) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:),  &
         intent(out) :: &
         trcr      ! ice tracers

      ! local variables

      integer (kind=int_kind) :: &
         n, it, itl, & ! loop indices
         ntr           ! tracer index

      real (kind=dbl_kind), dimension (:), allocatable :: &
         atrcr     ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind) :: &
         atrcrn    ! category value

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      aice0 = c1
      aice  = c0
      vice  = c0
      vsno  = c0

      allocate (atrcr(ntrcr))

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:) = c0

      do n = 1, ncat

            aice = aice + aicen(n)
            vice = vice + vicen(n)
            vsno = vsno + vsnon(n)

         do it = 1, ntrcr
            atrcrn = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                + trcr_base(it,2) * vicen(n) &
                                + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then  ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn = atrcrn * trcrn(ntr,n)
               enddo
            endif
            atrcr(it) = atrcr(it) + atrcrn
         enddo                  ! ntrcr
      enddo                     ! ncat

      ! Open water fraction
      aice0 = max (c1 - aice, c0)

      ! Tracers
      call colpkg_compute_tracers (ntrcr,     trcr_depend,   &
                                   atrcr(:),  aice,          &
                                   vice ,     vsno,          &
                                   trcr_base, n_trcr_strata, &
                                   nt_strata, trcr(:))

      deallocate (atrcr)

      end subroutine colpkg_aggregate

!=======================================================================

! Compute the strength of the ice pack, defined as the energy (J m-2)
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_ice_strength (ncat,               &
                                      aice,     vice,     &
                                      aice0,    aicen,    &
                                      vicen,    &
                                      strength)

      use ice_constants_colpkg, only: p333, c0, c1, c2, Cf, Cp, Pstar, Cstar, &
          rhoi, puny
      use ice_mechred, only: asum_ridging, ridge_itd

      integer (kind=int_kind), intent(in) :: & 
         ncat       ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         aice   , & ! concentration of ice
         vice   , & ! volume per unit area of ice  (m)
         aice0      ! concentration of open water

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen  , & ! concentration of ice
         vicen      ! volume per unit area of ice  (m)

      real (kind=dbl_kind), intent(inout) :: &
         strength   ! ice strength (N/m)

      ! local variables

      real (kind=dbl_kind) :: &
         asum   , & ! sum of ice and open water area
         aksum      ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic    ! participation function; fraction of ridging
                    ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin  , & ! minimum ridge thickness
         hrmax  , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp  , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg       ! mean ridge thickness/thickness of ridging ice

      integer (kind=int_kind) :: &
         n          ! thickness category index

      real (kind=dbl_kind) :: &
         hi     , & ! ice thickness (m)
         h2rdg  , & ! mean value of h^2 for new ridge
         dh2rdg     ! change in mean value of h^2 per unit area
                    ! consumed by ridging 

      if (kstrength == 1) then  ! Rothrock '75 formulation

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------

         call asum_ridging (ncat, aicen(:), aice0, asum)

         call ridge_itd (ncat,     aice0,      &
                         aicen(:), vicen(:),   &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         aksum,    apartic(:), &
                         hrmin(:), hrmax(:),   &
                         hrexp(:), krdg(:))

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler 1980 formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0)then
                  hi = vicen(n) / aicen(n)
                  h2rdg = p333 * (hrmax(n)**3 - hrmin(n)**3)  &
                               / (hrmax(n) - hrmin(n)) 
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif         ! aicen > puny
            enddo               ! n

         elseif (krdg_redist==1) then ! exponential formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0) then
                  hi = vicen(n) / aicen(n)
                  h2rdg =    hrmin(n)*hrmin(n) &
                        + c2*hrmin(n)*hrexp(n) &
                        + c2*hrexp(n)*hrexp(n)
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif
            enddo               ! n

         endif                  ! krdg_redist

         strength = Cf * Cp * strength / aksum
                       ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                       ! Cf accounts for frictional dissipation

      else                      ! kstrength /= 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------

         strength = Pstar*vice*exp(-Cstar*(c1-aice))

      endif                     ! kstrength

      end subroutine colpkg_ice_strength

!=======================================================================

      subroutine colpkg_atm_boundary(sfctype,                    &
                                     Tsf,         potT,          &
                                     uatm,        vatm,          &
                                     wind,        zlvl,          &
                                     Qa,          rhoa,          &
                                     strx,        stry,          &
                                     Tref,        Qref,          &
                                     delt,        delq,          &
                                     lhcoef,      shcoef,        &
                                     Cdn_atm,                    &
                                     Cdn_atm_ratio_n,            &
                                     uvel,        vvel,          &
                                     Uref)

      use ice_atmo, only: atmo_boundary_const, atmo_boundary_layer
      use ice_constants_colpkg, only: c0

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         Cdn_atm  , &    ! neutral drag coefficient
         Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind), &
         intent(inout) :: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(inout) :: &
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind), optional, intent(in) :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind), optional, intent(out) :: &
         Uref         ! reference height wind speed (m/s)

      real (kind=dbl_kind) :: &
         worku, workv, workr

      worku = c0
      workv = c0
      workr = c0
      if (present(uvel)) then
         worku = uvel
      endif
      if (present(uvel)) then
         worku = uvel
      endif

               if (trim(atmbndy) == 'constant') then
                  call atmo_boundary_const (sfctype,  calc_strair, &
                                            uatm,     vatm,     &
                                            wind,     rhoa,     &
                                            strx,     stry,     &
                                            Tsf,      potT,     &
                                            Qa,                 &
                                            delt,     delq,     &
                                            lhcoef,   shcoef,   &
                                            Cdn_atm)
               else ! default
                  call atmo_boundary_layer (sfctype,                 &
                                            calc_strair, formdrag,   &
                                            highfreq, natmiter,      &
                                            Tsf,      potT,          &
                                            uatm,     vatm,          &
                                            wind,     zlvl,          &
                                            Qa,       rhoa,          &
                                            strx,     stry,          &
                                            Tref,     Qref,          &
                                            delt,     delq,          &
                                            lhcoef,   shcoef,        &
                                            Cdn_atm,                 &
                                            Cdn_atm_ratio_n,         &
                                            worku,    workv,         &
                                            workr)
               endif ! atmbndy

      if (present(Uref)) then
         Uref = workr
      endif

      end subroutine colpkg_atm_boundary

!=======================================================================
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       ice_therm_vertical.

      subroutine colpkg_ocn_mixed_layer (alvdr_ocn, swvdr,      &
                                         alidr_ocn, swidr,      &
                                         alvdf_ocn, swvdf,      &
                                         alidf_ocn, swidf,      &
                                         sst,       flwout_ocn, &
                                         fsens_ocn, shcoef,     &
                                         flat_ocn,  lhcoef,     &
                                         evap_ocn,  flw,        &
                                         delt,      delq,       &
                                         aice,      fhocn,      &
                                         fswthru,   hmix,       &
                                         Tf,        qdp,        &
                                         frzmlt,    dt)

      use ice_constants_colpkg, only: c0, c1, c1000, &
          cp_ocn, Tffresh, stefan_boltzmann, Lvap, cprho

      real (kind=dbl_kind), intent(in) :: &
         alvdr_ocn , & ! visible, direct   (fraction)
         alidr_ocn , & ! near-ir, direct   (fraction)
         alvdf_ocn , & ! visible, diffuse  (fraction)
         alidf_ocn , & ! near-ir, diffuse  (fraction)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         flw       , & ! incoming longwave radiation (W/m^2)
         Tf        , & ! freezing temperature (C)
         hmix      , & ! mixed layer depth (m)
         delt      , & ! potential temperature difference   (K)
         delq      , & ! specific humidity difference   (kg/kg)
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef    , & ! transfer coefficient for latent heat
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fswthru   , & ! shortwave penetrating to ocean (W/m^2)
         aice      , & ! ice area fraction
         dt            ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         flwout_ocn, & ! outgoing longwave radiation (W/m^2)
         fsens_ocn , & ! sensible heat flux (W/m^2)
         flat_ocn  , & ! latent heat flux   (W/m^2)
         evap_ocn  , & ! evaporative water flux (kg/m^2/s)
         qdp       , & ! deep ocean heat flux (W/m^2), negative upward
         sst       , & ! sea surface temperature (C)
         frzmlt        ! freezing/melting potential (W/m^2)

      ! local variables

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      ! shortwave radiative flux
      swabs = (c1-alvdr_ocn) * swvdr + (c1-alidr_ocn) * swidr &
            + (c1-alvdf_ocn) * swvdf + (c1-alidf_ocn) * swidf 

      ! ocean surface temperature in Kelvin
      TsfK = sst + Tffresh

      ! longwave radiative flux
      flwout_ocn = -stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsens_ocn =  shcoef * delt
      flat_ocn  =  lhcoef * delq
      evap_ocn  = -flat_ocn / Lvap

      ! Compute sst change due to exchange with atm/ice above
      sst = sst + dt * ( &
            (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs) * (c1-aice) &
          + fhocn + fswthru)         &  ! these are *aice
          / (cprho*hmix)

      ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
      if (sst <= Tf .and. qdp > c0) qdp = c0

      ! computed T change due to exchange with deep layers:
      sst = sst - qdp*dt/(cprho*hmix)

      ! compute potential to freeze or melt ice
      frzmlt = (Tf-sst)*cprho*hmix/dt
      frzmlt = min(max(frzmlt,-frzmlt_max),frzmlt_max)

      ! if sst is below freezing, reset sst to Tf
      if (sst <= Tf) sst = Tf

      end subroutine colpkg_ocn_mixed_layer

!=======================================================================
! subroutine to set the column package internal parameters

      subroutine colpkg_init_parameters(&
           ktherm_in, &
           conduct_in, &
           fbot_xfer_type_in, &
           calc_Tsfc_in, &
           ustar_min_in, &
           a_rapid_mode_in, &
           Rac_rapid_mode_in, &
           aspect_rapid_mode_in, &
           dSdt_slow_mode_in, &
           phi_c_slow_mode_in, &
           phi_i_mushy_in, &
           shortwave_in, &
           albedo_type_in, &
           albicev_in, &
           albicei_in, &
           albsnowv_in, &
           albsnowi_in, &
           ahmax_in, &
           R_ice_in, &
           R_pnd_in, &
           R_snw_in, &
           dT_mlt_in, &
           rsnw_mlt_in, &
           kalg_in, &
           kstrength_in, &
           krdg_partic_in, &
           krdg_redist_in, &
           mu_rdg_in, &
           Cf_in, &
           atmbndy_in, &
           calc_strair_in, &
           formdrag_in, &
           highfreq_in, &
           natmiter_in, &
           oceanmixed_ice_in, &
           tfrz_option_in, &
           kitd_in, &
           kcatbound_in, &
           hs0_in, &
           frzpnd_in, &
           dpscale_in, &
           rfracmin_in, &
           rfracmax_in, &
           pndaspect_in, &
           hs1_in, &
           hp1_in)

        use ice_colpkg_shared, only: &
             ktherm, &
             conduct, &
             fbot_xfer_type, &
             calc_Tsfc, &
             ustar_min, &
             a_rapid_mode, &
             Rac_rapid_mode, &
             aspect_rapid_mode, &
             dSdt_slow_mode, &
             phi_c_slow_mode, &
             phi_i_mushy, &
             shortwave, &
             albedo_type, &
             albicev, &
             albicei, &
             albsnowv, &
             albsnowi, &
             ahmax, &
             R_ice, &
             R_pnd, &
             R_snw, &
             dT_mlt, &
             rsnw_mlt, &
             kalg, &
             kstrength, &
             krdg_partic, &
             krdg_redist, &
             mu_rdg, &
             Cf, &
             atmbndy, &
             calc_strair, &
             formdrag, &
             highfreq, &
             natmiter, &
             oceanmixed_ice, &
             tfrz_option, &
             kitd, &
             kcatbound, &
             hs0, &
             frzpnd, &
             dpscale, &
             rfracmin, &
             rfracmax, &
             pndaspect, &
             hs1, &
             hp1

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             ktherm_in          ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(in) :: &
             conduct_in, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(in) :: &
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(in) :: &
             ustar_min_in       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=dbl_kind), intent(in) :: &
             a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_in           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             shortwave_in, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
             albedo_type_in  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=dbl_kind), intent(in) :: &
             albicev_in  , & ! visible ice albedo for h > ahmax
             albicei_in  , & ! near-ir ice albedo for h > ahmax
             albsnowv_in , & ! cold snow albedo, visible
             albsnowi_in , & ! cold snow albedo, near IR
             ahmax_in        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=dbl_kind), intent(in) :: &
             R_ice_in    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
             R_pnd_in    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
             R_snw_in    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
             dT_mlt_in   , & ! change in temp for non-melt to melt snow grain 
                             ! radius change (C)
             rsnw_mlt_in , & ! maximum melting snow grain radius (10^-6 m)
             kalg_in         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: & ! defined in namelist 
             kstrength_in  , & ! 0 for simple Hibler (1979) formulation 
                               ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation 
                               ! 1 for exponential participation function 
             krdg_redist_in    ! 0 for Hibler (1980) formulation 
                               ! 1 for exponential redistribution function 
 
        real (kind=dbl_kind), intent(in) :: &  
             mu_rdg_in, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_in             ! ratio of ridging work to PE change in ridging (kstrength = 1)

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             atmbndy_in ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(in) :: &
             calc_strair_in, &  ! if true, calculate wind stress components
             formdrag_in,    &  ! if true, calculate form drag
             highfreq_in        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(in) :: &
             natmiter_in        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        logical (kind=log_kind), intent(in) :: &
             oceanmixed_ice_in           ! if true, use ocean mixed layer
        
        character(len=char_len), intent(in) :: &
             tfrz_option_in              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             kitd_in        , & ! type of itd conversions
                                !   0 = delta function
                                !   1 = linear remap
             kcatbound_in       !   0 = old category boundary formula
                                !   1 = new formula giving round numbers
                                !   2 = WMO standard
                                !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=dbl_kind), intent(in) :: &
             hs0_in             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(in) :: &
             frzpnd_in          ! pond refreezing parameterization
        
        real (kind=dbl_kind), intent(in) :: &
             dpscale_in, &      ! alter e-folding time scale for flushing 
             rfracmin_in, &     ! minimum retained fraction of meltwater
             rfracmax_in, &     ! maximum retained fraction of meltwater
             pndaspect_in, &    ! ratio of pond depth to pond fraction
             hs1_in             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=dbl_kind), intent(in) :: &
             hp1_in             ! critical parameter for pond ice thickness
        
        ktherm = ktherm_in
        conduct = conduct_in
        fbot_xfer_type = fbot_xfer_type_in
        calc_Tsfc = calc_Tsfc_in
        ustar_min = ustar_min_in
        a_rapid_mode = a_rapid_mode_in
        Rac_rapid_mode = Rac_rapid_mode_in
        aspect_rapid_mode = aspect_rapid_mode_in
        dSdt_slow_mode = dSdt_slow_mode_in
        phi_c_slow_mode = phi_c_slow_mode_in
        phi_i_mushy = phi_i_mushy_in
        shortwave = shortwave_in
        albedo_type = albedo_type_in
        albicev = albicev_in
        albicei = albicei_in
        albsnowv = albsnowv_in
        albsnowi = albsnowi_in
        ahmax = ahmax_in
        R_ice = R_ice_in
        R_pnd = R_pnd_in
        R_snw = R_snw_in
        dT_mlt = dT_mlt_in
        rsnw_mlt = rsnw_mlt_in
        kalg = kalg_in
        kstrength = kstrength_in
        krdg_partic = krdg_partic_in
        krdg_redist = krdg_redist_in
        mu_rdg = mu_rdg_in
        Cf = Cf_in
        atmbndy = atmbndy_in
        calc_strair = calc_strair_in
        formdrag = formdrag_in
        highfreq = highfreq_in
        natmiter = natmiter_in
        oceanmixed_ice = oceanmixed_ice_in
        tfrz_option = tfrz_option_in
        kitd = kitd_in
        kcatbound = kcatbound_in
        hs0 = hs0_in
        frzpnd = frzpnd_in
        dpscale = dpscale_in
        rfracmin = rfracmin_in
        rfracmax = rfracmax_in
        pndaspect = pndaspect_in
        hs1 = hs1_in
        hp1 = hp1_in

      end subroutine colpkg_init_parameters

!=======================================================================
! set tracer active flags

      subroutine colpkg_init_tracer_flags(&
           tr_iage_in      , & ! if .true., use age tracer
           tr_FY_in        , & ! if .true., use first-year area tracer
           tr_lvl_in       , & ! if .true., use level ice tracer
           tr_pond_in      , & ! if .true., use melt pond tracer
           tr_pond_cesm_in , & ! if .true., use cesm pond tracer
           tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
           tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
           tr_aero_in      , & ! if .true., use aerosol tracers
           tr_brine_in)        ! if .true., brine height differs from ice thickness

        use ice_colpkg_tracers, only: &
             tr_iage      , & ! if .true., use age tracer
             tr_FY        , & ! if .true., use first-year area tracer
             tr_lvl       , & ! if .true., use level ice tracer
             tr_pond      , & ! if .true., use melt pond tracer
             tr_pond_cesm , & ! if .true., use cesm pond tracer
             tr_pond_lvl  , & ! if .true., use level-ice pond tracer
             tr_pond_topo , & ! if .true., use explicit topography-based ponds
             tr_aero      , & ! if .true., use aerosol tracers
             tr_brine         ! if .true., brine height differs from ice thickness

        logical, intent(in) :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in         ! if .true., brine height differs from ice thickness

        tr_iage      = tr_iage_in
        tr_FY        = tr_FY_in
        tr_lvl       = tr_lvl_in
        tr_pond      = tr_pond_in
        tr_pond_cesm = tr_pond_cesm_in
        tr_pond_lvl  = tr_pond_lvl_in
        tr_pond_topo = tr_pond_topo_in
        tr_aero      = tr_aero_in
        tr_brine     = tr_brine_in

      end subroutine colpkg_init_tracer_flags

!=======================================================================

      subroutine colpkg_init_tracer_indices(&
           nt_Tsfc_in, & ! ice/snow temperature
           nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
           nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
           nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
           nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
           nt_iage_in, & ! volume-weighted ice age
           nt_FY_in, & ! area-weighted first-year ice area
           nt_alvl_in, & ! level ice area fraction
           nt_vlvl_in, & ! level ice volume fraction
           nt_apnd_in, & ! melt pond area fraction
           nt_hpnd_in, & ! melt pond depth
           nt_ipnd_in, & ! melt pond refrozen lid thickness
           nt_aero_in, & ! starting index for aerosols in ice
           nt_bgc_N_sk_in, & ! algae (skeletal layer)
           nt_bgc_C_sk_in, & ! 
           nt_bgc_chl_sk_in, & ! 
           nt_bgc_Nit_sk_in, & ! nutrients (skeletal layer) 
           nt_bgc_Am_sk_in, & ! 
           nt_bgc_Sil_sk_in, & !
           nt_bgc_DMSPp_sk_in, & ! trace gases (skeletal layer)
           nt_bgc_DMSPd_sk_in, & ! 
           nt_bgc_DMS_sk_in, & ! 
           nt_bgc_Nit_ml_in, & ! nutrients (ocean mixed layer) 
           nt_bgc_Am_ml_in, & ! 
           nt_bgc_Sil_ml_in, & !
           nt_bgc_DMSP_ml_in, & ! trace gases (ocean mixed layer)
           nt_bgc_DMS_ml_in)

        use ice_colpkg_tracers, only: &
             nt_Tsfc, & ! ice/snow temperature
             nt_qice, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno, & ! volume-weighted snow enthalpy (in layers)
             nt_sice, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage, & ! volume-weighted ice age
             nt_FY, & ! area-weighted first-year ice area
             nt_alvl, & ! level ice area fraction
             nt_vlvl, & ! level ice volume fraction
             nt_apnd, & ! melt pond area fraction
             nt_hpnd, & ! melt pond depth
             nt_ipnd, & ! melt pond refrozen lid thickness
             nt_aero, & ! starting index for aerosols in ice
             nt_bgc_N_sk, & ! algae (skeletal layer)
             nt_bgc_C_sk, & ! 
             nt_bgc_chl_sk, & ! 
             nt_bgc_Nit_sk, & ! nutrients (skeletal layer) 
             nt_bgc_Am_sk, & ! 
             nt_bgc_Sil_sk, & !
             nt_bgc_DMSPp_sk, & ! trace gases (skeletal layer)
             nt_bgc_DMSPd_sk, & ! 
             nt_bgc_DMS_sk, & ! 
             nt_bgc_Nit_ml, & ! nutrients (ocean mixed layer) 
             nt_bgc_Am_ml, & ! 
             nt_bgc_Sil_ml, & !
             nt_bgc_DMSP_ml, & ! trace gases (ocean mixed layer)
             nt_bgc_DMS_ml
        
        integer, intent(in) :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_aero_in, & ! starting index for aerosols in ice
             nt_bgc_N_sk_in, & ! algae (skeletal layer)
             nt_bgc_C_sk_in, & ! 
             nt_bgc_chl_sk_in, & ! 
             nt_bgc_Nit_sk_in, & ! nutrients (skeletal layer) 
             nt_bgc_Am_sk_in, & ! 
             nt_bgc_Sil_sk_in, & !
             nt_bgc_DMSPp_sk_in, & ! trace gases (skeletal layer)
             nt_bgc_DMSPd_sk_in, & ! 
             nt_bgc_DMS_sk_in, & ! 
             nt_bgc_Nit_ml_in, & ! nutrients (ocean mixed layer) 
             nt_bgc_Am_ml_in, & ! 
             nt_bgc_Sil_ml_in, & !
             nt_bgc_DMSP_ml_in, & ! trace gases (ocean mixed layer)
             nt_bgc_DMS_ml_in

        nt_Tsfc = nt_Tsfc_in
        nt_qice = nt_qice_in
        nt_qsno = nt_qsno_in
        nt_sice = nt_sice_in
        nt_fbri = nt_fbri_in
        nt_iage = nt_iage_in
        nt_FY = nt_FY_in
        nt_alvl = nt_alvl_in
        nt_vlvl = nt_vlvl_in
        nt_apnd = nt_apnd_in
        nt_hpnd = nt_hpnd_in
        nt_ipnd = nt_ipnd_in
        nt_aero = nt_aero_in
        nt_bgc_N_sk = nt_bgc_N_sk_in
        nt_bgc_C_sk = nt_bgc_C_sk_in
        nt_bgc_chl_sk = nt_bgc_chl_sk_in
        nt_bgc_Nit_sk = nt_bgc_Nit_sk_in
        nt_bgc_Am_sk = nt_bgc_Am_sk_in
        nt_bgc_Sil_sk = nt_bgc_Sil_sk_in
        nt_bgc_DMSPp_sk = nt_bgc_DMSPp_sk_in
        nt_bgc_DMSPd_sk = nt_bgc_DMSPd_sk_in
        nt_bgc_DMS_sk = nt_bgc_DMS_sk_in
        nt_bgc_Nit_ml = nt_bgc_Nit_ml_in
        nt_bgc_Am_ml = nt_bgc_Am_ml_in
        nt_bgc_Sil_ml = nt_bgc_Sil_ml_in
        nt_bgc_DMSP_ml = nt_bgc_DMSP_ml_in
        nt_bgc_DMS_ml = nt_bgc_DMS_ml_in

      end subroutine colpkg_init_tracer_indices

!=======================================================================
! set the number of column tracers

      subroutine colpkg_init_tracer_numbers(&
           ntrcr_in, &
           nbtrcr_in)

        use ice_colpkg_tracers, only: &
             ntrcr, &
             nbtrcr

        integer (kind=int_kind), intent(in) :: &
             ntrcr_in     ! number of tracers in use
        
        integer (kind=int_kind), intent(in) :: &
             nbtrcr_in    ! number of bgc tracers in use
        
        ntrcr = ntrcr_in
        nbtrcr = nbtrcr_in

      end subroutine colpkg_init_tracer_numbers

!=======================================================================

      end module ice_colpkg

!=======================================================================
