!  SVN:$Id: ice_flux_colpkg.F90 1175 2017-03-02 19:53:26Z akt $
!=======================================================================

! Flux manipulation routines for column package
!
! author Elizabeth C. Hunke, LANL
!
! 2014: Moved subroutines merge_fluxes, set_sfcflux from ice_flux.F90

      module ice_flux_colpkg

      use ice_kinds_mod
      use ice_constants_colpkg, only: c1, emissivity
      use ice_warnings, only: add_warning

      implicit none
      private
      public :: merge_fluxes, set_sfcflux               

!=======================================================================

      contains

!=======================================================================

! Aggregate flux information from all ice thickness categories
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_fluxes (aicen,                &    
                               flw,      coszn,      &
                               strairxn, strairyn,   &
                               Cdn_atm_ratio_n,      &
                               fsurfn,   fcondtopn,  &  
                               fsensn,   flatn,      & 
                               fswabsn,  flwoutn,    &
                               evapn,                &
                               Trefn,    Qrefn,      &
                               freshn,   fsaltn,     &
                               fhocnn,   fswthrun,   &
                               strairxT, strairyT,   &  
                               Cdn_atm_ratio,        &
                               fsurf,    fcondtop,   &
                               fsens,    flat,       & 
                               fswabs,   flwout,     &
                               evap,                 & 
                               Tref,     Qref,       &
                               fresh,    fsalt,      & 
                               fhocn,    fswthru,    &
                               melttn, meltsn, meltbn, congeln, snoicen, &
                               meltt,  melts,        &
                               meltb,                &
                               congel,  snoice,      &
                               Uref,     Urefn       )

      ! single category fluxes
      real (kind=dbl_kind), intent(in) :: &
          aicen   , & ! concentration of ice
          flw     , & ! downward longwave flux          (W/m**2)
          coszn   , & ! cosine of solar zenith angle 
          strairxn, & ! air/ice zonal  strss,           (N/m**2)
          strairyn, & ! air/ice merdnl strss,           (N/m**2)
          Cdn_atm_ratio_n, & ! ratio of total drag over neutral drag  
          fsurfn  , & ! net heat flux to top surface    (W/m**2)
          fcondtopn,& ! downward cond flux at top sfc   (W/m**2)
          fsensn  , & ! sensible heat flx               (W/m**2)
          flatn   , & ! latent   heat flx               (W/m**2)
          fswabsn , & ! shortwave absorbed heat flx     (W/m**2)
          flwoutn , & ! upwd lw emitted heat flx        (W/m**2)
          evapn   , & ! evaporation                     (kg/m2/s)
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
           
      real (kind=dbl_kind), optional, intent(in):: &
          Urefn       ! air speed reference level       (m/s)

      ! cumulative fluxes
      real (kind=dbl_kind), intent(inout) :: &
          strairxT, & ! air/ice zonal  strss,           (N/m**2)
          strairyT, & ! air/ice merdnl strss,           (N/m**2)
          Cdn_atm_ratio, & ! ratio of total drag over neutral drag
          fsurf   , & ! net heat flux to top surface    (W/m**2)
          fcondtop, & ! downward cond flux at top sfc   (W/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
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

      real (kind=dbl_kind), optional, intent(inout):: &
          Uref        ! air speed reference level       (m/s)

      !-----------------------------------------------------------------
      ! Merge fluxes
      ! NOTE: The albedo is aggregated only in cells where ice exists
      !       and (for the delta-Eddington scheme) where the sun is above
      !       the horizon. 
      !-----------------------------------------------------------------

      ! atmo fluxes

      strairxT   = strairxT + strairxn  * aicen
      strairyT   = strairyT + strairyn  * aicen
      Cdn_atm_ratio = Cdn_atm_ratio + &
                      Cdn_atm_ratio_n   * aicen
      fsurf      = fsurf    + fsurfn    * aicen
      fcondtop   = fcondtop + fcondtopn * aicen 
      fsens      = fsens    + fsensn    * aicen
      flat       = flat     + flatn     * aicen
      fswabs     = fswabs   + fswabsn   * aicen
      flwout     = flwout   &
           + (flwoutn - (c1-emissivity)*flw) * aicen
      evap       = evap     + evapn     * aicen
      Tref       = Tref     + Trefn     * aicen
      Qref       = Qref     + Qrefn     * aicen

      ! ocean fluxes
      if (present(Urefn) .and. present(Uref)) then
         Uref = Uref     + Urefn     * aicen
      endif

      fresh     = fresh     + freshn    * aicen
      fsalt     = fsalt     + fsaltn    * aicen
      fhocn     = fhocn     + fhocnn    * aicen
      fswthru   = fswthru   + fswthrun  * aicen

      ! ice/snow thickness

      meltt     = meltt     + melttn    * aicen
      meltb     = meltb     + meltbn    * aicen
      melts     = melts     + meltsn    * aicen
      congel    = congel    + congeln   * aicen
      snoice    = snoice    + snoicen   * aicen
      
      end subroutine merge_fluxes

!=======================================================================

! If model is not calculating surface temperature, set the surface
! flux values using values read in from forcing data or supplied via
! coupling (stored in ice_flux).
!
! If CICE is running in NEMO environment, convert fluxes from GBM values 
! to per unit ice area values. If model is not running in NEMO environment, 
! the forcing is supplied as per unit ice area values.
!
! authors Alison McLaren, Met Office

      subroutine set_sfcflux (aicen,               &
                              flatn_f,             &
                              fsensn_f,            &
                              fsurfn_f,            &
                              fcondtopn_f,         &
                              flatn,               &
                              fsensn,              &
                              fsurfn,              &
                              fcondtopn)

      ! ice state variables
      real (kind=dbl_kind), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         flatn_f     , & ! latent heat flux   (W/m^2) 
         fsensn_f    , & ! sensible heat flux (W/m^2) 
         fsurfn_f    , & ! net flux to top surface, not including fcondtopn
         fcondtopn_f     ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), intent(out):: &
         flatn       , & ! latent heat flux   (W/m^2) 
         fsensn      , & ! sensible heat flux   (W/m^2) 
         fsurfn      , & ! net flux to top surface, not including fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      ! local variables

      real (kind=dbl_kind)  :: &
         raicen          ! 1 or 1/aicen

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

      logical (kind=log_kind), parameter :: & 
         extreme_test=.false. ! test and write out extreme forcing data

      character(len=char_len_long) :: &
         warning ! warning message

      raicen        = c1

#ifdef CICE_IN_NEMO
!----------------------------------------------------------------------
! Convert fluxes from GBM values to per ice area values when 
! running in NEMO environment.  (When in standalone mode, fluxes
! are input as per ice area.)
!----------------------------------------------------------------------
      raicen        = c1 / aicen
#endif
      fsurfn   = fsurfn_f*raicen
      fcondtopn= fcondtopn_f*raicen
      flatn    = flatn_f*raicen
      fsensn   = fsensn_f*raicen

!----------------------------------------------------------------
! Flag up any extreme fluxes
!---------------------------------------------------------------

      if (extreme_test) then
         extreme_flag = .false.

         if (fcondtopn < -100.0_dbl_kind & 
              .or. fcondtopn > 20.0_dbl_kind) then
            extreme_flag = .true.
         endif
         
         if (fsurfn < -100.0_dbl_kind & 
              .or. fsurfn > 80.0_dbl_kind) then
            extreme_flag = .true.
         endif
         
         if (flatn < -20.0_dbl_kind & 
              .or. flatn > 20.0_dbl_kind) then
            extreme_flag = .true.
         endif

         if (extreme_flag) then

            if (fcondtopn < -100.0_dbl_kind & 
                 .or. fcondtopn > 20.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -100 > fcondtopn > 20'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,fcondtopn = ', & 
                    aicen,fcondtopn
               call add_warning(warning)
            endif
            
            if (fsurfn < -100.0_dbl_kind & 
                 .or. fsurfn > 80.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -100 > fsurfn > 40'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,fsurfn = ', & 
                    aicen,fsurfn
               call add_warning(warning)
            endif
            
            if (flatn < -20.0_dbl_kind & 
                 .or. flatn > 20.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -20 > flatn > 20'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,flatn = ', & 
                    aicen,flatn
               call add_warning(warning)
            endif
            
         endif  ! extreme_flag
      endif     ! extreme_test    
         
      end subroutine set_sfcflux 

!=======================================================================

      end module ice_flux_colpkg

!=======================================================================
