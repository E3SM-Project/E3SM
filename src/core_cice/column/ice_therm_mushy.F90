!  SVN:$Id: ice_therm_mushy.F90 1182 2017-03-16 19:29:26Z njeffery $
!=======================================================================

module ice_therm_mushy

  use ice_kinds_mod
  use ice_constants_colpkg, only: c0, c1, c2, c4, c8, c10, c1000, &
      p001, p01, p05, p1, p2, p5, pi, bignum, puny, ice_ref_salinity, &
      viscosity_dyn, rhow, rhoi, rhos, cp_ocn, cp_ice, Lfresh, gravit, &
      hs_min, ksno
  use ice_colpkg_shared, only: a_rapid_mode, Rac_rapid_mode, &
      aspect_rapid_mode, dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy

  use ice_therm_shared, only: ferrmax
  use ice_warnings, only: add_warning

  implicit none

  private
  public :: &
       temperature_changes_salinity, &
       permeability

  real(kind=dbl_kind), parameter :: &
       dTemp_errmax = 5.0e-4_dbl_kind ! max allowed change in temperature 
                                      ! between iterations

!=======================================================================

contains

!=======================================================================

  subroutine temperature_changes_salinity(dt,                 &
                                          nilyr,    nslyr,    &
                                          rhoa,     flw,      &
                                          potT,     Qa,       &
                                          shcoef,   lhcoef,   &
                                          fswsfc,   fswint,   &
                                          Sswabs,   Iswabs,   &
                                          hilyr,    hslyr,    &
                                          apond,    hpond,    &
                                          zqin,     zTin,     &
                                          zqsn,     zTsn,     &
                                          zSin,               &
                                          Tsf,      Tbot,     &
                                          sss,                &
                                          fsensn,   flatn,    &
                                          flwoutn,  fsurfn,   &
                                          fcondtop, fcondbot, &
                                          fadvheat, snoice,   &
                                          einit_old,          &
                                          lstop,    stop_label)

    ! solve the enthalpy and bulk salinity of the ice for a single column

    use ice_mushy_physics, only: &
         enthalpy_brine, &
         temperature_mush, &
         liquid_fraction, &
         temperature_snow, &
         temperature_mush_liquid_fraction, &
         liquidus_brine_salinity_mush, &
         conductivity_mush_array, &
         conductivity_snow_array

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers
    
    real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)
    
    real (kind=dbl_kind), intent(in) :: &
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         Tbot        , & ! ice bottom surfce temperature (deg C)
         sss             ! sea surface salinity (PSU)
         
    real (kind=dbl_kind), intent(inout) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint          ! SW absorbed in ice interior below surface (W m-2)
    
    real (kind=dbl_kind), intent(inout) :: &
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         apond       , & ! melt pond area fraction
         hpond           ! melt pond depth (m)

    real (kind=dbl_kind), intent(in) :: &
         einit_old       ! initial energy of melting (J m-2)
    
    real (kind=dbl_kind), dimension (:), intent(inout) :: &
         Sswabs      , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabs          ! SW radiation absorbed in ice layers (W m-2)
    
    real (kind=dbl_kind), intent(inout):: &
         fsurfn      , & ! net flux to top surface, excluding fcondtopn
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn         ! upward LW at surface (W m-2)
    
    real (kind=dbl_kind), intent(out):: &
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         fadvheat    , & ! flow of heat to ocean due to advection (W m-2)
         snoice          ! snow ice formation

    real (kind=dbl_kind), intent(inout):: &
         Tsf             ! ice/snow surface temperature (C)

    real (kind=dbl_kind), dimension (:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3)
         zTin        , & ! internal ice layer temperatures
         zSin        , & ! internal ice layer salinities
         zqsn        , & ! snow layer enthalpy (J m-3)
         zTsn            ! internal snow layer temperatures
    
    logical (kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag 

    character (len=*), intent(out) :: &
         stop_label   ! abort error message

    ! local variables
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         zqin0       , & ! ice layer enthalpy (J m-3) at start of timestep
         zTin0       , & ! internal ice layer temperatures (C) at start of timestep
         zSin0       , & ! internal ice layer salinities (ppt) at start of timestep
         phi         , & ! liquid fraction
         km          , & ! ice conductivity (W m-1 K-1)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(1:nilyr+1) :: &
         Sbr         , & ! brine salinity (ppt)
         qbr             ! brine enthalpy (J m-3)

    real(kind=dbl_kind), dimension(0:nilyr) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         zqsn0       , & ! snow layer enthalpy (J m-3) at start of timestep
         zTsn0       , & ! internal snow layer temperatures (C) at start of timestep
         ks              ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf0        , & ! ice/snow surface temperature (C) at start of timestep
         hin         , & ! ice thickness (m)
         hsn         , & ! snow thickness (m)
         hslyr_min   , & ! minimum snow layer thickness (m)
         w           , & ! vertical flushing Darcy velocity (m/s)
         qocn        , & ! ocean brine enthalpy (J m-3)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         Spond           ! melt pond salinity (ppt)

    integer(kind=int_kind) :: &
         k               ! ice/snow layer index

    logical(kind=log_kind) :: &
         lsnow           ! snow presence: T: has snow, F: no snow

    character(len=char_len_long) :: &
         warning         ! warning message
    
    lstop = .false.
    fadvheat   = c0
    snoice     = c0

    Tsf0  = Tsf
    zqsn0 = zqsn
    zqin0 = zqin
    zSin0 = zSin
    zTsn0 = zTsn
    zTin0 = zTin

    Spond = c0
    qpond = enthalpy_brine(c0)

    hslyr_min = hs_min / real(nslyr, dbl_kind)

    lsnow = (hslyr > hslyr_min)

    hin = hilyr * real(nilyr,dbl_kind)

    qocn = enthalpy_brine(Tbot)

    if (lsnow) then
       hsn = hslyr * real(nslyr,dbl_kind)
    else
       hsn = c0
    endif

    do k = 1, nilyr
       phi(k) = liquid_fraction(temperature_mush(zqin(k),zSin(k)),zSin(k))
    enddo ! k

    ! calculate vertical bulk darcy flow
    call flushing_velocity(zTin,   zSin,  &
                           phi,    nilyr, &
                           hin,    hsn,   &
                           hilyr,         &
                           hpond,  apond, & 
                           dt,     w)

    ! calculate quantities related to drainage
    call explicit_flow_velocities(nilyr,  zSin,   &
                                  zTin,   Tsf,    &
                                  Tbot,   q,      &
                                  dSdt,   Sbr,    &
                                  qbr,    dt,     &
                                  sss,    qocn,   &
                                  hilyr,  hin)
    
    ! calculate the conductivities
    call conductivity_mush_array(nilyr, zqin0, zSin0, km)

    if (lsnow) then
       ! case with snow

       ! calculate the snow conductivities
       call conductivity_snow_array(ks)

       ! run the two stage solver
       call two_stage_solver_snow(nilyr,       nslyr,      &
                                  Tsf,         Tsf0,       &
                                  zqsn,        zqsn0,      &
                                  zqin,        zqin0,      &
                                  zSin,        zSin0,      &
                                  zTsn,        zTsn0,      &
                                  zTin,        zTin0,      &
                                  phi,         Tbot,       &
                                  km,          ks,         &
                                  q,           dSdt,       &
                                  w,           dt,         &
                                  fswint,      fswsfc,     &
                                  rhoa,        flw,        &
                                  potT,        Qa,         &
                                  shcoef,      lhcoef,     &
                                  Iswabs,      Sswabs,     &
                                  qpond,       qocn,       &
                                  Spond,       sss,        &
                                  hilyr,       hslyr,      &
                                  fcondtop,    fcondbot,   &
                                  fadvheat,                &
                                  flwoutn,     fsensn,     &
                                  flatn,       fsurfn,     &
                                  lstop,       stop_label)

       if (lstop) then
          write(warning,*) "temperature_changes_salinity: Picard solver non-convergence (snow)"
          call add_warning(warning)
          return
       endif

       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nslyr
          zTsn(k) = temperature_snow(zqsn(k))
       enddo ! k

       do k = 1, nilyr
          zTin(k) = temperature_mush_liquid_fraction(zqin(k), phi(k))
          Sbr(k)  = liquidus_brine_salinity_mush(zTin(k)) 
          qbr(k)  = enthalpy_brine(zTin(k))
       enddo ! k
          
    else
       ! case without snow

       ! run the two stage solver
       call two_stage_solver_nosnow(nilyr,       nslyr,      &
                                    Tsf,         Tsf0,       &
                                    zqsn,        zqsn0,      &
                                    zqin,        zqin0,      &
                                    zSin,        zSin0,      &
                                    zTsn,        zTsn0,      &
                                    zTin,        zTin0,      &
                                    phi,         Tbot,       &
                                    km,          ks,         &
                                    q,           dSdt,       &
                                    w,           dt,         &
                                    fswint,      fswsfc,     &
                                    rhoa,        flw,        &
                                    potT,        Qa,         &
                                    shcoef,      lhcoef,     &
                                    Iswabs,      Sswabs,     &
                                    qpond,       qocn,       &
                                    Spond,       sss,        &
                                    hilyr,       hslyr,      &
                                    fcondtop,    fcondbot,   &
                                    fadvheat,                &
                                    flwoutn,     fsensn,     &
                                    flatn,       fsurfn,     &
                                    lstop,       stop_label)

       if (lstop) then
          write(warning,*) "temperature_changes_salinity: Picard solver non-convergence (no snow)"
          call add_warning(warning)
          return
       endif
          
       ! given the updated enthalpy and bulk salinity calculate other quantities
       do k = 1, nilyr
          zTin(k) = temperature_mush_liquid_fraction(zqin(k), phi(k))
          Sbr(k)  = liquidus_brine_salinity_mush(zTin(k)) 
          qbr(k)  = enthalpy_brine(zTin(k))
       enddo ! k
       
    endif

    if (lstop) then
       return
    end if

    ! drain ponds from flushing 
    call flush_pond(w, hin, hpond, apond, dt)

    ! flood snow ice
    call flood_ice(hsn,        hin,      &
                   nslyr,      nilyr,    & 
                   hslyr,      hilyr,    & 
                   zqsn,       zqin,     &
                   phi,        dt,       &
                   zSin,       Sbr,      &
                   sss,        qocn,     &
                   snoice,     fadvheat)

  end subroutine temperature_changes_salinity

!=======================================================================

  subroutine two_stage_solver_snow(nilyr,       nslyr,      &
                                   Tsf,         Tsf0,       &
                                   zqsn,        zqsn0,      &
                                   zqin,        zqin0,      &
                                   zSin,        zSin0,      &
                                   zTsn,        zTsn0,      &
                                   zTin,        zTin0,      &
                                   phi,         Tbot,       &
                                   km,          ks,         &
                                   q,           dSdt,       &
                                   w,           dt,         &
                                   fswint,      fswsfc,     &
                                   rhoa,        flw,        &
                                   potT,        Qa,         &
                                   shcoef,      lhcoef,     &
                                   Iswabs,      Sswabs,     &
                                   qpond,       qocn,       &
                                   Spond,       sss,        &
                                   hilyr,       hslyr,      &
                                   fcondtop,    fcondbot,   &
                                   fadvheat,                &
                                   flwoutn,     fsensn,     &
                                   flatn,       fsurfn,     &
                                   lstop,       stop_label)

    ! solve the vertical temperature and salt change for case with snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    integer (kind=int_kind), intent(in) :: &
         nilyr      , &  ! number of ice layers
         nslyr           ! number of snow layers

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! snow surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fadvheat        ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf0            ! snow surface temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zTsn            ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqsn0       , & ! snow layer enthalpy (J m-3) at beginning of timestep
         zTsn0       , & ! snow layer temperature (C) at beginning of timestep
         ks          , & ! snow conductivity (W m-1 K-1)
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3) 
         zSin        , & ! ice layer bulk salinity (ppt)
         zTin        , & ! ice layer temperature (C)
         phi             ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0       , & ! ice layer enthalpy (J m-3) at beginning of timestep
         zSin0       , & ! ice layer bulk salinity (ppt) at beginning of timestep
         zTin0       , & ! ice layer temperature (C) at beginning of timestep
         km          , & ! ice conductivity (W m-1 K-1)
         Iswabs      , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         dt          , & ! time step (s)
         Tbot        , & ! ice bottom surfce temperature (deg C)
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         w           , & ! vertical flushing Darcy velocity (m/s)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         qocn        , & ! ocean brine enthalpy (J m-3)
         Spond       , & ! melt pond salinity (ppt)
         sss             ! sea surface salinity (PSU)

    logical(kind=log_kind), intent(inout) :: &
         lstop           ! solver failure flag

    character(len=*), intent(out) :: &
         stop_label      ! fatal error message

    real(kind=dbl_kind) :: &
         fcondtop1   , & ! first stage downward cond flux at top surface (W m-2)
         fsurfn1     , & ! first stage net flux to top surface, excluding fcondtop
         Tsf1            ! first stage ice surface temperature (C)


    ! determine if surface is initially cold or melting
    if (Tsf < c0) then

       ! initially cold

       ! solve the system for cold and snow
       call picard_solver(nilyr,   nslyr,     &
                          .true.,  .true.,    &
                          Tsf,      zqsn,     &
                          zqin,     zSin,     &
                          zTin,     zTsn,     &
                          phi,      dt,       &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,                  &
                          lstop,    stop_label)

       ! halt if solver failed
       if (lstop) return

       ! check if solution is consistent - surface should still be cold
       if (Tsf < dTemp_errmax) then

          ! solution is consistent - have solution so finish
          return

       else
          
          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting
          Tsf1 = Tsf

          ! reset the solution to initial values
          Tsf  = c0
          zqsn = zqsn0
          zqin = zqin0
          zSin = zSin0

          ! solve the system for melting and snow             
          call picard_solver(nilyr,    nslyr,    &
                             .true.,   .false.,  &
                             Tsf,      zqsn,     &
                             zqin,     zSin,     &
                             zTin,     zTsn,     &
                             phi,      dt,       &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,                  &
                             lstop,    stop_label)

          ! halt if solver failed
          if (lstop) return

          ! check if solution is consistent 
          ! surface conductive heat flux should be less than 
          ! incoming surface heat flux
          if (fcondtop - fsurfn < ferrmax) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             call two_stage_inconsistency(1, Tsf1, c0, fcondtop, fsurfn)
             lstop = .true.
             stop_label = "two_stage_solver_snow: two_stage_inconsistency: cold"
             return

          endif ! surface flux consistency

       endif ! surface temperature consistency

    else

       ! initially melting
       Tsf = c0

       ! solve the system for melting and snow
       call picard_solver(nilyr,    nslyr,    &
                          .true.,   .false.,  &
                          Tsf,      zqsn,     &
                          zqin,     zSin,     &
                          zTin,     zTsn,     &
                          phi,      dt,       &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,                  &
                          lstop,    stop_label)

       ! halt if solver failed
       if (lstop) return
       
       ! check if solution is consistent 
       ! surface conductive heat flux should be less than 
       ! incoming surface heat flux
       if (fcondtop - fsurfn < ferrmax) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold
          fcondtop1 = fcondtop
          fsurfn1   = fsurfn

          ! reset the solution to initial values
          Tsf  = Tsf0
          zqsn = zqsn0
          zqin = zqin0
          zSin = zSin0

          ! solve the system for cold and snow      
          call picard_solver(nilyr,    nslyr,    &
                             .true.,   .true.,   &       
                             Tsf,      zqsn,     &
                             zqin,     zSin,     &
                             zTin,     zTsn,     &
                             phi,      dt,       &
                             hilyr,    hslyr,    & 
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,                  &
                             lstop,    stop_label)

          ! halt if solver failed
          if (lstop) return

          ! check if solution is consistent - surface should be cold
          if (Tsf < dTemp_errmax) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             ! failed to find a solution so need to refine solutions until consistency found
             call two_stage_inconsistency(2, Tsf, c0, fcondtop1, fsurfn1)
             lstop = .true.
             stop_label = "two_stage_solver_snow: two_stage_inconsistency: melting"
             return
 
          endif ! surface temperature consistency
          
       endif ! surface flux consistency

    endif

  end subroutine two_stage_solver_snow

!=======================================================================

  subroutine two_stage_solver_nosnow(nilyr,       nslyr,      &
                                     Tsf,         Tsf0,       &
                                     zqsn,        zqsn0,      &
                                     zqin,        zqin0,      &
                                     zSin,        zSin0,      &
                                     zTsn,        zTsn0,      &
                                     zTin,        zTin0,      &
                                     phi,         Tbot,       &
                                     km,          ks,         &
                                     q,           dSdt,       &
                                     w,           dt,         &
                                     fswint,      fswsfc,     &
                                     rhoa,        flw,        &
                                     potT,        Qa,         &
                                     shcoef,      lhcoef,     &
                                     Iswabs,      Sswabs,     &
                                     qpond,       qocn,       &
                                     Spond,       sss,        &
                                     hilyr,       hslyr,      &
                                     fcondtop,    fcondbot,   &
                                     fadvheat,                &
                                     flwoutn,     fsensn,     &
                                     flatn,       fsurfn,     &
                                     lstop,       stop_label)
    
    ! solve the vertical temperature and salt change for case with no snow
    ! 1) determine what type of surface condition existed previously - cold or melting
    ! 2) solve the system assuming this condition persists
    ! 3) check the consistency of the surface condition of the solution
    ! 4) If the surface condition is inconsistent resolve for the other surface condition
    ! 5) If neither solution is consistent the resolve the inconsistency

    use ice_mushy_physics, only: &
         liquidus_temperature_mush

    integer (kind=int_kind), intent(in) :: &
         nilyr      , &  ! number of ice layers
         nslyr           ! number of snow layers

    real(kind=dbl_kind), intent(inout) :: &
         Tsf             ! ice surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop    , & ! downward cond flux at top surface (W m-2)
         fcondbot    , & ! downward cond flux at bottom surface (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fadvheat        ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf0            ! ice surface temperature (C) at beginning of timestep

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zTsn            ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqsn0       , & ! snow layer enthalpy (J m-3) at beginning of timestep
         zTsn0       , & ! snow layer temperature (C) at beginning of timestep
         ks          , & ! snow conductivity (W m-1 K-1)
         Sswabs          ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqin        , & ! ice layer enthalpy (J m-3) 
         zSin        , & ! ice layer bulk salinity (ppt)
         zTin        , & ! ice layer temperature (C)
         phi             ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0       , & ! ice layer enthalpy (J m-3) at beginning of timestep
         zSin0       , & ! ice layer bulk salinity (ppt) at beginning of timestep
         zTin0       , & ! ice layer temperature (C) at beginning of timestep
         km          , & ! ice conductivity (W m-1 K-1)
         Iswabs      , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt            ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q               ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         dt          , & ! time step (s)
         hilyr       , & ! ice layer thickness (m)
         hslyr       , & ! snow layer thickness (m)
         Tbot        , & ! ice bottom surfce temperature (deg C)
         fswint      , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         potT        , & ! air potential temperature (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         w           , & ! vertical flushing Darcy velocity (m/s)
         qpond       , & ! melt pond brine enthalpy (J m-3)
         qocn        , & ! ocean brine enthalpy (J m-3)
         Spond       , & ! melt pond salinity (ppt)
         sss             ! sea surface salinity (PSU)

    logical, intent(inout) :: &
         lstop           ! solver failure flag

    character(len=*), intent(out) :: &
         stop_label      ! fatal error message

    real(kind=dbl_kind) :: &
         Tmlt        , & ! upper ice layer melting temperature (C)
         fcondtop1   , & ! first stage downward cond flux at top surface (W m-2)
         fsurfn1     , & ! first stage net flux to top surface, excluding fcondtop
         Tsf1            ! first stage ice surface temperature (C)

    ! initial surface melting temperature
    Tmlt = liquidus_temperature_mush(zSin0(1))

    ! determine if surface is initially cold or melting
    if (Tsf < Tmlt) then

       ! initially cold

       ! solve the system for cold and no snow          
       call picard_solver(nilyr,    nslyr,    &
                          .false.,  .true.,   &
                          Tsf,      zqsn,     &
                          zqin,     zSin,     &
                          zTin,     zTsn,     &
                          phi,      dt,       &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,                  &
                          lstop,    stop_label)

       ! halt if solver failed
       if (lstop) return

       ! check if solution is consistent - surface should still be cold
       if (Tsf < Tmlt + dTemp_errmax) then

          ! solution is consistent - have solution so finish
          return

       else
          ! solution is inconsistent - surface is warmer than melting 
          ! resolve assuming surface is melting
          Tsf1 = Tsf

          ! reset the solution to initial values
          Tsf  = liquidus_temperature_mush(zSin0(1))
          zqin = zqin0
          zSin = zSin0

          ! solve the system for melt and no snow
          call picard_solver(nilyr,    nslyr,    &
                             .false.,  .false.,  &
                             Tsf,      zqsn,     &
                             zqin,     zSin,     &
                             zTin,     zTsn,     &
                             phi,      dt,       &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      & 
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,                  &
                             lstop,    stop_label)

          ! halt if solver failed
          if (lstop) return

          ! check if solution is consistent 
          ! surface conductive heat flux should be less than 
          ! incoming surface heat flux
          if (fcondtop - fsurfn < ferrmax) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             call two_stage_inconsistency(3, Tsf1, Tmlt, fcondtop, fsurfn)
             lstop = .true.
             stop_label = "two_stage_solver_nosnow: two_stage_inconsistency: cold"
             return

          endif

       endif

    else
       ! initially melting

       ! solve the system for melt and no snow
       Tsf = Tmlt       

       call picard_solver(nilyr,    nslyr,    &
                          .false.,  .false.,  &
                          Tsf,      zqsn,     &
                          zqin,     zSin,     &
                          zTin,     zTsn,     &
                          phi,      dt,       &
                          hilyr,    hslyr,    &
                          km,       ks,       &
                          Iswabs,   Sswabs,   &
                          Tbot,               &
                          fswint,   fswsfc,   &
                          rhoa,     flw,      &
                          potT,     Qa,       &
                          shcoef,   lhcoef,   &
                          fcondtop, fcondbot, &
                          fadvheat,           &
                          flwoutn,  fsensn,   &
                          flatn,    fsurfn,   &
                          qpond,    qocn,     &
                          Spond,    sss,      &
                          q,        dSdt,     &
                          w,                  &
                          lstop,    stop_label)

       ! halt if solver failed
       if (lstop) return

       ! check if solution is consistent 
       ! surface conductive heat flux should be less than 
       ! incoming surface heat flux
       if (fcondtop - fsurfn < ferrmax) then

          ! solution is consistent - have solution so finish
          return

       else

          ! solution is inconsistent - resolve assuming other surface condition
          ! assume surface is cold
          fcondtop1 = fcondtop
          fsurfn1   = fsurfn

          ! reset the solution to initial values
          Tsf  = Tsf0
          zqin = zqin0
          zSin = zSin0

          ! solve the system for cold and no snow
          call picard_solver(nilyr,    nslyr,    &
                             .false.,  .true.,   &
                             Tsf,      zqsn,     &
                             zqin,     zSin,     &
                             zTin,     zTsn,     &
                             phi,      dt,       &
                             hilyr,    hslyr,    &
                             km,       ks,       &
                             Iswabs,   Sswabs,   &
                             Tbot,               &
                             fswint,   fswsfc,   &
                             rhoa,     flw,      &
                             potT,     Qa,       &
                             shcoef,   lhcoef,   &
                             fcondtop, fcondbot, &
                             fadvheat,           &
                             flwoutn,  fsensn,   &
                             flatn,    fsurfn,   &
                             qpond,    qocn,     &
                             Spond,    sss,      &
                             q,        dSdt,     &
                             w,                  &
                             lstop,    stop_label)

          ! halt if solver failed
          if (lstop) return

          ! check if solution is consistent - surface should be cold
          if (Tsf < Tmlt + dTemp_errmax) then

             ! solution is consistent - have solution so finish
             return

          else

             ! solution is inconsistent
             call two_stage_inconsistency(4, Tsf, Tmlt, fcondtop1, fsurfn1)
             lstop = .true.
             stop_label = "two_stage_solver_nosnow: two_stage_inconsistency: melting"
             return

          endif
          
       endif

    endif

  end subroutine two_stage_solver_nosnow

!=======================================================================

  subroutine two_stage_inconsistency(type, Tsf, Tmlt, fcondtop, fsurfn)

    integer (kind=int_kind), intent(in) :: &
         type

    real(kind=dbl_kind), intent(in) :: &
         Tsf, &
         Tmlt, &
         fcondtop, &
         fsurfn

    character(len=char_len_long) :: &
         warning ! warning message
    
    write(warning,*) "ice_therm_mushy: two stage inconsistency"
    call add_warning(warning)
    write(warning,*) "type:", type
    call add_warning(warning)

    if (type == 1) then

       write(warning,*) "First stage  : Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax"
       call add_warning(warning)
       write(warning,*) "             :", Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax
       call add_warning(warning)
       write(warning,*) "Second stage : fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax"
       call add_warning(warning)
       write(warning,*) "             :", fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax
       call add_warning(warning)

    else if (type == 2) then

       write(warning,*) "First stage  : Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax"
       call add_warning(warning)
       write(warning,*) "             :", Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax
       call add_warning(warning)
       write(warning,*) "Second stage : fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax"
       call add_warning(warning)
       write(warning,*) "             :", fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax
       call add_warning(warning)

    else if (type == 3) then

       write(warning,*) "First stage  : Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax"
       call add_warning(warning)
       write(warning,*) "             :", Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax
       call add_warning(warning)
       write(warning,*) "Second stage : fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax"
       call add_warning(warning)
       write(warning,*) "             :", fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax
       call add_warning(warning)

    else if (type == 4) then
       
       write(warning,*) "First stage  : fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax"
       call add_warning(warning)
       write(warning,*) "             :", fcondtop, fsurfn, ferrmax, fcondtop - fsurfn - ferrmax
       call add_warning(warning)
       write(warning,*) "Second stage : Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax"
       call add_warning(warning)
       write(warning,*) "             :", Tsf, Tmlt, dTemp_errmax, Tsf - Tmlt - dTemp_errmax
       call add_warning(warning)

    endif

  end subroutine two_stage_inconsistency

!=======================================================================
! Picard/TDMA based solver
!=======================================================================

  subroutine prep_picard(nilyr, nslyr,  &
                         lsnow, zqsn,   &
                         zqin,  zSin,   &
                         hilyr, hslyr,  &
                         km,    ks,     &
                         zTin,  zTsn,   &
                         Sbr,   phi,    &
                         dxp,   kcstar, &
                         einit)

    use ice_mushy_physics, only: &
         temperature_mush, &
         liquidus_brine_salinity_mush, &
         liquid_fraction, &
         temperature_snow

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    logical, intent(in) :: &
         lsnow      ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin   , & ! ice layer enthalpy (J m-3)
         zSin   , & ! ice layer bulk salinity (ppt)
         km     , & ! ice conductivity (W m-1 K-1)
         zqsn   , & ! snow layer enthalpy (J m-3)
         ks         ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr  , & ! ice layer thickness (m)
         hslyr      ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         zTin   , & ! ice layer temperature (C)
         Sbr    , & ! ice layer brine salinity (ppt)
         phi    , & ! ice layer liquid fraction
         zTsn   , & ! snow layer temperature (C)
         dxp    , & ! distances between grid points (m)
         kcstar     ! interface conductivities  (W m-1 K-1)

    real(kind=dbl_kind), intent(out) :: &
         einit      ! initial total energy (J)

    integer(kind=int_kind) :: k

    ! calculate initial ice temperatures
    do k = 1, nilyr
       zTin(k) = temperature_mush(zqin(k), zSin(k))
       Sbr(k)  = liquidus_brine_salinity_mush(zTin(k))
       phi(k)  = liquid_fraction(zTin(k), zSin(k))
    enddo ! k

    if (lsnow) then

       do k = 1, nslyr
          zTsn(k) = temperature_snow(zqsn(k))
       enddo ! k

    endif ! lsnow

    ! interface distances
    call calc_intercell_thickness(nilyr, nslyr, lsnow, hilyr, hslyr, dxp)

    ! interface conductivities
    call calc_intercell_conductivity(lsnow, nilyr, nslyr, &
                                     km, ks, hilyr, hslyr, kcstar)

    ! total energy content
    call total_energy_content(lsnow,        &
                              nilyr, nslyr, &
                              zqin,  zqsn,  &
                              hilyr, hslyr, &
                              einit)

  end subroutine prep_picard

!=======================================================================

  subroutine picard_solver(nilyr,    nslyr,    &
                           lsnow,    lcold,    &
                           Tsf,      zqsn,     &
                           zqin,     zSin,     &
                           zTin,     zTsn,     &
                           phi,      dt,       &
                           hilyr,    hslyr,    &
                           km,       ks,       &
                           Iswabs,   Sswabs,   &
                           Tbot,               &
                           fswint,   fswsfc,   &
                           rhoa,     flw,      &
                           potT,     Qa,       &
                           shcoef,   lhcoef,   &
                           fcondtop, fcondbot, &
                           fadvheat,           &
                           flwoutn,  fsensn,   &
                           flatn,    fsurfn,   &
                           qpond,    qocn,     &
                           Spond,    sss,      &
                           q,        dSdt,     &
                           w,                  &
                           lstop,    stop_label)

    use ice_therm_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    logical, intent(in) :: &
         lsnow         , & ! snow presence: T: has snow, F: no snow
         lcold             ! surface cold: T: surface is cold, F: surface is melting

    real(kind=dbl_kind), intent(inout) :: &
         Tsf               ! snow surface temperature (C)

    real(kind=dbl_kind), intent(out) :: &
         fcondtop      , & ! downward cond flux at top surface (W m-2)
         fcondbot      , & ! downward cond flux at bottom surface (W m-2)
         fadvheat          ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqin          , & ! ice layer enthalpy (J m-3)
         zSin          , & ! ice layer bulk salinity (ppt)
         zTin          , & ! ice layer temperature (C)
         phi           , & ! ice layer liquid fraction
         zqsn          , & ! snow layer enthalpy (J m-3)
         zTsn              ! snow layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         km            , & ! ice conductivity (W m-1 K-1)
         Iswabs        , & ! SW radiation absorbed in ice layers (W m-2)
         dSdt              ! gravity drainage desalination rate for slow mode (ppt s-1)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q                 ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         ks            , & ! snow conductivity (W m-1 K-1)
         Sswabs            ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(out) :: &
         flwoutn       , & ! upward LW at surface (W m-2)
         fsensn        , & ! surface downward sensible heat (W m-2)
         flatn         , & ! surface downward latent heat (W m-2)
         fsurfn            ! net flux to top surface, excluding fcondtop

    real(kind=dbl_kind), intent(in) :: &
         dt            , & ! time step (s)
         hilyr         , & ! ice layer thickness (m)
         hslyr         , & ! snow layer thickness (m)
         Tbot          , & ! ice bottom surfce temperature (deg C)
         fswint        , & ! SW absorbed in ice interior below surface (W m-2)
         fswsfc        , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa          , & ! air density (kg/m^3)
         flw           , & ! incoming longwave radiation (W/m^2)
         potT          , & ! air potential temperature (K)
         Qa            , & ! specific humidity (kg/kg)
         shcoef        , & ! transfer coefficient for sensible heat
         lhcoef        , & ! transfer coefficient for latent heat
         qpond         , & ! melt pond brine enthalpy (J m-3)
         qocn          , & ! ocean brine enthalpy (J m-3)
         Spond         , & ! melt pond salinity (ppt)
         sss           , & ! sea surface salinity (ppt)
         w                 ! vertical flushing Darcy velocity (m/s)
      
    logical(kind=log_kind), intent(inout) :: &
         lstop             ! solver failure flag 

    character(len=*), intent(out) :: &
         stop_label        ! fatal error message

    real(kind=dbl_kind), dimension(nilyr) :: &
         Sbr           , & ! ice layer brine salinity (ppt)
         qbr           , & ! ice layer brine enthalpy (J m-3)
         zTin0         , & ! ice layer temperature (C) at start of timestep
         zqin0         , & ! ice layer enthalpy (J m-3) at start of timestep
         zSin0         , & ! ice layer bulk salinity (ppt) at start of timestep
         zTin_prev         ! ice layer temperature at previous iteration
    
    real(kind=dbl_kind), dimension(nslyr) :: &
         zqsn0         , & ! snow layer enthalpy (J m-3) at start of timestep
         zTsn0         , & ! snow layer temperature (C) at start of timestep
         zTsn_prev         ! snow layer temperature at previous iteration

    real(kind=dbl_kind), dimension(nslyr+nilyr+1) :: &
         dxp           , & ! distances between grid points (m)
         kcstar            ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind) :: &
         Tsf0          , & ! snow surface temperature (C) at start of timestep
         dfsurfn_dTsf  , & ! derivative of net flux to top surface, excluding fcondtopn
         dflwoutn_dTsf , & ! derivative of longwave flux wrt surface temperature
         dfsensn_dTsf  , & ! derivative of sensible heat flux wrt surface temperature
         dflatn_dTsf   , & ! derivative of latent heat flux wrt surface temperature
         Tsf_prev      , & ! snow surface temperature at previous iteration
         einit         , & ! initial total energy (J)
         fadvheat_nit      ! heat to ocean due to advection (W m-2) during iteration

    logical :: &
         lconverged        ! has Picard solver converged?

    integer :: &
         nit               ! Picard iteration count

    integer, parameter :: &
         nit_max = 100     ! maximum number of Picard iterations

    lconverged = .false.

    ! prepare quantities for picard iteration
    call prep_picard(nilyr, nslyr,  &
                     lsnow, zqsn,   &
                     zqin,  zSin,   &
                     hilyr, hslyr,  &
                     km,    ks,     &
                     zTin,  zTsn,   &
                     Sbr,   phi,    &
                     dxp,   kcstar, &
                     einit)

    Tsf0  = Tsf
    zqin0 = zqin
    zqsn0 = zqsn
    zTin0 = zTin
    zTsn0 = zTsn
    zSin0 = zSin

    ! set prev variables
    Tsf_prev  = Tsf
    zTsn_prev = zTsn
    zTin_prev = zTin

    ! picard iteration
    picard: do nit = 1, nit_max

       ! surface heat flux
       call surface_heat_flux(Tsf,     fswsfc, &
                              rhoa,    flw,    &
                              potT,    Qa,     &
                              shcoef,  lhcoef, &
                              flwoutn, fsensn, &
                              flatn,   fsurfn)

       ! derivative of heat flux with respect to surface temperature
       call dsurface_heat_flux_dTsf(Tsf,          fswsfc,        &
                                    rhoa,         flw,           &
                                    potT,         Qa,            &
                                    shcoef,       lhcoef,        &
                                    dfsurfn_dTsf, dflwoutn_dTsf, &
                                    dfsensn_dTsf, dflatn_dTsf)

       ! tridiagonal solve of new temperatures
       call solve_heat_conduction(lsnow,     lcold,        &
                                  nilyr,     nslyr,        &
                                  Tsf,       Tbot,         &
                                  zqin0,     zqsn0,        &
                                  phi,       dt,           &
                                  qpond,     qocn,         &
                                  q,         w,            &
                                  hilyr,     hslyr,        &
                                  dxp,       kcstar,       &
                                  Iswabs,    Sswabs,       &
                                  fsurfn,    dfsurfn_dTsf, &
                                  zTin,      zTsn,nit)

       ! update brine enthalpy
       call picard_updates_enthalpy(nilyr, zTin, qbr)

       ! drainage fluxes
       call picard_drainage_fluxes(fadvheat_nit, q,    &
                                   qbr,          qocn, &
                                   hilyr,        nilyr)

       ! flushing fluxes
       call picard_flushing_fluxes(nilyr,           &
                                   fadvheat_nit, w, &
                                   qbr,             &
                                   qocn, qpond)

       ! perform convergence check
       call check_picard_convergence(nilyr,      nslyr,    &
                                     lsnow,                &
                                     lconverged, nit,      & 
                                     Tsf,        Tsf_prev, &
                                     zTin,       zTin_prev,&
                                     zTsn,       zTsn_prev,&
                                     phi,        Tbot,     &
                                     zqin,       zqsn,     &
                                     km,         ks,       &
                                     hilyr,      hslyr,    &
                                     fswint,               &
                                     einit,      dt,       &
                                     fcondtop,   fcondbot, &
                                     fadvheat_nit)

       if (lconverged) exit

       Tsf_prev  = Tsf
       zTsn_prev = zTsn
       zTin_prev = zTin

    enddo picard

    fadvheat = fadvheat_nit

    ! update the picard iterants
    call picard_updates(nilyr, zTin, &
                        Sbr, qbr)
    
    ! solve for the salinity
    call solve_salinity(zSin,  Sbr,   &
                        Spond, sss,   &
                        q,     dSdt,  &
                        w,     hilyr, &
                        dt,    nilyr)

    ! final surface heat flux
    call surface_heat_flux(Tsf,     fswsfc, &
                           rhoa,    flw,    &
                           potT,    Qa,     &
                           shcoef,  lhcoef, &
                           flwoutn, fsensn, &
                           flatn,   fsurfn)

    ! if not converged
    if (.not. lconverged) then

       call picard_nonconvergence(nilyr, nslyr,&
                                  Tsf0,  Tsf,  &
                                  zTsn0, zTsn, &
                                  zTin0, zTin, &
                                  zSin0, zSin, &
                                  zqsn0, zqsn, &
                                  zqin0, phi)
       lstop = .true.
       stop_label = "picard_solver: Picard solver non-convergence"

    endif

  end subroutine picard_solver

!=======================================================================

  subroutine picard_nonconvergence(nilyr, nslyr,&
                                   Tsf0,  Tsf,  &
                                   zTsn0, zTsn, &
                                   zTin0, zTin, &
                                   zSin0, zSin, &
                                   zqsn0, zqsn, &
                                   zqin0, phi)

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), intent(in) :: &
         Tsf0  , & ! snow surface temperature (C) at beginning of timestep
         Tsf       ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTsn0 , & ! snow layer temperature (C) at beginning of timestep
         zTsn  , & ! snow layer temperature (C)
         zqsn0 , &
         zqsn  , &
         zTin0 , & ! ice layer temperature (C)
         zTin  , & ! ice layer temperature (C)
         zSin0 , & ! ice layer bulk salinity (ppt)
         zSin  , & ! ice layer bulk salinity (ppt)
         phi   , & ! ice layer liquid fraction
         zqin0

    integer :: &
         k        ! vertical layer index

    character(len=char_len_long) :: &
         warning  ! warning message
    
    write(warning,*) "-------------------------------------"
    call add_warning(warning)

    write(warning,*) "picard convergence failed!"
    call add_warning(warning)
    write(warning,*) 0, Tsf0, Tsf
    call add_warning(warning)
    
    do k = 1, nslyr
       write(warning,*) k, zTsn0(k), zTsn(k), zqsn0(k)
       call add_warning(warning)
    enddo ! k          
    
    do k = 1, nilyr
       write(warning,*) k, zTin0(k), zTin(k), zSin0(k), zSin(k), phi(k), zqin0(k)
       call add_warning(warning)
    enddo ! k

    write(warning,*) "-------------------------------------"
    call add_warning(warning)

  end subroutine picard_nonconvergence

!=======================================================================
  
  subroutine check_picard_convergence(nilyr,      nslyr,    &
                                      lsnow,                &
                                      lconverged, nit,      & 
                                      Tsf,        Tsf_prev, &
                                      zTin,       zTin_prev,&
                                      zTsn,       zTsn_prev,&
                                      phi,        Tbot,     &
                                      zqin,       zqsn,     &
                                      km,         ks,       &
                                      hilyr,      hslyr,    &
                                      fswint,               &
                                      einit,      dt,       &
                                      fcondtop,   fcondbot, &
                                      fadvheat)

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    logical, intent(inout) :: &
         lconverged   ! has Picard solver converged?

    logical, intent(in) :: &
         lsnow        ! snow presence: T: has snow, F: no snow

    integer, intent(in) :: &
         nit          ! Picard iteration count

    real(kind=dbl_kind), intent(in) :: &
         dt       , & ! time step (s)
         Tsf      , & ! snow surface temperature (C)
         Tsf_prev , & ! snow surface temperature at previous iteration
         hilyr    , & ! ice layer thickness (m)
         hslyr    , & ! snow layer thickness (m)
         fswint   , & ! SW absorbed in ice interior below surface (W m-2)
         einit    , & ! initial total energy (J)
         Tbot     , & ! ice bottom surfce temperature (deg C)
         fadvheat     ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTin     , & ! ice layer temperature (C)
         zTin_prev, & ! ice layer temperature at previous iteration
         phi      , & ! ice layer liquid fraction
         km           ! ice conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         zqin         ! ice layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTsn     , & ! snow layer temperature (C)
         zTsn_prev, & ! snow layer temperature at previous iteration
         ks           ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         zqsn         ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), intent(out) :: &    
         fcondtop , & ! downward cond flux at top surface (W m-2)
         fcondbot     ! downward cond flux at bottom surface (W m-2)

    real(kind=dbl_kind) :: &
         ferr     , & ! energy flux error
         efinal   , & ! initial total energy (J) at iteration
         dzTsn    , & ! change in snow temperature (C) between iterations
         dzTin    , & ! change in ice temperature (C) between iterations
         dTsf         ! change in surface temperature (C) between iterations

    call picard_final(lsnow,        &
                      nilyr, nslyr, &
                      zqin,  zqsn,  &
                      zTin,  zTsn,  &
                      phi)

    call total_energy_content(lsnow,         &
                              nilyr,  nslyr, &
                              zqin,   zqsn,  &
                              hilyr,  hslyr, &
                              efinal)

    call maximum_variables_changes(lsnow,                  &
                                   Tsf,  Tsf_prev,  dTsf,  &
                                   zTsn, zTsn_prev, dzTsn, & 
                                   zTin, zTin_prev, dzTin)

    fcondbot = c2 * km(nilyr) * ((zTin(nilyr) - Tbot) / hilyr)

    if (lsnow) then
       fcondtop = c2 * ks(1) * ((Tsf - zTsn(1)) / hslyr)
    else
       fcondtop = c2 * km(1) * ((Tsf - zTin(1)) / hilyr)
    endif

    ferr = (efinal - einit) / dt - (fcondtop - fcondbot + fswint - fadvheat)

    lconverged = (dTsf  < dTemp_errmax .and. &
                  dzTsn < dTemp_errmax .and. &
                  dzTin < dTemp_errmax .and. &
                  abs(ferr) < 0.9_dbl_kind*ferrmax)

  end subroutine check_picard_convergence

!=======================================================================

  subroutine picard_drainage_fluxes(fadvheat, q,    &
                                    qbr,      qocn, &
                                    hilyr,    nilyr)

    integer (kind=int_kind), intent(in) :: &
         nilyr        ! number of ice layers

    real(kind=dbl_kind), intent(out) :: &
         fadvheat ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q        ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         qbr      ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         qocn , & ! ocean brine enthalpy (J m-3)
         hilyr    ! ice layer thickness (m)

    integer :: &
         k        ! vertical layer index

    fadvheat = c0

    ! calculate fluxes from base upwards
    do k = 1, nilyr-1

       fadvheat = fadvheat - q(k) * (qbr(k+1) - qbr(k))

    enddo ! k

    k = nilyr

    fadvheat = fadvheat - q(k) * (qocn - qbr(k))

  end subroutine picard_drainage_fluxes

!=======================================================================

  subroutine picard_flushing_fluxes(nilyr,         &
                                    fadvheat, w,   &
                                    qbr,           &
                                    qocn,     qpond)

    integer (kind=int_kind), intent(in) :: &
         nilyr        ! number of ice layers

    real(kind=dbl_kind), intent(inout) :: &
         fadvheat  ! flow of heat to ocean due to advection (W m-2)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         qbr       ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         w     , & ! vertical flushing Darcy velocity (m/s)
         qocn  , & ! ocean brine enthalpy (J m-3)
         qpond     ! melt pond brine enthalpy (J m-3)

    fadvheat = fadvheat + w * (qbr(nilyr) - qpond)

  end subroutine picard_flushing_fluxes

!=======================================================================

  subroutine maximum_variables_changes(lsnow,                  &
                                       Tsf,  Tsf_prev,  dTsf,  &
                                       zTsn, zTsn_prev, dzTsn, & 
                                       zTin, zTin_prev, dzTin)

    logical, intent(in) :: &
         lsnow        ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), intent(in) :: &
         Tsf      , & ! snow surface temperature (C)
         Tsf_prev     ! snow surface temperature at previous iteration

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTsn     , & ! snow layer temperature (C)
         zTsn_prev, & ! snow layer temperature at previous iteration
         zTin     , & ! ice layer temperature (C)
         zTin_prev    ! ice layer temperature at previous iteration

    real(kind=dbl_kind), intent(out) :: &
         dTsf     , & ! change in surface temperature (C) between iterations
         dzTsn    , & ! change in snow temperature (C) between iterations
         dzTin        ! change in surface temperature (C) between iterations

    dTsf = abs(Tsf - Tsf_prev)

    if (lsnow) then
       dzTsn = maxval(abs(zTsn - zTsn_prev))
    else ! lsnow
       dzTsn = c0 
    endif ! lsnow

    dzTin = maxval(abs(zTin - zTin_prev))

  end subroutine maximum_variables_changes

!=======================================================================

  subroutine total_energy_content(lsnow,         &
                                  nilyr,  nslyr, &
                                  zqin,   zqsn,  &
                                  hilyr,  hslyr, &
                                  energy)

    logical, intent(in) :: &
         lsnow     ! snow presence: T: has snow, F: no snow

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin  , & ! ice layer enthalpy (J m-3)
         zqsn      ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         hslyr     ! snow layer thickness (m)

    real(kind=dbl_kind), intent(out) :: &    
         energy    ! total energy of ice and snow

    integer :: &
         k         ! vertical layer index

    energy = c0
    
    if (lsnow) then
       
       do k = 1, nslyr

          energy = energy + hslyr * zqsn(k)

       enddo ! k

    endif ! lsnow

    do k = 1, nilyr

       energy = energy + hilyr * zqin(k)

    enddo ! k

  end subroutine total_energy_content

!=======================================================================

  subroutine picard_updates(nilyr, zTin, &
                            Sbr,   qbr)

    ! update brine salinity and liquid fraction based on new temperatures

    use ice_mushy_physics, only: &
         liquidus_brine_salinity_mush, &
         enthalpy_brine

    integer (kind=int_kind), intent(in) :: &
         nilyr   ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTin    ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         Sbr , & ! ice layer brine salinity (ppt)
         qbr     ! ice layer brine enthalpy (J m-3)

    integer :: &
         k       ! vertical layer index

    do k = 1, nilyr

       Sbr(k) = liquidus_brine_salinity_mush(zTin(k))
       qbr(k) = enthalpy_brine(zTin(k))

    enddo ! k

  end subroutine picard_updates

!=======================================================================

  subroutine picard_updates_enthalpy(nilyr, zTin, qbr)

    ! update brine salinity and liquid fraction based on new temperatures

    use ice_mushy_physics, only: &
         enthalpy_brine

    integer (kind=int_kind), intent(in) :: &
         nilyr   ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTin ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         qbr  ! ice layer brine enthalpy (J m-3)

    integer :: &
         k    ! vertical layer index

    do k = 1, nilyr

       qbr(k) = enthalpy_brine(zTin(k))

    enddo ! k

  end subroutine picard_updates_enthalpy

!=======================================================================

  subroutine picard_final(lsnow,        &
                          nilyr, nslyr, &
                          zqin,  zqsn,  &
                          zTin,  zTsn,  &
                          phi)

    use ice_mushy_physics, only: &
         enthalpy_mush_liquid_fraction, &
         enthalpy_snow

    integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nslyr    ! number of snow layers

    logical, intent(in) :: &
         lsnow   ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         zqin, & ! ice layer enthalpy (J m-3)
         zqsn    ! snow layer enthalpy (J m-3)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         phi , & ! ice layer liquid fraction
         zTsn    ! snow layer temperature (C)

    integer :: &
         k       ! vertical layer index

    do k = 1, nilyr
       zqin(k) = enthalpy_mush_liquid_fraction(zTin(k), phi(k))
    enddo ! k

    if (lsnow) then

       do k = 1, nslyr
          zqsn(k) = enthalpy_snow(zTsn(k))
       enddo ! k

    endif ! lsnow

  end subroutine picard_final

!=======================================================================

  subroutine calc_intercell_thickness(nilyr, nslyr, lsnow, hilyr, hslyr, dxp)

    integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nslyr    ! number of snow layers

    logical, intent(in) :: &
         lsnow     ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), intent(in) :: &
         hilyr , & ! ice layer thickness (m)
         hslyr     ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         dxp       ! distances between grid points (m)

    integer :: &
         l         ! vertical index

    if (lsnow) then

       dxp(1) = hslyr / c2
       
       do l = 2, nslyr

          dxp(l) = hslyr

       enddo ! l

       dxp(nslyr+1) = (hilyr + hslyr) / c2

       do l = nslyr+2, nilyr+nslyr

          dxp(l) = hilyr

       enddo ! l

       dxp(nilyr+nslyr+1) = hilyr / c2

    else ! lsnow

       dxp(1) = hilyr / c2

       do l = 2, nilyr

          dxp(l) = hilyr

       enddo ! l

       dxp(nilyr+1) = hilyr / c2

       do l = nilyr+2, nilyr+nslyr+1

          dxp(l) = c0

       enddo ! l

    endif ! lsnow

  end subroutine calc_intercell_thickness

!=======================================================================

  subroutine calc_intercell_conductivity(lsnow,        &
                                         nilyr, nslyr, &
                                         km,    ks,    &
                                         hilyr, hslyr, &
                                         kcstar)

    integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nslyr    ! number of snow layers

    logical, intent(in) :: &
         lsnow      ! snow presence: T: has snow, F: no snow

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         km     , & ! ice conductivity (W m-1 K-1)
         ks         ! snow conductivity (W m-1 K-1)

    real(kind=dbl_kind), intent(in) :: &
         hilyr  , & ! ice layer thickness (m)
         hslyr      ! snow layer thickness (m)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         kcstar     ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind) :: &
         fe         ! distance fraction at interface

    integer :: &
         k, &       ! vertical layer index
         l          ! vertical index

    if (lsnow) then

       kcstar(1) = ks(1)
       
       do l = 2, nslyr

          k = l
          kcstar(l) = (c2 * ks(k) * ks(k-1)) / (ks(k) + ks(k-1))

       enddo ! l

       fe = hilyr / (hilyr + hslyr)
       kcstar(nslyr+1) = c1 / ((c1 - fe) / ks(nslyr) + fe / km(1))

       do k = 2, nilyr

          l = k + nslyr
          kcstar(l) = (c2 * km(k) * km(k-1)) / (km(k) + km(k-1))

       enddo ! k

       kcstar(nilyr+nslyr+1) = km(nilyr)

    else ! lsnow

       kcstar(1) = km(1)

       do k = 2, nilyr
          
          l = k
          kcstar(l) = (c2 * km(k) * km(k-1)) / (km(k) + km(k-1))

       enddo ! k

       kcstar(nilyr+1) = km(nilyr)

       do l = nilyr+2, nilyr+nslyr+1

          kcstar(l) = c0

       enddo ! l

    endif ! lsnow

  end subroutine calc_intercell_conductivity

!=======================================================================
  
  subroutine solve_heat_conduction(lsnow,  lcold,        &
                                   nilyr,  nslyr,        &
                                   Tsf,    Tbot,         &
                                   zqin0,  zqsn0,        &
                                   phi,    dt,           &
                                   qpond,  qocn,         &
                                   q,      w,            &
                                   hilyr,  hslyr,        &
                                   dxp,    kcstar,       &
                                   Iswabs, Sswabs,       &
                                   fsurfn, dfsurfn_dTsf, &
                                   zTin,   zTsn,nit)

    logical, intent(in) :: &
         lsnow        , & ! snow presence: T: has snow, F: no snow
         lcold            ! surface cold: T: surface is cold, F: surface is melting

    integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nslyr    ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0        , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         zqsn0        , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)
    
    real(kind=dbl_kind), intent(inout) :: &
         Tsf              ! snow surface temperature (C)

    real(kind=dbl_kind), intent(in) :: &
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zTin             ! ice layer temperature (C)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         zTsn             ! snow layer temperature (C)

    integer, intent(in) :: &
         nit              ! Picard iteration count

    real(kind=dbl_kind), dimension(nilyr+nslyr+1) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b            , & ! right hand side of matrix solve
         T                ! ice and snow temperatures

    integer :: &
         nyn              ! matrix size

    ! set up matrix and right hand side - snow
    if (lsnow) then

       if (lcold) then

          call matrix_elements_snow_cold(Ap, As, An, b, nyn,   &
                                         nilyr,  nslyr,        &
                                         Tsf,    Tbot,         &
                                         zqin0,  zqsn0,        &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         w,                    &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)

       else ! lcold

          call matrix_elements_snow_melt(Ap, As, An, b, nyn,   &
                                         nilyr,  nslyr,        &
                                         Tsf,    Tbot,         &
                                         zqin0,  zqsn0,        &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         w,                    &
                                         hilyr,  hslyr,        &
                                         dxp,    kcstar,       &
                                         Iswabs, Sswabs,       &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)

       endif ! lcold

    else ! lsnow

       if (lcold) then

          call matrix_elements_nosnow_cold(Ap, As, An, b, nyn,   &
                                           nilyr,  nslyr,        &
                                           Tsf,    Tbot,         &
                                           zqin0,                &
                                           qpond,  qocn,         &
                                           phi,    q,            &
                                           w,                    &
                                           hilyr,                &
                                           dxp,    kcstar,       &
                                           Iswabs,               &
                                           fsurfn, dfsurfn_dTsf, &
                                           dt)

       else ! lcold

          call matrix_elements_nosnow_melt(Ap, As, An, b, nyn,   &
                                           nilyr,  nslyr,        &
                                           Tsf,    Tbot,         &
                                           zqin0,                &
                                           qpond,  qocn,         &
                                           phi,    q,            &
                                           w,                    &
                                           hilyr,                &
                                           dxp,    kcstar,       &
                                           Iswabs,               &
                                           fsurfn, dfsurfn_dTsf, &
                                           dt)

       endif ! lcold

    endif ! lsnow

    ! tridiag to get new temperatures
    call tdma_solve_sparse(nilyr, nslyr, &
              An(1:nyn), Ap(1:nyn), As(1:nyn), b(1:nyn), T(1:nyn), nyn)

    call update_temperatures(lsnow, lcold, &
                             nilyr, nslyr, &
                             T,     Tsf,   &
                             zTin,  zTsn)

  end subroutine solve_heat_conduction

!=======================================================================

  subroutine update_temperatures(lsnow, lcold, &
                                 nilyr, nslyr, &
                                 T,     Tsf,   &
                                 zTin,  zTsn)

    logical, intent(in) :: &
         lsnow , & ! snow presence: T: has snow, F: no snow
         lcold     ! surface cold: T: surface is cold, F: surface is melting

      integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         T         ! matrix solution vector

    real(kind=dbl_kind), intent(inout) :: &
         Tsf       ! snow surface temperature (C)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zTin  , & ! ice layer temperature (C)
         zTsn      ! snow layer temperature (C)

    integer :: &
         l     , & ! vertical index
         k         ! vertical layer index

    if (lsnow) then
       
       if (lcold) then

          Tsf = T(1)

          do k = 1, nslyr
             l = k + 1
             zTsn(k) = T(l)
          enddo ! k
          
          do k = 1, nilyr
             l = k + nslyr + 1
             zTin(k) = T(l)
          enddo ! k

       else ! lcold
          
          do k = 1, nslyr
             l = k
             zTsn(k) = T(l)
          enddo ! k
          
          do k = 1, nilyr
             l = k + nslyr
             zTin(k) = T(l)
          enddo ! k

       endif ! lcold

    else ! lsnow
       
       if (lcold) then

          Tsf = T(1)
          
          do k = 1, nilyr
             l = k + 1
             zTin(k) = T(l)
          enddo ! k

       else ! lcold

          do k = 1, nilyr
             l = k
             zTin(k) = T(l)
          enddo ! k

       endif ! lcold

    endif ! lsnow

  end subroutine update_temperatures

!=======================================================================
  
  subroutine matrix_elements_nosnow_melt(Ap, As, An, b, nyn,   &
                                         nilyr,  nslyr,        &
                                         Tsf,    Tbot,         &
                                         zqin0,                &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         w,                    &
                                         hilyr,                &
                                         dxp,    kcstar,       &
                                         Iswabs,               &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)
     
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0        , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi              ! ice layer liquid fraction

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index
    
    ! surface layer
    k = 1
    l = k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = -kcstar(k+1) / dxp(k+1) - &
             q(k) * cp_ocn * rhow
    An(l) = c0
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k) / dxp(k)) * Tsf + &
             w    * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = k
          
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(k+1) / dxp(k+1) + &
                kcstar(k)   / dxp(k) + &
                q(k) * cp_ocn * rhow + &
                w    * cp_ocn * rhow
       As(l) = -kcstar(k+1) / dxp(k+1) - &
                q(k) * cp_ocn * rhow
       An(l) = -kcstar(k)   / dxp(k) - &
                w    * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k
    
    ! bottom layer
    k = nilyr
    l = k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(k) / dxp(k) - &
            w     * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k+1) * Tbot) / dxp(k+1) + &
            q(k)  * qocn
    
    nyn = nilyr

  end subroutine matrix_elements_nosnow_melt

!=======================================================================

  subroutine matrix_elements_nosnow_cold(Ap, As, An, b, nyn,   &
                                         nilyr,  nslyr,        &
                                         Tsf,    Tbot,         &
                                         zqin0,                &
                                         qpond,  qocn,         &
                                         phi,    q,            &
                                         w,                    &
                                         hilyr,                &
                                         dxp,    kcstar,       &
                                         Iswabs,               &
                                         fsurfn, dfsurfn_dTsf, &
                                         dt)
     
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0        , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi              ! ice layer liquid fraction

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index

    ! surface temperature
    l = 1
    Ap(l) = dfsurfn_dTsf - kcstar(1) / dxp(1)
    As(l) = kcstar(1) / dxp(1)
    An(l) = c0
    b (l) = dfsurfn_dTsf * Tsf - fsurfn

    ! surface layer
    k = 1
    l = k + 1
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = -kcstar(k+1) / dxp(k+1) - &
             q(k) * cp_ocn * rhow
    An(l) = -kcstar(k)   / dxp(k)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
             w    * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = k + 1
          
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(k+1) / dxp(k+1) + &
                kcstar(k)   / dxp(k) + &
                q(k) * cp_ocn * rhow + &
                w    * cp_ocn * rhow
       As(l) = -kcstar(k+1) / dxp(k+1) - &
                q(k) * cp_ocn * rhow
       An(l) = -kcstar(k)   / dxp(k) - &
                w    * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k
    
    ! bottom layer
    k = nilyr
    l = k + 1
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(k+1) / dxp(k+1) + &
             kcstar(k)   / dxp(k) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(k) / dxp(k) - &
            w     * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(k+1) * Tbot) / dxp(k+1) + &
            q(k)  * qocn
    
    nyn = nilyr + 1

  end subroutine matrix_elements_nosnow_cold

!=======================================================================

  subroutine matrix_elements_snow_melt(Ap, As, An, b, nyn,   &
                                       nilyr,  nslyr,        &
                                       Tsf,    Tbot,         &
                                       zqin0,  zqsn0,        &
                                       qpond,  qocn,         &
                                       phi,    q,            &
                                       w,                    &
                                       hilyr,  hslyr,        &
                                       dxp,    kcstar,       &
                                       Iswabs, Sswabs,       &
                                       fsurfn, dfsurfn_dTsf, &
                                       dt)
     
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0        , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         zqsn0        , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l                ! vertical index

    ! surface layer
    k = 1
    l = k
    
    Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l)
    As(l) = -kcstar(l+1) / dxp(l+1)
    An(l) = c0
    b (l) = ((rhos * Lfresh + zqsn0(k)) / dt) * hslyr + Sswabs(k) + &
            (kcstar(l) * Tsf) / dxp(l)

    ! interior snow layers
    do k = 2, nslyr
       
       l = k

       Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
                kcstar(l+1) / dxp(l+1) + &
                kcstar(l)   / dxp(l)
       As(l) = -kcstar(l+1) / dxp(l+1)
       An(l) = -kcstar(l)   / dxp(l)
       b (l) = ((rhos * Lfresh + zqsn0(k)) / dt) * hslyr + Sswabs(k)
          
    enddo ! k
    
    ! top ice layer
    k = 1
    l = nslyr + k
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = -kcstar(l+1) / dxp(l+1) - &
             q(k) * cp_ocn * rhow
    An(l) = -kcstar(l)   / dxp(l)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
             w    * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = nslyr + k
       
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(l+1) / dxp(l+1) + &
                kcstar(l)   / dxp(l) + &
                q(k) * cp_ocn * rhow + &
                w    * cp_ocn * rhow
       As(l) = -kcstar(l+1) / dxp(l+1) - &
                q(k) * cp_ocn * rhow
       An(l) = -kcstar(l)   / dxp(l) - &
                w    * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k)
          
    enddo ! k

    ! bottom layer
    k = nilyr
    l = nilyr + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(l+1) / dxp(l+1) + &
             kcstar(l)   / dxp(l) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(l)   / dxp(l) - &
             w    * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(l+1) * Tbot) / dxp(l+1) + &
             q(k) * qocn
       
    nyn = nilyr + nslyr

  end subroutine matrix_elements_snow_melt

!=======================================================================

  subroutine matrix_elements_snow_cold(Ap, As, An, b, nyn,   &
                                       nilyr,  nslyr,        &
                                       Tsf,    Tbot,         &
                                       zqin0,  zqsn0,        &
                                       qpond,  qocn,         &
                                       phi,    q,            &
                                       w,                    &
                                       hilyr,  hslyr,        &
                                       dxp,    kcstar,       &
                                       Iswabs, Sswabs,       &
                                       fsurfn, dfsurfn_dTsf, &
                                       dt)
     
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Ap           , & ! diagonal of tridiagonal matrix
         As           , & ! lower off-diagonal of tridiagonal matrix
         An           , & ! upper off-diagonal of tridiagonal matrix
         b                ! right hand side of matrix solve

    integer, intent(out) :: &
         nyn              ! matrix size

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin0        , & ! ice layer enthalpy (J m-3) at beggining of timestep
         Iswabs       , & ! SW radiation absorbed in ice layers (W m-2)
         phi          , & ! ice layer liquid fraction
         zqsn0        , & ! snow layer enthalpy (J m-3) at start of timestep
         Sswabs           ! SW radiation absorbed in snow layers (W m-2)

    real(kind=dbl_kind), intent(in) :: &
         Tsf          , & ! snow surface temperature (C)
         dt           , & ! timestep (s)
         hilyr        , & ! ice layer thickness (m)
         hslyr        , & ! snow layer thickness (m)
         Tbot         , & ! ice bottom surfce temperature (deg C)
         qpond        , & ! melt pond brine enthalpy (J m-3)
         qocn         , & ! ocean brine enthalpy (J m-3)
         w            , & ! downwards vertical flushing Darcy velocity (m/s)
         fsurfn       , & ! net flux to top surface, excluding fcondtop
         dfsurfn_dTsf     ! derivative of net flux to top surface, excluding fcondtopn

    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &    
         q                ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         dxp          , & ! distances between grid points (m)
         kcstar           ! interface conductivities (W m-1 K-1)

    integer :: &
         k            , & ! vertical layer index
         l            , & ! matrix index
         m                ! vertical index

    ! surface temperature
    l = 1
    Ap(l) = dfsurfn_dTsf - kcstar(1) / dxp(1)
    As(l) = kcstar(1) / dxp(1)
    An(l) = c0
    b (l) = dfsurfn_dTsf * Tsf - fsurfn

    ! surface layer
    k = 1
    l = k + 1
    m = 1
    
    Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m)
    As(l) = -kcstar(m+1) / dxp(m+1)
    An(l) = -kcstar(m)   / dxp(m)
    b (l) = ((rhos * Lfresh + zqsn0(k)) / dt) * hslyr + Sswabs(k)

    ! interior snow layers
    do k = 2, nslyr
       
       l = k + 1
       m = k

       Ap(l) = ((rhos * cp_ice) / dt) * hslyr + &
                kcstar(m+1) / dxp(m+1) + &
                kcstar(m)   / dxp(m)
       As(l) = -kcstar(m+1) / dxp(m+1)
       An(l) = -kcstar(m)   / dxp(m)
       b (l) = ((rhos * Lfresh + zqsn0(k)) / dt) * hslyr + Sswabs(k)

    enddo ! k
    
    ! top ice layer
    k = 1
    l = nslyr + k + 1
    m = k + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = -kcstar(m+1) / dxp(m+1) - &
             q(k) * cp_ocn * rhow
    An(l) = -kcstar(m)   / dxp(m)
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
              w   * qpond

    ! interior ice layers
    do k = 2, nilyr-1
          
       l = nslyr + k + 1
       m = k + nslyr
       
       Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
                kcstar(m+1) / dxp(m+1) + &
                kcstar(m)   / dxp(m) + &
                q(k) * cp_ocn * rhow + &
                w    * cp_ocn * rhow
       As(l) = -kcstar(m+1) / dxp(m+1) - &
                q(k) * cp_ocn * rhow
       An(l) = -kcstar(m)   / dxp(m) - &
                w    * cp_ocn * rhow
       b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k)

    enddo ! k

    ! bottom layer
    k = nilyr
    l = nilyr + nslyr + 1
    m = k + nslyr
    
    Ap(l) = ((phi(k) * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice) / dt) * hilyr + &
             kcstar(m+1) / dxp(m+1) + &
             kcstar(m)   / dxp(m) + &
             q(k) * cp_ocn * rhow + &
             w    * cp_ocn * rhow
    As(l) = c0
    An(l) = -kcstar(m)   / dxp(m) - &
             w    * cp_ocn * rhow
    b (l) = (((c1 - phi(k)) * rhoi * Lfresh + zqin0(k)) / dt) * hilyr + Iswabs(k) + &
            (kcstar(m+1) * Tbot) / dxp(m+1) + &
             q(k) * qocn

    nyn = nilyr + nslyr + 1

  end subroutine matrix_elements_snow_cold

!=======================================================================

  subroutine solve_salinity(zSin,   Sbr,   &
                            Spond, sss,   &
                            q,     dSdt,  &
                            w,     hilyr, &
                            dt,    nilyr)

    integer (kind=int_kind), intent(in) :: &
         nilyr      ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zSin       ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         Sbr   , & ! ice layer brine salinity (ppt)
         dSdt      ! gravity drainage desalination rate for slow mode (ppt s-1)
         
    real(kind=dbl_kind), dimension(0:nilyr), intent(in) :: &
         q         ! upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), intent(in) :: &
         Spond , & ! melt pond salinity (ppt)
         sss   , & ! sea surface salinity (ppt)
         w     , & ! vertical flushing Darcy velocity (m/s)
         hilyr , & ! ice layer thickness (m)
         dt        ! timestep (s)

    integer :: &
         k         ! vertical layer index

    real(kind=dbl_kind), parameter :: &
         S_min = p01
    
    real(kind=dbl_kind), dimension(nilyr) :: &
         zSin0

    zSin0 = zSin

    k = 1
    zSin(k) = zSin(k) + max(S_min - zSin(k), &
         ((q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr + &
          dSdt(k)                               + &
          (w * (Spond    - Sbr(k))) / hilyr) * dt)

    do k = 2, nilyr-1

       zSin(k) = zSin(k) + max(S_min - zSin(k), &
            ((q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr + &
             dSdt(k)                               + &
             (w * (Sbr(k-1) - Sbr(k))) / hilyr) * dt)

    enddo ! k

    k = nilyr
    zSin(k) = zSin(k) + max(S_min - zSin(k), &
         ((q(k)  * (sss      - Sbr(k))) / hilyr + &
          dSdt(k)                               + &
          (w * (Sbr(k-1) - Sbr(k))) / hilyr) * dt)


    if (minval(zSin) < c0) then


         write(*,*) (q(k)  * (Sbr(k+1) - Sbr(k))) / hilyr, &
          dSdt(k)                               , &
          (w * (Spond    - Sbr(k))) / hilyr

       do k = 1, nilyr
       
          write(*,*) k, zSin(k), zSin0(k)

       enddo

       stop

    endif

  end subroutine solve_salinity

!=======================================================================

  subroutine tdma_solve_sparse(nilyr, nslyr, a, b, c, d, x, n)
    
    ! perform a tri-diagonal solve with TDMA using a sparse tridiagoinal matrix

    integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nslyr    ! number of snow layers

    integer(kind=int_kind), intent(in) :: &
         n      ! matrix size
    
    real(kind=dbl_kind), dimension(:), intent(in) :: &
         a  , & ! matrix lower off-diagonal
         b  , & ! matrix diagonal
         c  , & ! matrix upper off-diagonal
         d      ! right hand side vector
 
    real(kind=dbl_kind), dimension(:), intent(out) :: &
         x      ! solution vector

    real(kind=dbl_kind), dimension(nilyr+nslyr+1) :: &
         cp , & ! modified upper off-diagonal vector
         dp     ! modified right hand side vector

    integer(kind=int_kind) :: &
         i      ! vector index
    
    ! forward sweep
    cp(1) = c(1) / b(1)
    do i = 2, n-1
       cp(i) = c(i) / (b(i) - cp(i-1)*a(i))
    enddo
    
    dp(1) = d(1) / b(1)
    do i = 2, n
       dp(i) = (d(i) - dp(i-1)*a(i)) / (b(i) - cp(i-1)*a(i))
    enddo

    ! back substitution
    x(n) = dp(n)
    do i = n-1,1,-1
       x(i) = dp(i) - cp(i)*x(i+1)
    enddo

  end subroutine tdma_solve_sparse

!=======================================================================
! Effect of salinity
!=======================================================================

  function permeability(phi) result(perm)

    ! given the liquid fraction calculate the permeability
    ! See Golden et al. 2007

    real(kind=dbl_kind), intent(in) :: &
         phi                  ! liquid fraction

    real(kind=dbl_kind) :: &
         perm                 ! permeability (m2)
    
    real(kind=dbl_kind), parameter :: &
         phic = p05 ! critical liquid fraction for impermeability

    perm = 3.0e-8_dbl_kind * max(phi - phic, c0)**3

  end function permeability

!=======================================================================

  subroutine explicit_flow_velocities(nilyr, zSin, &
                                      zTin,  Tsf,  &
                                      Tbot,  q,    &
                                      dSdt,  Sbr,  &
                                      qbr,   dt,   &
                                      sss,   qocn, &
                                      hilyr, hin)

    ! calculate the rapid gravity drainage mode Darcy velocity and the
    ! slow mode drainage rate

    use ice_mushy_physics, only: &
         liquidus_brine_salinity_mush, &
         liquid_fraction, &
         enthalpy_brine, &
         density_brine

    integer (kind=int_kind), intent(in) :: &
         nilyr     ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zSin, &   ! ice layer bulk salinity (ppt)
         zTin      ! ice layer temperature (C)

    real(kind=dbl_kind), intent(in) :: &
         Tsf   , & ! ice/snow surface temperature (C)
         Tbot  , & ! ice bottom temperature (C)
         dt    , & ! time step (s)
         sss   , & ! sea surface salinty (ppt)
         qocn  , & ! ocean enthalpy (J m-3)
         hilyr , & ! ice layer thickness (m)
         hin       ! ice thickness (m)

    real(kind=dbl_kind), dimension(0:nilyr), intent(out) :: &
         q         ! rapid mode upward interface vertical Darcy flow (m s-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         dSdt      ! slow mode drainage rate (ppt s-1)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         Sbr , &   ! ice layer brine salinity (ppt)
         qbr       ! ice layer brine enthalpy (J m-3)

    real(kind=dbl_kind), parameter :: &
         kappal        = 8.824e-8_dbl_kind, & ! heat diffusivity of liquid
         ra_constants  = gravit / (viscosity_dyn * kappal), & ! for Rayleigh number
         fracmax       = p2               , & ! limiting advective layer fraction
         zSin_min      = p1               , & ! minimum bulk salinity (ppt)
         safety_factor = c10                  ! to prevent negative salinities

    real(kind=dbl_kind), dimension(1:nilyr) :: &
         phi           ! ice layer liquid fraction

    real(kind=dbl_kind), dimension(0:nilyr) :: &
         rho           ! ice layer brine density (kg m-3)

    real(kind=dbl_kind) :: &
         rho_ocn   , & ! ocean density (kg m-3)
         perm_min  , & ! minimum permeability from layer to ocean (m2)
         perm_harm , & ! harmonic mean of permeability from layer to ocean (m2)
         rho_sum   , & ! sum of the brine densities from layer to ocean (kg m-3)
         rho_pipe  , & ! density of the brine in the channel (kg m-3)
         z         , & ! distance to layer from top surface (m)
         perm      , & ! ice layer permeability (m2)
         drho      , & ! brine density difference between layer and ocean (kg m-3)
         Ra        , & ! local mush Rayleigh number
         rn        , & ! real value of number of layers considered
         L         , & ! thickness of the layers considered (m)
         dx        , & ! horizontal size of convective flow (m)
         dx2       , & ! square of the horizontal size of convective flow (m2)
         Am        , & ! A parameter for mush
         Bm        , & ! B parameter for mush
         Ap        , & ! A parameter for channel
         Bp        , & ! B parameter for channel
         qlimit    , & ! limit to vertical Darcy flow for numerical stability
         dS_guess  , & ! expected bulk salinity without limits
         alpha         ! desalination limiting factor

    integer(kind=int_kind) :: &
         k             ! ice layer index

    ! initial downward sweep - determine derived physical quantities
    do k = 1, nilyr
       
       Sbr(k) = liquidus_brine_salinity_mush(zTin(k))
       phi(k) = liquid_fraction(zTin(k), zSin(k))
       qbr(k) = enthalpy_brine(zTin(k))
       rho(k) = density_brine(Sbr(k))

    enddo ! k

    rho(0) = rho(1)

    ! ocean conditions
    Sbr(nilyr+1) = sss
    qbr(nilyr+1) = qocn
    rho_ocn = density_brine(sss)

    ! initialize accumulated quantities
    perm_min = bignum
    perm_harm = c0
    rho_sum = c0

    ! limit to q for numerical stability
    qlimit = (fracmax * hilyr) / dt

    ! no flow through ice top surface
    q(0) = c0

    ! first iterate over layers going up
    do k = nilyr, 1, -1

       ! vertical position from ice top surface
       z = ((real(k, dbl_kind) - p5) / real(nilyr, dbl_kind)) * hin

       ! permeabilities
       perm = permeability(phi(k))
       perm_min = min(perm_min,perm)
       perm_harm = perm_harm + (c1 / max(perm,1.0e-30_dbl_kind))

       ! densities
       rho_sum = rho_sum + rho(k)
       !rho_pipe = rho(k)
       rho_pipe = p5 * (rho(k) + rho(k-1))
       drho = max(rho(k) - rho_ocn, c0)

       ! mush Rayleigh number 
       Ra = drho * (hin-z) * perm_min * ra_constants

       ! height of mush layer to layer k
       rn = real(nilyr-k+1,dbl_kind)
       L = rn * hilyr

       ! horizontal size of convection
       dx = L * c2 * aspect_rapid_mode
       dx2 = dx**2

       ! determine vertical Darcy flow
       Am = (dx2 * rn) / (viscosity_dyn * perm_harm)
       Bm = (-gravit * rho_sum) / rn

       Ap = (pi * a_rapid_mode**4) / (c8 * viscosity_dyn)
       Bp = -rho_pipe * gravit

       q(k) = max((Am / dx2) * ((-Ap*Bp - Am*Bm) / (Am + Ap) + Bm), 1.0e-30_dbl_kind)

       ! modify by Rayleigh number and advection limit
       q(k) = min(q(k) * (max(Ra - Rac_rapid_mode, c0) / (Ra+puny)), qlimit)

       ! late stage drainage
       dSdt(k) = dSdt_slow_mode * (max((zSin(k) - phi_c_slow_mode*Sbr(k)), c0) &
                                *  max((Tbot - Tsf), c0)) / (hin + 0.001_dbl_kind)

       dSdt(k) = max(dSdt(k), (-zSin(k) * 0.5_dbl_kind) / dt)

       ! restrict flows to prevent too much salt loss
       dS_guess = (((q(k) * (Sbr(k+1) - Sbr(k))) / hilyr + dSdt(k)) * dt) * safety_factor

       if (abs(dS_guess) < puny) then
          alpha = c1
       else
          alpha = (zSin_min - zSin(k)) / dS_guess
       endif

       if (alpha < c0 .or. alpha > c1) alpha = c1

       q(k)    = q(k)    * alpha
       dSdt(k) = dSdt(k) * alpha

    enddo ! k

  end subroutine explicit_flow_velocities

!=======================================================================
! Flushing
!=======================================================================

  subroutine flushing_velocity(zTin,   zSin,  &
                               phi,    nilyr, &
                               hin,    hsn,   &
                               hilyr,         &
                               hpond,  apond, &
                               dt,     w)
   
    ! calculate the vertical flushing Darcy velocity (positive downward)

    use ice_mushy_physics, only: &
         density_brine, &
         liquidus_brine_salinity_mush

    use ice_colpkg_tracers, only: &
         tr_pond

    integer (kind=int_kind), intent(in) :: &
         nilyr         ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zTin      , & ! ice layer temperature (C)
         zSin      , & ! ice layer bulk salinity (ppt)
         phi           ! ice layer liquid fraction

    real(kind=dbl_kind), intent(in) :: &
         hilyr     , & ! ice layer thickness (m)
         hpond     , & ! melt pond thickness (m)
         apond     , & ! melt pond area (-)
         hsn       , & ! snow thickness (m)
         hin       , & ! ice thickness (m)
         dt            ! time step (s)

    real(kind=dbl_kind), intent(out) :: &
         w             ! vertical flushing Darcy flow rate (m s-1)

    real(kind=dbl_kind), parameter :: &
         advection_limit = 0.005_dbl_kind ! limit to fraction of brine in 
                                          ! any layer that can be advected 

    real(kind=dbl_kind) :: &
         perm       , & ! ice layer permeability (m2)
         ice_mass   , & ! mass of ice (kg m-2)
         perm_harm  , & ! harmonic mean of ice permeability (m2)
         hocn       , & ! ocean surface height above ice base (m)
         hbrine     , & ! brine surface height above ice base (m)
         w_down_max , & ! maximum downward flushing Darcy flow rate (m s-1) 
         phi_min    , & ! minimum porosity in the mush
         wlimit     , & ! limit to w to avoid advecting all brine in layer
         dhhead         ! hydraulic head (m)

    integer(kind=int_kind) :: &
         k              ! ice layer index

    ! initialize
    w = c0

    ! only flush if ponds are active
    if (tr_pond) then

       ice_mass  = c0
       perm_harm = c0
       phi_min   = c1

       do k = 1, nilyr

          ! liquid fraction
          !phi = liquid_fraction(zTin(k), zSin(k))
          phi_min = min(phi_min,phi(k))

          ! permeability
          perm = permeability(phi(k))

          ! ice mass
          ice_mass = ice_mass + phi(k)        * density_brine(liquidus_brine_salinity_mush(zTin(k))) + &
               (c1 - phi(k)) * rhoi

          ! permeability harmonic mean
          perm_harm = perm_harm + c1 / (perm + 1e-30_dbl_kind)

       enddo ! k

       ice_mass = ice_mass * hilyr

       perm_harm = real(nilyr,dbl_kind) / perm_harm 

       ! calculate ocean surface height above bottom of ice
       hocn = (ice_mass + hpond * apond * rhow + hsn * rhos) / rhow

       ! calculate brine height above bottom of ice
       hbrine = hin + hpond

       ! pressure head
       dhhead = max(hbrine - hocn,c0)

       ! darcy flow through ice
       w = (perm_harm * rhow * gravit * (dhhead / hin)) / viscosity_dyn

       ! maximum down flow to drain pond
       w_down_max = (hpond * apond) / dt

       ! limit flow
       w = min(w,w_down_max)

       ! limit amount of brine that can be advected out of any particular layer
       wlimit = (advection_limit * phi_min * hilyr) / dt

       if (abs(w) > puny) then
          w = w * max(min(abs(wlimit/w),c1),c0)
       else
          w = c0
       endif

       w = max(w, c0)

    endif

  end subroutine flushing_velocity

!=======================================================================

  subroutine flush_pond(w, hin, hpond, apond, dt)

    use ice_colpkg_tracers, only: &
         tr_pond

    ! given a flushing velocity drain the meltponds

    real(kind=dbl_kind), intent(in) :: &
         w     , & ! vertical flushing Darcy flow rate (m s-1)
         hin   , & ! ice thickness (m)
         apond , & ! melt pond area (-)
         dt        ! time step (s)

    real(kind=dbl_kind), intent(inout) :: &
         hpond     ! melt pond thickness (m)
    
    real(kind=dbl_kind), parameter :: &
         lambda_pond = c1 / (10.0_dbl_kind * 24.0_dbl_kind * 3600.0_dbl_kind), &
         hpond0 = 0.01_dbl_kind

    if (tr_pond) then
       if (apond > c0 .and. hpond > c0) then
          
          ! flush pond through mush
          hpond = hpond - w * dt / apond
          
          hpond = max(hpond, c0)
          
          ! exponential decay of pond
          hpond = hpond - lambda_pond * dt * (hpond + hpond0)
          
          hpond = max(hpond, c0)
          
       endif
    endif

  end subroutine flush_pond

 !=======================================================================

  subroutine flood_ice(hsn,    hin,      &
                       nslyr,  nilyr,    & 
                       hslyr,  hilyr,    & 
                       zqsn,   zqin,     &
                       phi,    dt,       &
                       zSin,   Sbr,      &
                       sss,    qocn,     &
                       snoice, fadvheat)

    ! given upwards flushing brine flow calculate amount of snow ice and
    ! convert snow to ice with appropriate properties

    use ice_mushy_physics, only: &
         density_brine

    integer (kind=int_kind), intent(in) :: &
         nilyr , & ! number of ice layers
         nslyr     ! number of snow layers

    real(kind=dbl_kind), intent(in) :: &
         dt                , & ! time step (s)
         hsn               , & ! snow thickness (m)
         hin               , & ! ice thickness (m)
         sss               , & ! sea surface salinity (ppt)
         qocn                  ! ocean brine enthalpy (J m-2)

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         zqsn              , & ! snow layer enthalpy (J m-2)
         zqin              , & ! ice layer enthalpy (J m-2)
         zSin              , & ! ice layer bulk salinity (ppt)
         phi                   ! ice liquid fraction

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         Sbr                   ! ice layer brine salinity (ppt)

    real(kind=dbl_kind), intent(inout) :: &
         hslyr             , & ! snow layer thickness (m)
         hilyr                 ! ice layer thickness (m)

    real(kind=dbl_kind), intent(out) :: &
         snoice                ! snow ice formation

   real(kind=dbl_kind), intent(inout) :: &
         fadvheat              ! advection heat flux to ocean

    real(kind=dbl_kind) :: &
         hin2              , & ! new ice thickness (m)
         hsn2              , & ! new snow thickness (m)
         hilyr2            , & ! new ice layer thickness (m)
         hslyr2            , & ! new snow layer thickness (m)
         dh                , & ! thickness of snowice formation (m)
         phi_snowice       , & ! liquid fraction of new snow ice
         rho_snowice       , & ! density of snowice (kg m-3)
         zSin_snowice      , & ! bulk salinity of new snowice (ppt)
         zqin_snowice      , & ! ice enthalpy of new snowice (J m-2)
         zqsn_snowice      , & ! snow enthalpy of snow thats becoming snowice (J m-2)
         freeboard_density , & ! negative of ice surface freeboard times the ocean density (kg m-2)
         ice_mass          , & ! mass of the ice (kg m-2)
         rho_ocn           , & ! density of the ocean (kg m-3)
         ice_density       , & ! density of ice layer (kg m-3)
         hadded            , & ! thickness rate of water used from ocean (m/s)
         wadded            , & ! mass rate of water used from ocean (kg/m^2/s)
         eadded            , & ! energy rate of water used from ocean (W/m^2) 
         sadded                ! salt rate of water used from ocean (kg/m^2/s)

    integer :: &
         k                     ! vertical index

    snoice = c0

    ! check we have snow
    if (hsn > puny) then
       
       rho_ocn = density_brine(sss)

       ! ice mass
       ice_mass = c0
       do k = 1, nilyr
          ice_density = min(phi(k) * density_brine(Sbr(k)) + (c1 - phi(k)) * rhoi,rho_ocn)
          ice_mass = ice_mass + ice_density
       enddo ! k
       ice_mass = ice_mass * hilyr

       ! negative freeboard times ocean density
       freeboard_density = max(ice_mass + hsn * rhos - hin * rho_ocn, c0)

       ! check if have flooded ice
       if (freeboard_density > c0) then

          ! sea ice fraction of newly formed snow ice
          phi_snowice = (c1 - rhos / rhoi)

          ! density of newly formed snowice
          rho_snowice = phi_snowice * rho_ocn + (c1 - phi_snowice) * rhoi

          ! calculate thickness of new ice added
          dh = freeboard_density / (rho_ocn - rho_snowice + rhos)
          dh = max(min(dh,hsn),c0)

          ! enthalpy of snow that becomes snowice
          call enthalpy_snow_snowice(nslyr, dh, hsn, zqsn, zqsn_snowice)

          ! change thicknesses
          hin2 = hin + dh
          hsn2 = hsn - dh

          hilyr2 = hin2 / real(nilyr,dbl_kind)
          hslyr2 = hsn2 / real(nslyr,dbl_kind)

          ! properties of new snow ice
          zSin_snowice = phi_snowice * sss
          zqin_snowice = phi_snowice * qocn + zqsn_snowice

          ! change snow properties
          call update_vertical_tracers_snow(nslyr, zqsn, hslyr, hslyr2)

          ! change ice properties
          call update_vertical_tracers_ice(nilyr, zqin, hilyr, hilyr2, &
               hin,  hin2,  zqin_snowice)
          call update_vertical_tracers_ice(nilyr, zSin, hilyr, hilyr2, &
               hin,  hin2,  zSin_snowice)
          call update_vertical_tracers_ice(nilyr, phi,  hilyr, hilyr2, &
               hin,  hin2,  phi_snowice)

          ! change thicknesses
          hilyr = hilyr2
          hslyr = hslyr2
          snoice = dh

          hadded = (dh * phi_snowice) / dt
          wadded = hadded * rhoi
          eadded = hadded * qocn
          sadded = wadded * ice_ref_salinity * p001

          ! conservation
          fadvheat = fadvheat - eadded

       endif

    endif

  end subroutine flood_ice

!=======================================================================

  subroutine enthalpy_snow_snowice(nslyr, dh, hsn, zqsn, zqsn_snowice)

    ! determine enthalpy of the snow being converted to snow ice
    
    integer (kind=int_kind), intent(in) :: &
         nslyr        ! number of snow layers

    real(kind=dbl_kind), intent(in) :: &
         dh       , & ! thickness of new snowice formation (m)
         hsn          ! initial snow thickness

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqsn         ! snow layer enthalpy (J m-2)

    real(kind=dbl_kind), intent(out) :: &
         zqsn_snowice ! enthalpy of snow becoming snowice (J m-2)

    real(kind=dbl_kind) :: &
         rnlyr        ! real value of number of snow layers turning to snowice

    integer(kind=int_kind) :: &
         nlyr     , & ! no of snow layers completely converted to snowice
         k            ! snow layer index

    zqsn_snowice = c0

    ! snow depth and snow layers affected by snowice formation
    if (hsn > puny) then
       rnlyr = (dh / hsn) * nslyr
       nlyr = min(floor(rnlyr),nslyr-1) ! nlyr=0 if nslyr=1

       ! loop over full snow layers affected
       ! not executed if nlyr=0
       do k = nslyr, nslyr-nlyr+1, -1
          zqsn_snowice = zqsn_snowice + zqsn(k) / rnlyr
       enddo ! k

       ! partially converted snow layer
       zqsn_snowice = zqsn_snowice + &
            ((rnlyr - real(nlyr,dbl_kind)) / rnlyr) * zqsn(nslyr-nlyr)
    endif

  end subroutine enthalpy_snow_snowice

!=======================================================================
  
  subroutine update_vertical_tracers_snow(nslyr, trc, hlyr1, hlyr2)

    ! given some snow ice formation regrid snow layers

    integer (kind=int_kind), intent(in) :: &
         nslyr       ! number of snow layers
    
    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         trc         ! vertical tracer

    real(kind=dbl_kind), intent(in) :: &
         hlyr1   , & ! old cell thickness
         hlyr2       ! new cell thickness

    real(kind=dbl_kind), dimension(1:nslyr) :: &
         trc2        ! temporary array for updated tracer
    
    ! vertical indexes for old and new grid
    integer(kind=int_kind) :: &
         k1      , & ! vertical index for old grid
         k2          ! vertical index for new grid
    
    real(kind=dbl_kind) :: &
         z1a     , & ! lower boundary of old cell
         z1b     , & ! upper boundary of old cell
         z2a     , & ! lower boundary of new cell
         z2b     , & ! upper boundary of new cell
         overlap     ! overlap between old and new cell
    
    ! loop over new grid cells
    do k2 = 1, nslyr
       
       ! initialize new tracer
       trc2(k2) = c0
       
       ! calculate upper and lower boundary of new cell
       z2a = (k2 - 1) * hlyr2
       z2b = k2       * hlyr2
       
       ! loop over old grid cells
       do k1 = 1, nslyr
          
          ! calculate upper and lower boundary of old cell
          z1a = (k1 - 1) * hlyr1
          z1b = k1       * hlyr1
          
          ! calculate overlap between old and new cell
          overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
          
          ! aggregate old grid cell contribution to new cell
          trc2(k2) = trc2(k2) + overlap * trc(k1)
          
       enddo ! k1

       ! renormalize new grid cell
       trc2(k2) = trc2(k2) / hlyr2
       
    enddo ! k2
    
    ! update vertical tracer array with the adjusted tracer
    trc = trc2

  end subroutine update_vertical_tracers_snow

!=======================================================================
  
  subroutine update_vertical_tracers_ice(nilyr, trc, hlyr1, hlyr2, &
                                         h1, h2, trc0)

    ! given some snow ice formation regrid ice layers

    integer (kind=int_kind), intent(in) :: &
         nilyr       ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(inout) :: &
         trc         ! vertical tracer
    
    real(kind=dbl_kind), intent(in) :: &
         hlyr1 , &   ! old cell thickness
         hlyr2 , &   ! new cell thickness
         h1    , &   ! old total thickness
         h2    , &   ! new total thickness
         trc0        ! tracer value of added snow ice on ice top
    
    real(kind=dbl_kind), dimension(1:nilyr) :: &
         trc2        ! temporary array for updated tracer
    
    integer(kind=int_kind) :: &
         k1 , &      ! vertical indexes for old grid
         k2          ! vertical indexes for new grid
    
    real(kind=dbl_kind) :: &
         z1a     , & ! lower boundary of old cell
         z1b     , & ! upper boundary of old cell
         z2a     , & ! lower boundary of new cell
         z2b     , & ! upper boundary of new cell
         overlap     ! overlap between old and new cell
    
    ! loop over new grid cells
    do k2 = 1, nilyr
       
       ! initialize new tracer
       trc2(k2) = c0
       
       ! calculate upper and lower boundary of new cell
       z2a = (k2 - 1) * hlyr2
       z2b = k2       * hlyr2

       ! calculate upper and lower boundary of added snow ice at top
       z1a = c0
       z1b = h2 - h1
       
       ! calculate overlap between added ice and new cell
       overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
       
       ! aggregate added ice contribution to new cell
       trc2(k2) = trc2(k2) + overlap * trc0

       ! loop over old grid cells
       do k1 = 1, nilyr
          
          ! calculate upper and lower boundary of old cell
          z1a = (k1 - 1) * hlyr1 + h2 - h1
          z1b = k1       * hlyr1 + h2 - h1
          
          ! calculate overlap between old and new cell
          overlap = max(min(z1b, z2b) - max(z1a, z2a), c0)
          
          ! aggregate old grid cell contribution to new cell
          trc2(k2) = trc2(k2) + overlap * trc(k1)
          
       enddo ! k1

       ! renormalize new grid cell
       trc2(k2) = trc2(k2) / hlyr2
       
    enddo ! k2
    
    ! update vertical tracer array with the adjusted tracer
    trc = trc2

  end subroutine update_vertical_tracers_ice

!=======================================================================

end module ice_therm_mushy

!=======================================================================
