module BiogeophysicsLakeMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BiogeophysicsLakeMod
!
! !DESCRIPTION:
! Calculates lake temperatures and surface fluxes.
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysicsLake

! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BiogeophysicsLake
!
! !INTERFACE:
  subroutine BiogeophysicsLake(lbc, ubc, lbp, ubp, num_lakec, filter_lakec, &
                               num_lakep, filter_lakep)
!
! !DESCRIPTION:
! Calculates lake temperatures and surface fluxes.
! Lake temperatures are determined from a one-dimensional thermal
! stratification model based on eddy diffusion concepts to
! represent vertical mixing of heat.
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz
!
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
! There are two types of lakes:
!   Deep lakes are 50 m.
!   Shallow lakes are 10 m deep.
!
!   For unfrozen deep lakes:    ke > 0 and    convective mixing
!   For unfrozen shallow lakes: ke = 0 and no convective mixing
!
! Use the Crank-Nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
!
! The solution conserves energy as:
!
! cw*([ts(      1)] n+1 - [ts(      1)] n)*dz(      1)/dt + ... +
! cw*([ts(nlevlak)] n+1 - [ts(nlevlak)] n)*dz(nlevlak)/dt = fin
!
! where:
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot
!            - hm + phi(1) + ... + phi(nlevlak)
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_time_manager       , only : get_step_size
    use clm_varpar         , only : nlevlak
    use clm_varcon         , only : hvap, hsub, hfus, cpair, cpliq, cpice, tkwat, tkice, &
                                    sb, vkc, grav, denh2o, tfrz, spval
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use TridiagonalMod     , only : Tridiagonal
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                ! column-index bounds
    integer, intent(in) :: lbp, ubp                ! pft-index bounds
    integer, intent(in) :: num_lakec               ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1) ! column filter for non-lake points
    integer, intent(in) :: num_lakep               ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_lakep(ubp-lbp+1) ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! Migrated to clm2.1 new data structures by Peter Thornton and M. Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: pcolumn(:)         ! pft's column index
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: cgridcell(:)       ! column's gridcell index
    real(r8), pointer :: forc_t(:)          ! atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_pbot(:)       ! atmospheric pressure (Pa)
    real(r8), pointer :: forc_hgt_u_pft(:)  ! observational height of wind at pft level [m]
    real(r8), pointer :: forc_hgt_t_pft(:)  ! observational height of temperature at pft level [m]
    real(r8), pointer :: forc_hgt_q_pft(:)  ! observational height of specific humidity at pft level [m]
    real(r8), pointer :: forc_th(:)         ! atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_q(:)          ! atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_u(:)          ! atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)          ! atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: forc_rho(:)        ! density (kg/m**3)
    real(r8), pointer :: forc_snow(:)       ! snow rate [mm/s]
    real(r8), pointer :: forc_rain(:)       ! rain rate [mm/s]
    real(r8), pointer :: t_grnd(:)          ! ground temperature (Kelvin)
    real(r8), pointer :: hc_soisno(:)       ! soil plus snow plus lake heat content (MJ/m2)
    real(r8), pointer :: h2osno(:)          ! snow water (mm H2O)
    real(r8), pointer :: snowdp(:)          ! snow height (m)
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: lat(:)             ! latitude (radians)
    real(r8), pointer :: dz(:,:)            ! layer thickness (m)
    real(r8), pointer :: z(:,:)             ! layer depth (m)
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: qflx_prec_grnd(:)  ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: qflx_snwcp_liq(:)  ! excess rainfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: qflx_snwcp_ice(:)  ! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_lh_grnd(:)    ! ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: t_veg(:)           ! vegetation temperature (Kelvin)
    real(r8), pointer :: t_ref2m(:)         ! 2 m height surface air temperature (Kelvin)
    real(r8), pointer :: q_ref2m(:)         ! 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: rh_ref2m(:)        ! 2 m height surface relative humidity (%)
    real(r8), pointer :: taux(:)            ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8), pointer :: tauy(:)            ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8), pointer :: qmelt(:)           ! snow melt [mm/s]
    real(r8), pointer :: ram1(:)            ! aerodynamical resistance (s/m)
    real(r8), pointer :: errsoi(:)          ! soil/lake energy conservation error (W/m**2)
    real(r8), pointer :: t_lake(:,:)        ! lake temperature (Kelvin)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer , parameter  :: idlak = 1     ! index of lake, 1 = deep lake, 2 = shallow lake
    integer , parameter  :: niters = 3    ! maximum number of iterations for surface temperature
    real(r8), parameter :: beta1 = 1._r8  ! coefficient of connective velocity (in computing W_*) [-]
    real(r8), parameter :: emg = 0.97_r8     ! ground emissivity (0.97 for snow)
    real(r8), parameter :: zii = 1000._r8 ! convective boundary height [m]
    real(r8), parameter :: p0 = 1._r8     ! neutral value of turbulent prandtl number
    integer  :: i,j,fc,fp,g,c,p         ! do loop or array index
    integer  :: fncopy                  ! number of values in pft filter copy
    integer  :: fnold                   ! previous number of pft filter values
    integer  :: fpcopy(num_lakep)       ! pft filter copy for iteration loop
    integer  :: num_unfrzc              ! number of values in unfrozen column filter
    integer  :: filter_unfrzc(ubc-lbc+1)! unfrozen column filter
    integer  :: iter                    ! iteration index
    integer  :: nmozsgn(lbp:ubp)        ! number of times moz changes sign
    integer  :: jtop(lbc:ubc)           ! number of levels for each column (all 1)
    real(r8) :: dtime                   ! land model time step (sec)
    real(r8) :: ax                      !
    real(r8) :: bx                      !
    real(r8) :: degdT                   ! d(eg)/dT
    real(r8) :: dqh(lbp:ubp)            ! diff of humidity between ref. height and surface
    real(r8) :: dth(lbp:ubp)            ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                    ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dzsur(lbc:ubc)          !
    real(r8) :: eg                      ! water vapor pressure at temperature T [pa]
    real(r8) :: hm                      ! energy residual [W/m2]
    real(r8) :: htvp(lbc:ubc)           ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8) :: obu(lbp:ubp)            ! monin-obukhov length (m)
    real(r8) :: obuold(lbp:ubp)         ! monin-obukhov length of previous iteration
    real(r8) :: qsatg(lbc:ubc)          ! saturated humidity [kg/kg]
    real(r8) :: qsatgdT(lbc:ubc)        ! d(qsatg)/dT
    real(r8) :: qstar                   ! moisture scaling parameter
    real(r8) :: ram(lbp:ubp)            ! aerodynamical resistance [s/m]
    real(r8) :: rah(lbp:ubp)            ! thermal resistance [s/m]
    real(r8) :: raw(lbp:ubp)            ! moisture resistance [s/m]
    real(r8) :: stftg3(lbp:ubp)         ! derivative of fluxes w.r.t ground temperature
    real(r8) :: temp1(lbp:ubp)          ! relation for potential temperature profile
    real(r8) :: temp12m(lbp:ubp)        ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(lbp:ubp)          ! relation for specific humidity profile
    real(r8) :: temp22m(lbp:ubp)        ! relation for specific humidity profile applied at 2-m
    real(r8) :: tgbef(lbc:ubc)          ! initial ground temperature
    real(r8) :: thm(lbp:ubp)            ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(r8) :: thv(lbc:ubc)            ! virtual potential temperature (kelvin)
    real(r8) :: thvstar                 ! virtual potential temperature scaling parameter
    real(r8) :: tksur                   ! thermal conductivity of snow/soil (w/m/kelvin)
    real(r8) :: tstar                   ! temperature scaling parameter
    real(r8) :: um(lbp:ubp)             ! wind speed including the stablity effect [m/s]
    real(r8) :: ur(lbp:ubp)             ! wind speed at reference height [m/s]
    real(r8) :: ustar(lbp:ubp)          ! friction velocity [m/s]
    real(r8) :: wc                      ! convective velocity [m/s]
    real(r8) :: zeta                    ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: zldis(lbp:ubp)          ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(lbp:ubp)         ! displacement (always zero) [m]
    real(r8) :: z0mg(lbp:ubp)           ! roughness length over ground, momentum [m]
    real(r8) :: z0hg(lbp:ubp)           ! roughness length over ground, sensible heat [m]
    real(r8) :: z0qg(lbp:ubp)           ! roughness length over ground, latent heat [m]
    real(r8) :: beta(2)                 ! fraction solar rad absorbed at surface: depends on lake type
    real(r8) :: za(2)                   ! base of surface absorption layer (m): depends on lake type
    real(r8) :: eta(2)                  ! light extinction coefficient (/m): depends on lake type
    real(r8) :: a(lbc:ubc,nlevlak)      ! "a" vector for tridiagonal matrix
    real(r8) :: b(lbc:ubc,nlevlak)      ! "b" vector for tridiagonal matrix
    real(r8) :: c1(lbc:ubc,nlevlak)     ! "c" vector for tridiagonal matrix
    real(r8) :: r(lbc:ubc,nlevlak)      ! "r" vector for tridiagonal solution
    real(r8) :: rhow(lbc:ubc,nlevlak)   ! density of water (kg/m**3)
    real(r8) :: phi(lbc:ubc,nlevlak)    ! solar radiation absorbed by layer (w/m**2)
    real(r8) :: kme(lbc:ubc,nlevlak)    ! molecular + eddy diffusion coefficient (m**2/s)
    real(r8) :: cwat                    ! specific heat capacity of water (j/m**3/kelvin)
    real(r8) :: ws(lbc:ubc)             ! surface friction velocity (m/s)
    real(r8) :: ks(lbc:ubc)             ! coefficient
    real(r8) :: in                      ! relative flux of solar radiation into layer
    real(r8) :: out                     ! relative flux of solar radiation out of layer
    real(r8) :: ri                      ! richardson number
    real(r8) :: fin(lbc:ubc)            ! heat flux into lake - flux out of lake (w/m**2)
    real(r8) :: ocvts(lbc:ubc)          ! (cwat*(t_lake[n  ])*dz
    real(r8) :: ncvts(lbc:ubc)          ! (cwat*(t_lake[n+1])*dz
    real(r8) :: m1                      ! intermediate variable for calculating r, a, b, c
    real(r8) :: m2                      ! intermediate variable for calculating r, a, b, c
    real(r8) :: m3                      ! intermediate variable for calculating r, a, b, c
    real(r8) :: ke                      ! eddy diffusion coefficient (m**2/s)
    real(r8) :: km                      ! molecular diffusion coefficient (m**2/s)
    real(r8) :: zin                     ! depth at top of layer (m)
    real(r8) :: zout                    ! depth at bottom of layer (m)
    real(r8) :: drhodz                  ! d [rhow] /dz (kg/m**4)
    real(r8) :: n2                      ! brunt-vaisala frequency (/s**2)
    real(r8) :: num                     ! used in calculating ri
    real(r8) :: den                     ! used in calculating ri
    real(r8) :: tav(lbc:ubc)            ! used in aver temp for convectively mixed layers
    real(r8) :: nav(lbc:ubc)            ! used in aver temp for convectively mixed layers
    real(r8) :: phidum                  ! temporary value of phi
    real(r8) :: u2m                     ! 2 m wind speed (m/s)
    real(r8) :: fm(lbp:ubp)             ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: e_ref2m                 ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                  ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m              ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT               ! derivative of 2 m height surface saturated specific humidity on t_ref2m
!
! Constants for lake temperature model
!
    data beta/0.4_r8, 0.4_r8/  ! (deep lake, shallow lake)
    data za  /0.6_r8, 0.5_r8/
    data eta /0.1_r8, 0.5_r8/
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    forc_t         => clm_a2l%forc_t
    forc_pbot      => clm_a2l%forc_pbot
    forc_th        => clm_a2l%forc_th
    forc_q         => clm_a2l%forc_q
    forc_u         => clm_a2l%forc_u
    forc_v         => clm_a2l%forc_v
    forc_rho       => clm_a2l%forc_rho
    forc_lwrad     => clm_a2l%forc_lwrad
    forc_snow      => clm_a2l%forc_snow
    forc_rain      => clm_a2l%forc_rain
    lat            => grc%lat

    ! Assign local pointers to derived type members (column-level)

    cgridcell      => col%gridcell
    dz             => cps%dz
    z              => cps%z
    t_lake         => ces%t_lake
    h2osno         => cws%h2osno
    snowdp         => cps%snowdp
    t_grnd         => ces%t_grnd
    hc_soisno      => ces%hc_soisno
    errsoi         => cebal%errsoi
    qmelt          => cwf%qmelt

    ! Assign local pointers to derived type members (pft-level)

    pcolumn        => pft%column
    pgridcell      => pft%gridcell
    sabg           => pef%sabg
    t_ref2m        => pes%t_ref2m
    q_ref2m        => pes%q_ref2m
    rh_ref2m       => pes%rh_ref2m
    t_veg          => pes%t_veg
    eflx_lwrad_out => pef%eflx_lwrad_out
    eflx_lwrad_net => pef%eflx_lwrad_net
    eflx_soil_grnd => pef%eflx_soil_grnd
    eflx_lh_tot    => pef%eflx_lh_tot
    eflx_lh_grnd   => pef%eflx_lh_grnd
    eflx_sh_grnd   => pef%eflx_sh_grnd
    eflx_sh_tot    => pef%eflx_sh_tot
    ram1           => pps%ram1
    taux           => pmf%taux
    tauy           => pmf%tauy
    qflx_prec_grnd => pwf%qflx_prec_grnd
    qflx_evap_soi  => pwf%qflx_evap_soi
    qflx_evap_tot  => pwf%qflx_evap_tot
    forc_hgt_u_pft => pps%forc_hgt_u_pft
    forc_hgt_t_pft => pps%forc_hgt_t_pft
    forc_hgt_q_pft => pps%forc_hgt_q_pft
    qflx_snwcp_ice => pwf%qflx_snwcp_ice    
    qflx_snwcp_liq => pwf%qflx_snwcp_liq    

    ! Determine step size

    dtime = get_step_size()

    ! Begin calculations

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       g = cgridcell(c)

       ! Initialize quantities computed below

       ocvts(c) = 0._r8
       ncvts(c) = 0._r8
       hc_soisno(c) = 0._r8

       ! Surface temperature and fluxes

       dzsur(c) = dz(c,1) + snowdp(c)

       ! Saturated vapor pressure, specific humidity and their derivatives
       ! at lake surface

       call QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height

       !zii = 1000.    ! m  (pbl height)
       thv(c) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))     ! virtual potential T
    end do

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       nmozsgn(p) = 0
       obuold(p) = 0._r8
       displa(p) = 0._r8
       thm(p) = forc_t(g) + 0.0098_r8*forc_hgt_t_pft(p)   ! intermediate variable

       ! Roughness lengths

       if (t_grnd(c) >= tfrz) then   ! for unfrozen lake
          z0mg(p) = 0.01_r8
       else                          ! for frozen lake
          z0mg(p) = 0.04_r8
       end if
       z0hg(p) = z0mg(p)
       z0qg(p) = z0mg(p)

       ! Latent heat

       if (forc_t(g) > tfrz) then
          htvp(c) = hvap
       else
          htvp(c) = hsub
       end if

       ! Initialize stability variables

       ur(p)    = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p)   = thm(p)-t_grnd(c)
       dqh(p)   = forc_q(g)-qsatg(c)
       dthv     = dth(p)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p) - 0._r8

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg(p), um(p), obu(p))

    end do

    iter = 1
    fncopy = num_lakep
    fpcopy(1:num_lakep) = filter_lakep(1:num_lakep)

    ! Begin stability iteration

    ITERATION : do while (iter <= niters .and. fncopy > 0)

       ! Determine friction velocity, and potential temperature and humidity
       ! profiles of the surface boundary layer

       call FrictionVelocity(lbp, ubp, fncopy, fpcopy, &
                             displa, z0mg, z0hg, z0qg, &
                             obu, iter, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, fm)

       do fp = 1, fncopy
          p = fpcopy(fp)
          c = pcolumn(p)
          g = pgridcell(p)

          tgbef(c) = t_grnd(c)
          if (t_grnd(c) > tfrz) then
             tksur = tkwat
          else
             tksur = tkice
          end if

          ! Determine aerodynamic resistances

          ram(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
          rah(p)  = 1._r8/(temp1(p)*ustar(p))
          raw(p)  = 1._r8/(temp2(p)*ustar(p))
          ram1(p) = ram(p)   !pass value to global variable

          ! Get derivative of fluxes with respect to ground temperature

          stftg3(p) = emg*sb*tgbef(c)*tgbef(c)*tgbef(c)

          ax  = sabg(p) + emg*forc_lwrad(g) + 3._r8*stftg3(p)*tgbef(c) &
               + forc_rho(g)*cpair/rah(p)*thm(p) &
               - htvp(c)*forc_rho(g)/raw(p)*(qsatg(c)-qsatgdT(c)*tgbef(c) - forc_q(g)) &
               + tksur*t_lake(c,1)/dzsur(c)

          bx  = 4._r8*stftg3(p) + forc_rho(g)*cpair/rah(p) &
               + htvp(c)*forc_rho(g)/raw(p)*qsatgdT(c) + tksur/dzsur(c)

          t_grnd(c) = ax/bx

          ! Surface fluxes of momentum, sensible and latent heat
          ! using ground temperatures from previous time step

          eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(p))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c))-forc_q(g))/raw(p)

          ! Re-calculate saturated vapor pressure, specific humidity and their
          ! derivatives at lake surface

          call QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

          dth(p)=thm(p)-t_grnd(c)
          dqh(p)=forc_q(g)-qsatg(c)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)

          !not used
          !dthv=dth(p)*(1.+0.61*forc_q(g))+0.61*forc_th(g)*dqh(p)
          thvstar=tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
          zeta=zldis(p)*vkc * grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0._r8) then     !stable
             zeta = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                     !unstable
             zeta = max(-100._r8,min(zeta,-0.01_r8))
             wc = beta1*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p)+wc*wc)
          end if
          obu(p) = zldis(p)/zeta

          if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1

          obuold(p) = obu(p)

       end do   ! end of filtered pft loop

       iter = iter + 1
       if (iter <= niters ) then
          ! Rebuild copy of pft filter for next pass through the ITERATION loop

          fnold = fncopy
          fncopy = 0
          do fp = 1, fnold
             p = fpcopy(fp)
             if (nmozsgn(p) < 3) then
                fncopy = fncopy + 1
                fpcopy(fncopy) = p
             end if
          end do   ! end of filtered pft loop
       end if

    end do ITERATION   ! end of stability iteration

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       ! initialize snow cap terms to zero for lake columns
       qflx_snwcp_ice(p) = 0._r8
       qflx_snwcp_liq(p) = 0._r8

       ! If there is snow on the ground and t_grnd > tfrz: reset t_grnd = tfrz.
       ! Re-evaluate ground fluxes. Energy imbalance used to melt snow.
       ! h2osno > 0.5 prevents spurious fluxes.
       ! note that qsatg and qsatgdT should be f(tgbef) (PET: not sure what this
       ! comment means)

       if (h2osno(c) > 0.5_r8 .AND. t_grnd(c) > tfrz) then
          t_grnd(c) = tfrz
          eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(p))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c)) - forc_q(g))/raw(p)
       end if

       ! Net longwave from ground to atmosphere

       eflx_lwrad_out(p) = (1._r8-emg)*forc_lwrad(g) + stftg3(p)*(-3._r8*tgbef(c)+4._r8*t_grnd(c))

       ! Ground heat flux

       eflx_soil_grnd(p) = sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) - &
            eflx_sh_grnd(p) - htvp(c)*qflx_evap_soi(p)

       taux(p) = -forc_rho(g)*forc_u(g)/ram(p)
       tauy(p) = -forc_rho(g)*forc_v(g)/ram(p)

       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)
       eflx_lh_tot(p)   = htvp(c)*qflx_evap_soi(p)
       eflx_lh_grnd(p)  = htvp(c)*qflx_evap_soi(p)

       ! 2 m height air temperature
       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity
       q_ref2m(p) = forc_q(g) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

       ! 2 m height relative humidity

       call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

       ! Energy residual used for melting snow
       if (h2osno(c) > 0._r8 .AND. t_grnd(c) >= tfrz) then
          hm = min(h2osno(c)*hfus/dtime, max(eflx_soil_grnd(p),0._r8))
       else
          hm = 0._r8
       end if
       qmelt(c) = hm/hfus             ! snow melt (mm/s)

       ! Prepare for lake layer temperature calculations below

       fin(c) = beta(idlak) * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
            eflx_sh_tot(p) + eflx_lh_tot(p) + hm)
       u2m = max(1.0_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))

       ws(c) = 1.2e-03_r8 * u2m
       ks(c) = 6.6_r8*sqrt(abs(sin(lat(g))))*(u2m**(-1.84_r8))

    end do

    ! Eddy diffusion +  molecular diffusion coefficient (constants):
    ! eddy diffusion coefficient used for unfrozen deep lakes only

    cwat = cpliq*denh2o ! a constant
    km = tkwat/cwat     ! a constant

    ! Lake density

    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          rhow(c,j) = 1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-277._r8))**1.68_r8 )
       end do
    end do

    do j = 1, nlevlak-1
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          drhodz = (rhow(c,j+1)-rhow(c,j)) / (z(c,j+1)-z(c,j))
          n2 = -grav / rhow(c,j) * drhodz
          num = 40._r8 * n2 * (vkc*z(c,j))**2
          den = max( (ws(c)**2) * exp(-2._r8*ks(c)*z(c,j)), 1.e-10_r8 )
          ri = ( -1._r8 + sqrt( max(1._r8+num/den, 0._r8) ) ) / 20._r8
          if (t_grnd(c) > tfrz) then
             ! valid for deep lake only (idlak == 1)
             ke = vkc*ws(c)*z(c,j)/p0 * exp(-ks(c)*z(c,j)) / (1._r8+37._r8*ri*ri)
          else
             ke = 0._r8
          end if
          kme(c,j) = km + ke
       end do
    end do

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       kme(c,nlevlak) = kme(c,nlevlak-1)
       ! set number of column levels for use by Tridiagonal below
       jtop(c) = 1
    end do

    ! Heat source term: unfrozen lakes only

    do j = 1, nlevlak
       do fp = 1, num_lakep
          p = filter_lakep(fp)
          c = pcolumn(p)

          zin  = z(c,j) - 0.5_r8*dz(c,j)
          zout = z(c,j) + 0.5_r8*dz(c,j)
          in  = exp( -eta(idlak)*max(  zin-za(idlak),0._r8 ) )
          out = exp( -eta(idlak)*max( zout-za(idlak),0._r8 ) )

          ! Assume solar absorption is only in the considered depth
          if (j == nlevlak) out = 0._r8
          if (t_grnd(c) > tfrz) then
             phidum = (in-out) * sabg(p) * (1._r8-beta(idlak))
          else if (j == 1) then
             phidum = sabg(p) * (1._r8-beta(idlak))
          else
             phidum = 0._r8
          end if
          phi(c,j) = phidum
       end do
    end do

    ! Sum cwat*t_lake*dz for energy check

    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          ocvts(c) = ocvts(c) + cwat*t_lake(c,j)*dz(c,j)
       end do
    end do

    ! Set up vector r and vectors a, b, c that define tridiagonal matrix

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       j = 1
       m2 = dz(c,j)/kme(c,j) + dz(c,j+1)/kme(c,j+1)
       m3 = dtime/dz(c,j)
       r(c,j) = t_lake(c,j) + (fin(c)+phi(c,j))*m3/cwat - (t_lake(c,j)-t_lake(c,j+1))*m3/m2
       a(c,j) = 0._r8
       b(c,j) = 1._r8 + m3/m2
       c1(c,j) = -m3/m2

       j = nlevlak
       m1 = dz(c,j-1)/kme(c,j-1) + dz(c,j)/kme(c,j)
       m3 = dtime/dz(c,j)
       r(c,j) = t_lake(c,j) + phi(c,j)*m3/cwat + (t_lake(c,j-1)-t_lake(c,j))*m3/m1
       a(c,j) = -m3/m1
       b(c,j) = 1._r8 + m3/m1
       c1(c,j) = 0._r8
    end do

    do j = 2, nlevlak-1
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          m1 = dz(c,j-1)/kme(c,j-1) + dz(c,j  )/kme(c,j  )
          m2 = dz(c,j  )/kme(c,j  ) + dz(c,j+1)/kme(c,j+1)
          m3 = dtime/dz(c,j)
          r(c,j) = t_lake(c,j) + phi(c,j)*m3/cwat + &
             (t_lake(c,j-1) - t_lake(c,j  ))*m3/m1 - &
             (t_lake(c,j  ) - t_lake(c,j+1))*m3/m2

          a(c,j) = -m3/m1
          b(c,j) = 1._r8 + m3/m1 + m3/m2
          c1(c,j) = -m3/m2
       end do
    end do

    ! Solve for t_lake: a, b, c, r, u

    call Tridiagonal(lbc, ubc, 1, nlevlak, jtop, num_lakec, filter_lakec, &
                     a, b, c1, r, t_lake(lbc:ubc,1:nlevlak))

    ! Convective mixing: make sure cwat*dz*ts is conserved.  Valid only for
    ! deep lakes (idlak == 1).

    num_unfrzc = 0
    do fc = 1, num_lakec
       c = filter_lakec(fc)
       if (t_grnd(c) > tfrz) then
          num_unfrzc = num_unfrzc + 1
          filter_unfrzc(num_unfrzc) = c
       end if
    end do

    do j = 1, nlevlak-1
       do fc = 1, num_unfrzc
          c = filter_unfrzc(fc)
          tav(c) = 0._r8
          nav(c) = 0._r8
       end do

       do i = 1, j+1
          do fc = 1, num_unfrzc
             c = filter_unfrzc(fc)
             if (rhow(c,j) > rhow(c,j+1)) then
                tav(c) = tav(c) + t_lake(c,i)*dz(c,i)
                nav(c) = nav(c) + dz(c,i)
             end if
          end do
       end do

       do fc = 1, num_unfrzc
          c = filter_unfrzc(fc)
          if (rhow(c,j) > rhow(c,j+1)) then
             tav(c) = tav(c)/nav(c)
          end if
       end do

       do i = 1, j+1
          do fc = 1, num_unfrzc
             c = filter_unfrzc(fc)
             if (nav(c) > 0._r8) then
                t_lake(c,i) = tav(c)
                rhow(c,i) = 1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,i)-277._r8))**1.68_r8 )
             end if
          end do
       end do
    end do

    ! Sum cwat*t_lake*dz and total energy into lake for energy check

    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          ncvts(c) = ncvts(c) + cwat*t_lake(c,j)*dz(c,j)
          hc_soisno(c) = hc_soisno(c) + cwat*t_lake(c,j)*dz(c,j) /1.e6_r8
          if (j == nlevlak) then 
             hc_soisno(c) = hc_soisno(c) +  &
                            cpice*h2osno(c)*t_grnd(c)*snowdp(c) /1.e6_r8
          endif
          fin(c) = fin(c) + phi(c,j)
       end do
    end do

    ! The following are needed for global average on history tape.

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)
       errsoi(c) = (ncvts(c)-ocvts(c)) / dtime - fin(c)
       t_veg(p) = forc_t(g)
       eflx_lwrad_net(p)  = eflx_lwrad_out(p) - forc_lwrad(g)
       qflx_prec_grnd(p) = forc_rain(g) + forc_snow(g)
    end do

  end subroutine BiogeophysicsLake

end module BiogeophysicsLakeMod
