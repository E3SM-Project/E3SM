module SLakeFluxesMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SLakeFluxesMod
!
! !DESCRIPTION:
! Calculates surface fluxes for lakes.
!
! !PUBLIC TYPES:
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeFluxes

! !REVISION HISTORY:
! Created by Zack Subin, 2009
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SLakeFluxes
!
! !INTERFACE:
  subroutine SLakeFluxes(lbc, ubc, lbp, ubp, num_lakec, filter_lakec, &
                               num_lakep, filter_lakep)
!
! !DESCRIPTION:
! Calculates lake temperatures and surface fluxes.
!
! Lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
! and soil layers with active temperature and gas diffusion below.
!
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_time_manager       , only : get_step_size
    use clm_varpar         , only : nlevlak
    use clm_varcon         , only : hvap, hsub, hfus, cpair, cpliq, tkwat, tkice, tkair, &
                                    sb, vkc, grav, denh2o, tfrz, spval, zsno
    use SLakeCon           , only : betavis, z0frzlake, tdmax, emg_lake, &
                                    minz0lake, cur0, cus, curm, fcrit
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use SLakeCon           , only : lake_use_old_fcrit_minz0
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
! subroutine driver
!
! !REVISION HISTORY:
! Author: Zack Subin, 2009
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: sabg_chk(:)        ! sum of soil/snow using current fsno, for balance check
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
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: lat(:)             ! latitude (radians)
    real(r8), pointer :: dz(:,:)            ! layer thickness for soil or snow (m)
    real(r8), pointer :: dz_lake(:,:)       ! layer thickness for lake (m)
    real(r8), pointer :: t_soisno(:,:)      ! soil (or snow) temperature (Kelvin)
    real(r8), pointer :: t_lake(:,:)        ! lake temperature (Kelvin)
    integer , pointer :: snl(:)             ! number of snow layers
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)    ! ice lens (kg/m2)
    real(r8), pointer :: savedtke1(:)       ! top level eddy conductivity from previous timestep (W/mK)
    real(r8), pointer :: lakedepth(:)       ! variable lake depth (m)
    real(r8), pointer :: lakefetch(:)       ! lake fetch from surface data (m)
! variables needed for SNICAR
    real(r8), pointer :: sabg_lyr(:,:)      ! absorbed solar radiation (pft,lyr) [W/m2]
!   Calculation of beta depending on NIR fraction of sabg
    real(r8), pointer :: fsds_nir_d(:)    ! incident direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsds_nir_i(:)    ! incident diffuse nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_d(:)     ! reflected direct beam nir solar radiation (W/m**2)
    real(r8), pointer :: fsr_nir_i(:)     ! reflected diffuse nir solar radiation (W/m**2)

!
! local pointers to implicit out arguments
!
    real(r8), pointer :: qflx_prec_grnd(:)  ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
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
    real(r8), pointer :: ram1(:)            ! aerodynamical resistance (s/m)
    real(r8), pointer :: ram1_lake(:)       ! aerodynamical resistance (s/m)
    real(r8), pointer :: ws(:)              ! surface friction velocity (m/s)
    real(r8), pointer :: ks(:)              ! coefficient passed to SLakeTemperature
                                            !     for calculation of decay of eddy diffusivity with depth
    real(r8), pointer :: eflx_gnet(:)       ! net heat flux into ground (W/m**2)
    real(r8), pointer :: ust_lake(:)        ! friction velocity (m/s)
    real(r8), pointer :: z0mg_col(:)        ! roughness length over ground, momentum [m]
    real(r8), pointer :: z0hg_col(:)        ! roughness length over ground, sensible heat [m]
    real(r8), pointer :: z0qg_col(:)        ! roughness length over ground, latent heat [m]
    real(r8), pointer :: qflx_snwcp_ice(:)  ! excess snowfall due to snow capping (mm H2O /s) [+]
    real(r8), pointer :: qflx_snwcp_liq(:)  ! excess rainfall due to snow capping (mm H2O /s) [+]

#ifdef LCH4
    real(r8), pointer :: lake_raw(:)        ! aerodynamic resistance for moisture (s/m)
#endif
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer , parameter  :: niters = 4    ! maximum number of iterations for surface temperature
    real(r8), parameter :: beta1 = 1._r8  ! coefficient of convective velocity (in computing W_*) [-]
    real(r8), parameter :: zii = 1000._r8 ! convective boundary height [m]
    integer  :: i,fc,fp,g,c,p           ! do loop or array index
    integer  :: fncopy                  ! number of values in pft filter copy
    integer  :: fnold                   ! previous number of pft filter values
    integer  :: fpcopy(num_lakep)       ! pft filter copy for iteration loop
    integer  :: iter                    ! iteration index
    integer  :: nmozsgn(lbp:ubp)        ! number of times moz changes sign
    integer  :: jtop(lbc:ubc)           ! top level for each column (no longer all 1)
    real(r8) :: dtime                   ! land model time step (sec)
    real(r8) :: ax                      ! used in iteration loop for calculating t_grnd (numerator of NR solution)
    real(r8) :: bx                      ! used in iteration loop for calculating t_grnd (denomin. of NR solution)
    real(r8) :: degdT                   ! d(eg)/dT
    real(r8) :: dqh(lbp:ubp)            ! diff of humidity between ref. height and surface
    real(r8) :: dth(lbp:ubp)            ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                    ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dzsur(lbc:ubc)          ! 1/2 the top layer thickness (m)
    real(r8) :: eg                      ! water vapor pressure at temperature T [pa]
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
    real(r8) :: tksur(lbc:ubc)          ! thermal conductivity of snow/soil (w/m/kelvin)
    real(r8) :: tsur(lbc:ubc)           ! top layer temperature
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
    real(r8) :: u2m                     ! 2 m wind speed (m/s)
    real(r8) :: fm(lbp:ubp)             ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: bw                      ! partial density of water (ice + liquid)
    real(r8) :: t_grnd_temp             ! Used in surface flux correction over frozen ground
    real(r8) :: betaprime(lbc:ubc)      ! Effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
    real(r8) :: e_ref2m                 ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                  ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m              ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT               ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: sabg_nir                ! NIR that is absorbed (W/m^2)

    ! For calculating roughness lengths
    real(r8) :: cur                     ! Charnock parameter (-)
    real(r8) :: fetch(lbc:ubc)          ! Fetch (m)
    real(r8) :: sqre0                   ! root of roughness Reynolds number
    real(r8), parameter :: kva0 = 1.51e-5_r8   ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
    real(r8) :: kva0temp                ! (K) temperature for kva0; will be set below
    real(r8), parameter :: kva0pres = 1.013e5_r8 ! (Pa) pressure for kva0
    real(r8) :: kva                     ! kinematic viscosity of air at ground temperature and forcing pressure
    real(r8), parameter :: prn = 0.713  ! Prandtl # for air at neutral stability
    real(r8), parameter :: sch = 0.66   ! Schmidt # for water in air at neutral stability

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
    lat            => clm3%g%lat

    ! Assign local pointers to derived type members (column-level)

    cgridcell      => clm3%g%l%c%gridcell
    dz             => clm3%g%l%c%cps%dz
    dz_lake        => clm3%g%l%c%cps%dz_lake
    t_soisno       => clm3%g%l%c%ces%t_soisno
    t_lake         => clm3%g%l%c%ces%t_lake
    snl            => clm3%g%l%c%cps%snl
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    t_grnd         => clm3%g%l%c%ces%t_grnd
    ws             => clm3%g%l%c%cps%ws
    ks             => clm3%g%l%c%cps%ks
    savedtke1      => clm3%g%l%c%cps%savedtke1
    ust_lake       => clm3%g%l%c%cps%ust_lake
    z0mg_col       => clm3%g%l%c%cps%z0mg
    z0hg_col       => clm3%g%l%c%cps%z0hg
    z0qg_col       => clm3%g%l%c%cps%z0qg
    lakedepth      => clm3%g%l%c%cps%lakedepth
    lakefetch      => clm3%g%l%c%cps%lakefetch
#ifdef LCH4
    lake_raw       => clm3%g%l%c%cch4%lake_raw
#endif
    sabg_chk       => clm3%g%l%c%p%pef%sabg_chk

    ! Assign local pointers to derived type members (pft-level)

    pcolumn        => clm3%g%l%c%p%column
    pgridcell      => clm3%g%l%c%p%gridcell
    sabg           => clm3%g%l%c%p%pef%sabg
    t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
    q_ref2m        => clm3%g%l%c%p%pes%q_ref2m
    t_veg          => clm3%g%l%c%p%pes%t_veg
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net
    eflx_soil_grnd => clm3%g%l%c%p%pef%eflx_soil_grnd
    eflx_lh_tot    => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_lh_grnd   => clm3%g%l%c%p%pef%eflx_lh_grnd
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot
    ram1           => clm3%g%l%c%p%pps%ram1
    ram1_lake      => clm3%g%l%c%p%pps%ram1_lake
    taux           => clm3%g%l%c%p%pmf%taux
    tauy           => clm3%g%l%c%p%pmf%tauy
    qflx_prec_grnd => clm3%g%l%c%p%pwf%qflx_prec_grnd
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_tot  => clm3%g%l%c%p%pwf%qflx_evap_tot
    eflx_gnet      => clm3%g%l%c%p%pef%eflx_gnet
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
    forc_hgt_t_pft => clm3%g%l%c%p%pps%forc_hgt_t_pft
    forc_hgt_q_pft => clm3%g%l%c%p%pps%forc_hgt_q_pft
    rh_ref2m       => clm3%g%l%c%p%pes%rh_ref2m
    sabg_lyr       => clm3%g%l%c%p%pef%sabg_lyr
    ! For calculation of NIR fraction of sabg
    fsds_nir_d    => clm3%g%l%c%p%pef%fsds_nir_d
    fsds_nir_i    => clm3%g%l%c%p%pef%fsds_nir_i
    fsr_nir_d     => clm3%g%l%c%p%pef%fsr_nir_d
    fsr_nir_i     => clm3%g%l%c%p%pef%fsr_nir_i
    qflx_snwcp_ice => clm3%g%l%c%p%pwf%qflx_snwcp_ice
    qflx_snwcp_liq => clm3%g%l%c%p%pwf%qflx_snwcp_liq


    ! Determine step size

    dtime = get_step_size()

    ! Begin calculations

    kva0temp = 20._r8 + tfrz

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = cgridcell(c)

       ! Set fetch for prognostic roughness length-- if not found in surface data.
       ! This is poorly constrained, and should eventually be based on global lake data
       ! For now, base on lake depth, assuming that small lakes are likely to be shallower
       ! The dependence will be weak, especially for large fetch
       ! http://www.chebucto.ns.ca/ccn/info/Science/SWCS/DATA/morphology.html#zr, based on
       ! Hutchinson, G.E. 1957 A treatise on limnology v.1. Geography, Physics and Chemistry,
       ! and Wetzel, R.G., and Likens, G.E.. 1991. Limnological Analyses, suggests lakes usually have
       ! depths less than 2% of their diameter.
       if (lakefetch(c) > 0._r8) then ! fetch available in surface data
          fetch(c) = lakefetch(c)
       else ! Estimate crudely based on lake depth
          if (lakedepth(c) < 4._r8) then
             fetch(c) = 100._r8 ! Roughly the smallest lakes resolveable in the GLWD
          else
             fetch(c) = 25._r8*lakedepth(c)
          end if
       end if

       ! Initialize roughness lengths

       if (t_grnd(c) > tfrz) then   ! for unfrozen lake
          z0mg(p) = z0mg_col(c)
          kva = kva0 * (t_grnd(c)/kva0temp)**1.5_r8 * kva0pres/forc_pbot(g) ! kinematic viscosity of air
          sqre0 = (max(z0mg(p)*ust_lake(c)/kva,0.1_r8))**0.5_r8   ! Square root of roughness Reynolds number
          z0hg(p) = z0mg(p) * exp( -vkc/prn*( 4._r8*sqre0 - 3.2_r8) ) ! SH roughness length
          z0qg(p) = z0mg(p) * exp( -vkc/sch*( 4._r8*sqre0 - 4.2_r8) ) ! LH roughness length
          z0qg(p) = max(z0qg(p), minz0lake)
          z0hg(p) = max(z0hg(p), minz0lake)
       else if (snl(c) == 0) then    ! frozen lake with ice
          z0mg(p) = z0frzlake
          z0hg(p) = z0mg(p)/exp(0.13_r8 * (ust_lake(c)*z0mg(p)/1.5e-5_r8)**0.45_r8) ! Consistent with BareGroundFluxes
          z0qg(p) = z0hg(p)
       else                          ! use roughness over snow as in Biogeophysics1
          z0mg(p) = zsno
          z0hg(p) = z0mg(p)/exp(0.13_r8 * (ust_lake(c)*z0mg(p)/1.5e-5_r8)**0.45_r8) ! Consistent with BareGroundFluxes
          z0qg(p) = z0hg(p)
       end if

       ! Surface temperature and fluxes

       ! Moved in from Biogeophysics1
       forc_hgt_u_pft(p) = forc_hgt_u_pft(p) + z0mg(p)
       forc_hgt_t_pft(p) = forc_hgt_t_pft(p) + z0mg(p)
       forc_hgt_q_pft(p) = forc_hgt_q_pft(p) + z0mg(p)

       ! Find top layer
       jtop(c) = snl(c) + 1

       if (snl(c) < 0) then
           betaprime(c) = sabg_lyr(p,jtop(c))/max(1.e-5_r8,sabg(p))  ! Assuming one pft
           dzsur(c) = dz(c,jtop(c))/2._r8
       else ! no snow layers
       ! Calculate the NIR fraction of absorbed solar.
           ! Total NIR absorbed:
           sabg_nir = fsds_nir_d(p) + fsds_nir_i(p) - fsr_nir_d(p) - fsr_nir_i(p)
           sabg_nir = min(sabg_nir, sabg(p))
           betaprime(c) = sabg_nir/max(1.e-5_r8,sabg(p))
           ! Some fraction of the "visible" may be absorbed in the surface layer.
           betaprime(c) = betaprime(c) + (1._r8-betaprime(c))*betavis
           dzsur(c) = dz_lake(c,1)/2._r8
       end if

       sabg_chk(p)  = sabg(p)
       ! Originally dzsur was 1*dz, but it should it be 1/2 dz.

       ! Saturated vapor pressure, specific humidity and their derivatives
       ! at lake surface

       call QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height

       thm(p) = forc_t(g) + 0.0098_r8*forc_hgt_t_pft(p)   ! intermediate variable
       thv(c) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))     ! virtual potential T
    end do



    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       nmozsgn(p) = 0
       obuold(p) = 0._r8
       displa(p) = 0._r8


       ! Latent heat

! Attn EK: This PERGRO code was commented on the basis of Gordon Bonan's suggestion, but this may deserve more look.
!#if (defined PERGRO)
!       htvp(c) = hvap
!#else
       if (t_grnd(c) > tfrz) then
          htvp(c) = hvap
       else
          htvp(c) = hsub
       end if
!#endif
       ! Zack Subin, 3/26/09: Changed to ground temperature rather than the air temperature above.

       ! Initialize stability variables

       ur(p)    = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p)   = thm(p)-t_grnd(c)
!#if (defined PERGRO)
!       dth(p)   = 0.0_r8
!#endif
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
          if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
             tksur(c) = savedtke1(c)
             ! Set this to the eddy conductivity from the last
             ! timestep, as the molecular conductivity will be orders of magnitude too small.
             ! It will be initialized in initSLakeMod to the molecular conductivity for the first timestep if arbinit.
             tsur(c) = t_lake(c,1)
          else if (snl(c) == 0) then  !frozen but no snow layers
             tksur(c) = tkice   ! This is an approximation because the whole layer may not be frozen, and it is not
                                ! accounting for the physical (but not nominal) expansion of the frozen layer.
             tsur(c) = t_lake(c,1)
          else
          !Need to calculate thermal conductivity of the top snow layer
             bw = (h2osoi_ice(c,jtop(c))+h2osoi_liq(c,jtop(c)))/dz(c,jtop(c))
             tksur(c) = tkair + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkair)
             tsur(c) = t_soisno(c,jtop(c))
          end if

          ! Determine aerodynamic resistances

          ram(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
          rah(p)  = 1._r8/(temp1(p)*ustar(p))
          raw(p)  = 1._r8/(temp2(p)*ustar(p))
#if (defined LCH4)
          lake_raw(c) = raw(p) ! Pass out for calculating ground ch4 conductance
#endif
          ram1(p) = ram(p)   !pass value to global variable
          ram1_lake(p) = ram1(p) ! for history

          ! Get derivative of fluxes with respect to ground temperature

          stftg3(p) = emg_lake*sb*tgbef(c)*tgbef(c)*tgbef(c)

          ! Changed surface temperature from t_lake(c,1) to tsur(c).
          ! Also adjusted so that if there are snow layers present, the top layer absorption
          ! from SNICAR is assigned to the surface skin.
          ax  = betaprime(c)*sabg(p) + emg_lake*forc_lwrad(g) + 3._r8*stftg3(p)*tgbef(c) &
               + forc_rho(g)*cpair/rah(p)*thm(p) &
               - htvp(c)*forc_rho(g)/raw(p)*(qsatg(c)-qsatgdT(c)*tgbef(c) - forc_q(g)) &
               + tksur(c)*tsur(c)/dzsur(c)
          !Changed sabg(p) to betaprime(c)*sabg(p).
          bx  = 4._r8*stftg3(p) + forc_rho(g)*cpair/rah(p) &
               + htvp(c)*forc_rho(g)/raw(p)*qsatgdT(c) + tksur(c)/dzsur(c)

          t_grnd(c) = ax/bx

          ! Update htvp
!#ifndef PERGRO
       if (t_grnd(c) > tfrz) then
          htvp(c) = hvap
       else
          htvp(c) = hsub
       end if
!#endif

          ! Surface fluxes of momentum, sensible and latent heat
          ! using ground temperatures from previous time step

          eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(p))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c))-forc_q(g))/raw(p)

          ! Re-calculate saturated vapor pressure, specific humidity and their
          ! derivatives at lake surface

          call QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

          dth(p)=thm(p)-t_grnd(c)
!#if (defined PERGRO)
!          dth(p)   = 0.0_r8
!#endif
          dqh(p)=forc_q(g)-qsatg(c)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)

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

          if (t_grnd(c) > tfrz .and. snl(c) == 0) then ! t_grnd hasn't been corrected yet if snow layers but above frz
             ! Update roughness lengths using approach in Subin et al. 2011
             ! Also allow wave development (phase speed) to be depth-limited as well as fetch-limited
             if (lake_use_old_fcrit_minz0) then
                ! Original formulation in Subin et al. 2011; converted Vickers & Mahrt 1997 to use u instead of u*
                ! assuming u = 0.1 u*.
                ! That probably slightly overestimates the dimensionless fetch as u* is often smaller than 0.1 u
                cur = cur0 + curm* exp( max( -(fetch(c)*grav/ur(p)/ur(p))**(1._r8/3._r8)/fcrit, &   ! Fetch-limited
                                             -(lakedepth(c)*grav/ur(p)/ur(p))**0.5_r8 ) )           ! depth-limited
                ! In this case fcrit is 22, not 100 in clm_varcon
             else
                ! Fetch relationship from Vickers & Mahrt 1997
                cur = cur0 + curm* exp( max( -(fetch(c)*grav/ustar(p)/ustar(p))**(1._r8/3._r8)/fcrit, &   ! Fetch-limited
                                          -(lakedepth(c)*grav/ur(p)/ur(p))**0.5_r8 ) )           ! depth-limited
             end if
             kva = kva0 * (t_grnd(c)/kva0temp)**1.5_r8 * kva0pres/forc_pbot(g) ! kinematic viscosity of air
             z0mg(p) = max(cus*kva/max(ustar(p),1.e-4_r8), cur*ustar(p)*ustar(p)/grav) ! momentum roughness length
                                      ! This lower limit on ustar is just to prevent floating point exceptions and
                                      ! should not be important
             z0mg(p) = max(z0mg(p), minz0lake) ! This limit is redundant with current values.
             sqre0 = (max(z0mg(p)*ustar(p)/kva,0.1_r8))**0.5_r8   ! Square root of roughness Reynolds number
             z0hg(p) = z0mg(p) * exp( -vkc/prn*( 4._r8*sqre0 - 3.2_r8) ) ! SH roughness length
             z0qg(p) = z0mg(p) * exp( -vkc/sch*( 4._r8*sqre0 - 4.2_r8) ) ! LH roughness length
             z0qg(p) = max(z0qg(p), minz0lake)
             z0hg(p) = max(z0hg(p), minz0lake)
          else if (snl(c) == 0) then
             ! in case it was above freezing and now below freezing
             z0mg(p) = z0frzlake
             z0hg(p) = z0mg(p)/exp(0.13_r8 * (ustar(p)*z0mg(p)/1.5e-5_r8)**0.45_r8) ! Consistent with BareGroundFluxes
             z0qg(p) = z0hg(p)
          else ! Snow layers
             ! z0mg won't have changed
             z0hg(p) = z0mg(p)/exp(0.13_r8 * (ustar(p)*z0mg(p)/1.5e-5_r8)**0.45_r8) ! Consistent with BareGroundFluxes
             z0qg(p) = z0hg(p)
          end if

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

       ! If there is snow on the ground or lake is frozen and t_grnd > tfrz: reset t_grnd = tfrz.
       ! Re-evaluate ground fluxes.
       ! [ZMS 1/7/11] Only for resolved snow layers, as unresolved snow does not have a temperature state and
       ! can accumulate on unfrozen lakes in SLakeHydrology; will be melted in SLakeTemperature or bring lake top
       ! to freezing.
       ! note that qsatg and qsatgdT should be f(tgbef) (PET: not sure what this
       ! comment means)
       ! Zack Subin, 3/27/09: Since they are now a function of whatever t_grnd was before cooling
       !    to freezing temperature, then this value should be used in the derivative correction term.
       ! Allow convection if ground temp is colder than lake but warmer than 4C, or warmer than 
       !    lake which is warmer than freezing but less than 4C.
       if ( (snl(c) < 0 .or. t_lake(c,1) <= tfrz) .and. t_grnd(c) > tfrz) then
          t_grnd_temp = t_grnd(c)
          t_grnd(c) = tfrz
          eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(p))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
       else if ( (t_lake(c,1) > t_grnd(c) .and. t_grnd(c) > tdmax) .or. &
                 (t_lake(c,1) < t_grnd(c) .and. t_lake(c,1) > tfrz .and. t_grnd(c) < tdmax) ) then
                 ! Convective mixing will occur at surface
          t_grnd_temp = t_grnd(c)
          t_grnd(c) = t_lake(c,1)
          eflx_sh_grnd(p) = forc_rho(g)*cpair*(t_grnd(c)-thm(p))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
       end if

          ! Update htvp
!#ifndef PERGRO
       if (t_grnd(c) > tfrz) then
          htvp(c) = hvap
       else
          htvp(c) = hsub
       end if
!#endif

       ! Net longwave from ground to atmosphere

       ! eflx_lwrad_out(p) = (1._r8-emg_lake)*forc_lwrad(g) + stftg3(p)*(-3._r8*tgbef(c)+4._r8*t_grnd(c))
       ! What is tgbef doing in this equation? Can't it be exact now? --Zack Subin, 4/14/09
       eflx_lwrad_out(p) = (1._r8-emg_lake)*forc_lwrad(g) + emg_lake*sb*t_grnd(c)**4._r8

       ! Ground heat flux

       eflx_soil_grnd(p) = sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) - &
            eflx_sh_grnd(p) - htvp(c)*qflx_evap_soi(p)
       ! The original code in BiogeophysicsLake had a bug that calculated incorrect fluxes but conserved energy.
       ! This is kept as the full sabg (not just that absorbed at surface) so that the energy balance check will be correct.
       !This is the effective energy flux into the ground including the lake [and now snow in CLM 4] solar absorption
       !below the surface.  This also keeps the output FGR similar to non-lakes by including the light & heat flux.
       ! The variable eflx_gnet will be used to pass the actual heat flux
       !from the ground interface into the lake.

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
       ! Effectively moved to SLakeTemp

       eflx_gnet(p) = betaprime(c) * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
            eflx_sh_tot(p) + eflx_lh_tot(p))
       ! This is the actual heat flux from the ground interface into the lake, not including
       ! the light that penetrates the surface.

       !u2m = max(1.0_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))
       ! u2 often goes below 1 m/s; it seems like the only reason for this minimum is to
       ! keep it from being zero in the ks equation below; 0.1 m/s is a better limit for
       ! stable conditions --ZS
       u2m = max(0.1_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))

       ws(c) = 1.2e-03_r8 * u2m
       ks(c) = 6.6_r8*sqrt(abs(sin(lat(g))))*(u2m**(-1.84_r8))

       ! Update column roughness lengths and friction velocity
       z0mg_col(c) = z0mg(p)
       z0hg_col(c) = z0hg(p)
       z0qg_col(c) = z0qg(p)
       ust_lake(c) = ustar(p)

    end do

    ! The following are needed for global average on history tape.

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)
       t_veg(p) = forc_t(g)
       eflx_lwrad_net(p)  = eflx_lwrad_out(p) - forc_lwrad(g)
       qflx_prec_grnd(p) = forc_rain(g) + forc_snow(g)

       ! Because they will be used in pft2col initialize here.
       ! This will be overwritten in SLakeHydrology
       qflx_snwcp_ice(p) = 0._r8
       qflx_snwcp_liq(p) = 0._r8
       
    end do

  end subroutine SLakeFluxes

end module SLakeFluxesMod
