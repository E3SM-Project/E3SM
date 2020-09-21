module LakeFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates surface fluxes and temperature for lakes.
  ! Created by Zack Subin, 2009
  !
  ! !USES
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use decompMod            , only : bounds_type
  use atm2lndType          , only : atm2lnd_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityType , only : frictionvel_type
  use LakeStateType        , only : lakestate_type
  use SolarAbsorbedType    , only : solarabs_type
  use TemperatureType      , only : temperature_type
  use WaterfluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use GridcellType         , only : grc_pp   
  use TopounitDataType     , only : top_as, top_af ! atmospheric state and flux variables  
  use ColumnType           , only : col_pp
  use ColumnDataType       , only : col_es, col_ws  
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf  
  !    
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LakeFluxes
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LakeFluxes(bounds, num_lakec, filter_lakec, num_lakep, filter_lakep, &
       atm2lnd_vars, solarabs_vars, frictionvel_vars, temperature_vars, &
       energyflux_vars, waterstate_vars, waterflux_vars, lakestate_vars) 
    !
    ! !DESCRIPTION:
    ! Calculates lake temperatures and surface fluxes.
    ! Lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
    ! and soil layers with active temperature and gas diffusion below.
    ! WARNING: This subroutine assumes lake columns have one and only one pft.
    !
    ! !USES:
    use elm_varpar          , only : nlevlak
    use elm_varcon          , only : hvap, hsub, hfus, cpair, cpliq, tkwat, tkice, tkair
    use elm_varcon          , only : sb, vkc, grav, denh2o, tfrz, spval, zsno
    use elm_varctl          , only : use_lch4
    use LakeCon             , only : betavis, z0frzlake, tdmax, emg_lake
    use LakeCon             , only : lake_use_old_fcrit_minz0
    use LakeCon             , only : minz0lake, cur0, cus, curm, fcrit
    use QSatMod             , only : QSat
    use FrictionVelocityMod , only : FrictionVelocity, MoninObukIni
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_lakec         ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_lakec(:)   ! column filter for non-lake points
    integer                , intent(in)    :: num_lakep         ! number of column non-lake points in pft filter
    integer                , intent(in)    :: filter_lakep(:)   ! patch filter for non-lake points
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    type(waterflux_type)   , intent(inout) :: waterflux_vars
    type(temperature_type) , intent(inout) :: temperature_vars
    type(lakestate_type)   , intent(inout) :: lakestate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: z0mg_col(:)               ! roughness length over ground, momentum [m]
    real(r8), pointer :: z0hg_col(:)               ! roughness length over ground, sensible heat [m]
    real(r8), pointer :: z0qg_col(:)               ! roughness length over ground, latent heat [m]
    integer , parameter  :: niters = 4             ! maximum number of iterations for surface temperature
    real(r8), parameter :: beta1 = 1._r8           ! coefficient of convective velocity (in computing W_*) [-]
    real(r8), parameter :: zii = 1000._r8          ! convective boundary height [m]
    integer  :: i,fc,fp,g,t,c,p                    ! do loop or array index
    integer  :: fncopy                             ! number of values in pft filter copy
    integer  :: fnold                              ! previous number of pft filter values
    integer  :: fpcopy(num_lakep)                  ! patch filter copy for iteration loop
    integer  :: iter                               ! iteration index
    integer  :: nmozsgn(bounds%begp:bounds%endp)   ! number of times moz changes sign
    integer  :: jtop(bounds%begc:bounds%endc)      ! top level for each column (no longer all 1)
    real(r8) :: ax                                 ! used in iteration loop for calculating t_grnd (numerator of NR solution)
    real(r8) :: bx                                 ! used in iteration loop for calculating t_grnd (denomin. of NR solution)
    real(r8) :: degdT                              ! d(eg)/dT
    real(r8) :: dqh(bounds%begp:bounds%endp)       ! diff of humidity between ref. height and surface
    real(r8) :: dth(bounds%begp:bounds%endp)       ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                               ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dzsur(bounds%begc:bounds%endc)     ! 1/2 the top layer thickness (m)
    real(r8) :: eg                                 ! water vapor pressure at temperature T [pa]
    real(r8) :: htvp(bounds%begc:bounds%endc)      ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8) :: obu(bounds%begp:bounds%endp)       ! monin-obukhov length (m)
    real(r8) :: obuold(bounds%begp:bounds%endp)    ! monin-obukhov length of previous iteration
    real(r8) :: qsatg(bounds%begc:bounds%endc)     ! saturated humidity [kg/kg]
    real(r8) :: qsatgdT(bounds%begc:bounds%endc)   ! d(qsatg)/dT
    real(r8) :: qstar                              ! moisture scaling parameter
    real(r8) :: ram(bounds%begp:bounds%endp)       ! aerodynamical resistance [s/m]
    real(r8) :: rah(bounds%begp:bounds%endp)       ! thermal resistance [s/m]
    real(r8) :: raw(bounds%begp:bounds%endp)       ! moisture resistance [s/m]
    real(r8) :: stftg3(bounds%begp:bounds%endp)    ! derivative of fluxes w.r.t ground temperature
    real(r8) :: temp1(bounds%begp:bounds%endp)     ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp)   ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)     ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp)   ! relation for specific humidity profile applied at 2-m
    real(r8) :: tgbef(bounds%begc:bounds%endc)     ! initial ground temperature
    real(r8) :: thm(bounds%begp:bounds%endp)       ! intermediate variable (forc_t+0.0098*forc_hgt_t_patch)
    real(r8) :: thv(bounds%begc:bounds%endc)       ! virtual potential temperature (kelvin)
    real(r8) :: thvstar                            ! virtual potential temperature scaling parameter
    real(r8) :: tksur(bounds%begc:bounds%endc)     ! thermal conductivity of snow/soil (w/m/kelvin)
    real(r8) :: tsur(bounds%begc:bounds%endc)      ! top layer temperature
    real(r8) :: tstar                              ! temperature scaling parameter
    real(r8) :: um(bounds%begp:bounds%endp)        ! wind speed including the stablity effect [m/s]
    real(r8) :: ur(bounds%begp:bounds%endp)        ! wind speed at reference height [m/s]
    real(r8) :: ustar(bounds%begp:bounds%endp)     ! friction velocity [m/s]
    real(r8) :: wc                                 ! convective velocity [m/s]
    real(r8) :: zeta                               ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: zldis(bounds%begp:bounds%endp)     ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(bounds%begp:bounds%endp)    ! displacement (always zero) [m]
    real(r8) :: z0mg(bounds%begp:bounds%endp)      ! roughness length over ground, momentum [m]
    real(r8) :: z0hg(bounds%begp:bounds%endp)      ! roughness length over ground, sensible heat [m]
    real(r8) :: z0qg(bounds%begp:bounds%endp)      ! roughness length over ground, latent heat [m]
    real(r8) :: u2m                                ! 2 m wind speed (m/s)
    real(r8) :: fm(bounds%begp:bounds%endp)        ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: bw                                 ! partial density of water (ice + liquid)
    real(r8) :: t_grnd_temp                        ! Used in surface flux correction over frozen ground
    real(r8) :: betaprime(bounds%begc:bounds%endc) ! Effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
    real(r8) :: e_ref2m                            ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                             ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m                         ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT                          ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: sabg_nir                           ! NIR that is absorbed (W/m^2)

    ! For calculating roughness lengths
    real(r8) :: cur                                ! Charnock parameter (-)
    real(r8) :: fetch(bounds%begc:bounds%endc)     ! Fetch (m)
    real(r8) :: sqre0                              ! root of roughness Reynolds number
    real(r8), parameter :: kva0 = 1.51e-5_r8       ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
    real(r8) :: kva0temp                           ! (K) temperature for kva0; will be set below
    real(r8), parameter :: kva0pres = 1.013e5_r8   ! (Pa) pressure for kva0
    real(r8) :: kva                                ! kinematic viscosity of air at ground temperature and forcing pressure
    real(r8), parameter :: prn = 0.713             ! Prandtl # for air at neutral stability
    real(r8), parameter :: sch = 0.66              ! Schmidt # for water in air at neutral stability
    !-----------------------------------------------------------------------

    associate(                                                           & 
         snl              =>    col_pp%snl                                , & ! Input:  [integer  (:)   ]  number of snow layers                              
         dz               =>    col_pp%dz                                 , & ! Input:  [real(r8) (:,:) ]  layer thickness for soil or snow (m)            
         dz_lake          =>    col_pp%dz_lake                            , & ! Input:  [real(r8) (:,:) ]  layer thickness for lake (m)                    
         lakedepth        =>    col_pp%lakedepth                          , & ! Input:  [real(r8) (:)   ]  variable lake depth (m)                           
         
         forc_t           =>    top_as%tbot                            , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)                  
         forc_th          =>    top_as%thbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)        
         forc_pbot        =>    top_as%pbot                            , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                         
         forc_q           =>    top_as%qbot                            , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)             
         forc_rho         =>    top_as%rhobot                          , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)                                 
         forc_lwrad       =>    top_af%lwrad                           , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)   
         forc_snow        =>    top_af%snow                            , & ! Input:  [real(r8) (:)   ]  snow rate (kg H2O/m**2/s, or mm liquid H2O/s)                                  
         forc_rain        =>    top_af%rain                            , & ! Input:  [real(r8) (:)   ]  rain rate (kg H2O/m**2/s, or mm liquid H2O/s)                                  
         forc_u           =>    top_as%ubot                            , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)    
         forc_v           =>    top_as%vbot                            , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)   
         
         fsds_nir_d       =>    solarabs_vars%fsds_nir_d_patch         , & ! Input:  [real(r8) (:)   ]  incident direct beam nir solar radiation (W/m**2) 
         fsds_nir_i       =>    solarabs_vars%fsds_nir_i_patch         , & ! Input:  [real(r8) (:)   ]  incident diffuse nir solar radiation (W/m**2)     
         fsr_nir_d        =>    solarabs_vars%fsr_nir_d_patch          , & ! Input:  [real(r8) (:)   ]  reflected direct beam nir solar radiation (W/m**2)
         fsr_nir_i        =>    solarabs_vars%fsr_nir_i_patch          , & ! Input:  [real(r8) (:)   ]  reflected diffuse nir solar radiation (W/m**2)    
         sabg_lyr         =>    solarabs_vars%sabg_lyr_patch           , & ! Input:  [real(r8) (:,:) ]  absorbed solar radiation (pft,lyr) [W/m2]       
         sabg_chk         =>    solarabs_vars%sabg_chk_patch           , & ! Output: [real(r8) (:)   ]  sum of soil/snow using current fsno, for balance check
         sabg             =>    solarabs_vars%sabg_patch               , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)       
         
         savedtke1        =>    lakestate_vars%savedtke1_col           , & ! Input:  [real(r8) (:)   ]  top level eddy conductivity from previous timestep (W/mK)
         lakefetch        =>    lakestate_vars%lakefetch_col           , & ! Input:  [real(r8) (:)   ]  lake fetch from surface data (m)                  
         
         h2osoi_liq       =>    col_ws%h2osoi_liq         , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osoi_ice       =>    col_ws%h2osoi_ice         , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                

         t_lake           =>    col_es%t_lake            , & ! Input:  [real(r8) (:,:) ]  lake temperature (Kelvin)                       
         t_soisno         =>    col_es%t_soisno          , & ! Input:  [real(r8) (:,:) ]  soil (or snow) temperature (Kelvin)             

         forc_hgt_u_patch =>    frictionvel_vars%forc_hgt_u_patch      , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft level [m]     
         forc_hgt_t_patch =>    frictionvel_vars%forc_hgt_t_patch      , & ! Input:  [real(r8) (:)   ]  observational height of temperature at pft level [m]
         forc_hgt_q_patch =>    frictionvel_vars%forc_hgt_q_patch      , & ! Input:  [real(r8) (:)   ]  observational height of specific humidity at pft level [m]

         q_ref2m          =>    veg_ws%q_ref2m          , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)      
         rh_ref2m         =>    veg_ws%rh_ref2m         , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)          
         qflx_evap_soi    =>    veg_wf%qflx_evap_soi     , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)          
         qflx_evap_tot    =>    veg_wf%qflx_evap_tot     , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     

         qflx_snwcp_ice   =>    veg_wf%qflx_snwcp_ice    , & ! Output: [real(r8) (:)   ]  excess snowfall due to snow capping (mm H2O /s) [+]
         qflx_snwcp_liq   =>    veg_wf%qflx_snwcp_liq    , & ! Output: [real(r8) (:)   ]  excess rainfall due to snow capping (mm H2O /s) [+] 
         qflx_prec_grnd   =>    veg_wf%qflx_prec_grnd    , & ! Output: [real(r8) (:)   ]  water onto ground including canopy runoff [kg/(m2 s)]
         qflx_dirct_rain  =>    veg_wf%qflx_dirct_rain   , & ! Output: [real(r8) (:)   ]  direct rain throughfall [mm H2O/s]
         qflx_leafdrip    =>    veg_wf%qflx_leafdrip     , & ! Output: [real(r8) (:)   ]  leaf rain drip [mm H2O/s]

         t_veg            =>    veg_es%t_veg             , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                   
         t_ref2m          =>    veg_es%t_ref2m         , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)       
         t_grnd           =>    col_es%t_grnd            , & ! Output: [real(r8) (:)   ]  ground temperature (Kelvin)                       

         ram1             =>    frictionvel_vars%ram1_patch            , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                    

         eflx_lwrad_out   =>    veg_ef%eflx_lwrad_out   , & ! Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)    
         eflx_lwrad_net   =>    veg_ef%eflx_lwrad_net   , & ! Output: [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm] 
         eflx_soil_grnd   =>    veg_ef%eflx_soil_grnd   , & ! Output: [real(r8) (:)   ]  soil heat flux (W/m**2) [+ = into soil]           
         eflx_lh_tot      =>    veg_ef%eflx_lh_tot      , & ! Output: [real(r8) (:)   ]  total latent heat flux (W/m8*2)  [+ to atm]       
         eflx_lh_grnd     =>    veg_ef%eflx_lh_grnd     , & ! Output: [real(r8) (:)   ]  ground evaporation heat flux (W/m**2) [+ to atm]  
         eflx_sh_grnd     =>    veg_ef%eflx_sh_grnd     , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot      =>    veg_ef%eflx_sh_tot      , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]      
         eflx_gnet        =>    veg_ef%eflx_gnet        , & ! Output: [real(r8) (:)   ]  net heat flux into ground (W/m**2)                
         taux             =>    veg_ef%taux             , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)              
         tauy             =>    veg_ef%tauy             , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)              

         ks               =>    lakestate_vars%ks_col                  , & ! Output: [real(r8) (:)   ]  coefficient passed to LakeTemperature            
         ws               =>    lakestate_vars%ws_col                  , & ! Output: [real(r8) (:)   ]  surface friction velocity (m/s)                   
         betaprime        =>    lakestate_vars%betaprime_col           , & ! Output: [real(r8) (:)   ]  fraction of solar rad absorbed at surface: equal to NIR fraction
         ram1_lake        =>    lakestate_vars%ram1_lake_patch         , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                    
         ust_lake         =>    lakestate_vars%ust_lake_col            , & ! Output: [real(r8) (:)   ]  friction velocity (m/s)                           
         lake_raw         =>    lakestate_vars%lake_raw_col            , & ! Output: [real(r8) (:)   ]  aerodynamic resistance for moisture (s/m)   
         
         begp             =>    bounds%begp                            , &
         endp             =>    bounds%endp                              &
         )

      ! the following cause a crash if they are set as associated
      z0mg_col => frictionvel_vars%z0mg_col
      z0hg_col => frictionvel_vars%z0hg_col
      z0qg_col => frictionvel_vars%z0qg_col

      kva0temp = 20._r8 + tfrz

      do fp = 1, num_lakep
         p = filter_lakep(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = col_pp%gridcell(c)

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
            kva = kva0 * (t_grnd(c)/kva0temp)**1.5_r8 * kva0pres/forc_pbot(t) ! kinematic viscosity of air
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

         forc_hgt_u_patch(p) = forc_hgt_u_patch(p) + z0mg(p)
         forc_hgt_t_patch(p) = forc_hgt_t_patch(p) + z0mg(p)
         forc_hgt_q_patch(p) = forc_hgt_q_patch(p) + z0mg(p)

         ! Find top layer
         jtop(c) = snl(c) + 1

         if (snl(c) < 0) then
            betaprime(c) = sabg_lyr(p,jtop(c))/max(1.e-5_r8,sabg(p))  ! Assuming one pft
            dzsur(c) = dz(c,jtop(c))/2._r8
         else ! no snow layers
            ! Calculate the NIR fraction of absorbed solar.
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

         call QSat(t_grnd(c), forc_pbot(t), eg, degdT, qsatg(c), qsatgdT(c))

         ! Potential, virtual potential temperature, and wind speed at the
         ! reference height

         thm(p) = forc_t(t) + 0.0098_r8*forc_hgt_t_patch(p)   ! intermediate variable
         thv(c) = forc_th(t)*(1._r8+0.61_r8*forc_q(t))     ! virtual potential T
      end do



      do fp = 1, num_lakep
         p = filter_lakep(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         nmozsgn(p) = 0
         obuold(p) = 0._r8
         displa(p) = 0._r8

         ! Latent heat

         if (t_grnd(c) > tfrz) then
            htvp(c) = hvap
         else
            htvp(c) = hsub
         end if
         ! Zack Subin, 3/26/09: Changed to ground temperature rather than the air temperature above.

         ! Initialize stability variables

         ur(p)    = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)))
         dth(p)   = thm(p)-t_grnd(c)
         dqh(p)   = forc_q(t)-qsatg(c)
         dthv     = dth(p)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p) - 0._r8

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

         call FrictionVelocity(begp, endp, fncopy, fpcopy, &
              displa(begp:endp), z0mg(begp:endp), z0hg(begp:endp), z0qg(begp:endp), &
              obu(begp:endp), iter, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
              temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp), &
              frictionvel_vars)

         do fp = 1, fncopy
            p = fpcopy(fp)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            g = veg_pp%gridcell(p)

            tgbef(c) = t_grnd(c)
            if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
               tksur(c) = savedtke1(c)
               ! Set this to the eddy conductivity from the last
               ! timestep, as the molecular conductivity will be orders of magnitude too small.
               ! It will be initialized in initLakeMod to the molecular conductivity for the first timestep if arbinit.
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
            if (use_lch4) then
               lake_raw(c) = raw(p) ! Pass out for calculating ground ch4 conductance
            end if
            ram1(p) = ram(p)       ! pass value to global variable
            ram1_lake(p) = ram1(p) ! for history

            ! Get derivative of fluxes with respect to ground temperature

            stftg3(p) = emg_lake*sb*tgbef(c)*tgbef(c)*tgbef(c)

            ! Changed surface temperature from t_lake(c,1) to tsur(c).
            ! Also adjusted so that if there are snow layers present, the top layer absorption
            ! from SNICAR is assigned to the surface skin.
            ax  = betaprime(c)*sabg(p) + emg_lake*forc_lwrad(t) + 3._r8*stftg3(p)*tgbef(c) &
                 + forc_rho(t)*cpair/rah(p)*thm(p) &
                 - htvp(c)*forc_rho(t)/raw(p)*(qsatg(c)-qsatgdT(c)*tgbef(c) - forc_q(t)) &
                 + tksur(c)*tsur(c)/dzsur(c)
            !Changed sabg(p) to betaprime(c)*sabg(p).
            bx  = 4._r8*stftg3(p) + forc_rho(t)*cpair/rah(p) &
                 + htvp(c)*forc_rho(t)/raw(p)*qsatgdT(c) + tksur(c)/dzsur(c)

            t_grnd(c) = ax/bx
            !prevent too large numerical step
            if(t_grnd(c)>tgbef(c)+20._r8)then
              t_grnd(c)=tgbef(c)+20._r8
            elseif(t_grnd(c)<tgbef(c)-20._r8)then
              t_grnd(c)=tgbef(c)-20._r8
            endif
            ! Update htvp
            if (t_grnd(c) > tfrz) then
               htvp(c) = hvap
            else
               htvp(c) = hsub
            end if

            ! Surface fluxes of momentum, sensible and latent heat
            ! using ground temperatures from previous time step

            eflx_sh_grnd(p) = forc_rho(t)*cpair*(t_grnd(c)-thm(p))/rah(p)
            qflx_evap_soi(p) = forc_rho(t)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c))-forc_q(t))/raw(p)

            ! Re-calculate saturated vapor pressure, specific humidity and their
            ! derivatives at lake surface

            call QSat(t_grnd(c), forc_pbot(t), eg, degdT, qsatg(c), qsatgdT(c))

            dth(p)=thm(p)-t_grnd(c)
            dqh(p)=forc_q(t)-qsatg(c)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)

            thvstar=tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
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
                  ! In this case fcrit is 22, not 100 in elm_varcon
               else
                  ! Fetch relationship from Vickers & Mahrt 1997
                  cur = cur0 + curm* exp( max( -(fetch(c)*grav/ustar(p)/ustar(p))**(1._r8/3._r8)/fcrit, &   ! Fetch-limited
                       -(lakedepth(c)*grav/ur(p)/ur(p))**0.5_r8 ) )           ! depth-limited
               end if


               kva = kva0 * (t_grnd(c)/kva0temp)**1.5_r8 * kva0pres/forc_pbot(t) ! kinematic viscosity of air
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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! If there is snow on the ground or lake is frozen and t_grnd > tfrz: reset t_grnd = tfrz.
         ! Re-evaluate ground fluxes.
         ! [ZMS 1/7/11] Only for resolved snow layers, as unresolved snow does not have a temperature state and
         ! can accumulate on unfrozen lakes in LakeHydrology; will be melted in LakeTemperature or bring lake top
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
            eflx_sh_grnd(p) = forc_rho(t)*cpair*(t_grnd(c)-thm(p))/rah(p)
            qflx_evap_soi(p) = forc_rho(t)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(t))/raw(p)
         else if ( (t_lake(c,1) > t_grnd(c) .and. t_grnd(c) > tdmax) .or. &
              (t_lake(c,1) < t_grnd(c) .and. t_lake(c,1) > tfrz .and. t_grnd(c) < tdmax) ) then
            ! Convective mixing will occur at surface
            t_grnd_temp = t_grnd(c)
            t_grnd(c) = t_lake(c,1)
            eflx_sh_grnd(p) = forc_rho(t)*cpair*(t_grnd(c)-thm(p))/rah(p)
            qflx_evap_soi(p) = forc_rho(t)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(t))/raw(p)
         end if

         ! Update htvp
         if (t_grnd(c) > tfrz) then
            htvp(c) = hvap
         else
            htvp(c) = hsub
         end if

         ! Net longwave from ground to atmosphere
         ! eflx_lwrad_out(p) = (1._r8-emg_lake)*forc_lwrad(c) + stftg3(p)*(-3._r8*tgbef(c)+4._r8*t_grnd(c))
         ! What is tgbef doing in this equation? Can't it be exact now? --Zack Subin, 4/14/09

         eflx_lwrad_out(p) = (1._r8-emg_lake)*forc_lwrad(t) + emg_lake*sb*t_grnd(c)**4._r8

         ! Ground heat flux

         eflx_soil_grnd(p) = sabg(p) + forc_lwrad(t) - eflx_lwrad_out(p) - &
              eflx_sh_grnd(p) - htvp(c)*qflx_evap_soi(p)
         ! The original code in Biogeophysiclake had a bug that calculated incorrect fluxes but conserved energy.
         ! This is kept as the full sabg (not just that absorbed at surface) so that the energy balance check will be correct.
         !This is the effective energy flux into the ground including the lake [and now snow in CLM 4] solar absorption
         !below the surface.  This also keeps the output FGR similar to non-lakes by including the light & heat flux.
         ! The variable eflx_gnet will be used to pass the actual heat flux
         !from the ground interface into the lake.

         taux(p) = -forc_rho(t)*forc_u(t)/ram(p)
         tauy(p) = -forc_rho(t)*forc_v(t)/ram(p)

         eflx_sh_tot(p)   = eflx_sh_grnd(p)
         qflx_evap_tot(p) = qflx_evap_soi(p)
         eflx_lh_tot(p)   = htvp(c)*qflx_evap_soi(p)
         eflx_lh_grnd(p)  = htvp(c)*qflx_evap_soi(p)

         ! 2 m height air temperature
         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

         ! 2 m height specific humidity
         q_ref2m(p) = forc_q(t) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)


         ! Energy residual used for melting snow
         ! Effectively moved to LakeTemp

         eflx_gnet(p) = betaprime(c) * sabg(p) + forc_lwrad(t) - (eflx_lwrad_out(p) + &
              eflx_sh_tot(p) + eflx_lh_tot(p))
         ! This is the actual heat flux from the ground interface into the lake, not including
         ! the light that penetrates the surface.

         !u2m = max(1.0_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))
         ! u2 often goes below 1 m/s; it seems like the only reason for this minimum is to
         ! keep it from being zero in the ks equation below; 0.1 m/s is a better limit for
         ! stable conditions --ZS
         u2m = max(0.1_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))

         ws(c) = 1.2e-03_r8 * u2m
         ks(c) = 6.6_r8*sqrt(abs(sin(grc_pp%lat(g))))*(u2m**(-1.84_r8))

         ! Update column roughness lengths and friction velocity
         z0mg_col(c) = z0mg(p)
         z0hg_col(c) = z0hg(p)
         z0qg_col(c) = z0qg(p)
         ust_lake(c) = ustar(p)

      end do

      ! The following are needed for global average on history tape.

      do fp = 1, num_lakep
         p = filter_lakep(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         
         t_veg(p) = forc_t(t)
         eflx_lwrad_net(p)  = eflx_lwrad_out(p) - forc_lwrad(t)
         qflx_prec_grnd(p) = forc_rain(t) + forc_snow(t)
         qflx_dirct_rain(p) = 0._r8
         qflx_leafdrip(p) = 0._r8

         ! Because they will be used in pft2col initialize here.
         ! This will be overwritten in LakeHydrology
         qflx_snwcp_ice(p) = 0._r8
         qflx_snwcp_liq(p) = 0._r8

      end do

    end associate

  end subroutine LakeFluxes

end module LakeFluxesMod
