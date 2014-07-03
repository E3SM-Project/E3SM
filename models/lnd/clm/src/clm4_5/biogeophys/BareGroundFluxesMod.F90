module BareGroundFluxesMod
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine BareGroundFluxes(bounds, num_nolakep, filter_nolakep)
    !
    ! !DESCRIPTION:
    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.
    !
    ! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l, a2l_downscaled_col
    use clm_varpar         , only : nlevgrnd
    use clm_varcon         , only : cpair, vkc, grav, denice, denh2o, istsoil
    use clm_varcon         , only : istcrop
    use shr_const_mod      , only : SHR_CONST_RGAS
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use QSatMod            , only : QSat
    use clm_varctl         , only : use_c13, use_c14, use_lch4
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_nolakep          ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(:)    ! pft filter for non-lake points
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: grnd_ch4_cond(:)  ! tracer conductance for boundary layer [m/s]
    integer, parameter  :: niters = 3      ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l                ! indices
    integer  :: filterp(bounds%endp-bounds%begp+1) ! pft filter for vegetated pfts
    integer  :: fn                               ! number of values in local pft filter
    integer  :: fp                               ! lake filter pft index
    integer  :: iter                             ! iteration index
    real(r8) :: zldis(bounds%begp:bounds%endp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(bounds%begp:bounds%endp)  ! displacement height [m]
    real(r8) :: zeta                             ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                               ! convective velocity [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                             ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(bounds%begp:bounds%endp)     ! Monin-Obukhov length (m)
    real(r8) :: ur(bounds%begp:bounds%endp)      ! wind speed at reference height [m/s]
    real(r8) :: um(bounds%begp:bounds%endp)      ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)   ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)   ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(bounds%begp:bounds%endp)   ! friction velocity [m/s]
    real(r8) :: tstar                            ! temperature scaling parameter
    real(r8) :: qstar                            ! moisture scaling parameter
    real(r8) :: thvstar                          ! virtual potential temperature scaling parameter
    real(r8) :: cf                               ! heat transfer coefficient from leaves [-]
    real(r8) :: ram                              ! aerodynamical resistance [s/m]
    real(r8) :: rah                              ! thermal resistance [s/m]
    real(r8) :: raw                              ! moisture resistance [s/m]
    real(r8) :: raih                             ! temporary variable [kg/m2/s]
    real(r8) :: raiw                             ! temporary variable [kg/m2/s]
    real(r8) :: fm(bounds%begp:bounds%endp)      ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_pft(bounds%begp:bounds%endp)
    real(r8) :: z0hg_pft(bounds%begp:bounds%endp)
    real(r8) :: z0qg_pft(bounds%begp:bounds%endp)
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(r8) :: www                    ! surface soil wetness [-]
    !------------------------------------------------------------------------------

   associate(& 
   t_soisno                       =>    ces%t_soisno                                     , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)                                           
   snl                            =>    cps%snl                                          , & ! Input:  [integer (:)]  number of snow layers                                                  
   t_h2osfc                       =>    ces%t_h2osfc                                     , & ! Input:  [real(r8) (:)]  surface water temperature                                             
   eflx_sh_snow                   =>    pef%eflx_sh_snow                                 , & ! Input:  [real(r8) (:)]  sensible heat flux from snow (W/m**2) [+ to atm]                      
   eflx_sh_soil                   =>    pef%eflx_sh_soil                                 , & ! Input:  [real(r8) (:)]  sensible heat flux from soil (W/m**2) [+ to atm]                      
   eflx_sh_h2osfc                 =>    pef%eflx_sh_h2osfc                               , & ! Input:  [real(r8) (:)]  sensible heat flux from soil (W/m**2) [+ to atm]                      
   qg_snow                        =>    cws%qg_snow                                      , & ! Input:  [real(r8) (:)]  specific humidity at snow surface [kg/kg]                             
   qg_soil                        =>    cws%qg_soil                                      , & ! Input:  [real(r8) (:)]  specific humidity at soil surface [kg/kg]                             
   qg_h2osfc                      =>    cws%qg_h2osfc                                    , & ! Input:  [real(r8) (:)]  specific humidity at h2osfc surface [kg/kg]                           
   qflx_ev_snow                   =>    pwf%qflx_ev_snow                                 , & ! Input:  [real(r8) (:)]  evaporation flux from snow (W/m**2) [+ to atm]                        
   qflx_ev_soil                   =>    pwf%qflx_ev_soil                                 , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]                        
   qflx_ev_h2osfc                 =>    pwf%qflx_ev_h2osfc                               , & ! Input:  [real(r8) (:)]  evaporation flux from h2osfc (W/m**2) [+ to atm]                      
   forc_u                         =>    clm_a2l%forc_u                                   , & ! Input:  [real(r8) (:)]  atmospheric wind speed in east direction (m/s)                        
   forc_v                         =>    clm_a2l%forc_v                                   , & ! Input:  [real(r8) (:)]  atmospheric wind speed in north direction (m/s)                       
   forc_th                        =>    a2l_downscaled_col%forc_th                       , & ! Input:  [real(r8) (:)]  atmospheric potential temperature (Kelvin)                            
   forc_t                         =>    a2l_downscaled_col%forc_t                        , & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)                                      
   forc_pbot                      =>    a2l_downscaled_col%forc_pbot                     , & ! Input:  [real(r8) (:)]  atmospheric pressure (Pa)                                             
   forc_rho                       =>    a2l_downscaled_col%forc_rho                      , & ! Input:  [real(r8) (:)]  density (kg/m**3)                                                     
   forc_q                         =>    a2l_downscaled_col%forc_q                        , & ! Input:  [real(r8) (:)]  atmospheric specific humidity (kg/kg)                                 
   frac_veg_nosno                 =>    pps%frac_veg_nosno                               , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1) [-]                
   dlrad                          =>    pef%dlrad                                        , & ! Output: [real(r8) (:)]  downward longwave radiation below the canopy [W/m2]                   
   ulrad                          =>    pef%ulrad                                        , & ! Output: [real(r8) (:)]  upward longwave radiation above the canopy [W/m2]                     
   t_grnd                         =>    ces%t_grnd                                       , & ! Input:  [real(r8) (:)]  ground surface temperature [K]                                        
   qg                             =>    cws%qg                                           , & ! Input:  [real(r8) (:)]  specific humidity at ground surface [kg/kg]                           
   z0mg_col                       =>    cps%z0mg                                         , & ! Input:  [real(r8) (:)]  roughness length, momentum [m]                                        
   z0hg_col                       =>    cps%z0hg                                         , & ! Input:  [real(r8) (:)]  roughness length, sensible heat [m]                                   
   z0qg_col                       =>    cps%z0qg                                         , & ! Input:  [real(r8) (:)]  roughness length, latent heat [m]                                     
   thv                            =>    ces%thv                                          , & ! Input:  [real(r8) (:)]  virtual potential temperature (kelvin)                                
   beta                           =>    cps%beta                                         , & ! Input:  [real(r8) (:)]  coefficient of conective velocity [-]                                 
   zii                            =>    cps%zii                                          , & ! Input:  [real(r8) (:)]  convective boundary height [m]                                        
   ram1                           =>    pps%ram1                                         , & ! Output: [real(r8) (:)]  aerodynamical resistance (s/m)                                        
   cgrnds                         =>    pef%cgrnds                                       , & ! Output: [real(r8) (:)]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]              
   cgrndl                         =>    pef%cgrndl                                       , & ! Output: [real(r8) (:)]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]               
   cgrnd                          =>    pef%cgrnd                                        , & ! Output: [real(r8) (:)]  deriv. of soil energy flux wrt to soil temp [w/m2/k]                  
   dqgdT                          =>    cws%dqgdT                                        , & ! Input:  [real(r8) (:)]  temperature derivative of "qg"                                        
   htvp                           =>    cps%htvp                                         , & ! Input:  [real(r8) (:)]  latent heat of evaporation (/sublimation) [J/kg]                      
   watsat                         =>    cps%watsat                                       , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)                      
   h2osoi_ice                     =>    cws%h2osoi_ice                                   , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                                    
   dz                             =>    cps%dz                                           , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                                     
   h2osoi_liq                     =>    cws%h2osoi_liq                                   , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                                                
   frac_sno                       =>    cps%frac_sno                                     , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)                           
   soilbeta                       =>    cws%soilbeta                                     , & ! Input:  [real(r8) (:)]  soil wetness relative to field capacity                               
   taux                           =>    pmf%taux                                         , & ! Output: [real(r8) (:)]  wind (shear) stress: e-w (kg/m/s**2)                                  
   tauy                           =>    pmf%tauy                                         , & ! Output: [real(r8) (:)]  wind (shear) stress: n-s (kg/m/s**2)                                  
   eflx_sh_grnd                   =>    pef%eflx_sh_grnd                                 , & ! Output: [real(r8) (:)]  sensible heat flux from ground (W/m**2) [+ to atm]                    
   eflx_sh_tot                    =>    pef%eflx_sh_tot                                  , & ! Output: [real(r8) (:)]  total sensible heat flux (W/m**2) [+ to atm]                          
   qflx_evap_soi                  =>    pwf%qflx_evap_soi                                , & ! Output: [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)                              
   qflx_evap_tot                  =>    pwf%qflx_evap_tot                                , & ! Output: [real(r8) (:)]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg                         
   t_ref2m                        =>    pes%t_ref2m                                      , & ! Output: [real(r8) (:)]  2 m height surface air temperature (Kelvin)                           
   q_ref2m                        =>    pes%q_ref2m                                      , & ! Output: [real(r8) (:)]  2 m height surface specific humidity (kg/kg)                          
   t_ref2m_r                      =>    pes%t_ref2m_r                                    , & ! Output: [real(r8) (:)]  Rural 2 m height surface air temperature (Kelvin)                     
   rh_ref2m_r                     =>    pes%rh_ref2m_r                                   , & ! Output: [real(r8) (:)]  Rural 2 m height surface relative humidity (%)                        
   plandunit                      =>    pft%landunit                                     , & ! Input:  [integer (:)]  pft's landunit index                                                   
   rh_ref2m                       =>    pes%rh_ref2m                                     , & ! Output: [real(r8) (:)]  2 m height surface relative humidity (%)                              
   t_veg                          =>    pes%t_veg                                        , & ! Output: [real(r8) (:)]  vegetation temperature (Kelvin)                                       
   thm                            =>    pes%thm                                          , & ! Input:  [real(r8) (:)]  intermediate variable (forc_t+0.0098*forc_hgt_t_pft)                  
   btran                          =>    pps%btran                                        , & ! Output: [real(r8) (:)]  transpiration wetness factor (0 to 1)                                 
   rssun                          =>    pps%rssun                                        , & ! Output: [real(r8) (:)]  sunlit stomatal resistance (s/m)                                      
   rssha                          =>    pps%rssha                                        , & ! Output: [real(r8) (:)]  shaded stomatal resistance (s/m)                                      
   rootr                          =>    pps%rootr                                        , & ! Output: [real(r8) (:,:)]  effective fraction of roots in each soil layer                      
   rresis                         =>    pps%rresis                                       , & ! Output: [real(r8) (:,:)]  root resistance by layer (0-1)  (nlevgrnd)                          
   psnsun                         =>    pcf%psnsun                                       , & ! Input:  [real(r8) (:)]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)                        
   psnsun_wc                      =>    pcf%psnsun_wc                                    , & ! Input:  [real(r8) (:)]  Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)        
   psnsun_wj                      =>    pcf%psnsun_wj                                    , & ! Input:  [real(r8) (:)]  RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)           
   psnsun_wp                      =>    pcf%psnsun_wp                                    , & ! Input:  [real(r8) (:)]  product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)        
   psnsha                         =>    pcf%psnsha                                       , & ! Input:  [real(r8) (:)]  shaded leaf photosynthesis (umol CO2 /m**2/ s)                        
   psnsha_wc                      =>    pcf%psnsha_wc                                    , & ! Input:  [real(r8) (:)]  Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)        
   psnsha_wj                      =>    pcf%psnsha_wj                                    , & ! Input:  [real(r8) (:)]  RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)           
   psnsha_wp                      =>    pcf%psnsha_wp                                    , & ! Input:  [real(r8) (:)]  product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)        
   fpsn                           =>    pcf%fpsn                                         , & ! Output: [real(r8) (:)]  photosynthesis (umol CO2 /m**2 /s)                                    
   fpsn_wc                        =>    pcf%fpsn_wc                                      , & ! Output: [real(r8) (:)]  Rubisco-limited photosynthesis (umol CO2 /m**2 /s)                    
   fpsn_wj                        =>    pcf%fpsn_wj                                      , & ! Output: [real(r8) (:)]  RuBP-limited photosynthesis (umol CO2 /m**2 /s)                       
   fpsn_wp                        =>    pcf%fpsn_wp                                      , & ! Output: [real(r8) (:)]  product-limited photosynthesis (umol CO2 /m**2 /s)                    
   forc_hgt_u_pft                 =>    pps%forc_hgt_u_pft                               , &
   begp                           =>    bounds%begp                                      , &
   endp                           =>    bounds%endp                                        &
   )

    grnd_ch4_cond  => pps%grnd_ch4_cond

    ! Filter pfts where frac_veg_nosno is zero

    fn = 0
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (frac_veg_nosno(p) == 0) then
          fn = fn + 1
          filterp(fn) = p
       end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1, fn
       p = filterp(f)
       c = pft%column(p)
       g = pft%gridcell(p)

       ! Initialization variables

       displa(p) = 0._r8
       dlrad(p)  = 0._r8
       ulrad(p)  = 0._r8

       ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p) = thm(p)-t_grnd(c)
       dqh(p) = forc_q(c) - qg(c)
       dthv = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p)

       ! Copy column roughness to local pft-level arrays

       z0mg_pft(p) = z0mg_col(c)
       z0hg_pft(p) = z0hg_col(c)
       z0qg_pft(p) = z0qg_col(c)

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_pft(p), um(p), obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1, niters

       call FrictionVelocity(begp, endp, fn, filterp, &
                             displa(begp:endp), z0mg_pft(begp:endp), z0hg_pft(begp:endp), z0qg_pft(begp:endp), &
                             obu(begp:endp), iter, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
                             temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp))

       do f = 1, fn
          p = filterp(f)
          c = pft%column(p)
          g = pft%gridcell(p)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)
          z0hg_pft(p) = z0mg_pft(p)/exp(0.13_r8 * (ustar(p)*z0mg_pft(p)/1.5e-5_r8)**0.45_r8)
          z0qg_pft(p) = z0hg_pft(p)
          thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
          zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0._r8) then                   !stable
             zeta = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                                      !unstable
             zeta = max(-100._r8,min(zeta,-0.01_r8))
             wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p) + wc*wc)
          end if
          obu(p) = zldis(p)/zeta
       end do

    end do ! end stability iteration

     do j = 1, nlevgrnd
       do f = 1, fn
          p = filterp(f)
          rootr(p,j) = 0._r8
          rresis(p,j) = 0._r8
        end do
     end do

    do f = 1, fn
       p = filterp(f)
       c = pft%column(p)
       g = pft%gridcell(p)
       l = pft%landunit(p)

       ! Determine aerodynamic resistances

       ram     = 1._r8/(ustar(p)*ustar(p)/um(p))
       rah     = 1._r8/(temp1(p)*ustar(p))
       raw     = 1._r8/(temp2(p)*ustar(p))
       raih    = forc_rho(c)*cpair/rah
       if (use_lch4) then
          grnd_ch4_cond(p) = 1._r8/raw
       end if

       ! Soil evaporation resistance
       www     = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
       www     = min(max(www,0.0_r8),1._r8)

       !changed by K.Sakaguchi. Soilbeta is used for evaporation
       if (dqh(p) .gt. 0._r8) then   !dew  (beta is not applied, just like rsoil used to be)
          raiw    = forc_rho(c)/(raw)
       else
       ! Lee and Pielke 1992 beta is applied
          raiw    = soilbeta(c)*forc_rho(c)/(raw)
       end if

       ram1(p) = ram  !pass value to global variable

       ! Output to pft-level data structures
       ! Derivative of fluxes with respect to ground temperature

       cgrnds(p) = raih
       cgrndl(p) = raiw*dqgdT(c)
       cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       taux(p)          = -forc_rho(c)*forc_u(g)/ram
       tauy(p)          = -forc_rho(c)*forc_v(g)/ram
       eflx_sh_grnd(p)  = -raih*dth(p)
       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       ! compute sensible heat fluxes individually
       eflx_sh_snow(p)  = -raih*(thm(p)-t_soisno(c,snl(c)+1))
       eflx_sh_soil(p)  = -raih*(thm(p)-t_soisno(c,1))
       eflx_sh_h2osfc(p)  = -raih*(thm(p)-t_h2osfc(c))
       qflx_evap_soi(p) = -raiw*dqh(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)
       ! compute latent heat fluxes individually
       qflx_ev_snow(p) = -raiw*(forc_q(c) - qg_snow(c))
       qflx_ev_soil(p) = -raiw*(forc_q(c) - qg_soil(c))
       qflx_ev_h2osfc(p) = -raiw*(forc_q(c) - qg_h2osfc(c))

       ! 2 m height air temperature

       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity

       q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

       ! 2 m height relative humidity
                                                                                
       call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

       rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
         rh_ref2m_r(p) = rh_ref2m(p)
         t_ref2m_r(p) = t_ref2m(p)
       end if

       ! Variables needed by history tape

       t_veg(p) = forc_t(c)
       btran(p) = 0._r8
       cf = forc_pbot(c)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
       rssun(p) = 1._r8/1.e15_r8 * cf
       rssha(p) = 1._r8/1.e15_r8 * cf

       ! Add the following to avoid NaN

       psnsun(p) = 0._r8
       psnsun_wc(p) = 0._r8
       psnsun_wj(p) = 0._r8
       psnsun_wp(p) = 0._r8
       psnsha(p) = 0._r8
       psnsha_wc(p) = 0._r8
       psnsha_wj(p) = 0._r8
       psnsha_wp(p) = 0._r8
       fpsn(p) = 0._r8
       fpsn_wc(p) = 0._r8
       fpsn_wj(p) = 0._r8
       fpsn_wp(p) = 0._r8
       ! adding code for isotopes, 8/17/05, PET

     if ( use_c13 ) then
        pps%alphapsnsun(p) = 0._r8
        pps%alphapsnsha(p) = 0._r8
        pepv%rc13_canair(p) = 0._r8
        pepv%rc13_psnsun(p) = 0._r8
        pepv%rc13_psnsha(p) = 0._r8
        pc13f%psnsun(p) = 0._r8
        pc13f%psnsha(p) = 0._r8
        pc13f%fpsn(p) = 0._r8
     endif

     if ( use_c14 ) then
        pepv%rc14_atm(p) = 0._r8
        ! pepv%rc14_canair(p) = 0._r8
        ! pepv%rc14_psnsun(p) = 0._r8
        ! pepv%rc14_psnsha(p) = 0._r8
        pc14f%psnsun(p) = 0._r8
        pc14f%psnsha(p) = 0._r8
        pc14f%fpsn(p) = 0._r8
     endif

    end do

    end associate 
   end subroutine BareGroundFluxes

end module BareGroundFluxesMod
