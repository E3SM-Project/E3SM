module  PhotosynthesisMod

#include "shr_assert.h"
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf photosynthesis and stomatal conductance calculation as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog, use_c13, use_c14, use_cn, use_cndv, use_ed
  use clm_varpar          , only : nlevcan
  use clm_varcon          , only : namep
  use decompMod           , only : bounds_type
  use QuadraticMod        , only : quadratic
  use EcophysConType      , only : ecophyscon
  use atm2lndType         , only : atm2lnd_type
  use CNStateType         , only : cnstate_type
  use CanopyStateType     , only : canopystate_type
  use TemperatureType     , only : temperature_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use PhotosynthesisType  , only : photosyns_type
  use PatchType           , only : pft
  use CNAllocationMod     , only : nu_com_leaf_physiology
  use PhosphorusStateType , only : phosphorusstate_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use clm_varctl          , only : cnallocate_carbon_only
  use clm_varctl          , only : cnallocate_carbonnitrogen_only
  use clm_varctl          , only : cnallocate_carbonphosphorus_only
  use clm_varctl          , only : iulog
  use pftvarcon           , only : noveg
    
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Photosynthesis       ! Leaf stomatal resistance and leaf photosynthesis
  public :: PhotosynthesisTotal  ! Determine of total photosynthesis
  public :: Fractionation        ! C13 fractionation during photosynthesis 

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: hybrid         ! hybrid solver for ci
  private :: ci_func        ! ci function
  private :: brent          ! brent solver for root of a single variable function
  private :: ft             ! photosynthesis temperature response
  private :: fth            ! photosynthesis temperature inhibition
  private :: fth25          ! scaling factor for photosynthesis temperature inhibition
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine Photosynthesis ( bounds, fn, filterp, &
       esat_tv, eair, oair, cair, rb, btran, &
       dayl_factor, atm2lnd_vars, temperature_vars, surfalb_vars, solarabs_vars, &
       canopystate_vars, photosyns_vars, nitrogenstate_vars,phosphorusstate_vars, phase)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    !
    ! !USES:
    use clm_varcon     , only : rgas, tfrz
    use clm_varctl     , only : cnallocate_carbon_only 
    use pftvarcon      , only : nbrdlf_dcd_tmp_shrub, nsoybean, nsoybeanirrig, npcropmin
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                         
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    real(r8)               , intent(in)    :: esat_tv( bounds%begp: )        ! saturation vapor pressure at t_veg (Pa) [pft]
    real(r8)               , intent(in)    :: eair( bounds%begp: )           ! vapor pressure of canopy air (Pa) [pft]
    real(r8)               , intent(in)    :: oair( bounds%begp: )           ! Atmospheric O2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: cair( bounds%begp: )           ! Atmospheric CO2 partial pressure (Pa) [pft]
    real(r8)               , intent(in)    :: rb( bounds%begp: )             ! boundary layer resistance (s/m) [pft]
    real(r8)               , intent(in)    :: btran( bounds%begp: )          ! transpiration wetness factor (0 to 1) [pft]
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: )    ! scalar (0-1) for daylength
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    type(surfalb_type)     , intent(in)    :: surfalb_vars
    type(solarabs_type)    , intent(in)    :: solarabs_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    type(nitrogenstate_type)  , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type), intent(inout) :: phosphorusstate_vars
    character(len=*)       , intent(in)    :: phase                          ! 'sun' or 'sha'
    
    !
    ! !LOCAL VARIABLES:
    !
    ! Leaf photosynthesis parameters
    real(r8) :: jmax_z(bounds%begp:bounds%endp,nlevcan)  ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: lnc(bounds%begp:bounds%endp)   ! leaf N concentration (gN leaf/m^2)
    real(r8) :: bbbopt(bounds%begp:bounds%endp)! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(r8) :: mbbopt(bounds%begp:bounds%endp)! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    real(r8) :: kn(bounds%begp:bounds%endp)    ! leaf nitrogen decay coefficient
    real(r8) :: vcmax25top     ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top      ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top       ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top       ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top        ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: vcmax25        ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25         ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25          ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25          ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25           ! leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kc25           ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25           ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: cp25           ! CO2 compensation point at 25C (Pa)

    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: tpuha          ! activation energy for tpu (J/mol)
    real(r8) :: lmrha          ! activation energy for lmr (J/mol)
    real(r8) :: kcha           ! activation energy for kc (J/mol)
    real(r8) :: koha           ! activation energy for ko (J/mol)
    real(r8) :: cpha           ! activation energy for cp (J/mol)

    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: tpuhd          ! deactivation energy for tpu (J/mol)
    real(r8) :: lmrhd          ! deactivation energy for lmr (J/mol)

    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: tpuse          ! entropy term for tpu (J/mol/K)
    real(r8) :: lmrse          ! entropy term for lmr (J/mol/K)

    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: tpuc           ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: lmrc           ! scaling factor for high temperature inhibition (25 C = 1.0)

    real(r8) :: fnps           ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii     ! empirical curvature parameter for electron transport rate

    real(r8) :: theta_ip          ! empirical curvature parameter for ap photosynthesis co-limitation

    ! Other
    integer  :: f,p,c,iv          ! indices
    real(r8) :: cf                ! s m**2/umol -> s/m
    real(r8) :: rsmax0            ! maximum stomatal resistance [s/m]
    real(r8) :: gb                ! leaf boundary layer conductance (m/s)
    real(r8) :: cs                ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: gs                ! leaf stomatal conductance (m/s)
    real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
    real(r8) :: sco               ! relative specificity of rubisco
    real(r8) :: ft                ! photosynthesis temperature response (statement function)
    real(r8) :: fth               ! photosynthesis temperature inhibition (statement function)
    real(r8) :: fth25             ! ccaling factor for photosynthesis temperature inhibition (statement function)
    real(r8) :: tl                ! leaf temperature in photosynthesis temperature function (K)
    real(r8) :: ha                ! activation energy in photosynthesis temperature function (J/mol)
    real(r8) :: hd                ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8) :: se                ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8) :: scaleFactor       ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: ciold             ! previous value of Ci for convergence check
    real(r8) :: gs_mol_err        ! gs_mol for error check
    real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    real(r8) :: ceair             ! vapor pressure of air, constrained (Pa)
    real(r8) :: fnr               ! (gRubisco/gN in Rubisco)
    real(r8) :: act25             ! (umol/mgRubisco/min) Rubisco activity at 25 C
    integer  :: niter             ! iteration loop index
    real(r8) :: nscaler           ! leaf nitrogen scaling coefficient

    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)

    real(r8) :: psn_wc_z(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to psn_z (umol CO2/m**2/s)

    real(r8) :: psncan            ! canopy sum of psn_z
    real(r8) :: psncan_wc         ! canopy sum of psn_wc_z
    real(r8) :: psncan_wj         ! canopy sum of psn_wj_z
    real(r8) :: psncan_wp         ! canopy sum of psn_wp_z
    real(r8) :: lmrcan            ! canopy sum of lmr_z
    real(r8) :: gscan             ! canopy sum of leaf conductance
    real(r8) :: laican            ! canopy sum of lai_z
    real(r8) :: rh_can

    real(r8) , pointer :: lai_z       (:,:)    
    real(r8) , pointer :: par_z       (:,:)    
    real(r8) , pointer :: vcmaxcint   (:)  
    real(r8) , pointer :: alphapsn    (:)   
    real(r8) , pointer :: psn         (:)        
    real(r8) , pointer :: psn_wc      (:)     
    real(r8) , pointer :: psn_wj      (:)     
    real(r8) , pointer :: psn_wp      (:)     
    real(r8) , pointer :: psn_z       (:,:)    
    real(r8) , pointer :: lmr         (:)        
    real(r8) , pointer :: lmr_z       (:,:)    
    real(r8) , pointer :: rs          (:)         
    real(r8) , pointer :: rs_z        (:,:)     
    real(r8) , pointer :: ci_z        (:,:)     
    real(r8) , pointer :: alphapsnsun (:)  
    real(r8) , pointer :: alphapsnsha (:)
    
    real(r8) :: lpc(bounds%begp:bounds%endp)   ! leaf P concentration (gP leaf/m^2)
    real(r8) :: sum_nscaler              
    real(r8) :: total_lai
    !------------------------------------------------------------------------------

    ! Temperature and soil water response functions

    ft(tl,ha) = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )
    fth(tl,hd,se,scaleFactor) = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )
    fth25(hd,se) = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    ! Enforce expected array sizes

    SHR_ASSERT_ALL((ubound(esat_tv)     == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(eair)        == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(oair)        == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(cair)        == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rb)          == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(btran)       == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dayl_factor) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))

    associate(                                                       & 
         c3psn         => ecophyscon%c3psn                         , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3                              
         leafcn        => ecophyscon%leafcn                        , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)                                                      
         flnr          => ecophyscon%flnr                          , & ! Input:  [real(r8) (:)   ]  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)       
         fnitr         => ecophyscon%fnitr                         , & ! Input:  [real(r8) (:)   ]  foliage nitrogen limitation factor (-)                                
         slatop        => ecophyscon%slatop                        , & ! Input:  [real(r8) (:)   ]  specific leaf area at top of canopy, projected area basis [m^2/gC]    

         forc_pbot     => atm2lnd_vars%forc_pbot_downscaled_col    , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             

         t_veg         => temperature_vars%t_veg_patch             , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t10           => temperature_vars%t_a10_patch             , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)                        
         tgcm          => temperature_vars%thm_patch               , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)                     

         nrad          => surfalb_vars%nrad_patch                  , & ! Input:  [integer  (:)   ]  pft number of canopy layers, above snow for radiative transfer  
         tlai_z        => surfalb_vars%tlai_z_patch                , & ! Input:  [real(r8) (:,:) ]  pft total leaf area index for canopy layer                              

         c3flag        => photosyns_vars%c3flag_patch              , & ! Output: [logical  (:)   ]  true if C3 and false if C4                                             
         ac            => photosyns_vars%ac_patch                  , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)              
         aj            => photosyns_vars%aj_patch                  , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)                 
         ap            => photosyns_vars%ap_patch                  , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag            => photosyns_vars%ag_patch                  , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)              
         an            => photosyns_vars%an_patch                  , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)                           
         gb_mol        => photosyns_vars%gb_mol_patch              , & ! Output: [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)                     
         gs_mol        => photosyns_vars%gs_mol_patch              , & ! Output: [real(r8) (:,:) ]  leaf stomatal conductance (umol H2O/m**2/s)                         
         vcmax_z       => photosyns_vars%vcmax_z_patch             , & ! Output: [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)                     
         cp            => photosyns_vars%cp_patch                  , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)                                           
         kc            => photosyns_vars%kc_patch                  , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)                                
         ko            => photosyns_vars%ko_patch                  , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)                                 
         qe            => photosyns_vars%qe_patch                  , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)          
         tpu_z         => photosyns_vars%tpu_z_patch               , & ! Output: [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)                 
         kp_z          => photosyns_vars%kp_z_patch                , & ! Output: [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)                     
         theta_cj      => photosyns_vars%theta_cj_patch            , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation 
         bbb           => photosyns_vars%bbb_patch                 , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)                 
         mbb           => photosyns_vars%mbb_patch                 , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship           
         rh_leaf       => photosyns_vars%rh_leaf_patch             , & ! Output: [real(r8) (:)   ]  fractional humidity at leaf surface (dimensionless)                   
         
         leafn         => nitrogenstate_vars%leafn_patch           , &
         leafn_storage => nitrogenstate_vars%leafn_storage_patch   , &
         leafn_xfer    => nitrogenstate_vars%leafn_xfer_patch      , &
         leafp         => phosphorusstate_vars%leafp_patch         , &
         leafp_storage => phosphorusstate_vars%leafp_storage_patch , &
         leafp_xfer    => phosphorusstate_vars%leafp_xfer_patch    , &
         i_vcmax       => ecophyscon%i_vc                          , &
         s_vcmax       => ecophyscon%s_vc                            &
         )
      
      if (phase == 'sun') then
         par_z     =>    solarabs_vars%parsun_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)                 
         lai_z     =>    canopystate_vars%laisun_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded                  
         vcmaxcint =>    surfalb_vars%vcmaxcintsun_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient                                     
         alphapsn  =>    photosyns_vars%alphapsnsun_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()                                   
         ci_z      =>    photosyns_vars%cisun_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)                                         
         rs        =>    photosyns_vars%rssun_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)                                        
         rs_z      =>    photosyns_vars%rssun_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)                        
         lmr       =>    photosyns_vars%lmrsun_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)                   
         lmr_z     =>    photosyns_vars%lmrsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)   
         psn       =>    photosyns_vars%psnsun_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]                 
         psn_z     =>    photosyns_vars%psnsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +] 
         psn_wc    =>    photosyns_vars%psnsun_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +] 
         psn_wj    =>    photosyns_vars%psnsun_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]    
         psn_wp    =>    photosyns_vars%psnsun_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +] 
      else if (phase == 'sha') then
         par_z     =>    solarabs_vars%parsha_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)                 
         lai_z     =>    canopystate_vars%laisha_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded                  
         vcmaxcint =>    surfalb_vars%vcmaxcintsha_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient                                    
         alphapsn  =>    photosyns_vars%alphapsnsha_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
         ci_z      =>    photosyns_vars%cisha_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)                                         
         rs        =>    photosyns_vars%rssha_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)                                        
         rs_z      =>    photosyns_vars%rssha_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)                        
         lmr       =>    photosyns_vars%lmrsha_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)                   
         lmr_z     =>    photosyns_vars%lmrsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)   
         psn       =>    photosyns_vars%psnsha_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]                 
         psn_z     =>    photosyns_vars%psnsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +] 
         psn_wc    =>    photosyns_vars%psnsha_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +] 
         psn_wj    =>    photosyns_vars%psnsha_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]    
         psn_wp    =>    photosyns_vars%psnsha_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +] 
      end if

      !==============================================================================!
      ! Photosynthesis and stomatal conductance parameters, from:
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      !==============================================================================!

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593

      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      do f = 1, fn
         p = filterp(f)
         c = pft%column(p)

         ! vcmax25 parameters, from CN

         fnr   = ecophyscon%fnr(pft%itype(p))   !7.16_r8
         act25 = ecophyscon%act25(pft%itype(p)) !3.6_r8   !umol/mgRubisco/min
         ! Convert rubisco activity units from umol/mgRubisco/min ->
         ! umol/gRubisco/s
         act25 = act25 * 1000.0_r8 / 60.0_r8

         ! Activation energy, from:
         ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !  Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
         ! except TPU from: Harley et al (1992) Plant, Cell and Environment
         ! 15:271-282

         kcha    = ecophyscon%kcha(pft%itype(p)) !79430._r8
         koha    = ecophyscon%koha(pft%itype(p)) !36380._r8
         cpha    = ecophyscon%cpha(pft%itype(p)) !37830._r8
         vcmaxha = ecophyscon%vcmaxha(pft%itype(p)) !72000._r8
         jmaxha  = ecophyscon%jmaxha(pft%itype(p)) !50000._r8
         tpuha   = ecophyscon%tpuha(pft%itype(p))  !72000._r8
         lmrha   = ecophyscon%lmrha(pft%itype(p))  !46390._r8

         ! High temperature deactivation, from:
         ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
         ! The factor "c" scales the deactivation to a value of 1.0 at 25C

         vcmaxhd = ecophyscon%vcmaxhd(pft%itype(p)) !200000._r8
         jmaxhd  = ecophyscon%jmaxhd(pft%itype(p))  !200000._r8
         tpuhd   = ecophyscon%tpuhd(pft%itype(p))   !200000._r8
         lmrhd   = ecophyscon%lmrhd(pft%itype(p))   !150650._r8
         lmrse   = ecophyscon%lmrse(pft%itype(p))   !490._r8
         lmrc    = fth25 (lmrhd, lmrse)

         ! C3 or C4 photosynthesis logical variable

         if (nint(c3psn(pft%itype(p))) == 1) then
            c3flag(p) = .true. 
         else if (nint(c3psn(pft%itype(p))) == 0) then
            c3flag(p) = .false.
         end if

         ! C3 and C4 dependent parameters

         if (c3flag(p)) then
            qe(p)       = ecophyscon%qe(pft%itype(p))       !0._r8
            theta_cj(p) = ecophyscon%theta_cj(pft%itype(p)) !0.98_r8
            bbbopt(p)   = ecophyscon%bbbopt(pft%itype(p))   !10000._r8
            mbbopt(p)   = ecophyscon%mbbopt(pft%itype(p))   !9._r8
         else
            qe(p)       = ecophyscon%qe(pft%itype(p))       !0.05_r8
            theta_cj(p) = ecophyscon%theta_cj(pft%itype(p)) !0.80_r8
            bbbopt(p)   = ecophyscon%bbbopt(pft%itype(p))   !40000._r8
            mbbopt(p)   = ecophyscon%mbbopt(pft%itype(p))   !4._r8
         end if

         ! Soil water stress applied to Ball-Berry parameters

         bbb(p) = max (bbbopt(p)*btran(p), 1._r8)
         mbb(p) = mbbopt(p)

         ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !
         !       kc25 = 404.9 umol/mol
         !       ko25 = 278.4 mmol/mol
         !       cp25 = 42.75 umol/mol
         !
         ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
         ! cp to account for variation in O2 using cp = 0.5 O2 / sco
         !

         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(c)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(c)
         sco  = 0.5_r8 * 0.209_r8 / (42.75_r8 / 1.e06_r8)
         cp25 = 0.5_r8 * oair(p) / sco

         kc(p) = kc25 * ft(t_veg(p), kcha)
         ko(p) = ko25 * ft(t_veg(p), koha)
         cp(p) = cp25 * ft(t_veg(p), cpha)

      end do

      ! Multi-layer parameters scaled by leaf nitrogen profile.
      ! Loop through each canopy layer to calculate nitrogen profile using
      ! cumulative lai at the midpoint of the layer
      
      do f = 1, fn
         p = filterp(f)
         if ( .not. nu_com_leaf_physiology) then
            ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
            lnc(p) = 1._r8 / (slatop(pft%itype(p)) * leafcn(pft%itype(p)))

            ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
            vcmax25top = lnc(p) * flnr(pft%itype(p)) * fnr * act25 * dayl_factor(p)
            if (.not. use_cn) then
               vcmax25top = vcmax25top * fnitr(pft%itype(p))
            else
               if ( CNAllocate_Carbon_only() ) vcmax25top = vcmax25top * fnitr(pft%itype(p))
            end if

            ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
            ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.
            jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

         else
         
            ! leaf level nutrient control on photosynthesis rate added by Q. Zhu Aug 2015
            
            if ( CNAllocate_Carbon_only() .or. cnallocate_carbonphosphorus_only()) then

               lnc(p) = 1._r8 / (slatop(pft%itype(p)) * leafcn(pft%itype(p)))
               vcmax25top = lnc(p) * flnr(pft%itype(p)) * fnr * act25 * dayl_factor(p)
               vcmax25top = vcmax25top * fnitr(pft%itype(p))
               jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

            else if ( cnallocate_carbonnitrogen_only() ) then ! only N control, from Kattge 2009 Global Change Biology 15 (4), 976-991

               ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
               sum_nscaler = 0.0_r8                                                       
               laican      = 0.0_r8
               total_lai   = 0.0_r8                                                      

               do iv = 1, nrad(p)
                  if (iv == 1) then
                     laican = 0.5_r8 * tlai_z(p,iv)
                  else
                     laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
                  end if
                  total_lai = total_lai + tlai_z(p,iv)
                  ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
                  ! profile. If sun/shade big leaf code, use canopy integrated factor.
                  if (nlevcan == 1) then                                               
                     nscaler = 1.0_r8                                                  
                  else if (nlevcan > 1) then                                           
                     nscaler = exp(-kn(p) * laican)                                    
                  end if
                  sum_nscaler = sum_nscaler + nscaler                                  
               end do

               if (total_lai > 0.0_r8 .and. sum_nscaler > 0.0_r8) then
                  ! dividing by LAI to convert total leaf nitrogen
                  ! from m2 ground to m2 leaf; dividing by sum_nscaler to
                  ! convert total leaf N to leaf N at canopy top
                  lnc(p) = leafn(p) / (total_lai * sum_nscaler)
                  lnc(p) = min(max(lnc(p),0.0_r8),3.0_r8) ! based on TRY database, doi: 10.1002/ece3.1173
               else                                                                    
                  lnc(p) = 0.0_r8                                                      
               end if

               vcmax25top = (i_vcmax(pft%itype(p)) + s_vcmax(pft%itype(p)) * lnc(p)) * dayl_factor(p) 
               jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

            else

               ! nu_com_leaf_physiology is true, vcmax25, jmax25 is derived from leafn, leafp concentration
               ! Anthony Walker 2014 DOI: 10.1002/ece3.1173

               if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
                  ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
                  sum_nscaler = 0.0_r8                                                       
                  laican      = 0.0_r8
                  total_lai   = 0.0_r8                                                      

                  do iv = 1, nrad(p)
                     if (iv == 1) then
                        laican = 0.5_r8 * tlai_z(p,iv)
                     else
                        laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
                     end if
                     total_lai = total_lai + tlai_z(p,iv)
                     ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
                     ! profile. If sun/shade big leaf code, use canopy integrated factor.
                     if (nlevcan == 1) then                                               
                        nscaler = 1.0_r8                                                  
                     else if (nlevcan > 1) then                                           
                        nscaler = exp(-kn(p) * laican)                                    
                     end if
                     sum_nscaler = sum_nscaler + nscaler                                  
                  end do

                  if (total_lai > 0.0_r8 .and. sum_nscaler > 0.0_r8) then
                     ! dividing by LAI to convert total leaf nitrogen
                     ! from m2 ground to m2 leaf; dividing by sum_nscaler to
                     ! convert total leaf N to leaf N at canopy top
                     lnc(p) = leafn(p) / (total_lai * sum_nscaler)
                     lpc(p) = leafp(p) / (total_lai * sum_nscaler)
                     lnc(p) = min(max(lnc(p),0.0_r8),3.0_r8) ! based on TRY database, doi: 10.1002/ece3.1173
                     lpc(p) = min(max(lpc(p),0.0_r8),0.5_r8) ! based on TRY database, doi: 10.1002/ece3.1173
                  else
                     lnc(p) = 0.0_r8
                     lpc(p) = 0.0_r8
                  end if

                  if (lnc(p) >= 0.1_r8 .and. lnc(p) <=3.0_r8 .and. lpc(p) >= 0.05_r8 .and. lpc(p) <= 0.5_r8) then
                     vcmax25top = exp(3.946_r8 + 0.921_r8*log(lnc(p)) + 0.121_r8*log(lpc(p)) + 0.282_r8*log(lnc(p))*log(lpc(p))) * dayl_factor(p)
                     jmax25top = exp(1.246_r8 + 0.886_r8*log(vcmax25top) + 0.089_r8*log(lpc(p))) * dayl_factor(p)
                  else if (lnc(p) < 0.1_r8 .or. lpc(p) < 0.05_r8) then
                     vcmax25top = 10.0_r8
                     jmax25top  = 10.0_r8
                  else
                     vcmax25top = 150.0_r8
                     jmax25top  = 250.0_r8
                  end if

               else
                  lnc(p)     = 0.0_r8
                  vcmax25top = 0.0_r8
                  jmax25top  = 0.0_r8
               end if
            end if
         end if
                 
         tpu25top  = 0.167_r8 * vcmax25top
         kp25top   = 20000._r8 * vcmax25top

         ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
         ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
         ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
         ! But not used as defined here if using sun/shade big leaf code. Instead,
         ! will use canopy integrated scaling factors from SurfaceAlbedo.

         if (dayl_factor(p) .eq. 0._r8) then
            kn(p) =  0._r8
         else
            kn(p) = exp(0.00963_r8 * vcmax25top/dayl_factor(p) - 2.43_r8)
         end if

         if (use_cn) then
            ! Leaf maintenance respiration to match the base rate used in CN
            ! but with the new temperature functions for C3 and C4 plants.
            !
            ! Base rate for maintenance respiration is from:
            ! M. Ryan, 1991. Effects of climate change on plant respiration.
            ! Ecological Applications, 1(2), 157-167.
            ! Original expression is br = 0.0106 molC/(molN h)
            ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
            !
            ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
            !
            ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
            ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
            !
            ! Then scale this value at the top of the canopy for canopy depth

            lmr25top = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
            lmr25top = lmr25top * lnc(p) / 12.e-06_r8

         else
            ! Leaf maintenance respiration in proportion to vcmax25top

            if (c3flag(p)) then
               lmr25top = vcmax25top * 0.015_r8
            else
               lmr25top = vcmax25top * 0.025_r8
            end if
         end if

         ! Loop through canopy layers (above snow). Respiration needs to be
         ! calculated every timestep. Others are calculated only if daytime

         laican = 0._r8
         do iv = 1, nrad(p)

            ! Cumulative lai at middle of layer

            if (iv == 1) then
               laican = 0.5_r8 * tlai_z(p,iv)
            else
               laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
            end if

            ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
            ! profile. If sun/shade big leaf code, use canopy integrated factor.

            if (nlevcan == 1) then
               nscaler = vcmaxcint(p)
               if (nu_com_leaf_physiology) nscaler = 1
            else if (nlevcan > 1) then
               nscaler = exp(-kn(p) * laican)
            end if

            ! Maintenance respiration

            lmr25 = lmr25top * nscaler
            if (c3flag(p)) then
               lmr_z(p,iv) = lmr25 * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
            else
               lmr_z(p,iv) = lmr25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z(p,iv) = lmr_z(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
            end if

            if (par_z(p,iv) <= 0._r8) then           ! night time

               vcmax_z(p,iv) = 0._r8
               jmax_z(p,iv) = 0._r8
               tpu_z(p,iv) = 0._r8
               kp_z(p,iv) = 0._r8

               if ( use_c13 ) then
                  alphapsn(p) = 1._r8
               end if

            else                                     ! day time

               vcmax25 = vcmax25top * nscaler
               jmax25 = jmax25top * nscaler
               tpu25 = tpu25top * nscaler
               kp25 = kp25top * nscaler

               ! Adjust for temperature

               vcmaxse = 668.39_r8 - 1.07_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               jmaxse  = 659.70_r8 - 0.75_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               tpuse = vcmaxse
               vcmaxc = fth25 (vcmaxhd, vcmaxse)
               jmaxc  = fth25 (jmaxhd, jmaxse)
               tpuc   = fth25 (tpuhd, tpuse)
               vcmax_z(p,iv) = vcmax25 * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,iv) = jmax25 * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,iv) = tpu25 * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

               if (.not. c3flag(p)) then
                  vcmax_z(p,iv) = vcmax25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
               end if

               kp_z(p,iv) = kp25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)

            end if

            ! Adjust for soil water

            vcmax_z(p,iv) = vcmax_z(p,iv) * btran(p)
            lmr_z(p,iv) = lmr_z(p,iv) * btran(p)

         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Leaf-level photosynthesis and stomatal conductance
      !==============================================================================!

      rsmax0 = 2.e4_r8

      do f = 1, fn
         p = filterp(f)
         c = pft%column(p)

         ! Leaf boundary layer conductance, umol/m**2/s

         cf = forc_pbot(c)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
         gb = 1._r8/rb(p)
         gb_mol(p) = gb * cf

         ! Loop through canopy layers (above snow). Only do calculations if daytime

         do iv = 1, nrad(p)

            if (par_z(p,iv) <= 0._r8) then           ! night time

               ac(p,iv) = 0._r8
               aj(p,iv) = 0._r8
               ap(p,iv) = 0._r8
               ag(p,iv) = 0._r8
               an(p,iv) = ag(p,iv) - lmr_z(p,iv)
               psn_z(p,iv) = 0._r8
               psn_wc_z(p,iv) = 0._r8
               psn_wj_z(p,iv) = 0._r8
               psn_wp_z(p,iv) = 0._r8
               rs_z(p,iv) = min(rsmax0, 1._r8/bbb(p) * cf)
               ci_z(p,iv) = 0._r8
               rh_leaf(p) = 0._r8

            else                                     ! day time

               !now the constraint is no longer needed, Jinyun Tang
               ceair = min( eair(p),  esat_tv(p) )
               rh_can = ceair / esat_tv(p)

               ! Electron transport rate for C3 plants. Convert par from W/m2 to 
               ! umol photons/m**2/s using the factor 4.6

               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,iv))
               cquad = qabs * jmax_z(p,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je = min(r1,r2)

               ! Iterative loop for ci beginning with initial guess

               if (c3flag(p)) then
                  ci_z(p,iv) = 0.7_r8 * cair(p)
               else
                  ci_z(p,iv) = 0.4_r8 * cair(p)
               end if

               niter = 0

               ! Increment iteration counter. Stop if too many iterations

               niter = niter + 1

               ! Save old ci

               ciold = ci_z(p,iv)

               !find ci and stomatal conductance
               call hybrid(ciold, p, iv, c, gb_mol(p), je, cair(p), oair(p), &
                    lmr_z(p,iv), par_z(p,iv), rh_can, gs_mol(p,iv), niter, &
                    atm2lnd_vars, photosyns_vars)

               ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

               if (an(p,iv) < 0._r8) gs_mol(p,iv) = bbb(p)

               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

               cs = cair(p) - 1.4_r8/gb_mol(p) * an(p,iv) * forc_pbot(c)
               cs = max(cs,1.e-06_r8)
               ci_z(p,iv) = cair(p) - an(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv))

               ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

               gs = gs_mol(p,iv) / cf
               rs_z(p,iv) = min(1._r8/gs, rsmax0)

               ! Photosynthesis. Save rate-limiting photosynthesis

               psn_z(p,iv) = ag(p,iv)

               psn_wc_z(p,iv) = 0._r8
               psn_wj_z(p,iv) = 0._r8
               psn_wp_z(p,iv) = 0._r8

               if (ac(p,iv) <= aj(p,iv) .and. ac(p,iv) <= ap(p,iv)) then
                  psn_wc_z(p,iv) =  psn_z(p,iv)
               else if (aj(p,iv) < ac(p,iv) .and. aj(p,iv) <= ap(p,iv)) then
                  psn_wj_z(p,iv) =  psn_z(p,iv)
               else if (ap(p,iv) < ac(p,iv) .and. ap(p,iv) < aj(p,iv)) then
                  psn_wp_z(p,iv) =  psn_z(p,iv)
               end if

               ! Make sure iterative solution is correct

               if (gs_mol(p,iv) < 0._r8) then
                  write (iulog,*)'Negative stomatal conductance:'
                  write (iulog,*)'p,iv,gs_mol= ',p,iv,gs_mol(p,iv)
                  call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
               end if

               ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b

               hs = (gb_mol(p)*ceair + gs_mol(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol(p,iv))*esat_tv(p))
               rh_leaf(p) = hs
               gs_mol_err = mbb(p)*max(an(p,iv), 0._r8)*hs/cs*forc_pbot(c) + bbb(p)

               if (abs(gs_mol(p,iv)-gs_mol_err) > 1.e-01_r8) then
                  write (iulog,*) 'Ball-Berry error check - stomatal conductance error:'
                  write (iulog,*) gs_mol(p,iv), gs_mol_err
               end if

            end if    ! night or day if branch
         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Canopy photosynthesis and stomatal conductance
      !==============================================================================!

      ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
      ! unit leaf area), which are used in other parts of the model. Here, laican
      ! sums to either laisun or laisha.

      do f = 1, fn
         p = filterp(f)

         psncan = 0._r8
         psncan_wc = 0._r8
         psncan_wj = 0._r8
         psncan_wp = 0._r8
         lmrcan = 0._r8
         gscan = 0._r8
         laican = 0._r8
         do iv = 1, nrad(p)
            psncan = psncan + psn_z(p,iv) * lai_z(p,iv)
            psncan_wc = psncan_wc + psn_wc_z(p,iv) * lai_z(p,iv)
            psncan_wj = psncan_wj + psn_wj_z(p,iv) * lai_z(p,iv)
            psncan_wp = psncan_wp + psn_wp_z(p,iv) * lai_z(p,iv)
            lmrcan = lmrcan + lmr_z(p,iv) * lai_z(p,iv)
            gscan = gscan + lai_z(p,iv) / (rb(p)+rs_z(p,iv))
            laican = laican + lai_z(p,iv)
         end do
         if (laican > 0._r8) then
            psn(p) = psncan / laican
            psn_wc(p) = psncan_wc / laican
            psn_wj(p) = psncan_wj / laican
            psn_wp(p) = psncan_wp / laican
            lmr(p) = lmrcan / laican
            rs(p) = laican / gscan - rb(p)
         else
            psn(p) =  0._r8
            psn_wc(p) =  0._r8
            psn_wj(p) =  0._r8
            psn_wp(p) =  0._r8
            lmr(p) = 0._r8
            rs(p) = 0._r8
         end if
      end do

    end associate

  end subroutine Photosynthesis

  !------------------------------------------------------------------------------
  subroutine PhotosynthesisTotal (fn, filterp, &
       atm2lnd_vars, cnstate_vars, canopystate_vars, photosyns_vars)
    !
    ! Determine total photosynthesis
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    !
    ! !LOCAL VARIABLES:
    integer :: f,fp,p,l,g               ! indices
    !-----------------------------------------------------------------------

    associate(                                             &
         forc_pco2   => atm2lnd_vars%forc_pco2_grc       , & ! Input:  [real(r8) (:) ]  partial pressure co2 (Pa)                                             
         forc_pc13o2 => atm2lnd_vars%forc_pc13o2_grc     , & ! Input:  [real(r8) (:) ]  partial pressure c13o2 (Pa)                                           
         forc_po2    => atm2lnd_vars%forc_po2_grc        , & ! Input:  [real(r8) (:) ]  partial pressure o2 (Pa)                                              

         rc14_atm    => cnstate_vars%rc14_atm_patch      , & ! Input : [real(r8) (:) ]  C14O2/C12O2 in atmosphere 

         laisun      => canopystate_vars%laisun_patch    , & ! Input:  [real(r8) (:) ]  sunlit leaf area                                                      
         laisha      => canopystate_vars%laisha_patch    , & ! Input:  [real(r8) (:) ]  shaded leaf area                                                      
         
         psnsun      => photosyns_vars%psnsun_patch      , & ! Input:  [real(r8) (:) ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)                        
         psnsha      => photosyns_vars%psnsha_patch      , & ! Input:  [real(r8) (:) ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)                        
         rc13_canair => photosyns_vars%rc13_canair_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in canopy air
         rc13_psnsun => photosyns_vars%rc13_psnsun_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in sunlit canopy psn flux             
         rc13_psnsha => photosyns_vars%rc13_psnsha_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in shaded canopy psn flux  
         alphapsnsun => photosyns_vars%alphapsnsun_patch , & ! Output: [real(r8) (:) ]  fractionation factor in sunlit canopy psn flux  
         alphapsnsha => photosyns_vars%alphapsnsha_patch , & ! Output: [real(r8) (:) ]  fractionation factor in shaded canopy psn flux   
         psnsun_wc   => photosyns_vars%psnsun_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)        
         psnsun_wj   => photosyns_vars%psnsun_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)           
         psnsun_wp   => photosyns_vars%psnsun_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)        
         psnsha_wc   => photosyns_vars%psnsha_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)        
         psnsha_wj   => photosyns_vars%psnsha_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)           
         psnsha_wp   => photosyns_vars%psnsha_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)        
         c13_psnsun  => photosyns_vars%c13_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 13CO2 /m**2/ s)  
         c13_psnsha  => photosyns_vars%c13_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 13CO2 /m**2/ s) 
         c14_psnsun  => photosyns_vars%c14_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 14CO2 /m**2/ s)  
         c14_psnsha  => photosyns_vars%c14_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 14CO2 /m**2/ s) 
         fpsn        => photosyns_vars%fpsn_patch        , & ! Output: [real(r8) (:) ]  photosynthesis (umol CO2 /m**2 /s)                                    
         fpsn_wc     => photosyns_vars%fpsn_wc_patch     , & ! Output: [real(r8) (:) ]  Rubisco-limited photosynthesis (umol CO2 /m**2 /s)                    
         fpsn_wj     => photosyns_vars%fpsn_wj_patch     , & ! Output: [real(r8) (:) ]  RuBP-limited photosynthesis (umol CO2 /m**2 /s)                       
         fpsn_wp     => photosyns_vars%fpsn_wp_patch       & ! Output: [real(r8) (:) ]  product-limited photosynthesis (umol CO2 /m**2 /s)                    
         )

      do f = 1, fn
         p = filterp(f)
         g = pft%gridcell(p)

         if (.not. use_ed) then
            fpsn(p)    = psnsun(p)   *laisun(p) + psnsha(p)   *laisha(p)
            fpsn_wc(p) = psnsun_wc(p)*laisun(p) + psnsha_wc(p)*laisha(p)
            fpsn_wj(p) = psnsun_wj(p)*laisun(p) + psnsha_wj(p)*laisha(p)
            fpsn_wp(p) = psnsun_wp(p)*laisun(p) + psnsha_wp(p)*laisha(p)
         end if

         if (use_cn) then
            if ( use_c13 ) then
               rc13_canair(p) = forc_pc13o2(g)/(forc_pco2(g) - forc_pc13o2(g))
               rc13_psnsun(p) = rc13_canair(p)/alphapsnsun(p)
               rc13_psnsha(p) = rc13_canair(p)/alphapsnsha(p)
               c13_psnsun(p)  = psnsun(p) * (rc13_psnsun(p)/(1._r8+rc13_psnsun(p)))
               c13_psnsha(p)  = psnsha(p) * (rc13_psnsha(p)/(1._r8+rc13_psnsha(p)))

               ! use fixed c13 ratio with del13C of -25 to test the overall c13 structure
               ! c13_psnsun(p) = 0.01095627 * psnsun(p)
               ! c13_psnsha(p) = 0.01095627 * psnsha(p)
            endif
            if ( use_c14 ) then
               c14_psnsun(p) = rc14_atm(p) * psnsun(p)
               c14_psnsha(p) = rc14_atm(p) * psnsha(p)
            endif
         end if

      end do

    end associate

  end subroutine PhotosynthesisTotal

  !------------------------------------------------------------------------------
  subroutine Fractionation(bounds, &
       fn, filterp, &
       atm2lnd_vars, canopystate_vars, cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, &
       phase)
    !
    ! !DESCRIPTION:
    ! C13 fractionation during photosynthesis is calculated here after the nitrogen
    ! limitation is taken into account in the CNAllocation module.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds               
    integer                , intent(in)    :: fn                   ! size of pft filter
    integer                , intent(in)    :: filterp(fn)          ! patch filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(solarabs_type)    , intent(in)    :: solarabs_vars
    type(surfalb_type)     , intent(in)    :: surfalb_vars
    type(photosyns_type)   , intent(in)    :: photosyns_vars
    character(len=*)       , intent(in)    :: phase                ! 'sun' or 'sha'
    !
    ! !LOCAL VARIABLES:
    real(r8) , pointer :: par_z (:,:)   ! needed for backwards compatiblity
    real(r8) , pointer :: alphapsn (:)  ! needed for backwards compatiblity
    integer  :: f,p,c,g,iv              ! indices
    real(r8) :: co2(bounds%begp:bounds%endp)  ! atmospheric co2 partial pressure (pa)
    real(r8) :: ci
    !------------------------------------------------------------------------------

    associate(                                                  & 
         forc_pbot   => atm2lnd_vars%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_pco2   => atm2lnd_vars%forc_pco2_grc            , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             

         c3psn       => ecophyscon%c3psn                      , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3                              

         nrad        => surfalb_vars%nrad_patch               , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer             

         downreg     => cnstate_vars%downreg_patch            , & ! Input:  [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)                 

         an          => photosyns_vars%an_patch               , & ! Input:  [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)                           
         gb_mol      => photosyns_vars%gb_mol_patch           , & ! Input:  [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)                     
         gs_mol      => photosyns_vars%gs_mol_patch             & ! Input:  [real(r8) (:,:) ]  leaf stomatal conductance (umol H2O/m**2/s)                         
         )

      if (phase == 'sun') then
         par_z    =>    solarabs_vars%parsun_z_patch     ! Input :  [real(r8) (:,:)]  par absorbed per unit lai for canopy layer (w/m**2)                 
         alphapsn =>    photosyns_vars%alphapsnsun_patch ! Output:  [real(r8) (:)]                                                                        
      else if (phase == 'sha') then
         par_z    =>    solarabs_vars%parsha_z_patch     ! Input :  [real(r8) (:,:)]  par absorbed per unit lai for canopy layer (w/m**2)                 
         alphapsn =>    photosyns_vars%alphapsnsha_patch ! Output:  [real(r8) (:)]                                                                        
      end if

      do f = 1, fn
         p = filterp(f)
         c= pft%column(p)
         g= pft%gridcell(p)

         co2(p) = forc_pco2(g)
         do iv = 1,nrad(p)
            if (par_z(p,iv) <= 0._r8) then           ! night time
               alphapsn(p) = 1._r8
            else                                     ! day time
               ci = co2(p) - ((an(p,iv) * (1._r8-downreg(p)) ) * &
                    forc_pbot(c) * &
                    (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv)))
               alphapsn(p) = 1._r8 + (((c3psn(pft%itype(p)) * &
                    (4.4_r8 + (22.6_r8*(ci/co2(p))))) + &
                    ((1._r8 - c3psn(pft%itype(p))) * 4.4_r8))/1000._r8)
            end if
         end do
      end do

    end associate

  end subroutine Fractionation

  !-------------------------------------------------------------------------------
  subroutine hybrid(x0, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z,&
       rh_can, gs_mol,iter, &
       atm2lnd_vars, photosyns_vars)
    !
    !! DESCRIPTION:
    ! use a hybrid solver to find the root of equation  
    ! f(x) = x- h(x),
    !we want to find x, s.t. f(x) = 0.
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarrantee convergence.

    !
    !! REVISION HISTORY:
    !Dec 14/2012: created by Jinyun Tang
    !
    !!USES:   
    !
    !! ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: x0              !initial guess and final value of the solution
    real(r8), intent(in) :: lmr_z              ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in) :: par_z              ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in) :: rh_can             ! canopy air relative humidity
    real(r8), intent(in) :: gb_mol             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: je                 ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: cair               ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in) :: oair               ! Atmospheric O2 partial pressure (Pa)
    integer,  intent(in) :: p, iv, c           ! pft, c3/c4, and column index
    real(r8), intent(out) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
    integer,  intent(out) :: iter              !number of iterations used, for record only   
    type(atm2lnd_type)  , intent(in)    :: atm2lnd_vars
    type(photosyns_type), intent(inout) :: photosyns_vars
    !
    !! LOCAL VARIABLES
    real(r8) :: a, b
    real(r8) :: fa, fb
    real(r8) :: x1, f0, f1
    real(r8) :: x, dx
    real(r8), parameter :: eps = 1.e-2_r8      !relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8
    integer,  parameter :: itmax = 40          !maximum number of iterations
    real(r8) :: tol,minx,minf

    call ci_func(x0, f0, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_vars, photosyns_vars)

    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8

    call ci_func(x1,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_vars, photosyns_vars)

    if(f1==0._r8)then
       x0 = x1
       return
    endif
    if(f1<minf)then
       minx=x1
       minf=f1
    endif

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
       iter = iter + 1
       dx = - f1 * (x1-x0)/(f1-f0)
       x = x1 + dx
       tol = abs(x) * eps
       if(abs(dx)<tol)then
          x0 = x
          exit
       endif
       x0 = x1
       f0 = f1
       x1 = x   

       call ci_func(x1,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
            atm2lnd_vars, photosyns_vars)

       if(f1<minf)then
          minx=x1
          minf=f1
       endif
       if(abs(f1)<=eps1)then
          x0 = x1
          exit
       endif

       !if a root zone is found, use the brent method for a robust backup strategy
       if(f1 * f0 < 0._r8)then

          call brent(x, x0,x1,f0,f1, tol, p, iv, c, gb_mol, je, cair, oair, &
               lmr_z, par_z, rh_can, gs_mol, &
               atm2lnd_vars, photosyns_vars)

          x0=x
          exit
       endif
       if(iter>itmax)then 
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.

          call ci_func(minx,f1, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
               atm2lnd_vars, photosyns_vars)

          exit
       endif
    enddo

  end subroutine hybrid

  !------------------------------------------------------------------------------
  subroutine brent(x, x1,x2,f1, f2, tol, ip, iv, ic, gb_mol, je, cair, oair,&
       lmr_z, par_z, rh_can, gs_mol, &
       atm2lnd_vars, photosyns_vars)
    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!ARGUMENTS:
    real(r8), intent(out) :: x                ! indepedent variable of the single value function ci_func(x)
    real(r8), intent(in) :: x1, x2, f1, f2    ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: tol               ! the error tolerance
    real(r8), intent(in) :: lmr_z             ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in) :: par_z             ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in) :: gb_mol            ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: je                ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: cair              ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in) :: oair              ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in) :: rh_can            ! inside canopy relative humidity 
    integer,  intent(in) :: ip, iv, ic        ! pft, c3/c4, and column index
    real(r8), intent(out) :: gs_mol           ! leaf stomatal conductance (umol H2O/m**2/s)
    type(atm2lnd_type)  , intent(in)    :: atm2lnd_vars
    type(photosyns_type), intent(inout) :: photosyns_vars
    !
    !!LOCAL VARIABLES:
    integer, parameter :: ITMAX=20            !maximum number of iterations
    real(r8), parameter :: EPS=1.e-2_r8       !relative error tolerance
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    !------------------------------------------------------------------------------

    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       write(iulog,*) 'root must be bracketed for brent'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==ITMAX)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*EPS*abs(b)+0.5_r8*tol  !Convergence check.   
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0.)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif

       call ci_func(b, fb, ip, iv, ic, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_vars, photosyns_vars)

       if(fb==0._r8)exit

    enddo

    if(iter==ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb
    x=b

    return
  end subroutine brent

  !-------------------------------------------------------------------------------
  function ft(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use clm_varcon  , only : rgas, tfrz   
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:   
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft

  !-------------------------------------------------------------------------------   
  function fth(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    use clm_varcon  , only : rgas, tfrz   
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:      
    real(r8) :: ans
    !-------------------------------------------------------------------------------   

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth

  !-------------------------------------------------------------------------------   
  function fth25(hd,se)result(ans)
    !
    !!DESCRIPTION:   
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use clm_varcon  , only : rgas, tfrz   
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temperature function (J/mol/K)
    !
    ! !LOCAL VARIABLES:   
    real(r8) :: ans
    !-------------------------------------------------------------------------------   

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25

  !------------------------------------------------------------------------------
  subroutine ci_func(ci, fval, p, iv, c, gb_mol, je, cair, oair, lmr_z, par_z,&
       rh_can, gs_mol, atm2lnd_vars, photosyns_vars)
    !
    !! DESCRIPTION:
    ! evaluate the function
    ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
    !
    ! remark:  I am attempting to maintain the original code structure, also
    ! considering one may be interested to output relevant variables for the
    ! photosynthesis model, I have decided to add these relevant variables to
    ! the relevant data types.
    !
    !!ARGUMENTS:
    real(r8)             , intent(in)    :: ci       ! intracellular leaf CO2 (Pa)
    real(r8)             , intent(in)    :: lmr_z    ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8)             , intent(in)    :: par_z    ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8)             , intent(in)    :: gb_mol   ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)             , intent(in)    :: je       ! electron transport rate (umol electrons/m**2/s)
    real(r8)             , intent(in)    :: cair     ! Atmospheric CO2 partial pressure (Pa)
    real(r8)             , intent(in)    :: oair     ! Atmospheric O2 partial pressure (Pa)
    real(r8)             , intent(in)    :: rh_can   ! canopy air realtive humidity
    integer              , intent(in)    :: p, iv, c ! pft, vegetation type and column indexes
    real(r8)             , intent(out)   :: fval     ! return function of the value f(ci)
    real(r8)             , intent(out)   :: gs_mol   ! leaf stomatal conductance (umol H2O/m**2/s)
    type(atm2lnd_type)   , intent(in)    :: atm2lnd_vars
    type(photosyns_type) , intent(inout) :: photosyns_vars
    !
    !local variables
    real(r8) :: ai                  ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: cs                  ! CO2 partial pressure at leaf surface (Pa)

    real(r8) :: aquad, bquad, cquad  ! terms for quadratic equations
    real(r8) :: r1, r2               ! roots of quadratic equation
    real(r8) :: fnps                 ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii           ! empirical curvature parameter for electron transport rate
    real(r8) :: theta_ip             ! empirical curvature parameter for ap photosynthesis co-limitation
    !------------------------------------------------------------------------------

    associate(& 
         forc_pbot  => atm2lnd_vars%forc_pbot_downscaled_col   , & ! Output: [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         c3flag     => photosyns_vars%c3flag_patch             , & ! Output: [logical  (:)   ]  true if C3 and false if C4                                             
         ac         => photosyns_vars%ac_patch                 , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)              
         aj         => photosyns_vars%aj_patch                 , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)                 
         ap         => photosyns_vars%ap_patch                 , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_vars%ag_patch                 , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)              
         an         => photosyns_vars%an_patch                 , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)                           
         vcmax_z    => photosyns_vars%vcmax_z_patch            , & ! Input:  [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)                     
         cp         => photosyns_vars%cp_patch                 , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)                                           
         kc         => photosyns_vars%kc_patch                 , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)                                
         ko         => photosyns_vars%ko_patch                 , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)                                 
         qe         => photosyns_vars%qe_patch                 , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)          
         tpu_z      => photosyns_vars%tpu_z_patch              , & ! Output: [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)                 
         kp_z       => photosyns_vars%kp_z_patch               , & ! Output: [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)                     
         theta_cj   => photosyns_vars%theta_cj_patch           , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation 
         bbb        => photosyns_vars%bbb_patch                , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)                 
         mbb        => photosyns_vars%mbb_patch                  & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship           
         )

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      if (c3flag(p)) then
         ! C3: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv) * max(ci-cp(p), 0._r8) / (ci+kc(p)*(1._r8+oair/ko(p)))

         ! C3: RuBP-limited photosynthesis
         aj(p,iv) = je * max(ci-cp(p), 0._r8) / (4._r8*ci+8._r8*cp(p))

         ! C3: Product-limited photosynthesis 
         ap(p,iv) = 3._r8 * tpu_z(p,iv)

      else

         ! C4: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv)

         ! C4: RuBP-limited photosynthesis
         aj(p,iv) = qe(p) * par_z * 4.6_r8

         ! C4: PEP carboxylase-limited (CO2-limited)
         ap(p,iv) = kp_z(p,iv) * max(ci, 0._r8) / forc_pbot(c)

      end if

      ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap

      aquad = theta_cj(p)
      bquad = -(ac(p,iv) + aj(p,iv))
      cquad = ac(p,iv) * aj(p,iv)
      call quadratic (aquad, bquad, cquad, r1, r2)
      ai = min(r1,r2)

      aquad = theta_ip
      bquad = -(ai + ap(p,iv))
      cquad = ai * ap(p,iv)
      call quadratic (aquad, bquad, cquad, r1, r2)
      ag(p,iv) = min(r1,r2)

      ! Net photosynthesis. Exit iteration if an < 0

      an(p,iv) = ag(p,iv) - lmr_z
      if (an(p,iv) < 0._r8) then
         fval = 0._r8
         return
      endif
      ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
      ! With an <= 0, then gs_mol = bbb

      cs = cair - 1.4_r8/gb_mol * an(p,iv) * forc_pbot(c)
      cs = max(cs,1.e-06_r8)
      aquad = cs
      bquad = cs*(gb_mol - bbb(p)) - mbb(p)*an(p,iv)*forc_pbot(c)
      cquad = -gb_mol*(cs*bbb(p) + mbb(p)*an(p,iv)*forc_pbot(c)*rh_can)
      call quadratic (aquad, bquad, cquad, r1, r2)
      gs_mol = max(r1,r2)

      ! Derive new estimate for ci

      fval =ci - cair + an(p,iv) * forc_pbot(c) * (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)

    end associate

  end subroutine ci_func

end module PhotosynthesisMod
