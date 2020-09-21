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
  use elm_varctl          , only : iulog, use_c13, use_c14, use_cn, use_fates
  use elm_varpar          , only : nlevcan
  use elm_varctl          , only : use_hydrstress
  use elm_varpar          , only : nvegwcs, mxpft
  use elm_varcon          , only : namep
  use decompMod           , only : bounds_type
  use QuadraticMod        , only : quadratic
  use VegetationPropertiesType      , only : veg_vp
  use atm2lndType         , only : atm2lnd_type
  use CNStateType         , only : cnstate_type
  use CanopyStateType     , only : canopystate_type
  use TemperatureType     , only : temperature_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use PhotosynthesisType  , only : photosyns_type
  use VegetationType           , only : veg_pp
  use AllocationMod     , only : nu_com_leaf_physiology
  use PhosphorusStateType , only : phosphorusstate_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use elm_varctl          , only : cnallocate_carbon_only
  use elm_varctl          , only : cnallocate_carbonnitrogen_only
  use elm_varctl          , only : cnallocate_carbonphosphorus_only
  use elm_varctl          , only : iulog
  use pftvarcon           , only : noveg
  use SharedParamsMod     , only : ParamsShareInst
  use TopounitDataType    , only : top_as
  use VegetationDataType  , only : veg_es, veg_ns, veg_ps  
  use VegetationDataType, only : veg_wf, veg_ws
  use ColumnDataType      , only : col_es, col_ws, col_wf 
  use SoilStateType              , only : soilstate_type
  use WaterFluxType              , only : waterflux_type
  use WaterStateType             , only : waterstate_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Photosynthesis       ! Leaf stomatal resistance and leaf photosynthesis
  public :: PhotosynthesisTotal  ! Determine of total photosynthesis
  public :: Fractionation        ! C13 fractionation during photosynthesis 
  public :: PhotosynthesisHydraulicStress ! Leaf stomatal resistance and leaf photosynthesis
                                          ! Simultaneous solution of
                                          ! sunlit/shaded per Pierre
                                          ! Gentine/Daniel Kennedy plant
                                          ! hydraulic stress method
  public :: plc                           ! Return value of vulnerability curve at x


  ! !PRIVATE MEMBER FUNCTIONS:
  private :: hybrid         ! hybrid solver for ci
  private :: ci_func        ! ci function
  private :: brent          ! brent solver for root of a single variable function
  private :: ft             ! photosynthesis temperature response
  private :: fth            ! photosynthesis temperature inhibition
  private :: fth25          ! scaling factor for photosynthesis temperature inhibition
  !------------------------------------------------------------------------
  ! For plant hydraulics approach
  private :: hybrid_PHS     ! hybrid solver for ci
  private :: ci_func_PHS    ! ci function
  private :: brent_PHS      ! brent solver for root of a single variable function
  private :: calcstress     ! compute the root water stress
  private :: getvegwp       ! calculate vegetation water potential (sun, sha, xylem, root)
  private :: getqflx        ! calculate sunlit and shaded transpiration
  private :: spacF          ! flux divergence across each vegetation segment
  private :: spacA          ! the inverse Jacobian matrix relating delta(vegwp) to f, d(vegwp)=A*f
  private :: d1plc          ! compute 1st deriv of conductance attenuation for each segment

  ! !PRIVATE DATA:
  integer, parameter, private :: sun=1     ! index for sunlit
  integer, parameter, private :: sha=2     ! index for shaded
  integer, parameter, private :: xyl=3     ! index for xylem
  integer, parameter, private :: root=4    ! index for root
  integer, parameter, private :: veg=0     ! index for vegetation
  integer, parameter, private :: soil=1    ! index for soil
  integer, parameter, private :: stomatalcond_mtd_bb1987     = 1   ! Ball-Berry 1987 method for photosynthesis
  integer, parameter, private :: stomatalcond_mtd_medlyn2011 = 2   ! Medlyn 2011 method for photosynthesis
  ! !PUBLIC VARIABLES:

  type :: photo_params_type
     real(r8), allocatable, public  :: krmax              (:)
     real(r8), allocatable, private :: kmax               (:,:)
     real(r8), allocatable, private :: psi50              (:,:)
     real(r8), allocatable, private :: ck                 (:,:)
     real(r8), allocatable, public  :: psi_soil_ref       (:)
     real(r8), allocatable, private :: lmr_intercept_atkin(:)
  contains
     procedure, private :: allocParams
     procedure, public :: readParams
  end type photo_params_type
  !
  type(photo_params_type), public, protected :: params_inst  ! params_inst is populated in readParamsMod
contains

  !------------------------------------------------------------------------------
  subroutine allocParams ( this )
    ! 
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)

    implicit none

    ! !ARGUMENTS:
    class(photo_params_type) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'allocParams'
    !-----------------------------------------------------------------------

    ! allocate parameters

    allocate( this%krmax       (0:mxpft) )          ; this%krmax(:)        = nan
    allocate( this%kmax        (0:mxpft,nvegwcs) )  ; this%kmax(:,:)       = nan
    allocate( this%psi50       (0:mxpft,nvegwcs) )  ; this%psi50(:,:)      = nan
    allocate( this%ck          (0:mxpft,nvegwcs) )  ; this%ck(:,:)         = nan
    allocate( this%psi_soil_ref(0:mxpft) )          ; this%psi_soil_ref(:) = nan

    if ( use_hydrstress .and. nvegwcs /= 4 )then
       call endrun(msg='Error:: the Plant Hydraulics Stress methodology is for the spacA function is hardcoded for nvegwcs==4' &
                   //errMsg(__FILE__, __LINE__))
    end if

  end subroutine allocParams

  !------------------------------------------------------------------------------
  subroutine readParams ( this, ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    implicit none

    ! !ARGUMENTS:
    !class(photosyns_type) :: this
    class(photo_params_type) :: this
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'readParams'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: temp1d(0:mxpft) ! temporary to read in parameter
    real(r8)           :: temp2d(0:mxpft,nvegwcs) ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters


    call params_inst%allocParams()

    tString = "krmax"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid,readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%krmax=temp1d
    tString = "psi_soil_ref"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid,readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%psi_soil_ref=temp1d
    tString = "lmr_intercept_atkin"
    call ncd_io(varname=trim(tString),data=temp1d, flag='read', ncid=ncid,readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%lmr_intercept_atkin=temp1d
    tString = "kmax"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid,readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%kmax=temp2d
    tString = "psi50"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%psi50=temp2d
    tString = "ck"
    call ncd_io(varname=trim(tString),data=temp2d, flag='read', ncid=ncid,readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%ck=temp2d

  end subroutine readParams



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
    ! Note: This subroutine is not called via FATES (RGK)
    !
    ! !USES:
    use elm_varcon     , only : rgas, tfrz
    use elm_varctl     , only : cnallocate_carbon_only 
    use pftvarcon      , only : nbrdlf_dcd_tmp_shrub, nsoybean, nsoybeanirrig, npcropmin
    use pftvarcon      , only : vcmax_np1, vcmax_np2, vcmax_np3, vcmax_np4, jmax_np1, jmax_np2, jmax_np3
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
    integer  :: f,p,c,t,iv        ! indices
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
         c3psn         => veg_vp%c3psn                         , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3                              
         leafcn        => veg_vp%leafcn                        , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)                                                      
         flnr          => veg_vp%flnr                          , & ! Input:  [real(r8) (:)   ]  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)       
         fnitr         => veg_vp%fnitr                         , & ! Input:  [real(r8) (:)   ]  foliage nitrogen limitation factor (-)                                
         slatop        => veg_vp%slatop                        , & ! Input:  [real(r8) (:)   ]  specific leaf area at top of canopy, projected area basis [m^2/gC]    

         forc_pbot     => top_as%pbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             

         t_veg         => veg_es%t_veg             , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t10           => veg_es%t_a10             , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)                        
         tgcm          => veg_es%thm               , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)                     

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
         
         leafn         => veg_ns%leafn           , &
         leafn_storage => veg_ns%leafn_storage   , &
         leafn_xfer    => veg_ns%leafn_xfer      , &
         leafp         => veg_ps%leafp         , &
         leafp_storage => veg_ps%leafp_storage , &
         leafp_xfer    => veg_ps%leafp_xfer    , &
         i_vcmax       => veg_vp%i_vc                          , &
         s_vcmax       => veg_vp%s_vc                            &
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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! vcmax25 parameters, from CN

         fnr   = veg_vp%fnr(veg_pp%itype(p))   !7.16_r8
         act25 = veg_vp%act25(veg_pp%itype(p)) !3.6_r8   !umol/mgRubisco/min
         ! Convert rubisco activity units from umol/mgRubisco/min ->
         ! umol/gRubisco/s
         act25 = act25 * 1000.0_r8 / 60.0_r8

         ! Activation energy, from:
         ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !  Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
         ! except TPU from: Harley et al (1992) Plant, Cell and Environment
         ! 15:271-282

         kcha    = veg_vp%kcha(veg_pp%itype(p)) !79430._r8
         koha    = veg_vp%koha(veg_pp%itype(p)) !36380._r8
         cpha    = veg_vp%cpha(veg_pp%itype(p)) !37830._r8
         vcmaxha = veg_vp%vcmaxha(veg_pp%itype(p)) !72000._r8
         jmaxha  = veg_vp%jmaxha(veg_pp%itype(p)) !50000._r8
         tpuha   = veg_vp%tpuha(veg_pp%itype(p))  !72000._r8
         lmrha   = veg_vp%lmrha(veg_pp%itype(p))  !46390._r8

         ! High temperature deactivation, from:
         ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
         ! The factor "c" scales the deactivation to a value of 1.0 at 25C

         vcmaxhd = veg_vp%vcmaxhd(veg_pp%itype(p)) !200000._r8
         jmaxhd  = veg_vp%jmaxhd(veg_pp%itype(p))  !200000._r8
         tpuhd   = veg_vp%tpuhd(veg_pp%itype(p))   !200000._r8
         lmrhd   = veg_vp%lmrhd(veg_pp%itype(p))   !150650._r8
         lmrse   = veg_vp%lmrse(veg_pp%itype(p))   !490._r8
         lmrc    = fth25 (lmrhd, lmrse)

         ! C3 or C4 photosynthesis logical variable

         if (nint(c3psn(veg_pp%itype(p))) == 1) then
            c3flag(p) = .true. 
         else if (nint(c3psn(veg_pp%itype(p))) == 0) then
            c3flag(p) = .false.
         end if

         ! C3 and C4 dependent parameters

         if (c3flag(p)) then
            qe(p)       = veg_vp%qe(veg_pp%itype(p))       !0._r8
            theta_cj(p) = veg_vp%theta_cj(veg_pp%itype(p)) !0.98_r8
            bbbopt(p)   = veg_vp%bbbopt(veg_pp%itype(p))   !10000._r8
            mbbopt(p)   = veg_vp%mbbopt(veg_pp%itype(p))   !9._r8
         else
            qe(p)       = veg_vp%qe(veg_pp%itype(p))       !0.05_r8
            theta_cj(p) = veg_vp%theta_cj(veg_pp%itype(p)) !0.80_r8
            bbbopt(p)   = veg_vp%bbbopt(veg_pp%itype(p))   !40000._r8
            mbbopt(p)   = veg_vp%mbbopt(veg_pp%itype(p))   !4._r8
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

         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(t)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(t)
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
            lnc(p) = 1._r8 / (slatop(veg_pp%itype(p)) * leafcn(veg_pp%itype(p)))

            ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
            vcmax25top = lnc(p) * flnr(veg_pp%itype(p)) * fnr * act25 * dayl_factor(p)
            if (.not. use_cn) then
               vcmax25top = vcmax25top * fnitr(veg_pp%itype(p))
            else
               if ( CNAllocate_Carbon_only() ) vcmax25top = vcmax25top * fnitr(veg_pp%itype(p))
            end if

            ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
            ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.
            jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

         else
         
            ! leaf level nutrient control on photosynthesis rate added by Q. Zhu Aug 2015
            
            if ( CNAllocate_Carbon_only() .or. cnallocate_carbonphosphorus_only()) then

               lnc(p) = 1._r8 / (slatop(veg_pp%itype(p)) * leafcn(veg_pp%itype(p)))
               vcmax25top = lnc(p) * flnr(veg_pp%itype(p)) * fnr * act25 * dayl_factor(p)
               vcmax25top = vcmax25top * fnitr(veg_pp%itype(p))
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
                  lnc(p) = min(max(lnc(p),0.25_r8),3.0_r8) ! based on doi: 10.1002/ece3.1173
               else                                                                    
                  lnc(p) = 0.0_r8                                                      
               end if

               vcmax25top = (i_vcmax(veg_pp%itype(p)) + s_vcmax(veg_pp%itype(p)) * lnc(p)) * dayl_factor(p)
               jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

            else

               ! nu_com_leaf_physiology is true, vcmax25, jmax25 is derived from leafn, leafp concentration
               ! Anthony Walker 2014 DOI: 10.1002/ece3.1173

               if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
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
                     lnc(p) = min(max(lnc(p),0.25_r8),3.0_r8) ! based on doi: 10.1002/ece3.1173
                     lpc(p) = min(max(lpc(p),0.014_r8),0.85_r8) ! based on doi: 10.1002/ece3.1173
                     vcmax25top = exp(vcmax_np1(veg_pp%itype(p)) + vcmax_np2(veg_pp%itype(p))*log(lnc(p)) + &
                          vcmax_np3(veg_pp%itype(p))*log(lpc(p)) + vcmax_np4(veg_pp%itype(p))*log(lnc(p))*log(lpc(p)))&
                          * dayl_factor(p)
                     jmax25top = exp(jmax_np1 + jmax_np2*log(vcmax25top) + jmax_np3*log(lpc(p))) * dayl_factor(p)
                     vcmax25top = min(max(vcmax25top, 10.0_r8), 150.0_r8)
                     jmax25top = min(max(jmax25top, 10.0_r8), 250.0_r8)
                  else
                     lnc(p) = 0.0_r8
                     lpc(p) = 0.0_r8
                     vcmax25top = 0.0_r8
                     jmax25top = 0.0_r8
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

            lmr25top = 2.525e-6_r8 * (ParamsShareInst%Q10_mr ** ((25._r8 - 20._r8)/10._r8))
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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! Leaf boundary layer conductance, umol/m**2/s

         cf = forc_pbot(t)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
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
               call hybrid(ciold, p, iv, c, t, gb_mol(p), je, cair(p), oair(p), &
                    lmr_z(p,iv), par_z(p,iv), rh_can, gs_mol(p,iv), niter, &
                    atm2lnd_vars, photosyns_vars)

               ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

               if (an(p,iv) < 0._r8) gs_mol(p,iv) = bbb(p)

               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

               cs = cair(p) - 1.4_r8/gb_mol(p) * an(p,iv) * forc_pbot(t)
               cs = max(cs,1.e-06_r8)
               ci_z(p,iv) = cair(p) - an(p,iv) * forc_pbot(t) * (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv))

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
               gs_mol_err = mbb(p)*max(an(p,iv), 0._r8)*hs/cs*forc_pbot(t) + bbb(p)

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

    ! Note: This subroutine is not called via FATES (RGK)

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
    integer :: f,fp,p,l,t,g               ! indices
    !-----------------------------------------------------------------------

    associate(                                             &
         forc_pco2   => top_as%pco2bot                   , & ! Input:  [real(r8) (:) ]  partial pressure co2 (Pa)                                             
         forc_pc13o2 => top_as%pc13o2bot                 , & ! Input:  [real(r8) (:) ]  partial pressure c13o2 (Pa)                                           
         forc_po2    => top_as%po2bot                    , & ! Input:  [real(r8) (:) ]  partial pressure o2 (Pa)                                              

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
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)

         if (.not.use_fates) then
            fpsn(p)    = psnsun(p)   *laisun(p) + psnsha(p)   *laisha(p)
            fpsn_wc(p) = psnsun_wc(p)*laisun(p) + psnsha_wc(p)*laisha(p)
            fpsn_wj(p) = psnsun_wj(p)*laisun(p) + psnsha_wj(p)*laisha(p)
            fpsn_wp(p) = psnsun_wp(p)*laisun(p) + psnsha_wp(p)*laisha(p)
         end if

         if (use_cn) then
            if ( use_c13 ) then
               rc13_canair(p) = forc_pc13o2(t)/(forc_pco2(t) - forc_pc13o2(t))
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
    integer  :: f,p,c,t,g,iv            ! indices
    real(r8) :: co2(bounds%begp:bounds%endp)  ! atmospheric co2 partial pressure (pa)
    real(r8) :: ci
    !------------------------------------------------------------------------------

    associate(                                                  & 
         forc_pbot   => top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_pco2   => top_as%pco2bot                        , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             

         c3psn       => veg_vp%c3psn                          , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3

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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         co2(p) = forc_pco2(t)
         do iv = 1,nrad(p)
            if (par_z(p,iv) <= 0._r8) then           ! night time
               alphapsn(p) = 1._r8
            else                                     ! day time
               ci = co2(p) - ((an(p,iv) * (1._r8-downreg(p)) ) * &
                    forc_pbot(t) * &
                    (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv)))
               alphapsn(p) = 1._r8 + (((c3psn(veg_pp%itype(p)) * &
                    (4.4_r8 + (22.6_r8*(ci/co2(p))))) + &
                    ((1._r8 - c3psn(veg_pp%itype(p))) * 4.4_r8))/1000._r8)
            end if
         end do
      end do

    end associate

  end subroutine Fractionation

  !-------------------------------------------------------------------------------
  subroutine hybrid(x0, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z,&
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
    integer,  intent(in) :: p, iv, c, t        ! pft, c3/c4, column, and topounit index
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

    call ci_func(x0, f0, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
         atm2lnd_vars, photosyns_vars)

    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8

    call ci_func(x1,f1, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
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

       call ci_func(x1,f1, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
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

          call brent(x, x0,x1,f0,f1, tol, p, iv, c, t, gb_mol, je, cair, oair, &
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

          call ci_func(minx,f1, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
               atm2lnd_vars, photosyns_vars)

          exit
       endif
    enddo

  end subroutine hybrid

  !------------------------------------------------------------------------------
  subroutine brent(x, x1,x2,f1, f2, tol, ip, iv, ic, it, gb_mol, je, cair, oair,&
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
    integer,  intent(in) :: ip, iv, ic, it    ! pft, c3/c4, column, and topounit index
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

       call ci_func(b, fb, ip, iv, ic, it, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
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
    use elm_varcon  , only : rgas, tfrz   
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
    use elm_varcon  , only : rgas, tfrz   
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
    use elm_varcon  , only : rgas, tfrz   
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
  subroutine ci_func(ci, fval, p, iv, c, t, gb_mol, je, cair, oair, lmr_z, par_z,&
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
    integer              , intent(in)    :: p,iv,c,t ! pft, vegetation type, column, and topounit indexes
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
         forc_pbot  => top_as%pbot                             , & ! Output: [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
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
         ap(p,iv) = kp_z(p,iv) * max(ci, 0._r8) / forc_pbot(t)

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

      cs = cair - 1.4_r8/gb_mol * an(p,iv) * forc_pbot(t)
      cs = max(cs,1.e-06_r8)
      aquad = cs
      bquad = cs*(gb_mol - bbb(p)) - mbb(p)*an(p,iv)*forc_pbot(t)
      cquad = -gb_mol*(cs*bbb(p) + mbb(p)*an(p,iv)*forc_pbot(t)*rh_can)
      call quadratic (aquad, bquad, cquad, r1, r2)
      gs_mol = max(r1,r2)

      ! Derive new estimate for ci

      fval =ci - cair + an(p,iv) * forc_pbot(t) * (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)

    end associate

  end subroutine ci_func
  !------------------------------------------------------------------------------

  subroutine PhotosynthesisHydraulicStress ( bounds, fn, filterp, &
       esat_tv, eair, oair, cair, rb, bsun, bsha, btran, dayl_factor,  &
       qsatl, qaf, &
       atm2lnd_inst, temperature_inst, soilstate_inst, waterstate_inst, &
       surfalb_inst, solarabs_inst, canopystate_inst, &
       photosyns_inst, waterflux_inst, nitrogenstate_vars,phosphorusstate_vars)

    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    ! Here, sunlit and shaded photosynthesis and stomatal conductance are solved 
    ! simultaneously per Pierre Gentine/Daniel Kennedy plant hydraulic stress
    ! method
    !
    ! !USES:
    use elm_varcon        , only : rgas, tfrz, rpi
    use elm_varctl        , only : cnallocate_carbon_only
    !use elm_varctl        , only : lnc_opt, reduce_dayl_factor, vcmax_opt    
    use elm_varpar        , only : nlevsoi
    use pftvarcon         , only : nbrdlf_dcd_tmp_shrub, npcropmin
    use pftvarcon         , only : vcmax_np1, vcmax_np2, vcmax_np3, vcmax_np4, jmax_np1, jmax_np2, jmax_np3
    use ColumnType        , only : col_pp        

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
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: )    ! scalar (0-1) for daylength
    real(r8)               , intent(in)    :: qsatl ( bounds%begp: )         ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)    :: qaf ( bounds%begp: )           ! humidity of canopy air [kg/kg]
    !real(r8)               , intent(in)    :: leafn( bounds%begp: )          ! leaf N (gN/m2)
    real(r8)               , intent(out)   :: bsun( bounds%begp: )           ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)   :: bsha( bounds%begp: )           ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)   :: btran( bounds%begp: )          ! transpiration wetness factor (0 to 1) [pft]

    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    type(surfalb_type)     , intent(in)    :: surfalb_inst
    type(solarabs_type)    , intent(in)    :: solarabs_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    !class(ozone_base_type) , intent(in)    :: ozone_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(nitrogenstate_type)  , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type), intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    !
    real(r8) :: froot_carbon( bounds%begp:bounds%endp )    ! fine root carbon (gC/m2) [pft]   
    real(r8) :: croot_carbon( bounds%begp:bounds%endp )    ! live coarse root carbon (gC/m2) [pft]   
    ! Leaf photosynthesis parameters
    real(r8) :: jmax_z(bounds%begp:bounds%endp,2,nlevcan) ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: bbbopt(bounds%begp:bounds%endp)           ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(r8) :: mbbopt(bounds%begp:bounds%endp)           ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    real(r8) :: kn(bounds%begp:bounds%endp)               ! leaf nitrogen decay coefficient
    real(r8) :: vcmax25top     ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top      ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top       ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top       ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top        ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: vcmax25_sun    ! sunlit leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: vcmax25_sha    ! shaded leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25_sun     ! sunlit leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: jmax25_sha     ! shaded leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25_sun      ! sunlit leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: tpu25_sha      ! shaded leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25_sun      ! sunlit leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25_sha      ! shaded leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25_sun       ! sunlit leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kp25_sha       ! shaded leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
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

    real(r8) :: theta_ip       ! empirical curvature parameter for ap photosynthesis co-limitation

    ! Other
    integer  :: f,p,c,t,iv          ! indices
    integer, parameter :: sun = 1 ! index for sunlit leaves
    integer, parameter :: sha = 2 ! index for shaded leaves

    real(r8) :: cf                ! s m**2/umol -> s/m
    real(r8) :: rsmax0            ! maximum stomatal resistance [s/m]
    real(r8) :: gb                ! leaf boundary layer conductance (m/s)
    real(r8) :: cs_sun            ! CO2 partial pressure at sunlit leaf surface (Pa)
    real(r8) :: cs_sha            ! CO2 partial pressure at shaded leaf surface (Pa)
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
    real(r8) :: je_sun            ! sunlit leaf electron transport rate (umol electrons/m**2/s)
    real(r8) :: je_sha            ! shaded leaf electron transport rate (umol electrons/m**2/s)
    real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
    real(r8) :: r1,r2             ! roots of quadratic equation
    real(r8) :: ceair             ! vapor pressure of air, constrained (Pa)
    real(r8) :: fnr               ! (gRubisco/gN in Rubisco)
    real(r8) :: act25             ! (umol/mgRubisco/min) Rubisco activity at 25 C
    integer  :: iter1             ! number of iterations used, for record only
    integer  :: iter2             ! number of iterations used, for record only 
    real(r8) :: nscaler           ! leaf nitrogen scaling coefficient
    real(r8) :: nscaler_sun       ! sunlit leaf nitrogen scaling coefficient
    real(r8) :: nscaler_sha       ! shaded leaf nitrogen scaling coefficient

    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)

    real(r8) :: psn_wc_z_sun(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z_sun(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z_sun(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to sunlit psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wc_z_sha(bounds%begp:bounds%endp,nlevcan) ! Rubisco-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z_sha(bounds%begp:bounds%endp,nlevcan) ! RuBP-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z_sha(bounds%begp:bounds%endp,nlevcan) ! product-limited contribution to shaded psn_z (umol CO2/m**2/s)
    real(r8) :: rh_leaf_sun(bounds%begp:bounds%endp)          ! fractional humidity at sunlit leaf surface (dimensionless)
    real(r8) :: rh_leaf_sha(bounds%begp:bounds%endp)          ! fractional humidity at shaded leaf surface (dimensionless)

    real(r8) :: psncan_sun            ! canopy sum of sunlit psn_z
    real(r8) :: psncan_wc_sun         ! canopy sum of sunlit psn_wc_z
    real(r8) :: psncan_wj_sun         ! canopy sum of sunlit psn_wj_z
    real(r8) :: psncan_wp_sun         ! canopy sum of sunlit psn_wp_z
    real(r8) :: lmrcan_sun            ! canopy sum of sunlit lmr_z
    real(r8) :: gscan_sun             ! canopy sum of sunlit leaf conductance
    real(r8) :: laican_sun            ! canopy sum of sunlit lai_z
    real(r8) :: psncan_sha            ! canopy sum of shaded psn_z
    real(r8) :: psncan_wc_sha         ! canopy sum of shaded psn_wc_z
    real(r8) :: psncan_wj_sha         ! canopy sum of shaded psn_wj_z
    real(r8) :: psncan_wp_sha         ! canopy sum of shaded psn_wp_z
    real(r8) :: lmrcan_sha            ! canopy sum of shaded lmr_z
    real(r8) :: gscan_sha             ! canopy sum of shaded leaf conductance
    real(r8) :: laican_sha            ! canopy sum of shaded lai_z
    real(r8) :: laican                ! canopy sum of lai_z
    real(r8) :: rh_can                ! canopy air relative humidity

    real(r8) , pointer :: an_sun          (:,:) ! net sunlit leaf photosynthesis (umol CO2/m**2/s)
    real(r8) , pointer :: an_sha          (:,:) ! net shaded leaf photosynthesis (umol CO2/m**2/s)
    real(r8) , pointer :: lai_z_sun       (:,:) ! leaf area index for canopy layer, sunlit
    real(r8) , pointer :: par_z_sun       (:,:) ! par absorbed per unit lai for canopy layer, sunlit (w/m**2)
    real(r8) , pointer :: vcmaxcint_sun   (:)   ! leaf to canopy scaling coefficient, sunlit
    real(r8) , pointer :: alphapsn_sun    (:)   ! 13C fractionation factor for PSN, sunlit ()
    real(r8) , pointer :: psn_sun         (:)   ! foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wc_sun      (:)   ! Rubisco-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wj_sun      (:)   ! RuBP-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +] 
    real(r8) , pointer :: psn_wp_sun      (:)   ! product-limited foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_z_sun       (:,:) ! canopy layer: foliage photosynthesis, sunlit (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: lmr_sun         (:)   ! leaf maintenance respiration rate, sunlit (umol CO2/m**2/s)
    real(r8) , pointer :: lmr_z_sun       (:,:) ! canopy layer: leaf maintenance respiration rate, sunlit (umol CO2/m**2/s)
    real(r8) , pointer :: rs_sun          (:)   ! leaf stomatal resistance, sunlit (s/m)
    real(r8) , pointer :: rs_z_sun        (:,:) ! canopy layer: leaf stomatal resistance, sunlit (s/m)
    real(r8) , pointer :: ci_z_sun        (:,:) ! intracellular leaf CO2, sunlit (Pa)
    real(r8) , pointer :: o3coefv_sun     (:)   ! o3 coefficient used in photo calculation, sunlit
    real(r8) , pointer :: o3coefg_sun     (:)   ! o3 coefficient used in rs calculation, sunlit
    real(r8) , pointer :: gs_mol_sun      (:,:) ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8) , pointer :: lai_z_sha       (:,:) ! leaf area index for canopy layer, shaded
    real(r8) , pointer :: par_z_sha       (:,:) ! par absorbed per unit lai for canopy layer, shaded (w/m**2)
    real(r8) , pointer :: vcmaxcint_sha   (:)   ! leaf to canopy scaling coefficient, shaded
    real(r8) , pointer :: alphapsn_sha    (:)   ! 13C fractionation factor for PSN, shaded ()
    real(r8) , pointer :: psn_sha         (:)   ! foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wc_sha      (:)   ! Rubisco-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_wj_sha      (:)   ! RuBP-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +] 
    real(r8) , pointer :: psn_wp_sha      (:)   ! product-limited foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: psn_z_sha       (:,:) ! canopy layer: foliage photosynthesis, shaded (umol co2 /m**2/ s) [always +]
    real(r8) , pointer :: lmr_sha         (:)   ! leaf maintenance respiration rate, shaded (umol CO2/m**2/s)
    real(r8) , pointer :: lmr_z_sha       (:,:) ! canopy layer: leaf maintenance respiration rate, shaded (umol CO2/m**2/s)
    real(r8) , pointer :: rs_sha          (:)   ! leaf stomatal resistance, shaded (s/m)
    real(r8) , pointer :: rs_z_sha        (:,:) ! canopy layer: leaf stomatal resistance, shaded (s/m)
    real(r8) , pointer :: ci_z_sha        (:,:) ! intracellular leaf CO2, shaded (Pa)
    real(r8) , pointer :: o3coefv_sha     (:)   ! o3 coefficient used in photo calculation, shaded
    real(r8) , pointer :: o3coefg_sha     (:)   ! o3 coefficient used in rs calculation, shaded
    real(r8) , pointer :: gs_mol_sha      (:,:) ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    real(r8) , pointer :: vegwp           (:,:) ! vegetation water matric potential (mm)

    real(r8) :: sum_nscaler
    real(r8) :: total_lai                
    integer  :: nptreemax                
!scs
    integer  :: j                       ! index
    real(r8) :: rs_resis                ! combined soil-root resistance [s]
    real(r8) :: r_soil                  ! root spacing [m]
    real(r8) :: root_biomass_density    ! root biomass density [g/m3]
    real(r8) :: root_cross_sec_area     ! root cross sectional area [m2]
    real(r8) :: root_length_density     ! root length density [m/m3]
    real(r8) :: froot_average_length    ! average coarse root length [m]
    real(r8) :: croot_average_length    ! average coarse root length [m]
    real(r8) :: soil_conductance        ! soil to root hydraulic conductance [1/s]
    real(r8) :: root_conductance        ! root hydraulic conductance [1/s]
    real(r8) :: rai(nlevsoi)            ! root area index [m2/m2]
    real(r8) :: rai_n(nlevsoi)            ! root area index [m2/m2]
    real(r8) :: fs(nlevsoi)             ! root conductance scale factor (reduction in conductance due to decreasing (more negative) root water potential)
    real(r8) :: fs_n(nlevsoi)             ! root conductance scale factor (reduction in conductance due to decreasing (more negative) root water potential)
    real(r8) :: gsminsun                ! Minimum stomatal conductance sunlit
    real(r8) :: gsminsha                ! Minimum stomatal conductance shaded
    real(r8) :: gs_slope_sun            ! Slope stomatal conductance sunlit
    real(r8) :: gs_slope_sha            ! Slope stomatal conductance shaded
    real(r8), parameter :: croot_lateral_length = 0.25_r8   ! specified lateral coarse root length [m]
    real(r8), parameter :: c_to_b = 2.0_r8           !(g biomass /g C)
    real(r8) :: lnc(bounds%begp:bounds%endp)   ! leaf N concentration (gN leaf/m^2)
    real(r8) :: lpc(bounds%begp:bounds%endp)   ! leaf N concentration (gN leaf/m^2)
    real(r8) :: grav2(nlevsoi)        ! soil layer gravitational potential relative to surface (mm H2O) 
!Note that root density is for dry biomass not carbon. CLM provides root biomass
!as carbon. The conversion is 0.5 g C / g biomass

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
    SHR_ASSERT_ALL((ubound(bsun)        == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(bsha)        == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(btran)       == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dayl_factor) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(qsatl)       == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(qaf)         == (/bounds%endp/)), errMsg(__FILE__, __LINE__))

    associate(                                                 &
            qflx_rootsoi_col    => col_wf%qflx_rootsoi    , & ! Output: [real(r8) (:,:) ]  
         k_soil_root  => soilstate_inst%k_soil_root_patch    , & ! Input: [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         hk_l         =>    soilstate_inst%hk_l_col          , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity (mm/s) 
         hksat        => soilstate_inst%hksat_col            , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         smp          => soilstate_inst%smp_l_col            , & ! Input: [real(r8) (:,:) ]  soil matrix potential [mm]

         root_conductance_patch => soilstate_inst%root_conductance_patch , & ! Output:   [real(r8) (:,:)] root conductance
         soil_conductance_patch => soilstate_inst%soil_conductance_patch , & ! Output:   [real(r8) (:,:)] soil conductance
         rootfr       => soilstate_inst%rootfr_patch         , & ! Input: [real(r8) (:,:)]
         dz           => col_pp%dz                              , & ! Input: [real(r8) (:,:) ]  layer thickness (m)
         z            => col_pp%z                               , & ! Input: [real(r8) (:,:) ]  layer depth (m)

         c3psn      => veg_vp%c3psn                          , & ! Input:  photosynthetic pathway: 0. = c4, 1. = c3
         leafcn     => veg_vp%leafcn                         , & ! Input:  leaf C:N (gC/gN)
         flnr       => veg_vp%flnr                           , & ! Input:  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
         fnitr      => veg_vp%fnitr                          , & ! Input:  foliage nitrogen limitation factor (-)
         slatop     => veg_vp%slatop                         , & ! Input:  specific leaf area at top of canopy, projected area basis [m^2/gC]
         stem_leaf     => veg_vp%stem_leaf                         , & ! allocation parameter: new stem c per new leaf C (gC/gC)
         froot_leaf     => veg_vp%froot_leaf                         , & ! allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem     => veg_vp%croot_stem                         , & ! allocation parameter: new coarse root C per new stem C (gC/gC)
         forc_pbot  => top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)

         t_veg         => veg_es%t_veg             , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t10           => veg_es%t_a10             , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)                        
         tgcm          => veg_es%thm               , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)                     
         nrad       => surfalb_inst%nrad_patch               , & ! Input:  [integer  (:)   ]  pft number of canopy layers, above snow for radiative transfer
         tlai_z     => surfalb_inst%tlai_z_patch             , & ! Input:  [real(r8) (:,:) ]  pft total leaf area index for canopy layer
         tlai       => canopystate_inst%tlai_patch           , & ! Input:  [real(r8)(:)    ]  one-sided leaf area index, no burying by snow  
         tsai       => canopystate_inst%tsai_patch           , & ! Input:  [real(r8)(:)    ]  one-sided leaf area index, no burying by snow  
         c3flag     => photosyns_inst%c3flag_patch           , & ! Output: [logical  (:)   ]  true if C3 and false if C4
         ac         => photosyns_inst%ac_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_phs_patch      , & ! Output: [real(r8) (:,:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         tpu_z      => photosyns_inst%tpu_z_phs_patch        , & ! Output: [real(r8) (:,:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_phs_patch         , & ! Output: [real(r8) (:,:,:) ]  initial slope of CO2 response curve (C4 plants)
         gb_mol     => photosyns_inst%gb_mol_patch           , & ! Output: [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         cp         => photosyns_inst%cp_patch               , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch               , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch               , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         theta_cj   => photosyns_inst%theta_cj_patch         , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         rh_leaf    => photosyns_inst%rh_leaf_patch          , & ! Output: [real(r8) (:)   ]  fractional humidity at leaf surface (dimensionless)
         leafn         => veg_ns%leafn           , &
         leafn_storage => veg_ns%leafn_storage   , &
         leafn_xfer    => veg_ns%leafn_xfer      , &
         leafp         => veg_ps%leafp         , &
         leafp_storage => veg_ps%leafp_storage , &
         leafp_xfer    => veg_ps%leafp_xfer    , &
         i_vcmax       => veg_vp%i_vc                          , &
         s_vcmax       => veg_vp%s_vc                          , &
         bsw           => soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         sucsat        => soilstate_inst%sucsat_col             ,  & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         ivt           => veg_pp%itype                             & ! Input:  [integer  (:)   ]  patch vegetation type
      )
      an_sun        =>    photosyns_inst%an_sun_patch         ! Output: [real(r8) (:,:) ]  net sunlit leaf photosynthesis (umol CO2/m**2/s)
      an_sha        =>    photosyns_inst%an_sha_patch         ! Output: [real(r8) (:,:) ]  net shaded leaf photosynthesis (umol CO2/m**2/s)
      par_z_sun     =>    solarabs_inst%parsun_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
      lai_z_sun     =>    canopystate_inst%laisun_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
      vcmaxcint_sun =>    surfalb_inst%vcmaxcintsun_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
      alphapsn_sun  =>    photosyns_inst%alphapsnsun_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
      ci_z_sun      =>    photosyns_inst%cisun_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
      rs_sun        =>    photosyns_inst%rssun_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
      rs_z_sun      =>    photosyns_inst%rssun_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
      lmr_sun       =>    photosyns_inst%lmrsun_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
      lmr_z_sun     =>    photosyns_inst%lmrsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
      psn_sun       =>    photosyns_inst%psnsun_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_z_sun     =>    photosyns_inst%psnsun_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wc_sun    =>    photosyns_inst%psnsun_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wj_sun    =>    photosyns_inst%psnsun_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wp_sun    =>    photosyns_inst%psnsun_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      gs_mol_sun    =>    photosyns_inst%gs_mol_sun_patch     ! Output: [real(r8) (:,:) ]  sunlit leaf stomatal conductance (umol H2O/m**2/s)
      par_z_sha     =>    solarabs_inst%parsha_z_patch        ! Input:  [real(r8) (:,:) ]  par absorbed per unit lai for canopy layer (w/m**2)
      lai_z_sha     =>    canopystate_inst%laisha_z_patch     ! Input:  [real(r8) (:,:) ]  leaf area index for canopy layer, sunlit or shaded
      vcmaxcint_sha =>    surfalb_inst%vcmaxcintsha_patch     ! Input:  [real(r8) (:)   ]  leaf to canopy scaling coefficient
      alphapsn_sha  =>    photosyns_inst%alphapsnsha_patch    ! Input:  [real(r8) (:)   ]  13C fractionation factor for PSN ()
      ci_z_sha      =>    photosyns_inst%cisha_z_patch        ! Output: [real(r8) (:,:) ]  intracellular leaf CO2 (Pa)
      rs_sha        =>    photosyns_inst%rssha_patch          ! Output: [real(r8) (:)   ]  leaf stomatal resistance (s/m)
      rs_z_sha      =>    photosyns_inst%rssha_z_patch        ! Output: [real(r8) (:,:) ]  canopy layer: leaf stomatal resistance (s/m)
      lmr_sha       =>    photosyns_inst%lmrsha_patch         ! Output: [real(r8) (:)   ]  leaf maintenance respiration rate (umol CO2/m**2/s)
      lmr_z_sha     =>    photosyns_inst%lmrsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
      psn_sha       =>    photosyns_inst%psnsha_patch         ! Output: [real(r8) (:)   ]  foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_z_sha     =>    photosyns_inst%psnsha_z_patch       ! Output: [real(r8) (:,:) ]  canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wc_sha    =>    photosyns_inst%psnsha_wc_patch      ! Output: [real(r8) (:)   ]  Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wj_sha    =>    photosyns_inst%psnsha_wj_patch      ! Output: [real(r8) (:)   ]  RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      psn_wp_sha    =>    photosyns_inst%psnsha_wp_patch      ! Output: [real(r8) (:)   ]  product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
      gs_mol_sha    =>    photosyns_inst%gs_mol_sha_patch     ! Output: [real(r8) (:,:) ]  shaded leaf stomatal conductance (umol H2O/m**2/s)
      vegwp         => canopystate_inst%vegwp_patch           ! Input/Output: [real(r8) (:,:) ]  vegetation water matric potential (mm)
      !==============================================================================!
      ! Photosynthesis and stomatal conductance parameters, from:
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      !==============================================================================!
      ! calculate root-soil interface conductance 
      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         froot_carbon(p) = tlai(p) / slatop(ivt(p))*froot_leaf(ivt(p))
         croot_carbon(p) = tlai(p)/slatop(ivt(p))*stem_leaf(ivt(p)) &
                  *croot_stem(ivt(p))

         do j = 1,nlevsoi

            ! calculate conversion from conductivity to conductance
            root_biomass_density = c_to_b * froot_carbon(p) * rootfr(p,j) / dz(c,j)
            ! ensure minimum root biomass (using 1gC/m2)
            root_biomass_density = max(c_to_b*1._r8,root_biomass_density)

            ! Root length density: m root per m3 soil
            root_cross_sec_area = rpi*veg_vp%root_radius(veg_pp%itype(p))**2
            root_length_density = root_biomass_density / (veg_vp%root_density(veg_pp%itype(p)) * root_cross_sec_area)

            ! Root-area index (RAI)
            rai(j) = (tsai(p)+tlai(p)) * veg_vp%froot_leaf(veg_pp%itype(p)) * rootfr(p,j)

            ! fix coarse root_average_length to specified length
            croot_average_length = croot_lateral_length

            ! calculate r_soil using Gardner/spa equation (Bonan, GMD, 2014)
            r_soil = sqrt(1./(rpi*root_length_density))

            ! length scale approach
            soil_conductance = min(hksat(c,j),hk_l(c,j))/(1.e3*r_soil)

            ! use vegetation plc function to adjust root conductance
            fs(j)=  plc(smp(c,j),p,c,root,veg)

            ! krmax is root conductance per area per length
            root_conductance = (fs(j)*rai(j)*params_inst%krmax(veg_pp%itype(p)))/(croot_average_length + z(c,j))

            soil_conductance = max(soil_conductance, 1.e-16_r8)
            root_conductance = max(root_conductance, 1.e-16_r8)

            root_conductance_patch(p,j) = root_conductance
            soil_conductance_patch(p,j) = soil_conductance

            ! sum resistances in soil and root
            rs_resis = 1._r8/soil_conductance + 1._r8/root_conductance

            ! conductance is inverse resistance
            ! explicitly set conductance to zero for top soil layer
            if(rai(j)*rootfr(p,j) > 0._r8 .and. j > 1) then
               k_soil_root(p,j) =  1._r8/rs_resis
            else
               k_soil_root(p,j) =  0.
            endif
         end do
      enddo

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116,
      ! doi:10.1029/2010JG001593

      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! vcmax25 parameters, from CN

         fnr   = veg_vp%fnr(veg_pp%itype(p))   !7.16_r8
         act25 = veg_vp%act25(veg_pp%itype(p)) !3.6_r8   !umol/mgRubisco/min
         ! Convert rubisco activity units from umol/mgRubisco/min ->
         ! umol/gRubisco/s
         act25 = act25 * 1000.0_r8 / 60.0_r8

         ! Activation energy, from:
         ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !  Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
         ! except TPU from: Harley et al (1992) Plant, Cell and Environment
         ! 15:271-282

         kcha    = veg_vp%kcha(veg_pp%itype(p)) !79430._r8
         koha    = veg_vp%koha(veg_pp%itype(p)) !36380._r8
         cpha    = veg_vp%cpha(veg_pp%itype(p)) !37830._r8
         vcmaxha = veg_vp%vcmaxha(veg_pp%itype(p)) !72000._r8
         jmaxha  = veg_vp%jmaxha(veg_pp%itype(p)) !50000._r8
         tpuha   = veg_vp%tpuha(veg_pp%itype(p))  !72000._r8
         lmrha   = veg_vp%lmrha(veg_pp%itype(p))  !46390._r8

         ! High temperature deactivation, from:
         ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
         ! The factor "c" scales the deactivation to a value of 1.0 at 25C

         vcmaxhd = veg_vp%vcmaxhd(veg_pp%itype(p)) !200000._r8
         jmaxhd  = veg_vp%jmaxhd(veg_pp%itype(p))  !200000._r8
         tpuhd   = veg_vp%tpuhd(veg_pp%itype(p))   !200000._r8
         lmrhd   = veg_vp%lmrhd(veg_pp%itype(p))   !150650._r8
         lmrse   = veg_vp%lmrse(veg_pp%itype(p))   !490._r8
         lmrc    = fth25 (lmrhd, lmrse)

         ! C3 or C4 photosynthesis logical variable

         if (nint(c3psn(veg_pp%itype(p))) == 1) then
            c3flag(p) = .true.
         else if (nint(c3psn(veg_pp%itype(p))) == 0) then
            c3flag(p) = .false.
         end if

         ! C3 and C4 dependent parameters

         if (c3flag(p)) then
            qe(p)       = veg_vp%qe(veg_pp%itype(p))       !0._r8
            theta_cj(p) = veg_vp%theta_cj(veg_pp%itype(p)) !0.98_r8
            bbbopt(p)   = veg_vp%bbbopt(veg_pp%itype(p))   !10000._r8
            mbbopt(p)   = veg_vp%mbbopt(veg_pp%itype(p))   !9._r8
         else
            qe(p)       = veg_vp%qe(veg_pp%itype(p))       !0.05_r8
            theta_cj(p) = veg_vp%theta_cj(veg_pp%itype(p)) !0.80_r8
            bbbopt(p)   = veg_vp%bbbopt(veg_pp%itype(p))   !40000._r8
            mbbopt(p)   = veg_vp%mbbopt(veg_pp%itype(p))   !4._r8
         end if

         ! Soil water stress applied to Ball-Berry parameters

         bbb(p) = bbbopt(p)
         mbb(p) = mbbopt(p)

         ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment
         ! 24:253-259
         !
         !       ko25 = 278.4 mmol/mol
         !       cp25 = 42.75 umol/mol
         !
         ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and
         ! re-calculate
         ! cp to account for variation in O2 using cp = 0.5 O2 / sco
         !

         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(t)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(t)
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
            ! Leaf nitrogen concentration at the top of the canopy (g N leaf /
            ! m**2 leaf)
            lnc(p) = 1._r8 / (slatop(veg_pp%itype(p)) * leafcn(veg_pp%itype(p)))

            ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
            vcmax25top = lnc(p) * flnr(veg_pp%itype(p)) * fnr * act25 *dayl_factor(p)
            if (.not. use_cn) then
               vcmax25top = vcmax25top * fnitr(veg_pp%itype(p))
            else
               if ( CNAllocate_Carbon_only() ) vcmax25top = vcmax25top *fnitr(veg_pp%itype(p))
            end if

            ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116,
            ! doi:10.1029/2010JG001593
            ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of
            ! Experimental Botany 44:907-920.
            jmax25top = (2.59_r8 -0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

         else

            ! leaf level nutrient control on photosynthesis rate added by Q. Zhu
            ! Aug 2015

            if ( CNAllocate_Carbon_only() .or.cnallocate_carbonphosphorus_only()) then

               lnc(p) = 1._r8 / (slatop(veg_pp%itype(p)) * leafcn(veg_pp%itype(p)))
               vcmax25top = lnc(p) * flnr(veg_pp%itype(p)) * fnr * act25 *dayl_factor(p)
               vcmax25top = vcmax25top * fnitr(veg_pp%itype(p))
               jmax25top = (2.59_r8 -0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

            else if ( cnallocate_carbonnitrogen_only() ) then ! only N control,from Kattge 2009 Global Change Biology 15 (4), 976-991

               ! Leaf nitrogen concentration at the top of the canopy (g N leaf
               ! / m**2 leaf)
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
                  ! Scale for leaf nitrogen profile. If multi-layer code, use
                  ! explicit
                  ! profile. If sun/shade big leaf code, use canopy integrated
                  ! factor.
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
                  lnc(p) = min(max(lnc(p),0.25_r8),3.0_r8) ! based on doi: 10.1002/ece3.1173
               else
                  lnc(p) = 0.0_r8
               end if

               vcmax25top = (i_vcmax(veg_pp%itype(p)) + s_vcmax(veg_pp%itype(p)) *lnc(p)) * dayl_factor(p)
               jmax25top = (2.59_r8 -0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

            else

               ! nu_com_leaf_physiology is true, vcmax25, jmax25 is derived from
               ! leafn, leafp concentration
               ! Anthony Walker 2014 DOI: 10.1002/ece3.1173

               if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                  ! Leaf nitrogen concentration at the top of the canopy (g N
                  ! leaf / m**2 leaf)
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
                     ! Scale for leaf nitrogen profile. If multi-layer code, use
                     ! explicit
                     ! profile. If sun/shade big leaf code, use canopy
                     ! integrated factor.
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
                     lnc(p) = min(max(lnc(p),0.25_r8),3.0_r8) ! based on doi:10.1002/ece3.1173
                     lpc(p) = min(max(lpc(p),0.014_r8),0.85_r8) ! based on doi:10.1002/ece3.1173
                     vcmax25top = exp(vcmax_np1(veg_pp%itype(p)) + vcmax_np2(veg_pp%itype(p))*log(lnc(p)) + &
                          vcmax_np3(veg_pp%itype(p))*log(lpc(p)) + vcmax_np4(veg_pp%itype(p))*log(lnc(p))*log(lpc(p)))&
                          * dayl_factor(p)
                     jmax25top = exp(jmax_np1 + jmax_np2*log(vcmax25top) + jmax_np3*log(lpc(p))) * dayl_factor(p)
                     vcmax25top = min(max(vcmax25top, 10.0_r8), 150.0_r8)
                     jmax25top = min(max(jmax25top, 10.0_r8), 250.0_r8)
                  else
                     lnc(p) = 0.0_r8
                     lpc(p) = 0.0_r8
                     vcmax25top = 0.0_r8
                     jmax25top = 0.0_r8
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

         ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116,
         ! doi:10.1029/2010JG001593 used
         ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010)
         ! Biogeosciences, 7, 1833-1859
         ! Remove daylength factor from vcmax25 so that kn is based on maximum
         ! vcmax25
         ! But not used as defined here if using sun/shade big leaf code.
         ! Instead,
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
            ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN
            ! s)
            !
            ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
            !
            ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
            ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
            !
            ! Then scale this value at the top of the canopy for canopy depth

            lmr25top = 2.525e-6_r8 * (ParamsShareInst%Q10_mr ** ((25._r8 - 20._r8)/10._r8))
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

         !KO What to do about lmr25 (nscaler?)
         !KO  Is the multi-layer canopy option (nlevcan > 1) still correct here?. 
         !KO  Daniel has just defined nscaler_sun(sha)  = vcmaxcint_sun(sha) and 
         !KO  lmr25_sun(sha) = lmr25top*nscaler_sun(sha)
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
               nscaler_sun = vcmaxcint_sun(p)
               nscaler_sha = vcmaxcint_sha(p)
            else if (nlevcan > 1) then
               nscaler_sun = exp(-kn(p) * laican)
               nscaler_sha = exp(-kn(p) * laican)
            end if

            ! Maintenance respiration

            lmr25_sun = lmr25top * nscaler_sun
            lmr25_sha = lmr25top * nscaler_sha

            if (c3flag(p)) then
               lmr_z_sun(p,iv) = lmr25_sun * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
               lmr_z_sha(p,iv) = lmr25_sha * ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
            else
               lmr_z_sun(p,iv) = lmr25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z_sun(p,iv) = lmr_z_sun(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
               lmr_z_sha(p,iv) = lmr25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z_sha(p,iv) = lmr_z_sha(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
            end if

            if (par_z_sun(p,iv) <= 0._r8) then        ! night time

               vcmax_z(p,sun,iv) = 0._r8
               jmax_z(p,sun,iv) = 0._r8
               tpu_z(p,sun,iv) = 0._r8
               kp_z(p,sun,iv) = 0._r8

               vcmax_z(p,sha,iv) = 0._r8
               jmax_z(p,sha,iv) = 0._r8
               tpu_z(p,sha,iv) = 0._r8
               kp_z(p,sha,iv) = 0._r8

               if ( use_c13 ) then
                  alphapsn_sun(p) = 1._r8
                  alphapsn_sha(p) = 1._r8
               end if

            else                                     ! day time

               vcmax25_sun = vcmax25top * nscaler_sun
               jmax25_sun = jmax25top * nscaler_sun
               tpu25_sun = tpu25top * nscaler_sun        
               vcmax25_sha = vcmax25top * nscaler_sha
               jmax25_sha = jmax25top * nscaler_sha
               tpu25_sha = tpu25top * nscaler_sha       
 
               kp25_sun = kp25top * nscaler_sun
               kp25_sha = kp25top * nscaler_sha

               ! Adjust for temperature

               vcmaxse = 668.39_r8 - 1.07_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               jmaxse  = 659.70_r8 - 0.75_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               tpuse = vcmaxse
               vcmaxc = fth25 (vcmaxhd, vcmaxse)
               jmaxc  = fth25 (jmaxhd, jmaxse)
               tpuc   = fth25 (tpuhd, tpuse)
               vcmax_z(p,sun,iv) = vcmax25_sun * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,sun,iv) = jmax25_sun * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,sun,iv) = tpu25_sun * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)
               vcmax_z(p,sha,iv) = vcmax25_sha * ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
               jmax_z(p,sha,iv) = jmax25_sha * ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
               tpu_z(p,sha,iv) = tpu25_sha * ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

               if (.not. c3flag(p)) then
                  vcmax_z(p,sun,iv) = vcmax25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,sun,iv) = vcmax_z(p,sun,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,sun,iv) = vcmax_z(p,sun,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
                  vcmax_z(p,sha,iv) = vcmax25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,sha,iv) = vcmax_z(p,sha,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,sha,iv) = vcmax_z(p,sha,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
               end if

               kp_z(p,sun,iv) = kp25_sun * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               kp_z(p,sha,iv) = kp25_sha * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)

            end if


         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Leaf-level photosynthesis and stomatal conductance
      !==============================================================================!

      rsmax0 = 2.e4_r8

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)

         ! Leaf boundary layer conductance, umol/m**2/s

         cf = forc_pbot(t)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
         gb = 1._r8/rb(p)
         gb_mol(p) = gb * cf

         ! Loop through canopy layers (above snow). Only do calculations if daytime

         do iv = 1, nrad(p)

            if (par_z_sun(p,iv) <= 0._r8) then        ! night time

               !zqz temporary signal for night time
               vegwp(p,1)=1._r8
               gsminsun = bbb(p)
               gsminsha = bbb(p)


               call calcstress(p,c,t,vegwp(p,:),bsun(p),bsha(p),gb_mol(p),bbb(p),bbb(p), &
                    qsatl(p),qaf(p), atm2lnd_inst,canopystate_inst,waterstate_inst, &
                    soilstate_inst,temperature_inst, waterflux_inst)

               ac(p,sun,iv) = 0._r8
               aj(p,sun,iv) = 0._r8
               ap(p,sun,iv) = 0._r8
               ag(p,sun,iv) = 0._r8
               an_sun(p,iv) = ag(p,sun,iv) - bsun(p) * lmr_z_sun(p,iv)
               psn_z_sun(p,iv) = 0._r8
               psn_wc_z_sun(p,iv) = 0._r8
               psn_wj_z_sun(p,iv) = 0._r8
               psn_wp_z_sun(p,iv) = 0._r8
               !KO  Follow CLM photosynthesis to limit bbb
               rs_z_sun(p,iv) = min(rsmax0, 1._r8/(max( bsun(p)*bbb(p), 1._r8 )) * cf)
               ci_z_sun(p,iv) = 0._r8
               rh_leaf_sun(p) = 0._r8

               ac(p,sha,iv) = 0._r8
               aj(p,sha,iv) = 0._r8
               ap(p,sha,iv) = 0._r8
               ag(p,sha,iv) = 0._r8
               an_sha(p,iv) = ag(p,sha,iv) - bsha(p) * lmr_z_sha(p,iv)
               psn_z_sha(p,iv) = 0._r8
               psn_wc_z_sha(p,iv) = 0._r8
               psn_wj_z_sha(p,iv) = 0._r8
               psn_wp_z_sha(p,iv) = 0._r8
               !KO  Follow CLM photosynthesis to limit bbb
               rs_z_sha(p,iv) = min(rsmax0, 1._r8/(max( bsha(p)*bbb(p), 1._r8 )) * cf)
               ci_z_sha(p,iv) = 0._r8
               rh_leaf_sha(p) = 0._r8

            else                                     ! day time

               !now the constraint is no longer needed, Jinyun Tang
               ceair = min( eair(p),  esat_tv(p) )
               rh_can = ceair / esat_tv(p)


               ! Electron transport rate for C3 plants. Convert par from W/m2 to
               ! umol photons/m**2/s using the factor 4.6

               ! sun
               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z_sun(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,sun,iv))
               cquad = qabs * jmax_z(p,sun,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je_sun = min(r1,r2)

               ! sha
               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z_sha(p,iv) * 4.6_r8
               aquad = theta_psii
               bquad = -(qabs + jmax_z(p,sha,iv))
               cquad = qabs * jmax_z(p,sha,iv)
               call quadratic (aquad, bquad, cquad, r1, r2)
               je_sha = min(r1,r2)

               ! Iterative loop for ci beginning with initial guess

               if (c3flag(p)) then
                  ci_z_sun(p,iv) = 0.7_r8 * cair(p)
                  ci_z_sha(p,iv) = 0.7_r8 * cair(p)
               else
                  ci_z_sun(p,iv) = 0.4_r8 * cair(p)
                  ci_z_sha(p,iv) = 0.4_r8 * cair(p)
               end if

               !find ci and stomatal conductance
               call hybrid_PHS(ci_z_sun(p,iv), ci_z_sha(p,iv), p, iv, c, t, gb_mol(p), bsun(p),bsha(p), je_sun, &
                               je_sha, cair(p), oair(p), lmr_z_sun(p,iv), lmr_z_sha(p,iv), &
                               par_z_sun(p,iv), par_z_sha(p,iv), rh_can, gs_mol_sun(p,iv), gs_mol_sha(p,iv), &
                               qsatl(p), qaf(p), iter1, iter2, atm2lnd_inst, photosyns_inst, &
                               canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst, waterflux_inst)

               gsminsun     = bbb(p)
               gsminsha     = bbb(p)
               gs_slope_sun = mbb(p)
               gs_slope_sha = mbb(p)

               ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

               if (an_sun(p,iv) < 0._r8) gs_mol_sun(p,iv) = max(bsun(p)*gsminsun, 1._r8 )
               if (an_sha(p,iv) < 0._r8) gs_mol_sha(p,iv) = max(bsha(p)*gsminsha, 1._r8 )


               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

               cs_sun = cair(p) - 1.4_r8/gb_mol(p) * an_sun(p,iv) * forc_pbot(t)
               cs_sun = max(cs_sun,1.e-06_r8)
               ci_z_sun(p,iv) = cair(p) - an_sun(p,iv) * forc_pbot(t) * &
                                (1.4_r8*gs_mol_sun(p,iv)+1.6_r8*gb_mol(p)) / &
                                (gb_mol(p)*gs_mol_sun(p,iv))

               cs_sha = cair(p) - 1.4_r8/gb_mol(p) * an_sha(p,iv) * forc_pbot(t)
               cs_sha = max(cs_sha,1.e-06_r8)
               ci_z_sha(p,iv) = cair(p) - an_sha(p,iv) * forc_pbot(t) * &
                                (1.4_r8*gs_mol_sha(p,iv)+1.6_r8*gb_mol(p)) / &
                                (gb_mol(p)*gs_mol_sha(p,iv))

               ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

               gs = gs_mol_sun(p,iv) / cf
               rs_z_sun(p,iv) = min(1._r8/gs, rsmax0)
               gs = gs_mol_sha(p,iv) / cf
               rs_z_sha(p,iv) = min(1._r8/gs, rsmax0)

               ! Photosynthesis. Save rate-limiting photosynthesis

               psn_z_sun(p,iv) = ag(p,sun,iv)

               psn_wc_z_sun(p,iv) = 0._r8
               psn_wj_z_sun(p,iv) = 0._r8
               psn_wp_z_sun(p,iv) = 0._r8

               if (ac(p,sun,iv) <= aj(p,sun,iv) .and. ac(p,sun,iv) <= ap(p,sun,iv)) then
                  psn_wc_z_sun(p,iv) =  psn_z_sun(p,iv)
               else if (aj(p,sun,iv) < ac(p,sun,iv) .and. aj(p,sun,iv) <= ap(p,sun,iv)) then
                  psn_wj_z_sun(p,iv) =  psn_z_sun(p,iv)
               else if (ap(p,sun,iv) < ac(p,sun,iv) .and. ap(p,sun,iv) < aj(p,sun,iv)) then
                  psn_wp_z_sun(p,iv) =  psn_z_sun(p,iv)
               end if

               psn_z_sha(p,iv) = ag(p,sha,iv)
               !psn_z_sha(p,iv) = psn_z_sha(p,iv) * o3coefv_sha(p)

               psn_wc_z_sha(p,iv) = 0._r8
               psn_wj_z_sha(p,iv) = 0._r8
               psn_wp_z_sha(p,iv) = 0._r8

               if (ac(p,sha,iv) <= aj(p,sha,iv) .and. ac(p,sha,iv) <= ap(p,sha,iv)) then
                  psn_wc_z_sha(p,iv) =  psn_z_sha(p,iv)
               else if (aj(p,sha,iv) < ac(p,sha,iv) .and. aj(p,sha,iv) <= ap(p,sha,iv)) then
                  psn_wj_z_sha(p,iv) =  psn_z_sha(p,iv)
               else if (ap(p,sha,iv) < ac(p,sha,iv) .and. ap(p,sha,iv) < aj(p,sha,iv)) then
                  psn_wp_z_sha(p,iv) =  psn_z_sha(p,iv)
               end if

               ! Make sure iterative solution is correct

               if (gs_mol_sun(p,iv) < 0._r8 .or. gs_mol_sha(p,iv) < 0._r8) then
                  write (iulog,*)'Negative stomatal conductance:'
                  write (iulog,*)'p,iv,gs_mol_sun,gs_mol_sha= ',p,iv,gs_mol_sun(p,iv),gs_mol_sha(p,iv)
                  call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
               end if

               ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b

               hs = (gb_mol(p)*ceair + gs_mol_sun(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol_sun(p,iv))*esat_tv(p))
               rh_leaf_sun(p) = hs
               !KO  Follow CLM photosynthesis to limit bbb
               gs_mol_err = mbb(p)*max(an_sun(p,iv), 0._r8)*hs/cs_sun*forc_pbot(t) + max( bsun(p)*bbb(p), 1._r8 )

               if (abs(gs_mol_sun(p,iv)-gs_mol_err) > 1.e-01_r8) then
                  write (iulog,*) 'Ball-Berry error check - sunlit stomatal conductance error:'
                  write (iulog,*) gs_mol_sun(p,iv), gs_mol_err
               end if

               hs = (gb_mol(p)*ceair + gs_mol_sha(p,iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol_sha(p,iv))*esat_tv(p))
               rh_leaf_sha(p) = hs
               !KO  Follow CLM photosynthesis to limit bbb
               gs_mol_err = mbb(p)*max(an_sha(p,iv), 0._r8)*hs/cs_sha*forc_pbot(t) + max( bsha(p)*bbb(p), 1._r8)

               if (abs(gs_mol_sha(p,iv)-gs_mol_err) > 1.e-01_r8) then
                  write (iulog,*) 'Ball-Berry error check - shaded stomatal conductance error:'
                  write (iulog,*) gs_mol_sha(p,iv), gs_mol_err
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

    grav2(1:nlevsoi) = z(c,1:nlevsoi) * 1000._r8
      do f = 1, fn
         p = filterp(f)

         psncan_sun = 0._r8
         psncan_wc_sun = 0._r8
         psncan_wj_sun = 0._r8
         psncan_wp_sun = 0._r8
         lmrcan_sun = 0._r8
         gscan_sun = 0._r8
         laican_sun = 0._r8
         do iv = 1, nrad(p)
            psncan_sun = psncan_sun + psn_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wc_sun = psncan_wc_sun + psn_wc_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wj_sun = psncan_wj_sun + psn_wj_z_sun(p,iv) * lai_z_sun(p,iv)
            psncan_wp_sun = psncan_wp_sun + psn_wp_z_sun(p,iv) * lai_z_sun(p,iv)
            lmrcan_sun = lmrcan_sun + lmr_z_sun(p,iv) * lai_z_sun(p,iv)
            gscan_sun = gscan_sun + lai_z_sun(p,iv) / (rb(p)+rs_z_sun(p,iv))
            laican_sun = laican_sun + lai_z_sun(p,iv)
         end do
         if (laican_sun > 0._r8) then
            psn_sun(p) = psncan_sun / laican_sun
            psn_wc_sun(p) = psncan_wc_sun / laican_sun
            psn_wj_sun(p) = psncan_wj_sun / laican_sun
            psn_wp_sun(p) = psncan_wp_sun / laican_sun
            lmr_sun(p) = lmrcan_sun / laican_sun
            rs_sun(p) = laican_sun / gscan_sun - rb(p)
         else
            psn_sun(p) =  0._r8
            psn_wc_sun(p) =  0._r8
            psn_wj_sun(p) =  0._r8
            psn_wp_sun(p) =  0._r8
            lmr_sun(p) = 0._r8
            rs_sun(p) = 0._r8
         end if
         psncan_sha = 0._r8
         psncan_wc_sha = 0._r8
         psncan_wj_sha = 0._r8
         psncan_wp_sha = 0._r8
         lmrcan_sha = 0._r8
         gscan_sha = 0._r8
         laican_sha = 0._r8
         do iv = 1, nrad(p)
            psncan_sha = psncan_sha + psn_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wc_sha = psncan_wc_sha + psn_wc_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wj_sha = psncan_wj_sha + psn_wj_z_sha(p,iv) * lai_z_sha(p,iv)
            psncan_wp_sha = psncan_wp_sha + psn_wp_z_sha(p,iv) * lai_z_sha(p,iv)
            lmrcan_sha = lmrcan_sha + lmr_z_sha(p,iv) * lai_z_sha(p,iv)
            gscan_sha = gscan_sha + lai_z_sha(p,iv) / (rb(p)+rs_z_sha(p,iv))
            laican_sha = laican_sha + lai_z_sha(p,iv)
         end do
         if (laican_sha > 0._r8) then
            psn_sha(p) = psncan_sha / laican_sha
            psn_wc_sha(p) = psncan_wc_sha / laican_sha
            psn_wj_sha(p) = psncan_wj_sha / laican_sha
            psn_wp_sha(p) = psncan_wp_sha / laican_sha
            lmr_sha(p) = lmrcan_sha / laican_sha
            rs_sha(p) = laican_sha / gscan_sha - rb(p)
         else
            psn_sha(p) =  0._r8
            psn_wc_sha(p) =  0._r8
            psn_wj_sha(p) =  0._r8
            psn_wp_sha(p) =  0._r8
            lmr_sha(p) = 0._r8
            rs_sha(p) = 0._r8
         end if
         
         !KO  Here's how I'm combining bsun and bsha to get btran
         !KO  But this is not really an indication of soil moisture stress that can be
         !KO  used for, e.g., irrigation?
         if ( laican_sha+laican_sun > 0._r8 ) then
            btran(p) = bsun(p) * (laican_sun / (laican_sun + laican_sha)) + &
                       bsha(p) * (laican_sha / (laican_sun + laican_sha))         
         else
            !KO  Btran has a valid value even if there is no exposed lai (elai=0).  
            !KO  In this case, bsun and bsha should have the same value and btran 
            !KO  can be set to either bsun or bsha.  But this needs to be checked.
            btran(p) = bsun(p)
         end if

      end do

    end associate

  end subroutine PhotosynthesisHydraulicStress
  !------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  subroutine hybrid_PHS(x0sun, x0sha, p, iv, c, t, gb_mol, bsun, bsha, jesun, jesha, &
       cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
       gs_mol_sun, gs_mol_sha, qsatl, qaf, iter1, iter2, atm2lnd_inst, photosyns_inst, &
       canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst, waterflux_inst)
    !
    !! DESCRIPTION:
    !use a hybrid solver to find the root of the ci_func equation for sunlit and shaded leaves
    ! f(x) = x- h(x)                                                                                                                                               
    !we want to find x, s.t. f(x) = 0.
    !outside loop iterates for bsun/bsha, which are functions of stomatal conductance
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarantee convergence.
    !
    !! REVISION HISTORY:
    !
    !
    !!USES:
    !
    !! ARGUMENTS:
    implicit none
    real(r8), intent(inout) :: x0sun,x0sha              ! initial guess and final value of the solution for cisun/cisha
    integer , intent(in)    :: p                        ! pft index
    integer , intent(in)    :: iv                       ! radiation canopy layer index
    integer , intent(in)    :: c                        ! column index
    integer , intent(in)    :: t                        ! topounit index
    real(r8), intent(in)    :: gb_mol                   ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: bsun                     ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), intent(out)   :: bsha                     ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), intent(in)    :: jesun                    ! sunlit leaf electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: jesha                    ! shaded leaf electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: cair                     ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in)    :: oair                     ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in)    :: lmr_z_sun                ! sunlit canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: lmr_z_sha                ! shaded canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: par_z_sun                ! par absorbed per unit lai for sunlit canopy layer (w/m**2)
    real(r8), intent(in)    :: par_z_sha                ! par absorbed per unit lai for shaded canopy layer (w/m**2)
    real(r8), intent(in)    :: rh_can                   ! canopy air relative humidity
    real(r8), intent(out)   :: gs_mol_sun               ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: gs_mol_sha               ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in)    :: qsatl                    ! leaf specific humidity [kg/kg]
    real(r8), intent(in)    :: qaf                      ! humidity of canopy air [kg/kg]
    integer,  intent(out)   :: iter1                    ! number of iterations used to find appropriate bsun/bsha
    integer,  intent(out)   :: iter2                    ! number of iterations used to find cisun/cisha
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    !! LOCAL VARIABLES
    real(r8) :: x(nvegwcs) ! working copy of vegwp(p,:)
    real(r8) :: gs0sun   ! unstressed sunlit stomatal conductance
    real(r8) :: gs0sha   ! unstressed shaded stomatal conductance
    logical  :: havegs   ! signals direction of calculation gs->qflx or qflx->gs
    real(r8) :: soilflux ! total soil column transpiration [mm/s] 
    real(r8) :: x1sun    ! second guess for cisun
    real(r8) :: f0sun    ! error of cifunc(x0sun)
    real(r8) :: f1sun    ! error of cifunc(x1sun)
    real(r8) :: xsun     ! open variable for brent to return cisun solution
    real(r8) :: dxsun    ! delta cisun from iter_i to iter_i+1
    real(r8) :: x1sha    ! second guess for cisha
    real(r8) :: f0sha    ! error of cifunc(x0sha)
    real(r8) :: f1sha    ! error of cifunc(x1sha)
    real(r8) :: xsha     ! open variable for brent to return cisha solution
    real(r8) :: dxsha    ! delta cisha from iter_i to iter_i+1
    real(r8) :: b0sun    ! bsun from previous iter
    real(r8) :: b0sha    ! bsha from previous iter
    real(r8) :: dbsun    ! delta(bsun) from iter_i to iter_i+1
    real(r8) :: dbsha    ! delta(bsun) from iter_i to iter_i+1
    logical  :: bflag    ! signals to call calcstress to recalc bsun/bsha (or not)
    real(r8) :: tolsun   ! error tolerance for cisun solution [Pa]
    real(r8) :: tolsha   ! error tolerance for cisun solution [Pa]
    real(r8) :: minf     ! storage spot for best cisun/cisha solution
    real(r8) :: minxsun  ! cisun associated with minf
    real(r8) :: minxsha  ! cisha associated with minf
    real(r8), parameter :: toldb = 1.e-2_r8  ! tolerance for satisfactory bsun/bsha solution
    real(r8), parameter :: eps = 1.e-2_r8    ! relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8    ! absolute accuracy threshold for fsun/fsha
    integer , parameter :: itmax = 3         ! maximum number of iterations zqz (increase later)
    !------------------------------------------------------------------------------
    
    associate(                                                    &
         qflx_tran_veg => veg_wf%qflx_tran_veg    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         vegwp         => canopystate_inst%vegwp_patch            & ! Input/Output: [real(r8) (:,:) ]  vegetation water matric potential (mm)
    )

    
    x1sun = x0sun
    x1sha = x0sha
    bflag = .false.
    b0sun = -1._r8
    b0sha = -1._r8
    bsun  = 1._r8
    bsha  = 1._r8
    iter1 = 0
    
    do                       !outer loop updates bsun/bsha and makes two ci_func calls for interpolation
       x=vegwp(p,:)
       iter1=iter1+1
       iter2=0
       x0sun=max(0.1_r8,x1sun)  !need to make sure x0 .neq. x1
       x1sun=0.99_r8*x1sun
       x0sha=max(0.1_r8,x1sha)
       x1sha=0.99_r8*x1sha
       tolsun = abs(x1sun) * eps
       tolsha = abs(x1sha) * eps
       
       ! this ci_func_PHS call updates bsun/bsha (except on first iter)
       call ci_func_PHS(x,x0sun, x0sha, f0sun, f0sha, p, iv, c, t, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
            temperature_inst, waterflux_inst)
       
       ! update bsun/bsha convergence vars
       dbsun=b0sun-bsun
       dbsha=b0sha-bsha
       b0sun=bsun
       b0sha=bsha
       bflag=.false.
       
       ! this ci_func_PHS call creates second point for ci interpolation
       call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, t, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
            temperature_inst, waterflux_inst)
       
       do                !inner loop finds ci
          if ( (abs(f0sun) < eps1) .and. (abs(f0sha) < eps1) ) then
             x1sun=x0sun
             x1sha=x0sha
             exit
          endif
          if ( (abs(f1sun) < eps1) .and. (abs(f1sha) < eps1) ) then
             exit
          endif
          iter2=iter2+1
          
          if ( (f1sun - f0sun) == 0._r8) then
             !makes next x1sun the midpt between current x1 & x0
             dxsun = 0.5_r8*(x1sun+x0sun)-x1sun
          else
             dxsun=-f1sun*(x1sun-x0sun)/(f1sun-f0sun)
          end if
          if ( (f1sha - f0sha) == 0._r8) then
             dxsha = 0.5_r8*(x1sha+x0sha)-x1sha
          else
             dxsha=-f1sha*(x1sha-x0sha)/(f1sha-f0sha)
          end if
          x0sun=x1sun
          x1sun=x1sun+dxsun
          x0sha=x1sha
          x1sha=x1sha+dxsha
          
          call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, t, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
               gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
               qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
               temperature_inst, waterflux_inst)

          if ( (abs(dxsun) < tolsun ) .and. (abs(dxsha) <tolsha) ) then
             x0sun=x1sun
             x0sha=x1sha
             exit
          endif
          if (iter2 .eq. 1) then
             !initialize best ci vars
             minf=abs(f1sun+f1sha)
             minxsun=x1sun
             minxsha=x1sha
          else
             if (abs(f1sun+f1sha)<minf) then
                !update best ci vars
                minf=abs(f1sun+f1sha)
                minxsun=x1sun
                minxsha=x1sha
             endif
          endif
          
          if ( (abs(f1sun) < eps1) .and. (abs(f1sha) < eps1) ) then
             exit
          endif
          
          if ( (f1sun*f0sun < 0._r8) .and. (f1sha*f0sha < 0._r8) ) then
             
             call brent_PHS(xsun, x0sun, x1sun, f0sun, f1sun, xsha, x0sha, x1sha, f0sha, f1sha, &
                  tolsun, p, iv, c, t, gb_mol, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,&
                  rh_can, gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf, atm2lnd_inst, photosyns_inst, &
                  canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst, waterflux_inst)
             x0sun=xsun
             x0sha=xsha
             exit
          endif
          
          if (iter2 > itmax) then
             x1sun=minxsun
             x1sha=minxsha
             call ci_func_PHS(x,x1sun, x1sha, f1sun, f1sha, p, iv, c, t, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
                  gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
                  qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
                  temperature_inst, waterflux_inst)
             exit
          endif
          
       enddo
       
       !update unstressed stomatal conductance
       if (bsun>0.01_r8) then
          gs0sun=gs_mol_sun/bsun
       endif
       if (bsha>0.01_r8) then
          gs0sha=gs_mol_sha/bsha
       endif
       
       bflag=.true.
       
       if ( (abs(dbsun) < toldb) .and. (abs(dbsha) < toldb) ) then
          exit
       endif
       
       if (iter1 > itmax) then
          exit
       endif
    
    enddo
    x0sun=x1sun
    x0sha=x1sha
    
    !set vegwp for the final gs_mol solution
    call getvegwp(p, c, t, x, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf, soilflux, &
         atm2lnd_inst, canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst)
    vegwp(p,:)=x
    if (soilflux<0._r8) soilflux = 0._r8
    qflx_tran_veg(p) = soilflux
    
    end associate
    
  end subroutine hybrid_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine brent_PHS(xsun, x1sun, x2sun, f1sun, f2sun, xsha, x1sha, x2sha, f1sha, f2sha, &
       tol, ip, iv, ic, it, gb_mol, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha,&
       rh_can, gs_mol_sun, gs_mol_sha, bsun, bsha, qsatl, qaf, atm2lnd_inst, photosyns_inst, &
       canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst, waterflux_inst)
    !------------------------------------------------------------------------------
    implicit none
    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol. Performed for cisun and cisha.
    !
    !!REVISION HISTORY:
    !
    !!ARGUMENTS:
    real(r8), intent(out)   :: xsun                 ! independent variable of the single value function ci_func(x)
    real(r8), intent(in)    :: x1sun, x2sun         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: f1sun, f2sun         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(out)   :: xsha                 ! independent variable of the single value function ci_func(x)
    real(r8), intent(in)    :: x1sha, x2sha         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: f1sha, f2sha         ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in)    :: tol                  ! the error tolerance
    integer , intent(in)    :: ip, iv, ic, it       ! pft, c3/c4, column index, topounit index
    real(r8), intent(in)    :: gb_mol               ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(in)    :: jesun,jesha          ! electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in)    :: cair                 ! Atmospheric CO2 partial pressure (Pa)
    real(r8), intent(in)    :: oair                 ! Atmospheric O2 partial pressure (Pa)
    real(r8), intent(in)    :: lmr_z_sun, lmr_z_sha ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(in)    :: par_z_sun, par_z_sha ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in)    :: rh_can               ! inside canopy relative humidity
    real(r8), intent(out)   :: gs_mol_sun           ! sunlit leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(out)   :: gs_mol_sha           ! shaded leaf stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(inout) :: bsun                 ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), intent(inout) :: bsha                 ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), intent(in)    :: qsatl                ! leaf specific humidity [kg/kg]
    real(r8), intent(in)    :: qaf                  ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)       :: atm2lnd_inst
    type(photosyns_type)   , intent(inout)    :: photosyns_inst
    type(canopystate_type) , intent(inout)    :: canopystate_inst
    type(waterstate_type)  , intent(inout)    :: waterstate_inst
    type(waterflux_type)   , intent(inout)    :: waterflux_inst
    type(soilstate_type)   , intent(inout)    :: soilstate_inst
    type(temperature_type) , intent(in)       :: temperature_inst
    !------------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer                 :: phase                ! sun==1, sha==2
    integer , parameter     :: nphs = 2             ! number of phases for sun/shade
    integer , parameter     :: itmax = 20           ! maximum number of iterations
    real(r8), parameter     :: eps = 1.e-4_r8       ! relative error tolerance
    integer                 :: iter                 !
    real(r8)                :: a(nphs),b(nphs),c(nphs),d(nphs),e(nphs),fa(nphs),fb(nphs),fc(nphs)
    real(r8)                :: p(nphs),q(nphs),r(nphs),s(nphs),tol1(nphs),xm(nphs)
    real(r8)                :: x(nvegwcs)           !dummy variable passed to cifunc
    logical , parameter     :: bflag = .false.      !indicates the cifunc should not call calcstress
    !------------------------------------------------------------------------------
    
    a(:)=(/x1sun,x1sha/)
    b(:)=(/x2sun,x2sha/)
    fa(:)=(/f1sun,f1sha/)
    fb(:)=(/f2sun,f2sha/)
    
    do phase=1, nphs
       if ( (fa(phase) > 0._r8 .and. fb(phase) > 0._r8) .or. (fa(phase) < 0._r8 .and. fb(phase) < 0._r8) ) then
          write(iulog,*) 'root must be bracketed for brent'
          call endrun(msg=errmsg(__FILE__, __LINE__))
       endif
    enddo
    
    c=b
    fc=fb
    iter = 0
    do
       if( iter == itmax ) exit
       iter=iter+1
       
       do phase=1, nphs
          if( (fb(phase) > 0._r8 .and. fc(phase) > 0._r8) .or. (fb(phase) < 0._r8 .and. fc(phase) < 0._r8)) then
             c(phase)=a(phase)   !Rename a, b, c and adjust bounding interval d.
             fc(phase)=fa(phase)
             d(phase)=b(phase)-a(phase)
             e(phase)=d(phase)
          endif
          if( abs(fc(phase)) < abs(fb(phase)) ) then
             a(phase)=b(phase)
             b(phase)=c(phase)
             c(phase)=a(phase)
             fa(phase)=fb(phase)
             fb(phase)=fc(phase)
             fc(phase)=fa(phase)
          endif
       enddo
       tol1=2._r8*eps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       
       if( abs(xm(sun)) <= tol1(sun) .or. fb(sun) == 0._r8 ) then
          if( abs(xm(sha)) <= tol1(sha) .or. fb(sha) == 0._r8 ) then
             xsun=b(sun)
             xsha=b(sha)
             return
          endif
       endif
       
       do phase=1, nphs
          if( abs(e(phase)) >= tol1(phase) .and. abs(fa(phase)) > abs(fb(phase)) ) then
             s(phase)=fb(phase)/fa(phase) !Attempt inverse quadratic interpolation.
             if(a(phase) == c(phase)) then
                p(phase)=2._r8*xm(phase)*s(phase)
                q(phase)=1._r8-s(phase)
             else
                q(phase)=fa(phase)/fc(phase)
                r(phase)=fb(phase)/fc(phase)
                p(phase)=s(phase)*(2._r8*xm(phase)*q(phase)*(q(phase)-r(phase))-(b(phase)-a(phase))*(r(phase)-1._r8))
                q(phase)=(q(phase)-1._r8)*(r(phase)-1._r8)*(s(phase)-1._r8)
             endif
             if( p(phase) > 0._r8 ) q(phase)=-q(phase) !Check whether in bounds.
             p(phase)=abs(p(phase))
             if( 2._r8*p(phase) < min(3._r8*xm(phase)*q(phase)-abs(tol1(phase)*q(phase)),abs(e(phase)*q(phase))) ) then
                e(phase)=d(phase) !Accept interpolation.
                d(phase)=p(phase)/q(phase)
             else
                d(phase)=xm(phase)  !Interpolation failed, use bisection.
                e(phase)=d(phase)
             endif
          else !Bounds decreasing too slowly, use bisection.
             d(phase)=xm(phase)
             e(phase)=d(phase)
          endif
          a(phase)=b(phase) !Move last best guess to a.
          fa(phase)=fb(phase)
          if( abs(d(phase)) > tol1(phase) ) then !Evaluate new trial root.
             b(phase)=b(phase)+d(phase)
          else
             b(phase)=b(phase)+sign(tol1(phase),xm(phase))
          endif
       enddo
       
       call ci_func_PHS(x,b(sun), b(sha), fb(sun), fb(sha), ip, iv, ic, it, bsun, bsha, bflag, gb_mol, gs_mol_sun, gs_mol_sha,&
            gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
            qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
            temperature_inst, waterflux_inst)
       
       if( (fb(sun) == 0._r8) .and. (fb(sha) == 0._r8) ) exit
    enddo
    if( iter == itmax) write(iulog,*) 'brent exceeding maximum iterations', b, fb
    xsun=b(sun)
    xsha=b(sha)
    
    return
    
  end subroutine brent_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine ci_func_PHS(x,cisun, cisha, fvalsun, fvalsha, p, iv, c, t, bsun, bsha, bflag, gb_mol, gs0sun, gs0sha,&
       gs_mol_sun, gs_mol_sha, jesun, jesha, cair, oair, lmr_z_sun, lmr_z_sha, par_z_sun, par_z_sha, rh_can, &
       qsatl, qaf, atm2lnd_inst, photosyns_inst, canopystate_inst, waterstate_inst, soilstate_inst, &
       temperature_inst, waterflux_inst)
    !------------------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! evaluate the function
    ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an for sunlit and shaded leaves
    !
    ! !REVISION HISTORY:
    !
    !
    ! !USES:
    use elm_varpar        , only : nlevsoi
    implicit none
    !
    ! !ARGUMENTS:
    real(r8)               , intent(inout) :: x(nvegwcs)         ! working copy of vegwp(p,:) 
    real(r8)               , intent(in)    :: cisun,cisha        ! intracellular leaf CO2 (Pa)
    real(r8)               , intent(out)   :: fvalsun,fvalsha    ! return function of the value f(ci)
    integer                , intent(in)    :: p,c,iv,t           ! pft, column, radiation indexes, and topounit index
    real(r8)               , intent(inout) :: bsun               ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(inout) :: bsha               ! shaded canopy transpiration wetness factor (0 to 1)
    logical                , intent(in)    :: bflag              ! signals to call calcstress to recalc bsun/bsha (or not)
    real(r8)               , intent(in)    :: gb_mol             ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)    :: gs0sun,gs0sha      ! local gs_mol copies
    real(r8)               , intent(inout) :: gs_mol_sun,gs_mol_sha !leaf stomatal conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)    :: jesun, jesha       ! electron transport rate (umol electrons/m**2/s)
    real(r8)               , intent(in)    :: cair               ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)    :: oair               ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)    :: lmr_z_sun, lmr_z_sha ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8)               , intent(in)    :: par_z_sun, par_z_sha ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8)               , intent(in)    :: rh_can             ! canopy air relative humidity
    real(r8)               , intent(in)    :: qsatl              ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)    :: qaf                ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(canopystate_type) , intent(in)    :: canopystate_inst
    type(waterstate_type)  , intent(in)    :: waterstate_inst
    type(waterflux_type)   , intent(in)    :: waterflux_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(temperature_type) , intent(in)    :: temperature_inst

    ! !LOCAL VARIABLES:
    real(r8) :: ai                   ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: cs_sun,cs_sha        ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: aquad, bquad, cquad  ! terms for quadratic equations
    real(r8) :: r1, r2               ! roots of quadratic equation
    real(r8) :: fnps                 ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii           ! empirical curvature parameter for electron transport rate
    real(r8) :: theta_ip             ! empirical curvature parameter for ap photosynthesis co-limitation
    real(r8) :: term                 ! intermediate in Medlyn stomatal model
    !
    !------------------------------------------------------------------------------
    
    associate(                                                 &
         forc_pbot     => top_as%pbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         c3flag     => photosyns_inst%c3flag_patch           , & ! Input:  [logical  (:)   ]    true if C3 and false if C4
         ac         => photosyns_inst%ac_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_inst%aj_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_inst%ap_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_inst%ag_phs_patch           , & ! Output: [real(r8) (:,:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_inst%vcmax_z_phs_patch      , & ! Input:  [real(r8) (:,:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         cp         => photosyns_inst%cp_patch               , & ! Output: [real(r8) (:)   ]    CO2 compensation point (Pa)
         kc         => photosyns_inst%kc_patch               , & ! Output: [real(r8) (:)   ]    Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_inst%ko_patch               , & ! Output: [real(r8) (:)   ]    Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_inst%qe_patch               , & ! Output: [real(r8) (:)   ]    quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z      => photosyns_inst%tpu_z_phs_patch        , & ! Output: [real(r8) (:,:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_inst%kp_z_phs_patch         , & ! Output: [real(r8) (:,:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj   => photosyns_inst%theta_cj_patch         , & ! Output: [real(r8) (:)   ]    empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_inst%bbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_inst%mbb_patch              , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         an_sun     => photosyns_inst%an_sun_patch           , & ! Output: [real(r8) (:,:) ]  net sunlit leaf photosynthesis (umol CO2/m**2/s)
         an_sha     => photosyns_inst%an_sha_patch             & ! Output: [real(r8) (:,:) ]  net shaded leaf photosynthesis (umol CO2/m**2/s)
         )
    
    !------------------------------------------------------------------------------
    ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    fnps = 0.15_r8
    theta_psii = 0.7_r8
    theta_ip = 0.95_r8
    
    if (bflag) then   !zqz what if bsun==0 ... doesn't break... but follow up

       call calcstress(p,c,t,x,bsun,bsha,gb_mol,gs0sun,gs0sha,qsatl,qaf, &
            atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst, &
            temperature_inst, waterflux_inst)
    endif
    
    if (c3flag(p)) then
       ! C3: Rubisco-limited photosynthesis
       ac(p,sun,iv) = bsun * vcmax_z(p,sun,iv) * max(cisun-cp(p), 0._r8) / (cisun+kc(p)*(1._r8+oair/ko(p)))
       ac(p,sha,iv) = bsha * vcmax_z(p,sha,iv) * max(cisha-cp(p), 0._r8) / (cisha+kc(p)*(1._r8+oair/ko(p)))
       
       ! C3: RuBP-limited photosynthesis
       aj(p,sun,iv) = jesun * max(cisun-cp(p), 0._r8) / (4._r8*cisun+8._r8*cp(p))
       aj(p,sha,iv) = jesha * max(cisha-cp(p), 0._r8) / (4._r8*cisha+8._r8*cp(p))
       
       ! C3: Product-limited photosynthesis
       ap(p,sun,iv) = 3._r8 * tpu_z(p,sun,iv)
       ap(p,sha,iv) = 3._r8 * tpu_z(p,sha,iv)
       
    else
       ! C4: Rubisco-limited photosynthesis
       ac(p,sun,iv) = bsun * vcmax_z(p,sun,iv)
       ac(p,sha,iv) = bsha * vcmax_z(p,sha,iv)
       
       ! C4: RuBP-limited photosynthesis
       aj(p,sun,iv) = qe(p) * par_z_sun * 4.6_r8
       aj(p,sha,iv) = qe(p) * par_z_sha * 4.6_r8
       
       ! C4: PEP carboxylase-limited (CO2-limited)
       ap(p,sun,iv) = kp_z(p,sun,iv) * max(cisun, 0._r8) / forc_pbot(t)
       ap(p,sha,iv) = kp_z(p,sha,iv) * max(cisha, 0._r8) / forc_pbot(t)
       
    end if
    
    ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap
    
    ! Sunlit
    aquad = theta_cj(p)
    bquad = -(ac(p,sun,iv) + aj(p,sun,iv))
    cquad = ac(p,sun,iv) * aj(p,sun,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ai = min(r1,r2)
    
    aquad = theta_ip
    bquad = -(ai + ap(p,sun,iv))
    cquad = ai * ap(p,sun,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ag(p,sun,iv) = max(0._r8,min(r1,r2))
    
    ! Shaded
    aquad = theta_cj(p)
    bquad = -(ac(p,sha,iv) + aj(p,sha,iv))
    cquad = ac(p,sha,iv) * aj(p,sha,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ai = min(r1,r2)
    
    aquad = theta_ip
    bquad = -(ai + ap(p,sha,iv))
    cquad = ai * ap(p,sha,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ag(p,sha,iv) = max(0._r8,min(r1,r2))
    
    ! Net photosynthesis. Exit iteration if an < 0
    an_sun(p,iv) = ag(p,sun,iv) - bsun * lmr_z_sun
    an_sha(p,iv) = ag(p,sha,iv) - bsha * lmr_z_sha
    
    if (an_sun(p,iv) < 0._r8) then
       gs_mol_sun = bbb(p)
       gs_mol_sun = max( bsun*gs_mol_sun, 1._r8)
       fvalsun = 0._r8  ! really tho? zqz
    endif
    if (an_sha(p,iv) < 0._r8) then
       gs_mol_sha = bbb(p)
       gs_mol_sha = max( bsha*gs_mol_sha, 1._r8)
       fvalsha = 0._r8
    endif
    if ((an_sun(p,iv) < 0._r8) .AND. (an_sha(p,iv) < 0._r8)) then
       return
    endif
    
    ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
    ! With an <= 0, then gs_mol = bbb
    
    ! Sunlit
    cs_sun = cair - 1.4_r8/gb_mol * an_sun(p,iv) * forc_pbot(t)
    cs_sun = max(cs_sun,10.e-06_r8)

    aquad = cs_sun
    bquad = cs_sun*(gb_mol - max(bsun*bbb(p),1._r8)) - mbb(p)*an_sun(p,iv)*forc_pbot(t)
    cquad = -gb_mol*(cs_sun*max(bsun*bbb(p),1._r8) + mbb(p)*an_sun(p,iv)*forc_pbot(t)*rh_can)
    call quadratic (aquad, bquad, cquad, r1, r2)
    gs_mol_sun = max(r1,r2)
    
    ! Shaded
    cs_sha = cair - 1.4_r8/gb_mol * an_sha(p,iv) * forc_pbot(t)
    cs_sha = max(cs_sha,10.e-06_r8)
    
    aquad = cs_sha
    bquad = cs_sha*(gb_mol - max(bsha*bbb(p),1._r8)) - mbb(p)*an_sha(p,iv)*forc_pbot(t)
    cquad = -gb_mol*(cs_sha*max(bsha*bbb(p),1._r8) + mbb(p)*an_sha(p,iv)*forc_pbot(t)*rh_can)
    call quadratic (aquad, bquad, cquad, r1, r2)
    gs_mol_sha = max(r1,r2)
    
    ! Derive new estimate for cisun,cisha
    if (an_sun(p,iv) >= 0._r8) then
       if (gs_mol_sun > 0._r8) then
          fvalsun =cisun - cair + an_sun(p,iv) * forc_pbot(t) * (1.4_r8*gs_mol_sun+1.6_r8*gb_mol) / (gb_mol*gs_mol_sun)
       else
          fvalsun =cisun - cair
       endif
    endif
    if (an_sha(p,iv) >= 0._r8) then
       if (gs_mol_sha > 0._r8) then
          fvalsha =cisha - cair + an_sha(p,iv) * forc_pbot(t) * (1.4_r8*gs_mol_sha+1.6_r8*gb_mol) / (gb_mol*gs_mol_sha)
       else
          fvalsha =cisha - cair
       endif
    endif
    end associate
  end subroutine ci_func_PHS
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine calcstress(p,c,t,x,bsun,bsha,gb_mol,gs_mol_sun,gs_mol_sha,qsatl,qaf, &
       atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst, &
       temperature_inst, waterflux_inst)
    !
    ! DESCRIPTIONS
    ! compute the transpiration stress using a plant hydraulics approach
    ! calls spacF, spacA, and getvegwp
    !
    ! USES
    use elm_varpar        , only : nlevsoi
    use elm_varcon        , only : rgas
    !!
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    integer                , intent(in)  :: t               ! topounit index
    real(r8)               , intent(inout)  :: x(nvegwcs)   ! working copy of vegwp(p,:)
    real(r8)               , intent(out) :: bsun            ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out) :: bsha            ! shaded sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(in)  :: gb_mol          ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: gs_mol_sun      ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: gs_mol_sha      ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: qsatl           ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)  :: qaf             ! humidity of canopy air [kg/kg]
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterstate_type)  , intent(in)  :: waterstate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterflux_type)   , intent(in)  :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: A(nvegwcs,nvegwcs)    ! matrix relating d(vegwp) and f: d(vegwp)=A*f 
    real(r8) :: f(nvegwcs)            ! flux divergence (mm/s)
    real(r8) :: dx(nvegwcs)           ! change in vegwp from one iter to the next [mm]
    real(r8) :: efpot                 ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry_sun            ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha            ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: qflx_sun              ! [kg/m2/s]
    real(r8) :: qflx_sha              ! [kg/m2/s]
    real(r8) :: gs0sun,gs0sha         ! local gs_mol copies
    real(r8) :: qsun,qsha             ! attenuated transpiration fluxes
    integer  :: j                     ! index
    real(r8) :: cf                    ! s m**2/umol -> s/m
    integer  :: iter                  ! newton's method iteration number
    logical  :: flag                  ! signal that matrix was not invertible
    logical  :: night                 ! signal to store vegwp within this routine, b/c it is night-time and full suite won't be called
    integer, parameter  :: itmax=50   ! exit newton's method if iters>itmax
    real(r8), parameter :: tolf=1.e-6,toldx=1.e-9 !tolerances for a satisfactory solution
    logical  :: havegs                ! signals direction of calculation gs->qflx or qflx->gs 
    real(r8) :: soilflux              ! total soil column transpiration [mm/s] 
    real(r8), parameter :: tol_lai=.001_r8 ! minimum lai where transpiration is calc'd 
    !------------------------------------------------------------------------------
    
    associate(                                                    &
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         elai          => canopystate_inst%elai_patch           , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai          => canopystate_inst%esai_patch           , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         fdry          => veg_ws%fdry            , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]
         forc_pbot     => top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_rho      => top_as%rhobot                         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
!         forc_rho      => atm2lnd_inst%forc_rho_downscaled_col  , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
!         forc_pbot     => atm2lnd_inst%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         tgcm          => veg_es%thm            , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)
         bsw           => soilstate_inst%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         qflx_tran_veg => veg_wf%qflx_tran_veg    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         sucsat        => soilstate_inst%sucsat_col               & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         )

    !temporary flag for night time vegwp(sun)>0  
    if (x(sun)>0._r8) then
       night=.TRUE.
       x(sun)=x(sha)
    else
       night=.FALSE.
    endif
    
    !copy to avoid rewriting gs_mol_sun
    gs0sun=gs_mol_sun
    gs0sha=gs_mol_sha
    
    !compute transpiration demand
    havegs=.true.
    call getqflx(p,c,t,gb_mol,gs0sun,gs0sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterstate_inst, temperature_inst)
    
    if ((laisun(p)>tol_lai .or. laisha(p)>tol_lai).and.&
         (qflx_sun>0._r8 .or. qflx_sha>0._r8))then

    !newton's method solves for matching fluxes through the spac
    iter=0
    do
       
       iter=iter+1

       call spacF(p,c,x,f,qflx_sun,qflx_sha, &
            atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst,temperature_inst,waterflux_inst)
          
       if ( sqrt(sum(f*f)) < tolf*(qflx_sun+qflx_sha) ) then  !fluxes balanced -> exit
          flag = .false.
          exit
       end if
       if ( iter>itmax ) then                                 !exceeds max iters -> exit
          flag = .false.
          exit
       end if
       
       call spacA(p,c,x,A,qflx_sun,qflx_sha,flag, &
            atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst,temperature_inst,waterflux_inst)

       if (flag) then
          ! cannot invert the matrix, solve for x algebraically assuming no flux                            
          exit
       end if

       if (laisun(p)>tol_lai.and.laisha(p)>tol_lai)then
          dx = matmul(A,f)
       else
          !reduces to 3x3 system
          !in this case, dx is not always [sun,sha,xyl,root]
          !sun and sha flip depending on which is lai==0
          dx(sun)=0._r8
          dx(sha:root)=matmul(A(sha:root,sha:root),f(sha:root))
       endif
       
       
       if ( maxval(abs(dx)) > 50000._r8) then
          dx = 50000._r8 * dx / maxval(abs(dx))  !rescale step to max of 50000
       end if


       if (laisun(p)>tol_lai.and.laisha(p)>tol_lai)then
          x=x+dx
       elseif (laisha(p)>tol_lai) then
          x=x+dx
          x(sun)=x(xyl) ! psi_sun = psi_xyl because laisun==0
       else
          x(xyl:root)=x(xyl:root)+dx(xyl:root)
          x(sun)=x(sun)+dx(sha)  ! implementation ugly bit, chose to flip dx(sun) and dx(sha) for laisha==0 case
          x(sha)=x(xyl) ! psi_sha = psi_xyl because laisha==0
         
       endif


       if ( sqrt(sum(dx*dx)) < toldx) then
          !step in vegwp small -> exit
          exit
       end if
       
       ! this is a catch to force spac gradient to atmosphere
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(sun) > x(xyl) )  x(sun) = x(xyl)
       if ( x(sha) > x(xyl) )  x(sha) = x(xyl)
       
    end do

    else
       ! both qflxsun and qflxsha==0
	flag=.true.
    end if

    if (flag) then
       ! solve algebraically
       call getvegwp(p, c, t, x, gb_mol, gs0sun, gs0sha, qsatl, qaf, soilflux, &
               atm2lnd_inst, canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst)
       bsun = plc(x(sun),p,c,sun,veg)
       bsha = plc(x(sha),p,c,sha,veg)
    else     
    ! compute attenuated flux
    qsun=qflx_sun*plc(x(sun),p,c,sun,veg)
    qsha=qflx_sha*plc(x(sha),p,c,sha,veg)
    
    ! retrieve stressed stomatal conductance
    havegs=.FALSE.
    call getqflx(p,c,t,gb_mol,gs0sun,gs0sha,qsun,qsha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterstate_inst, temperature_inst)
    
    ! compute water stress
    ! .. generally -> B= gs_stressed / gs_unstressed
    ! .. when gs=0 -> B= plc( x )
    if (qflx_sun>0._r8) then
       bsun = gs0sun/gs_mol_sun
    else
       bsun = plc(x(sun),p,c,sun,veg)
    endif
    if (qflx_sha>0._r8) then
       bsha = gs0sha/gs_mol_sha
    else
       bsha = plc(x(sha),p,c,sha,veg)
    endif
    endif
    if ( bsun < 0.01_r8 ) bsun = 0._r8
    if ( bsha < 0.01_r8 ) bsha = 0._r8

    !zqz is this the best place to do this?
    ! was looking like qflx_tran_veg/vegwp was not being set at night time
    ! set vegwp for the final gs_mol solution
    if (night) then
       gs0sun=bsun*gs_mol_sun
       gs0sha=bsha*gs_mol_sha
       call getvegwp(p, c, t, x, gb_mol, gs0sun, gs0sha, qsatl, qaf, soilflux, &
            atm2lnd_inst, canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst)
       if (soilflux<0._r8) soilflux = 0._r8
       qflx_tran_veg(p) = soilflux
    endif
    
    
    end associate
  
  end subroutine calcstress
   
   !------------------------------------------------------------------------------
   
  !------------------------------------------------------------------------------
  subroutine spacA(p,c,x,invA,qflx_sun,qflx_sha,flag, &
       atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst, &
       temperature_inst, waterflux_inst)
    
    !
    ! DESCRIPTION
    !  Returns invA, the inverse matrix relating delta(vegwp) to f
    !   d(vegwp)=invA*f
    !   evaluated at vegwp(p)
    !
    ! The methodology is currently hardcoded for linear algebra assuming the
    ! number of vegetation segments is four. Thus the matrix A and it's inverse
    ! invA are both 4x4 matrices. A more general method could be done using for
    ! example a LINPACK linear algebra solver.
    !
    ! USES
    use elm_varpar        , only : nlevsoi
    use elm_varcon        , only : rgas
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O] 
    real(r8)               , intent(out) :: invA(nvegwcs,nvegwcs)   ! matrix relating d(vegwp) and f: d(vegwp)=invA*f
    real(r8)               , intent(in)  :: qflx_sun        ! Sunlit leaf transpiration [kg/m2/s] 
    real(r8)               , intent(in)  :: qflx_sha        ! Shaded leaf transpiration [kg/m2/s]
    logical                , intent(out) :: flag            ! tells calling function that the matrix is not invertible
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterstate_type)  , intent(in)  :: waterstate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterflux_type)   , intent(in)  :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto1                 ! sunlit transpiration reduction function [-]
    real(r8) :: fsto2                 ! shaded transpiration reduction function [-] 
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-] 
    real(r8) :: dfsto1                ! 1st derivative of fsto1 w.r.t. change in vegwp
    real(r8) :: dfsto2                ! 1st derivative of fsto2 w.r.t. change in vegwp
    real(r8) :: dfx                   ! 1st derivative of fx w.r.t. change in vegwp
    real(r8) :: dfr                   ! 1st derivative of fr w.r.t. change in vegwp
    real(r8) :: A(nvegwcs,nvegwcs)    ! matrix relating vegwp to flux divergence f=A*d(vegwp)
    real(r8) :: leading               ! inverse of determiniant
    real(r8) :: determ                ! determinant of matrix
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: invfactor             ! 
    real(r8), parameter :: tol_lai=.001_r8 ! minimum lai where transpiration is calc'd
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
#ifndef NDEBUG
    ! Only execute this code if DEBUG=TRUE
    if ( nvegwcs /= 4 )then
       call endrun(msg='Error:: this function is hardcoded for 4x4 matrices with nvegwcs==4'//errMsg(__FILE__, __LINE__))
    end if
#endif
    
    associate(                                                    &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         ivt           => veg_pp%itype                             & ! Input:  [integer  (:)   ]  patch vegetation type
         )
    
    ! initialize all elements to zero
    A = 0._r8
    invA = 0._r8

    grav1 = htop(p)*1000._r8
    
    !compute conductance attentuation for each segment
    fsto1=  plc(x(sun),p,c,sun,veg)
    fsto2=  plc(x(sha),p,c,sha,veg)
    fx=     plc(x(xyl),p,c,xyl,veg)
    fr=     plc(x(root),p,c,root,veg)
    
    !compute 1st deriv of conductance attenuation for each segment
    dfsto1=  d1plc(x(sun),p,c,sun,veg)
    dfsto2=  d1plc(x(sha),p,c,sha,veg)
    dfx=     d1plc(x(xyl),p,c,xyl,veg)
    dfr=     d1plc(x(root),p,c,root,veg)
    
    !A - f=A*d(vegwp)
    A(1,1)= - laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx&
         - qflx_sun * dfsto1
    A(1,3)= laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * dfx * (x(xyl)-x(sun))&
         + laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx
    A(2,2)= - laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx&
         - qflx_sha * dfsto2
    A(2,3)= laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * dfx * (x(xyl)-x(sha))&
         + laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx
    A(3,1)= laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx
    A(3,2)= laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx
    A(3,3)= - laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * dfx * (x(xyl)-x(sun))&
         - laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx&
         - laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * dfx * (x(xyl)-x(sha))&
         - laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx&
         - tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr
    A(3,4)= tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * dfr * (x(root)-x(xyl)-grav1)&
         + tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr
    A(4,3)= tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr
    A(4,4)= - tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr&
         - tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * dfr * (x(root)-x(xyl)-grav1)&
         - sum(k_soil_root(p,1:nlevsoi))

    invfactor=1._r8
    A=invfactor*A

    !matrix inversion
    if (laisun(p)>tol_lai .and. laisha(p)>tol_lai) then
       ! general case

       determ=A(4,4)*A(2,2)*A(3,3)*A(1,1) - A(4,4)*A(2,2)*A(3,1)*A(1,3)&
            - A(4,4)*A(3,2)*A(2,3)*A(1,1) - A(4,3)*A(1,1)*A(2,2)*A(3,4)
       if ( abs(determ) <= 1.e-50_r8 ) then
          flag = .true.  !tells calling function that the matrix is not invertible
          return
       else
          flag = .false.
       end if       
    
       leading = 1._r8/determ

       !algebraic inversion of the matrix
       invA(1,1)=leading*A(4,4)*A(2,2)*A(3,3) - leading*A(4,4)*A(3,2)*A(2,3) - leading*A(4,3)*A(2,2)*A(3,4)
       invA(2,1)=leading*A(2,3)*A(4,4)*A(3,1)
       invA(3,1)=-leading*A(4,4)*A(2,2)*A(3,1)
       invA(4,1)=leading*A(4,3)*A(2,2)*A(3,1)
       invA(1,2)=leading*A(1,3)*A(4,4)*A(3,2)
       invA(2,2)=leading*A(4,4)*A(3,3)*A(1,1)-leading*A(4,4)*A(3,1)*A(1,3)-leading*A(4,3)*A(1,1)*A(3,4)
       invA(3,2)=-leading*A(1,1)*A(4,4)*A(3,2)
       invA(4,2)=leading*A(4,3)*A(1,1)*A(3,2)
       invA(1,3)=-leading*A(1,3)*A(2,2)*A(4,4)
       invA(2,3)=-leading*A(2,3)*A(1,1)*A(4,4)
       invA(3,3)=leading*A(2,2)*A(1,1)*A(4,4)
       invA(4,3)=-leading*A(4,3)*A(1,1)*A(2,2)
       invA(1,4)=leading*A(1,3)*A(3,4)*A(2,2)
       invA(2,4)=leading*A(2,3)*A(3,4)*A(1,1)
       invA(3,4)=-leading*A(3,4)*A(1,1)*A(2,2)
       invA(4,4)=leading*A(2,2)*A(3,3)*A(1,1)-leading*A(2,2)*A(3,1)*A(1,3)-leading*A(3,2)*A(2,3)*A(1,1)
       invA=invfactor*invA !undo inversion scaling
    else
       ! if laisun or laisha ==0 invert the corresponding 3x3 matrix
       ! if both are zero, this routine is not called
       if (laisha(p)<=tol_lai) then
          ! shift nonzero matrix values so that both 3x3 cases can be inverted with the same code
          A(2,2)=A(1,1)
          A(3,2)=A(3,1)
          A(2,3)=A(1,3)
       endif
       determ=A(2,2)*A(3,3)*A(4,4)-A(3,4)*A(2,2)*A(4,3)-A(2,3)*A(3,2)*A(4,4)
       if ( abs(determ) <= 1.e-50_r8 ) then
          flag = .true.  !tells calling function that the matrix is not invertible
          return
       else
          flag = .false.
       end if
       
       !algebraic inversion of the 3x3 matrix stored in A(2:4,2:4)
       invA(2,2)=A(3,3)*A(4,4)-A(3,4)*A(4,3)
       invA(2,3)=-A(2,3)*A(4,4)
       invA(2,4)=A(3,4)*A(2,3)
       invA(3,2)=-A(3,2)*A(4,4)
       invA(3,3)=A(2,2)*A(4,4)
       invA(3,4)=-A(3,4)*A(2,2)
       invA(4,2)=A(3,2)*A(4,3)
       invA(4,3)=-A(2,2)*A(4,3)
       invA(4,4)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
       invA=1._r8/determ*invA
       
    endif

    end associate
    
  end subroutine spacA
  
  !--------------------------------------------------------------------------------
  
  !------------------------------------------------------------------------------
  subroutine spacF(p,c,x,f,qflx_sun,qflx_sha, &
       atm2lnd_inst,canopystate_inst,waterstate_inst,soilstate_inst, &
       temperature_inst, waterflux_inst)
    !
    ! DESCRIPTION
    ! Returns f, the flux divergence across each vegetation segment
    !  calculated for vegwp(p,:) as passed in via x
    !
    ! USES
    use elm_varpar        , only : nlevsoi
    use elm_varcon        , only : rgas
    use ColumnType        , only : col_pp
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p               ! pft index
    integer                , intent(in)  :: c               ! column index
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O]
    real(r8)               , intent(out) :: f(nvegwcs)      ! water flux divergence [mm/s]
    real(r8)               , intent(in)  :: qflx_sun        ! Sunlit leaf transpiration [kg/m2/s] 
    real(r8)               , intent(in)  :: qflx_sha        ! Shaded leaf transpiration [kg/m2/s] 
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterstate_type)  , intent(in)  :: waterstate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    type(waterflux_type)   , intent(in)  :: waterflux_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto1                 ! sunlit transpiration reduction function [-]
    real(r8) :: fsto2                 ! shaded transpiration reduction function [-]
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-]
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O) 
    real(r8) :: grav2(nlevsoi)        ! soil layer gravitational potential relative to surface (mm H2O) 
    real(r8) :: temp                  ! used to copy f(sun) to f(sha) for special case
    real(r8), parameter :: tol_lai=.001_r8  ! needs to be the same as in calcstress and spacA (poor form, refactor)<
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
    
    associate(                                              &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input:  [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input:  [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         smp           => soilstate_inst%smp_l_col              , & ! Input: [real(r8) (:,:) ]  soil matrix potential [mm]
         qflx_tran_veg => veg_wf%qflx_tran_veg    , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         z             => col_pp%z                                   & ! Input:  [real(r8) (:,:) ]  layer node depth (m)
         )
    
    grav1 = htop(p) * 1000._r8
    grav2(1:nlevsoi) = z(c,1:nlevsoi) * 1000._r8
    
    fsto1=  plc(x(sun),p,c,sun,veg)
    fsto2=  plc(x(sha),p,c,sha,veg)
    fx=     plc(x(xyl),p,c,xyl,veg)
    fr=     plc(x(root),p,c,root,veg)
    
    !compute flux divergence across each plant segment
    f(sun)= qflx_sun * fsto1 - laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx * (x(xyl)-x(sun))
    f(sha)= qflx_sha * fsto2 - laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx * (x(xyl)-x(sha))
    f(xyl)= laisun(p) * params_inst%kmax(veg_pp%itype(p),sun) * fx * (x(xyl)-x(sun))&
         + laisha(p) * params_inst%kmax(veg_pp%itype(p),sha) * fx * (x(xyl)-x(sha)) &
         - tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr * (x(root)-x(xyl)-grav1)
    f(root)= tsai(p) * params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr * (x(root)-x(xyl)-grav1) &
         + sum( k_soil_root(p,1:nlevsoi) * (x(root)+grav2(1:nlevsoi)) ) &
         - sum( k_soil_root(p,1:nlevsoi) * smp(c,1:nlevsoi) )

    waterflux_inst%sapflow_patch = & 
            tsai(p)*params_inst%kmax(veg_pp%itype(p),xyl) / htop(p) * fr * (x(root)-x(xyl)-grav1)

    if (laisha(p)<tol_lai) then
       ! special case for laisha ~ 0
       ! flip sunlit and shade fluxes to match special case handling in spacA
       temp=f(sun)
       f(sun)=f(sha)
       f(sha)=temp
    endif

    end associate

  end subroutine spacF
  
  !--------------------------------------------------------------------------------
  subroutine getvegwp(p, c, t, x, gb_mol, gs_mol_sun, gs_mol_sha, qsatl, qaf, soilflux, &
       atm2lnd_inst, canopystate_inst, waterstate_inst, soilstate_inst, temperature_inst)
    ! !DESCRIPTION:
    !  Calculates transpiration and returns corresponding vegwp in x
    !
    ! !USES:
    ! calls getqflx
    use elm_varpar  , only : nlevsoi
    use ColumnType  , only : col_pp
    implicit none
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: p                ! pft index
    integer                , intent(in)  :: c                ! column index
    integer                , intent(in)  :: t                ! topounit index
    real(r8)               , intent(out) :: x(nvegwcs)       ! working copy of veg water potential for patch p
    real(r8)               , intent(in)  :: gb_mol           ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(inout)  :: gs_mol_sun    ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(inout)  :: gs_mol_sha    ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)  :: qsatl            ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)  :: qaf              ! humidity of canopy air [kg/kg]
    real(r8)               , intent(out) :: soilflux         ! total soil column transpiration [mm/s]
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterstate_type)  , intent(in)  :: waterstate_inst
    type(soilstate_type)   , intent(in)  :: soilstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: qflx_sun                 ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) :: qflx_sha                 ! Shaded leaf transpiration [kg/m2/s] 
    real(r8) :: fx                       ! fraction of maximum conductance, xylem-to-leaf [-]  
    real(r8) :: fr                       ! fraction of maximum conductance, root-to-xylem [-]  
    real(r8) :: grav1                    ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: grav2(nlevsoi)           ! soil layer gravitational potential relative to surface (mm H2O) 
    integer  :: j                        ! index
    logical  :: havegs                   ! signals direction of calculation gs->qflx or qflx->gs 
    !----------------------------------------------------------------------
    associate(                                                    &
         k_soil_root  => soilstate_inst%k_soil_root_patch       , & ! Input:  [real(r8) (:,:) ]  soil-root interface conductance (mm/s)
         laisun        => canopystate_inst%laisun_patch         , & ! Input: [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input: [real(r8) (:)   ]  shaded leaf area
         htop          => canopystate_inst%htop_patch           , & ! Input: [real(r8) (:)   ]  patch canopy top (m)
         tsai          => canopystate_inst%tsai_patch           , & ! Input:  [real(r8) (:)   ]  patch canopy one-sided stem area index, no burying by snow
         smp           => soilstate_inst%smp_l_col              , & ! Input: [real(r8) (:,:) ]  soil matrix potential [mm]
         rootfr        => soilstate_inst%rootfr_patch           , & ! Input: [real(r8) (:,:) ]  fraction of roots in each soil layer
         bsw           => soilstate_inst%bsw_col                , & ! Input: [real(r8) (:,:) ]  Clapp and Hornberger "b"
         !ivt           => patch%itype                           , & ! Input: [integer  (:)   ]  patch vegetation type
         hk_l          => soilstate_inst%hk_l_col               , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity (mm/s)
         hksat         => soilstate_inst%hksat_col              , & ! Input: [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         sucsat        => soilstate_inst%sucsat_col             , & ! Input: [real(r8) (:,:) ]  minimum soil suction (mm)
         z             => col_pp%z                                   & ! Input: [real(r8) (:,:) ]  layer node depth (m)
         )
    
    grav1 = 1000._r8 *htop(p)
    grav2(1:nlevsoi) = 1000._r8 * z(c,1:nlevsoi)
    
    !compute transpiration demand
    havegs=.true.
    call getqflx(p,c,t,gb_mol,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
         atm2lnd_inst, canopystate_inst, waterstate_inst, temperature_inst)
    
    !calculate root water potential
    if ( abs(sum(k_soil_root(p,1:nlevsoi))) == 0._r8 ) then
       x(root) = sum(smp(c,1:nlevsoi) - grav2)/nlevsoi
    else
       x(root) = (sum(k_soil_root(p,1:nlevsoi)*(smp(c,1:nlevsoi)-grav2))-qflx_sun-qflx_sha) &
                  /sum(k_soil_root(p,1:nlevsoi))
    endif
    
    !calculate xylem water potential
    fr = plc(x(root),p,c,root,veg)
    if ( (tsai(p) > 0._r8) .and. (fr > 0._r8) ) then
       x(xyl) = x(root) - grav1 - (qflx_sun+qflx_sha)/(fr*params_inst%kmax(veg_pp%itype(p),root)/htop(p)*tsai(p))!removed htop conversion
    else
       x(xyl) = x(root) - grav1
    endif
    
    !calculate sun/sha leaf water potential
    fx = plc(x(xyl),p,c,xyl,veg)
    if ( (laisha(p) > 0._r8) .and. (fx > 0._r8) ) then
       x(sha) = x(xyl) - (qflx_sha/(fx*params_inst%kmax(veg_pp%itype(p),xyl)*laisha(p)))
    else
       x(sha) = x(xyl)
    endif
    if ( (laisun(p) > 0._r8) .and. (fx > 0._r8) ) then
       x(sun) = x(xyl) - (qflx_sun/(fx*params_inst%kmax(veg_pp%itype(p),xyl)*laisun(p)))
    else
       x(sun) = x(xyl)
    endif

    !calculate soil flux
    soilflux = 0._r8
    do j = 1,nlevsoi
       soilflux = soilflux + k_soil_root(p,j)*(smp(c,j)-x(root)-grav2(j))
    enddo

    end associate

  end subroutine getvegwp
  
  !--------------------------------------------------------------------------------
  subroutine getqflx(p,c,t,gb_mol,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatl,qaf,havegs, &
       atm2lnd_inst, canopystate_inst, waterstate_inst, temperature_inst)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    use elm_varcon        , only : rgas
    implicit none
    !
    ! !ARGUMENTS:
    integer  , intent(in)     :: p          ! pft index
    integer  , intent(in)     :: c          ! column index
    integer  , intent(in)     :: t          ! topunit index
    real(r8) , intent(in)     :: gb_mol     ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: gs_mol_sun ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: gs_mol_sha ! Ball-Berry leaf conductance (umol H2O/m**2/s)
    real(r8) , intent(inout)  :: qflx_sun   ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) , intent(inout)  :: qflx_sha   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatl      ! leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    logical  , intent(in)     :: havegs     ! signals direction of calculation gs->qflx or qflx->gs
    type(atm2lnd_type)     , intent(in)  :: atm2lnd_inst
    type(canopystate_type) , intent(in)  :: canopystate_inst
    type(waterstate_type)  , intent(in)  :: waterstate_inst
    type(temperature_type) , intent(in)  :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                      ! heat conductance for leaf [m/s]
    real(r8) :: efpot                    ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry_sun               ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cf                       ! s m**2/umol -> s/m
    !----------------------------------------------------------------------
    
    associate(                                                    &
         laisun        => canopystate_inst%laisun_patch         , & ! Input: [real(r8) (:)   ]  sunlit leaf area
         laisha        => canopystate_inst%laisha_patch         , & ! Input: [real(r8) (:)   ]  shaded leaf area
         elai          => canopystate_inst%elai_patch           , & ! Input: [real(r8) (:)   ]  one-sided leaf area index with burying by snow
         esai          => canopystate_inst%esai_patch           , & ! Input: [real(r8) (:)   ]  one-sided stem area index with burying by snow
         fdry          => veg_ws%fdry            , & ! Input: [real(r8) (:)   ]  fraction of foliage that is green and dry [-]
         forc_pbot     => top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_rho      => top_as%rhobot                         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)
         tgcm          => veg_es%thm              & ! Input: [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)
         )
    
    
    cf       = forc_pbot(t)/(rgas*1.e-3_r8*tgcm(p))*1.e6_r8  ! gb->gbmol conversion factor
    wtl      = (elai(p)+esai(p))*gb_mol
    efpot    = forc_rho(t)*wtl*(qsatl-qaf)
    if (havegs) then

       if ( (efpot > 0._r8) .and. (elai(p) > 0._r8) ) then
          if ( gs_mol_sun > 0._r8 ) then
             rppdry_sun = fdry(p)/gb_mol*(laisun(p)/(1._r8/gb_mol+1._r8/gs_mol_sun))/elai(p)
             qflx_sun   = efpot*rppdry_sun/cf
          else
             qflx_sun   = 0._r8
          end if
          if ( gs_mol_sha > 0._r8 ) then
             rppdry_sha = fdry(p)/gb_mol*(laisha(p)/(1._r8/gb_mol+1._r8/gs_mol_sha))/elai(p)
             qflx_sha   = efpot*rppdry_sha/cf
          else
             qflx_sha   = 0._r8
          end if
       else
          qflx_sun      = 0._r8
          qflx_sha      = 0._r8
       end if
       
    else
       if (qflx_sun > 0._r8) then
          gs_mol_sun=gb_mol*qflx_sun*cf*elai(p)/(efpot*fdry(p)*laisun(p)-qflx_sun*cf*elai(p))
       else
          gs_mol_sun=0._r8
       endif
       if (qflx_sha > 0._r8) then
          gs_mol_sha=gb_mol*qflx_sha*cf*elai(p)/(efpot*fdry(p)*laisha(p)-qflx_sha*cf*elai(p))
       else
          gs_mol_sha=0._r8
       endif
       
    endif

    end associate

  end subroutine getqflx

  !--------------------------------------------------------------------------------
  function plc(x,p,c,level,plc_method)
    ! !DESCRIPTION
    ! Return value of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in)  :: x             ! water potential input
    integer  , intent(in)  :: p             ! index for pft
    integer  , intent(in)  :: c             ! index for column
    integer  , intent(in)  :: level         ! veg segment lvl (1:nvegwcs) 
    integer  , intent(in)  :: plc_method    !
    real(r8)               :: plc           ! attenuated conductance [0:1] 0=no flow
    !
    ! !PARAMETERS
    integer , parameter :: vegetation_weibull=0  ! case number
    !------------------------------------------------------------------------------
    associate(                                                    &
         ivt  => veg_pp%itype                             & ! Input: [integer  (:)   ]  patch vegetation type
             )
    
    select case (plc_method)
       !possible to add other methods later
    case (vegetation_weibull)
       if( x >= 0.0_r8 ) then
         plc = 1._r8
       else
         plc=2._r8**(-(x/params_inst%psi50(veg_pp%itype(p),level))**params_inst%ck(veg_pp%itype(p),level))
       end if
       if ( plc < 0.005_r8) plc = 0._r8
    case default
       print *,'must choose plc method'
    end select

    end associate
    
  end function plc
  !--------------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------------
  function d1plc(x,p,c,level,plc_method)
    ! !DESCRIPTION
    ! Return 1st derivative of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in) :: x                ! water potential input
    integer  , intent(in) :: p                ! index for pft
    integer  , intent(in) :: c                ! index for column
    integer  , intent(in) :: level            ! veg segment lvl (1:nvegwcs)
    integer  , intent(in) :: plc_method       ! 0 for vegetation, 1 for soil
    real(r8)              :: d1plc            ! first deriv of plc curve at x
    !
    ! !PARAMETERS
    integer , parameter :: vegetation_weibull=0  ! case number
    !------------------------------------------------------------------------------
    associate(                                                    &
         ivt           => veg_pp%itype                             & ! Input: [integer  (:)   ]  patch vegetation type
             )

    select case (plc_method)
       !possible to add other methods later
    case (vegetation_weibull)
       if( x >= 0.d0 ) then
         d1plc = 0._r8
       else
         d1plc= -params_inst%ck(ivt(p),level) * log(2._r8) * (2._r8**(-(x/params_inst%psi50(ivt(p),level)) &
              **params_inst%ck(ivt(p),level))) &
              * ((x/params_inst%psi50(ivt(p),level))**params_inst%ck(ivt(p),level)) / x
       end if 
    case default
       print *,'must choose plc method'
    end select

    end associate

  end function d1plc

end module PhotosynthesisMod
