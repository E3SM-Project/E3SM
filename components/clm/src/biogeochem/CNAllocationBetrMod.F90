module CNAllocationBetrMod

#include "shr_assert.h"
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in allocation model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : use_c13, use_c14, use_nitrif_denitrif
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use subgridAveMod       , only : p2c
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use PhotosynthesisType  , only : photosyns_type
  use CropType            , only : crop_type
  use EcophysConType      , only : ecophyscon
  use LandunitType        , only : lun                
  use ColumnType          , only : col                
  use PatchType           , only : pft                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readCNAllocBetrParams
  public :: CNAllocationBetrInit         ! Initialization
  public :: calc_plant_nutrient_demand
  public :: plantCNAlloc
  type :: CNAllocParamsType
     real(r8) :: bdnr              ! bulk denitrification rate (1/s)
     real(r8) :: dayscrecover      ! number of days to recover negative cpool
     real(r8) :: compet_plant_no3  ! (unitless) relative compettiveness of plants for NO3
     real(r8) :: compet_plant_nh4  ! (unitless) relative compettiveness of plants for NH4
     real(r8) :: compet_decomp_no3 ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8) :: compet_decomp_nh4 ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8) :: compet_denit      ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8) :: compet_nit        ! (unitless) relative competitiveness of nitrifiers for NH4
  end type CNAllocParamsType
  !
  ! CNAllocParamsInst is populated in readCNAllocParams which is called in 
  type(CNAllocParamsType),protected ::  CNAllocParamsInst
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=*), parameter, public :: suplnAll='ALL'       ! Supplemental Nitrogen for all PFT's
  character(len=*), parameter, public :: suplnNon='NONE'      ! No supplemental Nitrogen
  character(len=15)          , public :: suplnitro = suplnNon ! Supplemental Nitrogen mode
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8)              :: dt                   !decomp timestep (seconds)
  real(r8)              :: bdnr                 !bulk denitrification rate (1/s)
  real(r8)              :: dayscrecover         !number of days to recover negative cpool
  real(r8), allocatable :: arepr(:)             !reproduction allocation coefficient
  real(r8), allocatable :: aroot(:)             !root allocation coefficient
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNAllocBetrParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNAllocParamsType'                  !
    character(len=100) :: errCode = '-Error reading in parameters file:' !
    logical            :: readv                                          ! has variable been read in or not
    real(r8)           :: tempr                                          ! temporary to read in parameter
    character(len=100) :: tString                                        ! temp. var for reading
    !-----------------------------------------------------------------------

    ! read in parameters

    tString='bdnr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%bdnr=tempr

    tString='dayscrecover'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%dayscrecover=tempr

    tString='compet_plant_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_plant_no3=tempr

    tString='compet_plant_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_plant_nh4=tempr
   
    tString='compet_decomp_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_decomp_no3=tempr

    tString='compet_decomp_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_decomp_nh4=tempr
   
    tString='compet_denit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_denit=tempr

    tString='compet_nit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNAllocParamsInst%compet_nit=tempr   

  end subroutine readCNAllocBetrParams

  !-----------------------------------------------------------------------
  subroutine CNAllocationBetrInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size
    use clm_varpar      , only: crop_prog
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'CNAllocationInit'
    logical           :: carbon_only
    !-----------------------------------------------------------------------

    if ( crop_prog )then
       allocate(arepr(bounds%begp:bounds%endp)); arepr(bounds%begp : bounds%endp) = nan
       allocate(aroot(bounds%begp:bounds%endp)); aroot(bounds%begp : bounds%endp) = nan
    end if


    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    bdnr         = CNAllocParamsInst%bdnr * (dt/secspday)
    dayscrecover = CNAllocParamsInst%dayscrecover

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
       Carbon_only = .false.
    case(suplnAll)
       Carbon_only = .true.
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(__FILE__, __LINE__))
    end select

  end subroutine CNAllocationBetrInit

  
!-------------------------------------------------------------------------------

 subroutine plantCNAlloc(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars,  cnstate_vars,  carbonstate_vars, carbonflux_vars, &
       c13_carbonflux_vars, c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars)
  !
  ! DESCRIPTION
  !
  ! do plant productivity downregulation after considering nutrient limitation
 

  !
  ! !USES:
  use shr_sys_mod      , only: shr_sys_flush
  use clm_varctl       , only: iulog, cnallocate_carbon_only
  use pftvarcon        , only: npcropmin, declfact, bfact, aleaff, arootf, astemf
  use pftvarcon        , only: arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
  use clm_varpar       , only: nlevsoi, nlevdecomp
  use clm_varcon       , only: nitrif_n2o_loss_frac, secspday
  use landunit_varcon  , only: istsoil, istcrop
  use clm_time_manager , only: get_step_size
  use clm_varctl       , only: use_c13, use_c14
 implicit none
  !
  ! !ARGUMENTS:
  type(bounds_type)        , intent(in)    :: bounds
  integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
  integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
  integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
  integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
  type(photosyns_type)     , intent(in)    :: photosyns_vars
  type(cnstate_type)       , intent(inout) :: cnstate_vars
  type(carbonstate_type)   , intent(in)    :: carbonstate_vars
  type(carbonflux_type)    , intent(inout) :: carbonflux_vars
  type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
  type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
  type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
  type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars

    !
    ! !LOCAL VARIABLES:

    !
  integer :: c,p,l,pi,j                                              !indices
  integer :: fp                                                      !lake filter pft index
  integer :: fc                                                      !lake filter column index
  real(r8):: f1,f2,f3,f4,g1,g2                                       !allocation parameters
  real(r8):: cnl,cnfr,cnlw,cndw                                      !C:N ratios for leaf, fine root, and wood
  real(r8):: fcur                                                    !fraction of current psn displayed as growth
  real(r8):: gresp_storage                                           !temporary variable for growth resp to storage
  real(r8):: nlc                                                     !temporary variable for total new leaf carbon allocation
  real(r8):: f5                                                      !grain allocation parameter
  real(r8):: cng                                                     !C:N ratio for grain (= cnlw for now; slevis)   
    !-----------------------------------------------------------------------

  associate(                                                                               &
     ivt                          => pft%itype                                           , & ! Input:  [integer  (:) ]  pft vegetation type                                
         
     woody                        => ecophyscon%woody                                    , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
     froot_leaf                   => ecophyscon%froot_leaf                               , & ! Input:  [real(r8) (:)   ]  allocation parameter: new fine root C per new leaf C (gC/gC)
     croot_stem                   => ecophyscon%croot_stem                               , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
     stem_leaf                    => ecophyscon%stem_leaf                                , & ! Input:  [real(r8) (:)   ]  allocation parameter: new stem c per new leaf C (gC/gC)
     flivewd                      => ecophyscon%flivewd                                  , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
     leafcn                       => ecophyscon%leafcn                                   , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)                        
     frootcn                      => ecophyscon%frootcn                                  , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)                   
     livewdcn                     => ecophyscon%livewdcn                                 , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
     deadwdcn                     => ecophyscon%deadwdcn                                 , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
     fcur2                        => ecophyscon%fcur                                     , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
     graincn                      => ecophyscon%graincn                                  , & ! Input:  [real(r8) (:)   ]  grain C:N (gC/gN)                       
     psnsun                       => photosyns_vars%psnsun_patch                         , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
     psnsha                       => photosyns_vars%psnsha_patch                         , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)         
     leafc                        => carbonstate_vars%leafc_patch                        , & ! Input:  [real(r8) (:)   ]                                          
     frootc                       => carbonstate_vars%frootc_patch                       , & ! Input:  [real(r8) (:)   ]                                          
     livestemc                    => carbonstate_vars%livestemc_patch                    , & ! Input:  [real(r8) (:)   ]                                                   
     croplive                     => cnstate_vars%croplive_patch                         , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested     
     aleaf                        => cnstate_vars%aleaf_patch                            , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient             
     astem                        => cnstate_vars%astem_patch                            , & ! Output: [real(r8) (:)   ]  stem allocation coefficient             
     fpg                          => cnstate_vars%fpg_col                                , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)    
     c_allometry                  => cnstate_vars%c_allometry_patch                      , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)                
     n_allometry                  => cnstate_vars%n_allometry_patch                      , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)                =
     downreg                      => cnstate_vars%downreg_patch                          , & ! Output: [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)

     annsum_npp                   => carbonflux_vars%annsum_npp_patch                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation  
     gpp                          => carbonflux_vars%gpp_before_downreg_patch            , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
     availc                       => carbonflux_vars%availc_patch                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
     excess_cflux                 => carbonflux_vars%excess_cflux_patch                  , & ! Output: [real(r8) (:)   ]  C flux not allocated due to downregulation (gC/m2/s)
     plant_calloc                 => carbonflux_vars%plant_calloc_patch                  , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)        
     psnsun_to_cpool              => carbonflux_vars%psnsun_to_cpool_patch               , & ! Output: [real(r8) (:)   ]
     psnshade_to_cpool            => carbonflux_vars%psnshade_to_cpool_patch             , & ! Output: [real(r8) (:)   ]

     cpool_to_leafc               => carbonflux_vars%cpool_to_leafc_patch                , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_leafc_storage       => carbonflux_vars%cpool_to_leafc_storage_patch        , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_frootc              => carbonflux_vars%cpool_to_frootc_patch               , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_frootc_storage      => carbonflux_vars%cpool_to_frootc_storage_patch       , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_livestemc           => carbonflux_vars%cpool_to_livestemc_patch            , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_livestemc_storage   => carbonflux_vars%cpool_to_livestemc_storage_patch    , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_deadstemc           => carbonflux_vars%cpool_to_deadstemc_patch            , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_deadstemc_storage   => carbonflux_vars%cpool_to_deadstemc_storage_patch    , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_livecrootc          => carbonflux_vars%cpool_to_livecrootc_patch           , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_livecrootc_storage  => carbonflux_vars%cpool_to_livecrootc_storage_patch   , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_deadcrootc          => carbonflux_vars%cpool_to_deadcrootc_patch           , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_deadcrootc_storage  => carbonflux_vars%cpool_to_deadcrootc_storage_patch   , & ! Output: [real(r8) (:)   ]                                          
     cpool_to_gresp_storage       => carbonflux_vars%cpool_to_gresp_storage_patch        , & ! Output: [real(r8) (:)   ]  allocation to growth respiration storage (gC/m2/s)
     cpool_to_grainc              => carbonflux_vars%cpool_to_grainc_patch               , & ! Output: [real(r8) (:)   ]  allocation to grain C (gC/m2/s)         
     cpool_to_grainc_storage      => carbonflux_vars%cpool_to_grainc_storage_patch       , & ! Output: [real(r8) (:)   ]  allocation to grain C storage (gC/m2/s) 
       
     retransn                     => nitrogenstate_vars%retransn_patch                   , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N  

     plant_ndemand                => nitrogenflux_vars%plant_ndemand_patch               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
     plant_nalloc                 => nitrogenflux_vars%plant_nalloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)        
     npool_to_grainn              => nitrogenflux_vars%npool_to_grainn_patch             , & ! Output: [real(r8) (:)   ]  allocation to grain N (gN/m2/s)         
     npool_to_grainn_storage      => nitrogenflux_vars%npool_to_grainn_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to grain N storage (gN/m2/s) 
     retransn_to_npool            => nitrogenflux_vars%retransn_to_npool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
     sminn_to_npool               => nitrogenflux_vars%sminn_to_npool_patch              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
     npool_to_leafn               => nitrogenflux_vars%npool_to_leafn_patch              , & ! Output: [real(r8) (:)   ]  allocation to leaf N (gN/m2/s)          
     npool_to_leafn_storage       => nitrogenflux_vars%npool_to_leafn_storage_patch      , & ! Output: [real(r8) (:)   ]  allocation to leaf N storage (gN/m2/s)  
     npool_to_frootn              => nitrogenflux_vars%npool_to_frootn_patch             , & ! Output: [real(r8) (:)   ]  allocation to fine root N (gN/m2/s)     
     npool_to_frootn_storage      => nitrogenflux_vars%npool_to_frootn_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to fine root N storage (gN/m2/s)
     npool_to_livestemn           => nitrogenflux_vars%npool_to_livestemn_patch          , & ! Output: [real(r8) (:)   ]                                          
     npool_to_livestemn_storage   => nitrogenflux_vars%npool_to_livestemn_storage_patch  , & ! Output: [real(r8) (:)   ]                                          
     npool_to_deadstemn           => nitrogenflux_vars%npool_to_deadstemn_patch          , & ! Output: [real(r8) (:)   ]                                          
     npool_to_deadstemn_storage   => nitrogenflux_vars%npool_to_deadstemn_storage_patch  , & ! Output: [real(r8) (:)   ]                                          
     npool_to_livecrootn          => nitrogenflux_vars%npool_to_livecrootn_patch         , & ! Output: [real(r8) (:)   ]                                          
     npool_to_livecrootn_storage  => nitrogenflux_vars%npool_to_livecrootn_storage_patch , & ! Output: [real(r8) (:)   ]                                          
     npool_to_deadcrootn          => nitrogenflux_vars%npool_to_deadcrootn_patch         , & ! Output: [real(r8) (:)   ]                                          
     npool_to_deadcrootn_storage  => nitrogenflux_vars%npool_to_deadcrootn_storage_patch , & ! Output: [real(r8) (:)   ]                                          
     frootn_to_retransn           => nitrogenflux_vars%frootn_to_retransn_patch          , & ! Output: [real(r8) (:)   ]                                          
     sminn_to_plant               => nitrogenflux_vars%sminn_to_plant_col                , & ! Output: [real(r8) (:)   ]                                          

     c13cf => c13_carbonflux_vars, &
     c14cf => c14_carbonflux_vars  &
   )
 
 
      ! start new pft loop to distribute the available N between the
      ! competing patches on the basis of relative demand, and allocate C and N to
      ! new growth and storage

   do fp=1,num_soilp
      p = filter_soilp(fp)
      c = pft%column(p)

         ! set some local allocation variables
      f1 = froot_leaf(ivt(p))
      f2 = croot_stem(ivt(p))

      ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
      ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
      ! There was an error in this formula in previous version, where the coefficient
      ! was 0.004 instead of 0.0025.
      ! This variable allocation is only for trees. Shrubs have a constant
      ! allocation as specified in the pft-physiology file.  The value is also used
      ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
      if (stem_leaf(ivt(p)) == -1._r8) then
         f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
      else
         f3 = stem_leaf(ivt(p))
      end if

      f4 = flivewd(ivt(p))
      g1 = grperc(ivt(p))
      g2 = grpnow(ivt(p))
      cnl = leafcn(ivt(p))
      cnfr = frootcn(ivt(p))
      cnlw = livewdcn(ivt(p))
      cndw = deadwdcn(ivt(p))
      fcur = fcur2(ivt(p))

      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         if (croplive(p)) then
            f1 = aroot(p) / aleaf(p)
            f3 = astem(p) / aleaf(p)
            f5 = arepr(p) / aleaf(p)
            g1 = 0.25_r8
         else
            f1 = 0._r8
            f3 = 0._r8
            f5 = 0._r8
            g1 = 0.25_r8
         end if
      end if

      ! increase fcur linearly with ndays_active, until fcur reaches 1.0 at
      ! ndays_active = days/year.  This prevents the continued storage of C and N.
      ! turning off this correction (PET, 12/11/03), instead using bgtr in
      ! phenology algorithm.
      !fcur = fcur + (1._r8 - fcur)*lgsf(p)
      sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
      plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)


      ! calculate the associated carbon allocation, and the excess
      ! carbon flux that must be accounted for through downregulation
      plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
      excess_cflux(p) = availc(p) - plant_calloc(p)

      ! reduce gpp fluxes due to N limitation
      if (gpp(p) > 0.0_r8) then
         downreg(p) = excess_cflux(p)/gpp(p)
         psnsun_to_cpool(p)   = psnsun_to_cpool(p)  *(1._r8 - downreg(p))
         psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1._r8 - downreg(p))
         if ( use_c13 ) then
            c13cf%psnsun_to_cpool_patch(p)   = c13cf%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
            c13cf%psnshade_to_cpool_patch(p) = c13cf%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
         endif

         if ( use_c14 ) then
            c14cf%psnsun_to_cpool_patch(p)   = c14cf%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
            c14cf%psnshade_to_cpool_patch(p) = c14cf%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
         endif
      end if

      ! calculate the amount of new leaf C dictated by these allocation
      ! decisions, and calculate the daily fluxes of C and N to current
      ! growth and storage pools

      ! fcur is the proportion of this day's growth that is displayed now,
      ! the remainder going into storage for display next year through the
      ! transfer pools

      nlc = plant_calloc(p) / c_allometry(p)

      cpool_to_leafc(p)          = nlc * fcur
      cpool_to_leafc_storage(p)  = nlc * (1._r8 - fcur)
      cpool_to_frootc(p)         = nlc * f1 * fcur
      cpool_to_frootc_storage(p) = nlc * f1 * (1._r8 - fcur)
      if (woody(ivt(p)) == 1._r8) then
         cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
         cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
         cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
         cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
         cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
         cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
         cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
         cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
         cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_grainc(p)             = nlc * f5 * fcur
         cpool_to_grainc_storage(p)     = nlc * f5 * (1._r8 -fcur)
      end if

      ! corresponding N fluxes
      npool_to_leafn(p)          = (nlc / cnl) * fcur
      npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur)
      npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
      npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
      if (woody(ivt(p)) == 1._r8) then
         npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
         npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
         npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         cng = graincn(ivt(p))
         npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
         npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
         npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_grainn(p)             = (nlc * f5 / cng) * fcur
         npool_to_grainn_storage(p)     = (nlc * f5 / cng) * (1._r8 -fcur)
      end if

      ! Calculate the amount of carbon that needs to go into growth
      ! respiration storage to satisfy all of the storage growth demands.
      ! Allows for the fraction of growth respiration that is released at the
      ! time of fixation, versus the remaining fraction that is stored for
      ! release at the time of display. Note that all the growth respiration
      ! fluxes that get released on a given timestep are calculated in growth_resp(),
      ! but that the storage of C for growth resp during display of transferred
      ! growth is assigned here.

      gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p)
      if (woody(ivt(p)) == 1._r8) then
         gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_grainc_storage(p)
      end if
      cpool_to_gresp_storage(p) = gresp_storage * g1 * (1._r8 - g2)

   end do ! end pft loop

  end associate 
 end subroutine plantCNAlloc


 !-----------------------------------------------------------------------------
 subroutine calc_plant_nutrient_demand(bounds, num_soilc, filter_soilc,  num_soilp, filter_soilp,&
       photosyns_vars, crop_vars, canopystate_vars,                           &
       cnstate_vars, carbonstate_vars, carbonflux_vars,                       &
       c13_carbonflux_vars, c14_carbonflux_vars,                              &
       nitrogenstate_vars,  nitrogenflux_vars)!, plantsoilnutrientflux_vars )

  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CanopyStateType     , only : canopystate_type
  use CanopyStateType     , only : canopystate_type
  use PhotosynthesisType  , only : photosyns_type
!  use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type

  implicit none
  ! !ARGUMENTS:
  type(bounds_type)        , intent(in)           :: bounds                     !
  integer                  , intent(in)           :: num_soilc                  ! number of soil columns in filter
  integer                  , intent(in)           :: filter_soilc(:)            ! filter for soil columns
  integer                  , intent(in)           :: num_soilp                  ! number of soil patches in filter
  integer                  , intent(in)           :: filter_soilp(:)            ! filter for soil patches
  type(photosyns_type)     , intent(in)           :: photosyns_vars             !
  type(crop_type)          , intent(in)           :: crop_vars                  !
  type(canopystate_type)   , intent(in)           :: canopystate_vars           !
  type(carbonstate_type)   , intent(in)           :: carbonstate_vars           !
  type(cnstate_type)       , intent(inout)        :: cnstate_vars               !
  type(carbonflux_type)    , intent(inout)        :: carbonflux_vars            !
  type(carbonflux_type)    , intent(inout)        :: c13_carbonflux_vars        !
  type(carbonflux_type)    , intent(inout)        :: c14_carbonflux_vars        !
  type(nitrogenstate_type) , intent(inout)        :: nitrogenstate_vars         !
  type(nitrogenflux_type)  , intent(inout)        :: nitrogenflux_vars          !
!  type(plantsoilnutrientflux_type), intent(inout) :: plantsoilnutrientflux_vars !
  
  call calc_plant_nitrogen_demand(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
    photosyns_vars, canopystate_vars, crop_vars, carbonstate_vars, &
    cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
    c13_carbonflux_vars, c14_carbonflux_vars) !,                             &
!    plantsoilnutrientflux_vars%plant_totn_demand_flx_col(bounds%begc:bounds%endc))
  
  !this can used to plug in phosphorus?   
 end subroutine calc_plant_nutrient_demand

 !-----------------------------------------------------------------------------
 subroutine calc_plant_nitrogen_demand(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
   photosyns_vars, canopystate_vars, crop_vars, carbonstate_vars, &
    cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
    c13_carbonflux_vars, c14_carbonflux_vars) !, plant_totn_demand_flx_col)
  !
  ! DESCRIPTION
  ! compute plant nitrogen demand
  !

  ! !USES:
  use pftvarcon           , only : npcropmin, declfact, bfact, aleaff, arootf, astemf
  use pftvarcon           , only : arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
  use clm_varcon          , only : secspday
  use clm_varctl          , only : use_c13, use_c14 
  use clm_time_manager    , only : get_step_size
  use subgridAveMod       , only : p2c    
  implicit none
 
  ! !ARGUMENTS:
  type(bounds_type)        , intent(in)    :: bounds                                             !
  integer                  , intent(in)    :: num_soilc                                          ! number of soil columns in filter
  integer                  , intent(in)    :: filter_soilc(:)                                    ! filter for soil columns
  integer                  , intent(in)    :: num_soilp                                          ! number of soil patches in filter
  integer                  , intent(in)    :: filter_soilp(:)                                    ! filter for soil patches
  type(photosyns_type)     , intent(in)    :: photosyns_vars                                     !
  type(crop_type)          , intent(in)    :: crop_vars                                          !
  type(canopystate_type)   , intent(in)    :: canopystate_vars                                   !
  type(carbonstate_type)   , intent(in)    :: carbonstate_vars                                   !
  type(cnstate_type)       , intent(inout) :: cnstate_vars                                       !
  type(carbonflux_type)    , intent(inout) :: carbonflux_vars                                    !
  type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars                                !
  type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars                                !
  type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars                                 !
  type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars                                  !
!  real(r8)                 , intent(inout) :: plant_totn_demand_flx_col(bounds%begc:bounds%endc) !
  real(r8)                  :: plant_totn_demand_flx_col(bounds%begc:bounds%endc) !
  integer                                  :: c,p,l,pi,j                                         ! indices
  integer                                  :: fp                                                 ! lake filter pft index
  integer                                  :: fc                                                 ! lake filter column index
  real(r8)                                 :: mr                                                 ! maintenance respiration (gC/m2/s)
  real(r8)                                 :: f1,f2,f3,f4,g1,g2                                  ! allocation parameters
  real(r8)                                 :: cnl,cnfr,cnlw,cndw                                 ! C:N ratios for leaf, fine root, and wood
  real(r8)                                 :: curmr, curmr_ratio                                 ! xsmrpool temporary variables
  real(r8)                                 :: f5                                                 ! grain allocation parameter
  real(r8)                                 :: cng                                                ! C:N ratio for grain (= cnlw for now; slevis)
  real(r8)                                 :: fleaf                                              ! fraction allocated to leaf
  real(r8)                                 :: t1                                                 ! temporary variable
  real(r8)                                 :: dt                                                 ! model time step
  real(r8)                                 :: dayscrecover                                       !
  
  associate(                                                                              &
    ivt                          => pft%itype                                           , & ! Input:  [integer  (:) ]  pft vegetation type                                       
    woody                        => ecophyscon%woody                                    , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
    froot_leaf                   => ecophyscon%froot_leaf                               , & ! Input:  [real(r8) (:)   ]  allocation parameter: new fine root C per new leaf C (gC/gC)
    croot_stem                   => ecophyscon%croot_stem                               , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
    stem_leaf                    => ecophyscon%stem_leaf                                , & ! Input:  [real(r8) (:)   ]  allocation parameter: new stem c per new leaf C (gC/gC)
    flivewd                      => ecophyscon%flivewd                                  , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
    leafcn                       => ecophyscon%leafcn                                   , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)                        
    frootcn                      => ecophyscon%frootcn                                  , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)                   
    livewdcn                     => ecophyscon%livewdcn                                 , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
    deadwdcn                     => ecophyscon%deadwdcn                                 , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
    graincn                      => ecophyscon%graincn                                  , & ! Input:  [real(r8) (:)   ]  grain C:N (gC/gN)                       
    fleafcn                      => ecophyscon%fleafcn                                  , & ! Input:  [real(r8) (:)   ]  leaf c:n during organ fill              
    ffrootcn                     => ecophyscon%ffrootcn                                 , & ! Input:  [real(r8) (:)   ]  froot c:n during organ fill             
    fstemcn                      => ecophyscon%fstemcn                                  , & ! Input:  [real(r8) (:)   ]  stem c:n during organ fill              

    psnsun                       => photosyns_vars%psnsun_patch                         , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
    psnsha                       => photosyns_vars%psnsha_patch                         , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
    c13_psnsun                   => photosyns_vars%c13_psnsun_patch                     , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
    c13_psnsha                   => photosyns_vars%c13_psnsha_patch                     , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
    c14_psnsun                   => photosyns_vars%c14_psnsun_patch                     , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
    c14_psnsha                   => photosyns_vars%c14_psnsha_patch                     , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         
    laisun                       => canopystate_vars%laisun_patch                       , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index        
    laisha                       => canopystate_vars%laisha_patch                       , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index        

    hui                          => crop_vars%gddplant_patch                            , & ! Input:  [real(r8) (:)   ]  =gdd since planting (gddplant)          
    leafout                      => crop_vars%gddtsoi_patch                             , & ! Input:  [real(r8) (:)   ]  =gdd from top soil layer temperature    

    xsmrpool                     => carbonstate_vars%xsmrpool_patch                     , & ! Input:  [real(r8) (:)   ]  (gC/m2) temporary photosynthate C pool  
    leafc                        => carbonstate_vars%leafc_patch                        , & ! Input:  [real(r8) (:)   ]                                          
    frootc                       => carbonstate_vars%frootc_patch                       , & ! Input:  [real(r8) (:)   ]                                          
    livestemc                    => carbonstate_vars%livestemc_patch                    , & ! Input:  [real(r8) (:)   ]                                          
         
    gddmaturity                  => cnstate_vars%gddmaturity_patch                      , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest                   
    huileaf                      => cnstate_vars%huileaf_patch                          , & ! Input:  [real(r8) (:)   ]  heat unit index needed from planting to leaf emergence
    huigrain                     => cnstate_vars%huigrain_patch                         , & ! Input:  [real(r8) (:)   ]  same to reach vegetative maturity       
    croplive                     => cnstate_vars%croplive_patch                         , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested     
    peaklai                      => cnstate_vars%peaklai_patch                          , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max        
    aleafi                       => cnstate_vars%aleafi_patch                           , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
    astemi                       => cnstate_vars%astemi_patch                           , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
    aleaf                        => cnstate_vars%aleaf_patch                            , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient             
    astem                        => cnstate_vars%astem_patch                            , & ! Output: [real(r8) (:)   ]  stem allocation coefficient             
    grain_flag                   => cnstate_vars%grain_flag_patch                       , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not             
    c_allometry                  => cnstate_vars%c_allometry_patch                      , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)                
    n_allometry                  => cnstate_vars%n_allometry_patch                      , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)                
    tempsum_potential_gpp        => cnstate_vars%tempsum_potential_gpp_patch            , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP   
    tempmax_retransn             => cnstate_vars%tempmax_retransn_patch                 , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)
    annsum_potential_gpp         => cnstate_vars%annsum_potential_gpp_patch             , & ! Output: [real(r8) (:)   ]  annual sum of potential GPP             
    annmax_retransn              => cnstate_vars%annmax_retransn_patch                  , & ! Output: [real(r8) (:)   ]  annual max of retranslocated N pool     

    leaf_mr                      => carbonflux_vars%leaf_mr_patch                       , & ! Input:  [real(r8) (:)   ]                                          
    froot_mr                     => carbonflux_vars%froot_mr_patch                      , & ! Input:  [real(r8) (:)   ]                                          
    livestem_mr                  => carbonflux_vars%livestem_mr_patch                   , & ! Input:  [real(r8) (:)   ]                                          
    livecroot_mr                 => carbonflux_vars%livecroot_mr_patch                  , & ! Input:  [real(r8) (:)   ]                                          
    grain_mr                     => carbonflux_vars%grain_mr_patch                      , & ! Input:  [real(r8) (:)   ]                                          
    annsum_npp                   => carbonflux_vars%annsum_npp_patch                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation  
    gpp                          => carbonflux_vars%gpp_before_downreg_patch            , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
    availc                       => carbonflux_vars%availc_patch                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
    xsmrpool_recover             => carbonflux_vars%xsmrpool_recover_patch              , & ! Output: [real(r8) (:)   ]  C flux assigned to recovery of negative cpool (gC/m2/s)
    psnsun_to_cpool              => carbonflux_vars%psnsun_to_cpool_patch               , & ! Output: [real(r8) (:)   ]
    psnshade_to_cpool            => carbonflux_vars%psnshade_to_cpool_patch             , & ! Output: [real(r8) (:)   ]

    leaf_curmr                   => carbonflux_vars%leaf_curmr_patch                    , &
    froot_curmr                  => carbonflux_vars%froot_curmr_patch                   , & ! Output: [real(r8) (:)   ]                                          
    livestem_curmr               => carbonflux_vars%livestem_curmr_patch                , & ! Output: [real(r8) (:)   ]                                          
    livecroot_curmr              => carbonflux_vars%livecroot_curmr_patch               , & ! Output: [real(r8) (:)   ]                                          
    grain_curmr                  => carbonflux_vars%grain_curmr_patch                   , & ! Output: [real(r8) (:)   ]                                          
    leaf_xsmr                    => carbonflux_vars%leaf_xsmr_patch                     , & ! Output: [real(r8) (:)   ]                                          
    froot_xsmr                   => carbonflux_vars%froot_xsmr_patch                    , & ! Output: [real(r8) (:)   ]                                          
    livestem_xsmr                => carbonflux_vars%livestem_xsmr_patch                 , & ! Output: [real(r8) (:)   ]                                          
    livecroot_xsmr               => carbonflux_vars%livecroot_xsmr_patch                , & ! Output: [real(r8) (:)   ]                                          
    grain_xsmr                   => carbonflux_vars%grain_xsmr_patch                    , & ! Output: [real(r8) (:)   ]                                          
    cpool_to_xsmrpool            => carbonflux_vars%cpool_to_xsmrpool_patch             , & ! Output: [real(r8) (:)   ]                                                  
    retransn                     => nitrogenstate_vars%retransn_patch                   , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N  
    plant_ndemand                => nitrogenflux_vars%plant_ndemand_patch               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
    avail_retransn               => nitrogenflux_vars%avail_retransn_patch              , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
    retransn_to_npool            => nitrogenflux_vars%retransn_to_npool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
    leafn_to_retransn            => nitrogenflux_vars%leafn_to_retransn_patch           , & ! Output: [real(r8) (:)   ]                                          
    frootn_to_retransn           => nitrogenflux_vars%frootn_to_retransn_patch          , & ! Output: [real(r8) (:)   ]                                          
    livestemn_to_retransn        => nitrogenflux_vars%livestemn_to_retransn_patch       , & ! Output: [real(r8) (:)   ]                                                          
    c13cf => c13_carbonflux_vars, &
    c14cf => c14_carbonflux_vars  &
  )

  ! set time steps
  dt = real( get_step_size(), r8 )

  dayscrecover = CNAllocParamsInst%dayscrecover
  ! loop over patches to assess the total plant N demand
  do fp=1,num_soilp
    p = filter_soilp(fp)

    ! get the time step total gross photosynthesis
    ! this is coming from the canopy fluxes code, and is the
    ! gpp that is used to control stomatal conductance.
    ! For the nitrogen downregulation code, this is assumed
    ! to be the potential gpp, and the actual gpp will be
    ! reduced due to N limitation. 

    ! Convert psn from umol/m2/s -> gC/m2/s

    ! The input psn (psnsun and psnsha) are expressed per unit LAI
    ! in the sunlit and shaded canopy, respectively. These need to be
    ! scaled by laisun and laisha to get the total gpp for allocation

    ! Note that no associate statement is used for the isotope carbon fluxes below 
    ! since they are not always allocated AND nag compiler will complain if you try to
    ! to have an associate statement with unallocated memory

    psnsun_to_cpool(p)   = psnsun(p) * laisun(p) * 12.011e-6_r8
    psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6_r8

    if ( use_c13 ) then
      c13cf%psnsun_to_cpool_patch(p)   = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
      c13cf%psnshade_to_cpool_patch(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
    endif

    if ( use_c14 ) then
      c14cf%psnsun_to_cpool_patch(p)   = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
      c14cf%psnshade_to_cpool_patch(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
    endif

    gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

    ! get the time step total maintenance respiration
    ! These fluxes should already be in gC/m2/s

    mr = leaf_mr(p) + froot_mr(p)
    if (woody(ivt(p)) == 1.0_r8) then
      mr = mr + livestem_mr(p) + livecroot_mr(p)
    else if (ivt(p) >= npcropmin) then
      if (croplive(p)) mr = mr + livestem_mr(p) + grain_mr(p)
    end if

    ! carbon flux available for allocation
    availc(p) = gpp(p) - mr

    ! new code added for isotope calculations, 7/1/05, PET
    ! If mr > gpp, then some mr comes from gpp, the rest comes from
    ! cpool (xsmr)
    if (mr > 0._r8 .and. availc(p) < 0._r8) then
      curmr = gpp(p)
      curmr_ratio = curmr / mr
    else
      curmr_ratio = 1._r8
    end if
    leaf_curmr(p) = leaf_mr(p) * curmr_ratio
    leaf_xsmr(p) = leaf_mr(p) - leaf_curmr(p)
    froot_curmr(p) = froot_mr(p) * curmr_ratio
    froot_xsmr(p) = froot_mr(p) - froot_curmr(p)
    livestem_curmr(p) = livestem_mr(p) * curmr_ratio
    livestem_xsmr(p) = livestem_mr(p) - livestem_curmr(p)
    livecroot_curmr(p) = livecroot_mr(p) * curmr_ratio
    livecroot_xsmr(p) = livecroot_mr(p) - livecroot_curmr(p)
    grain_curmr(p) = grain_mr(p) * curmr_ratio
    grain_xsmr(p) = grain_mr(p) - grain_curmr(p)

    ! no allocation when available c is negative
    availc(p) = max(availc(p),0.0_r8)

   ! test for an xsmrpool deficit
   if (xsmrpool(p) < 0.0_r8) then
     ! Running a deficit in the xsmrpool, so the first priority is to let
     ! some availc from this timestep accumulate in xsmrpool.
     ! Determine rate of recovery for xsmrpool deficit

     xsmrpool_recover(p) = -xsmrpool(p)/(dayscrecover*secspday)
     if (xsmrpool_recover(p) < availc(p)) then
       ! available carbon reduced by amount for xsmrpool recovery
       availc(p) = availc(p) - xsmrpool_recover(p)
     else
       ! all of the available carbon goes to xsmrpool recovery
       xsmrpool_recover(p) = availc(p)
       availc(p) = 0.0_r8
     end if
   cpool_to_xsmrpool(p) = xsmrpool_recover(p)
  end if

  f1 = froot_leaf(ivt(p))
  f2 = croot_stem(ivt(p))

  ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
  ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
  ! This variable allocation is only for trees. Shrubs have a constant
  ! allocation as specified in the pft-physiology file.  The value is also used
  ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

  if (stem_leaf(ivt(p)) == -1._r8) then
    f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
  else
    f3 = stem_leaf(ivt(p))
  end if

  f4   = flivewd(ivt(p))
  g1   = grperc(ivt(p))
  g2   = grpnow(ivt(p))
  cnl  = leafcn(ivt(p))
  cnfr = frootcn(ivt(p))
  cnlw = livewdcn(ivt(p))
  cndw = deadwdcn(ivt(p))

  ! calculate f1 to f5 for prog crops following AgroIBIS subr phenocrop

  f5 = 0._r8 ! continued intializations from above

  if (ivt(p) >= npcropmin) then ! skip 2 generic crops

    if (croplive(p)) then
      ! same phases appear in subroutine CropPhenology

      ! Phase 1 completed:
      ! ==================
      ! if hui is less than the number of gdd needed for filling of grain
      ! leaf emergence also has to have taken place for lai changes to occur
      ! and carbon assimilation
      ! Next phase: leaf emergence to start of leaf decline

      if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p)) then

        ! allocation rules for crops based on maturity and linear decrease
        ! of amount allocated to roots over course of the growing season

        if (peaklai(p) == 1) then ! lai at maximum allowed
          arepr(p) = 0._r8
          aleaf(p) = 1.e-5_r8
          astem(p) = 0._r8
          aroot(p) = 1._r8 - arepr(p) - aleaf(p) - astem(p)
        else
          arepr(p) = 0._r8
          aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -   &
              (arooti(ivt(p)) - arootf(ivt(p))) *  &
              min(1._r8, hui(p)/gddmaturity(p))))
          fleaf = fleafi(ivt(p)) * (exp(-bfact(ivt(p))) -         &
              exp(-bfact(ivt(p))*hui(p)/huigrain(p))) / &
                (exp(-bfact(ivt(p)))-1) ! fraction alloc to leaf (from J Norman alloc curve)
          aleaf(p) = max(1.e-5_r8, (1._r8 - aroot(p)) * fleaf)
          astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
        end if

        ! AgroIBIS included here an immediate adjustment to aleaf & astem if the 
        ! predicted lai from the above allocation coefficients exceeded laimx.
        ! We have decided to live with lais slightly higher than laimx by
        ! enforcing the cap in the following tstep through the peaklai logic above.

        astemi(p) = astem(p) ! save for use by equations after shift
        aleafi(p) = aleaf(p) ! to reproductive phenology stage begins
        grain_flag(p) = 0._r8 ! setting to 0 while in phase 2

        ! Phase 2 completed:
        ! ==================
        ! shift allocation either when enough gdd are accumulated or maximum number
        ! of days has elapsed since planting

      else if (hui(p) >= huigrain(p)) then

        aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
                       (arooti(ivt(p)) - arootf(ivt(p))) * min(1._r8, hui(p)/gddmaturity(p))))
        if (astemi(p) > astemf(ivt(p))) then
          astem(p) = max(0._r8, max(astemf(ivt(p)), astem(p) * &
                   (1._r8 - min((hui(p)-                 &
                   huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                    huigrain(p)),1._r8)**allconss(ivt(p)) )))
        end if
        if (aleafi(p) > aleaff(ivt(p))) then
          aleaf(p) = max(1.e-5_r8, max(aleaff(ivt(p)), aleaf(p) * &
                   (1._r8 - min((hui(p)-                    &
                     huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                     huigrain(p)),1._r8)**allconsl(ivt(p)) )))
        end if

        !Beth's retranslocation of leafn, stemn, rootn to organ
        !Filter excess plant N to retransn pool for organ N
        !Only do one time then hold grain_flag till onset next season

        ! slevis: Will astem ever = astemf exactly?
        ! Beth's response: ...looks like astem can equal astemf under the right circumstances. 
        !It might be worth a rewrite to capture what I was trying to do, but the retranslocation for 
        !corn and wheat begins at the beginning of the grain fill stage, but for soybean I was holding it 
        !until after the leaf and stem decline were complete. Looking at how astem is calculated, once the 
        !stem decline is near complete, astem should (usually) be set to astemf. The reason for holding off 
        !on soybean is that the retranslocation scheme begins at the beginning of the grain phase, when the 
        !leaf and stem are still growing, but declining. Since carbon is still getting allocated and now 
        !there is more nitrogen available, the nitrogen can be diverted from grain. For corn and wheat 
        !the impact was probably enough to boost productivity, but for soybean the nitrogen was better off 
        !fulfilling the grain fill. It seems that if the peak lai is reached for soybean though that this 
        !would be bypassed altogether, not the intended outcome. I checked several of my output files and 
        !they all seemed to be going through the retranslocation loop for soybean - good news.

        if (ivt(p) /= nsoybean .or. astem(p) == astemf(ivt(p))) then
          if (grain_flag(p) == 0._r8) then
            t1 = 1 / dt
            leafn_to_retransn(p) = t1 * ((leafc(p) / leafcn(ivt(p))) - (leafc(p) / &
                             fleafcn(ivt(p))))
            livestemn_to_retransn(p) = t1 * ((livestemc(p) / livewdcn(ivt(p))) - (livestemc(p) / &
                             fstemcn(ivt(p))))
            frootn_to_retransn(p) = 0._r8
            if (ffrootcn(ivt(p)) > 0._r8) then
              frootn_to_retransn(p) = t1 * ((frootc(p) / frootcn(ivt(p))) - (frootc(p) / &
                                ffrootcn(ivt(p))))
            end if
            grain_flag(p) = 1._r8
          end if
        end if

        arepr(p) = 1._r8 - aroot(p) - astem(p) - aleaf(p)

      else                   ! pre emergence
        aleaf(p) = 1.e-5_r8 ! allocation coefficients should be irrelevant
        astem(p) = 0._r8    ! because crops have no live carbon pools;
        aroot(p) = 0._r8    ! this applies to this "else" and to the "else"
        arepr(p) = 0._r8    ! a few lines down
      end if

      f1 = aroot(p) / aleaf(p)
      f3 = astem(p) / aleaf(p)
      f5 = arepr(p) / aleaf(p)
      g1 = 0.25_r8

    else   ! .not croplive
      f1 = 0._r8
      f3 = 0._r8
      f5 = 0._r8
      g1 = 0.25_r8
    end if
  end if

  ! based on available C, use constant allometric relationships to
  ! determine N requirements

  if (woody(ivt(p)) == 1.0_r8) then
    c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
    n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
  else if (ivt(p) >= npcropmin) then ! skip generic crops
    cng = graincn(ivt(p))
    c_allometry(p) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
    n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
  else
    c_allometry(p) = 1._r8+g1+f1+f1*g1
    n_allometry(p) = 1._r8/cnl + f1/cnfr
  end if
  plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))

  ! retranslocated N deployment depends on seasonal cycle of potential GPP
  ! (requires one year run to accumulate demand)

  tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

  ! Adding the following line to carry max retransn info to CN Annual Update
  tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))

  ! Beth's code: crops pull from retransn pool only during grain fill;
  !              retransn pool has N from leaves, stems, and roots for
  !              retranslocation

  if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
    avail_retransn(p) = plant_ndemand(p)
  else if (ivt(p) < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
    avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
  else
    avail_retransn(p) = 0.0_r8
  end if

  ! make sure available retrans N doesn't exceed storage
  avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)

  ! modify plant N demand according to the availability of
  ! retranslocated N
  ! take from retransn pool at most the flux required to meet
  ! plant ndemand

  if (plant_ndemand(p) > avail_retransn(p)) then
    retransn_to_npool(p) = avail_retransn(p)
  else
    retransn_to_npool(p) = plant_ndemand(p)
  end if
  plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)

  end do ! end pft loop

  ! now use the p2c routine to get the column-averaged plant_ndemand
  call p2c(bounds, num_soilc, filter_soilc, &
           plant_ndemand(bounds%begp:bounds%endp), &
           plant_totn_demand_flx_col(bounds%begc:bounds%endc))
           
  ! obtain the nutrient uptake potential based on fine root profile
           
  end associate 
 end subroutine calc_plant_nitrogen_demand
   
end module CNAllocationBetrMod
