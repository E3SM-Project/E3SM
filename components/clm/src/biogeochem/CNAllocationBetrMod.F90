module CNAllocationBeTRMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in allocation model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : use_c13, use_c14, use_nitrif_denitrif, spinup_state
  use clm_varctl          , only : nyears_ad_carbon_only
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use subgridAveMod       , only : p2c
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  !!! add phosphorus
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use CNStateType         , only : cnstate_type
  use PhotosynthesisType  , only : photosyns_type
  use CropType            , only : crop_type
  use EcophysConType      , only : ecophyscon
  use LandunitType        , only : lun
  use ColumnType          , only : col
  use PatchType           , only : pft
  ! bgc interface & pflotran module switches
  use clm_varctl          , only: use_bgc_interface,use_clm_bgc, use_pflotran, pf_cmode
  use clm_varctl          , only : nu_com
  use SoilStatetype       , only : soilstate_type
  use WaterStateType      , only : waterstate_type
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  !
  implicit none
  save

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNAllocationBeTRInit         ! Initialization
  private:: calc_plantN_kineticpar

  !!-----------------------------------------------------------------------------------------------------
  !! CNAllocation is divided into 3 subroutines/phases:
  private :: CNAllocation1_PlantNPDemand     !!Plant N/P Demand;       called in CNEcosystemDynNoLeaching1
  public :: CNAllocation3_PlantCNPAlloc     !!Plant C/N/P Allocation; called in CNDecompAlloc2
  !!-----------------------------------------------------------------------------------------------------
  private :: dynamic_plant_alloc        ! dynamic plant carbon allocation based on different nutrient stress
  public :: update_plant_stoichiometry ! prognostic plant C:N:P stoichiometry

  !
  !
  ! !PUBLIC DATA MEMBERS:
  character(len=*), parameter, public :: suplnAll='ALL'  ! Supplemental Nitrogen for all PFT's
  character(len=*), parameter, public :: suplnNon='NONE' ! No supplemental Nitrogen
  character(len=15), public :: suplnitro = suplnNon      ! Supplemental Nitrogen mode
  !! add phosphorus  - X. YANG
  character(len=*), parameter, public :: suplpAll='ALL'  ! Supplemental Phosphorus for all PFT's
  character(len=*), parameter, public :: suplpNon='NONE' ! No supplemental Phosphorus
  character(len=15), public :: suplphos = suplpAll    ! Supplemental Phosphorus mode
  !! add competition, - Q. Zhu
  logical,          public :: nu_com_leaf_physiology = .false.
  logical,          public :: nu_com_root_kinetics = .false.
  logical,          public :: nu_com_phosphatase = .false.
  logical,          public :: nu_com_nfix = .false.

  !
  ! !PRIVATE DATA MEMBERS:
  real(r8)              :: dt                   !decomp timestep (seconds)
  real(r8)              :: dayscrecover         !number of days to recover negative cpool
  real(r8), allocatable :: arepr(:)             !reproduction allocation coefficient
  real(r8), allocatable :: aroot(:)             !root allocation coefficient
  real(r8), allocatable :: col_plant_ndemand(:) !column-level plant N demand
  real(r8), allocatable :: col_plant_pdemand(:) !column-level plant P demand

  logical :: crop_supln  = .false.             !Prognostic crop receives supplemental Nitrogen

  real(r8), allocatable :: decompmicc(:,:)                ! column-level soil microbial decomposer biomass gC/m3

  real(r8), parameter   :: E_plant_scalar  = 0.0000125_r8 ! scaling factor for plant fine root biomass to calculate nutrient carrier enzyme abundance
  real(r8), parameter   :: E_decomp_scalar = 0.05_r8      ! scaling factor for plant fine root biomass to calculate nutrient carrier enzyme abundance

  real(r8)              :: e_km_nh4                       ! temp variable of sum(E/KM) for NH4 competition BGC mode
  real(r8)              :: e_km_no3                       ! temp variable of sum(E/KM) for NO3 competition BGC mode
  real(r8)              :: e_km_p                         ! temp variable of sum(E/KM) for P competition
  real(r8)              :: e_km_n                         ! temp variable of sum(E/KM) for N competition CN mode

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNAllocationBeTRInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varcon      , only: secspday
    use clm_time_manager, only: get_step_size, get_curr_date
    use clm_varpar      , only: crop_prog
    use clm_varctl      , only: iulog, cnallocate_carbon_only_set
    use clm_varctl      , only: cnallocate_carbonnitrogen_only_set
    use clm_varctl      , only: cnallocate_carbonphosphorus_only_set
    use shr_infnan_mod  , only: nan => shr_infnan_nan, assignment(=)
    use clm_varpar      , only: nlevdecomp
    use CNAllocationMod , only: CNAllocParamsInst
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'CNAllocationBeTRInit'
    integer :: yr, mon, day, sec
    logical :: carbon_only
    logical :: carbonnitrogen_only
    logical :: carbonphosphorus_only
    !-----------------------------------------------------------------------

    if ( crop_prog )then
       allocate(arepr(bounds%begp:bounds%endp)); arepr(bounds%begp : bounds%endp) = nan
       allocate(aroot(bounds%begp:bounds%endp)); aroot(bounds%begp : bounds%endp) = nan
    end if
    allocate(col_plant_ndemand(bounds%begc:bounds%endc)); col_plant_ndemand(bounds%begc : bounds%endc) = nan
    allocate(col_plant_pdemand(bounds%begc:bounds%endc)); col_plant_pdemand(bounds%begc : bounds%endc) = nan
    allocate(decompmicc(bounds%begc:bounds%endc,1:nlevdecomp)); decompmicc(bounds%begc:bounds%endc,1:nlevdecomp) = nan

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    dayscrecover = CNAllocParamsInst%dayscrecover

    ! Change namelist settings into private logical variables
    select case(suplnitro)
    case(suplnNon)
        select case (suplphos)
          case(suplpNon)
             Carbon_only = .false.
             CarbonNitrogen_only = .false.
             CarbonPhosphorus_only=.false.
             crop_supln  = .false.
          case(suplpAll)
             Carbon_only = .false.
             CarbonNitrogen_only = .true.
             CarbonPhosphorus_only=.false.
             crop_supln  = .false.
        end select
    case(suplnAll)
        select case (suplphos)
          case(suplpNon)
             Carbon_only = .false.
             CarbonNitrogen_only = .false.
             CarbonPhosphorus_only=.true.
             crop_supln  = .false.
          case(suplpAll)
             Carbon_only = .true.
             CarbonNitrogen_only = .false.
             CarbonPhosphorus_only=.false.
             crop_supln  = .false.
        end select
    case default
       write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
            suplnNon, ' or ', suplnAll
       call endrun(msg='ERROR: supplemental Nitrogen flag is not correct'//&
            errMsg(__FILE__, __LINE__))
    end select

    ! phosphorus conditions of plants are needed, in order to use new fixation and phosphatase
    ! activity subroutines, under carbon only or carbon nitrogen only mode, fixation and phosphatase
    ! activity are set to false
    if (carbon_only .or. carbonnitrogen_only) then
        nu_com_nfix = .false.
        nu_com_phosphatase = .false.
    end if

    call get_curr_date(yr, mon, day, sec)
    if (spinup_state == 1 .and. yr .le. nyears_ad_carbon_only) then
      Carbon_only = .true.
     end if

    call cnallocate_carbon_only_set(carbon_only)
    call cnallocate_carbonnitrogen_only_set(carbonnitrogen_only)
    call cnallocate_carbonphosphorus_only_set(carbonphosphorus_only)

  end subroutine CNAllocationBeTRInit
!!-------------------------------------------------------------------------------------------------
  subroutine SetPlantMicNPDemand(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars, crop_vars, canopystate_vars, cnstate_vars,             &
       carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,  &
       nitrogenstate_vars, nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars, PlantMicKinetics_vars)
  implicit none

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars

   !calculate the plant nutrient demand
   call CNAllocation1_PlantNPDemand(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars, crop_vars, canopystate_vars, cnstate_vars,             &
       carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,  &
       nitrogenstate_vars, nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)

   !extract the kinetic parameters
   call calc_plantN_kineticpar(bounds, num_soilc, filter_soilc               , &
                            num_soilp, filter_soilp                         , &
                            cnstate_vars                                    , &
                            carbonstate_vars                                , &
                            nitrogenstate_vars                              , &
                            phosphorusstate_vars                            , &
                            carbonflux_vars                                 , &
                            PlantMicKinetics_vars                             )

  end subroutine SetPlantMicNPDemand
!!-------------------------------------------------------------------------------------------------
  subroutine CNAllocation1_PlantNPDemand (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars, crop_vars, canopystate_vars, cnstate_vars,             &
       carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,  &
       nitrogenstate_vars, nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
    !! PHASE-1 of CNAllocation: loop over patches to assess the total plant N demand and P demand
    ! !USES:
    use shr_sys_mod      , only: shr_sys_flush
    use clm_varctl       , only: iulog,cnallocate_carbon_only,cnallocate_carbonnitrogen_only,&
                                 cnallocate_carbonphosphorus_only
    use pftvarcon        , only: npcropmin, declfact, bfact, aleaff, arootf, astemf, noveg
    use pftvarcon        , only: arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
    use clm_varpar       , only: nlevdecomp
    use clm_varcon       , only: nitrif_n2o_loss_frac, secspday
    use clm_varctl       , only: cnallocate_carbon_only_set
!    use landunit_varcon  , only: istsoil, istcrop
    use clm_time_manager , only: get_step_size, get_curr_date
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!     !!  add phosphorus  -X.YANG
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:

    real(r8) :: fpi_no3_vr(bounds%begc:bounds%endc,1:nlevdecomp) ! fraction of potential immobilization supplied by no3(no units)
    real(r8) :: fpi_nh4_vr(bounds%begc:bounds%endc,1:nlevdecomp) ! fraction of potential immobilization supplied by nh4 (no units)
    real(r8) :: sum_nh4_demand_vr(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_nh4_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_vr(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)

    real(r8) :: sum_pdemand_scaled(bounds%begc:bounds%endc,1:nlevdecomp)  ! sum of total P demand, scaled with relative competitiveness
    real(r8) :: excess_immob_nh4_vr(bounds%begc:bounds%endc,1:nlevdecomp) ! nh4 excess flux, if soil microbes are more P limited
    real(r8) :: excess_immob_no3_vr(bounds%begc:bounds%endc,1:nlevdecomp) ! no3 excess flux, if soil microbes are more P limited
    real(r8) :: excess_immob_p_vr(bounds%begc:bounds%endc,1:nlevdecomp)   ! P excess flux, if soil microbes are more N limited
    real(r8) :: compet_plant_no3(bounds%begp:bounds%endp)                 ! (unitless) relative compettiveness of plants for NO3 BGC mode
    real(r8) :: compet_plant_nh4(bounds%begp:bounds%endp)                 ! (unitless) relative compettiveness of plants for NH4 BGC mode
    real(r8) :: compet_plant_n(bounds%begp:bounds%endp)                   ! (unitless) relative compettiveness of plants for N CN mode
    real(r8) :: compet_plant_p(bounds%begp:bounds%endp)                   ! (unitless) relative competitiveness of plant for P
    real(r8) :: compet_leach_no3                                          ! (unitless) relative competitiveness of leaching for NO3
    real(r8) :: compet_decomp_p                                           ! (unitless) relative competitiveness of immobilizer for P
    real(r8) :: compet_minsurf_p                                          ! (unitless) relative competitiveness of mineral surface for P
    real(r8) :: compet_leach_p                                            ! (unitless) relative competitiveness of leaching for P

    !
    integer :: c,p,l,j                                               !indices
    integer :: fp                                                    !lake filter pft index
    integer :: fc                                                    !lake filter column index
    real(r8):: mr                                                    !maintenance respiration (gC/m2/s)
    real(r8):: f1,f2,f3,f4,g1,g2                                     !allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw                                    !C:N ratios for leaf, fine root, and wood

    real(r8):: curmr, curmr_ratio                                    !xsmrpool temporary variables
!    real(r8):: sum_ndemand_vr(bounds%begc:bounds%endc, 1:nlevdecomp) !total column N demand (gN/m3/s) at a given level
!    real(r8):: sminn_tot(bounds%begc:bounds%endc)
    real(r8):: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)

    real(r8) f5                                                      !grain allocation parameter
    real(r8) cng                                                     !C:N ratio for grain (= cnlw for now; slevis)
    real(r8) fleaf                                                   !fraction allocated to leaf
    real(r8) t1                                                      !temporary variable
    integer :: yr, mon, day, sec

    !! Local P variables
    real(r8):: cpl,cpfr,cplw,cpdw,cpg                                    !C:N ratios for leaf, fine root, and wood
    real(r8):: sum_pdemand_vr(bounds%begc:bounds%endc, 1:nlevdecomp)     !total column P demand (gN/m3/s) at a given level
    real(r8):: puptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8):: solutionp_tot(bounds%begc:bounds%endc)
    integer :: plimit(bounds%begc:bounds%endc,0:nlevdecomp)              !flag for P limitation
    real(r8):: residual_sminp_vr(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8):: residual_sminp(bounds%begc:bounds%endc)
    real(r8):: residual_plant_pdemand(bounds%begc:bounds%endc)
    real(r8):: sum_ndemand_scaled(bounds%begc:bounds%endc, 1:nlevdecomp) !total column N demand (gN/m3/s) at a given level

    real(r8), pointer :: desorb_to_solutionp_vr       (:,:)
    real(r8), pointer :: primp_to_labilep_vr_col      (:,:)
    real(r8), pointer :: biochem_pmin_vr_col          (:,:)
    real(r8), pointer :: secondp_to_labilep_vr_col    (:,:)
    real(r8), pointer :: labilep_to_secondp_vr_col    (:,:)
    real(r8), pointer :: adsorb_to_labilep_vr         (:,:)
    real(r8), pointer :: plant_pdemand_vr_patch       (:,:)
    real(r8), pointer :: plant_n_uptake_flux          (:)

    real(r8), pointer :: labilep_vr                   (:,:)
    real(r8), pointer :: secondp_vr                   (:,:)
    real(r8), pointer :: actual_leafcp                (:)
    real(r8), pointer :: actual_frootcp               (:)
    real(r8), pointer :: actual_livewdcp              (:)
    real(r8), pointer :: actual_deadwdcp              (:)
    real(r8), pointer :: leafp                        (:)
    real(r8), pointer :: plant_p_uptake_flux          (:)

    real(r8), pointer :: col_plant_ndemand_vr         (:,:)
    real(r8), pointer :: col_plant_nh4demand_vr       (:,:)
    real(r8), pointer :: col_plant_no3demand_vr       (:,:)
    real(r8), pointer :: col_plant_pdemand_vr         (:,:)
    real(r8), pointer :: plant_nh4demand_vr_patch     (:,:)
    real(r8), pointer :: plant_no3demand_vr_patch     (:,:)
    real(r8), pointer :: plant_ndemand_vr_patch       (:,:)
    real(r8), pointer :: actual_immob_no3             (:)
    real(r8), pointer :: actual_immob_nh4             (:)
    real(r8), pointer :: benefit_pgpp_pleafc          (:)

    !-----------------------------------------------------------------------

    associate(                                                                                   &
         ivt                          => pft%itype                                             , & ! Input:  [integer  (:) ]  pft vegetation type

         woody                        => ecophyscon%woody                                      , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf                   => ecophyscon%froot_leaf                                 , & ! Input:  [real(r8) (:)   ]  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem                   => ecophyscon%croot_stem                                 , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf                    => ecophyscon%stem_leaf                                  , & ! Input:  [real(r8) (:)   ]  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd                      => ecophyscon%flivewd                                    , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn                       => ecophyscon%leafcn                                     , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
         frootcn                      => ecophyscon%frootcn                                    , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)
         livewdcn                     => ecophyscon%livewdcn                                   , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn                     => ecophyscon%deadwdcn                                   , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
         fcur2                        => ecophyscon%fcur                                       , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
         graincn                      => ecophyscon%graincn                                    , & ! Input:  [real(r8) (:)   ]  grain C:N (gC/gN)
         fleafcn                      => ecophyscon%fleafcn                                    , & ! Input:  [real(r8) (:)   ]  leaf c:n during organ fill
         ffrootcn                     => ecophyscon%ffrootcn                                   , & ! Input:  [real(r8) (:)   ]  froot c:n during organ fill
         fstemcn                      => ecophyscon%fstemcn                                    , & ! Input:  [real(r8) (:)   ]  stem c:n during organ fill

         psnsun                       => photosyns_vars%psnsun_patch                           , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         psnsha                       => photosyns_vars%psnsha_patch                           , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun                   => photosyns_vars%c13_psnsun_patch                       , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsha                   => photosyns_vars%c13_psnsha_patch                       , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsun                   => photosyns_vars%c14_psnsun_patch                       , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsha                   => photosyns_vars%c14_psnsha_patch                       , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

         laisun                       => canopystate_vars%laisun_patch                         , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                       => canopystate_vars%laisha_patch                         , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index

         hui                          => crop_vars%gddplant_patch                              , & ! Input:  [real(r8) (:)   ]  =gdd since planting (gddplant)
         leafout                      => crop_vars%gddtsoi_patch                               , & ! Input:  [real(r8) (:)   ]  =gdd from top soil layer temperature

         xsmrpool                     => carbonstate_vars%xsmrpool_patch                       , & ! Input:  [real(r8) (:)   ]  (gC/m2) temporary photosynthate C pool
         leafc                        => carbonstate_vars%leafc_patch                          , & ! Input:  [real(r8) (:)   ]
         frootc                       => carbonstate_vars%frootc_patch                         , & ! Input:  [real(r8) (:)   ]
         livestemc                    => carbonstate_vars%livestemc_patch                      , & ! Input:  [real(r8) (:)   ]
         plant_ndemand_col            => nitrogenflux_vars%plant_ndemand_col                 , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_col            => phosphorusflux_vars%plant_pdemand_col               , & ! Output:  [real(r8) (:,:) ]
         plant_ndemand_vr_col         => nitrogenflux_vars%plant_ndemand_vr_col              , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_vr_col         => phosphorusflux_vars%plant_pdemand_vr_col            , & ! Output:  [real(r8) (:,:) ]

         gddmaturity                  => cnstate_vars%gddmaturity_patch                        , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest
         huileaf                      => cnstate_vars%huileaf_patch                            , & ! Input:  [real(r8) (:)   ]  heat unit index needed from planting to leaf emergence
         huigrain                     => cnstate_vars%huigrain_patch                           , & ! Input:  [real(r8) (:)   ]  same to reach vegetative maturity
         croplive                     => cnstate_vars%croplive_patch                           , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
         peaklai                      => cnstate_vars%peaklai_patch                            , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         !lgsf                        => cnstate_vars%lgsf_patch                               , & ! Input:  [real(r8) (:)   ]  long growing season factor [0-1]
         aleafi                       => cnstate_vars%aleafi_patch                             , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         astemi                       => cnstate_vars%astemi_patch                             , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         aleaf                        => cnstate_vars%aleaf_patch                              , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnstate_vars%astem_patch                              , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         fpg                          => cnstate_vars%fpg_col                                  , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi                          => cnstate_vars%fpi_col                                  , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => cnstate_vars%fpi_vr_col                               , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

         !!! add phosphorus
         leafcp                       => ecophyscon%leafcp                                     , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
         frootcp                      => ecophyscon%frootcp                                    , & ! Input:  [real(r8) (:)   ]  fine root C:P (gC/gP)
         livewdcp                     => ecophyscon%livewdcp                                   , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
         deadwdcp                     => ecophyscon%deadwdcp                                   , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:P (gC/gP)
         graincp                      => ecophyscon%graincp                                    , & ! Input:  [real(r8) (:)   ]  grain C:P (gC/gP)
         fpg_p                        => cnstate_vars%fpg_p_col                                , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi_p                        => cnstate_vars%fpi_p_col                                , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                             , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

         nfixation_prof               => cnstate_vars%nfixation_prof_col                       , & ! Output: [real(r8) (:,:) ]
         grain_flag                   => cnstate_vars%grain_flag_patch                         , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not
         c_allometry                  => cnstate_vars%c_allometry_patch                        , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry                  => cnstate_vars%n_allometry_patch                        , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         tempsum_potential_gpp        => cnstate_vars%tempsum_potential_gpp_patch              , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP
         tempmax_retransn             => cnstate_vars%tempmax_retransn_patch                   , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)
         annsum_potential_gpp         => cnstate_vars%annsum_potential_gpp_patch               , & ! Output: [real(r8) (:)   ]  annual sum of potential GPP
         annmax_retransn              => cnstate_vars%annmax_retransn_patch                    , & ! Output: [real(r8) (:)   ]  annual max of retranslocated N pool
         downreg                      => cnstate_vars%downreg_patch                            , & ! Output: [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)

         leaf_mr                      => carbonflux_vars%leaf_mr_patch                         , & ! Input:  [real(r8) (:)   ]
         froot_mr                     => carbonflux_vars%froot_mr_patch                        , & ! Input:  [real(r8) (:)   ]
         livestem_mr                  => carbonflux_vars%livestem_mr_patch                     , & ! Input:  [real(r8) (:)   ]
         livecroot_mr                 => carbonflux_vars%livecroot_mr_patch                    , & ! Input:  [real(r8) (:)   ]
         grain_mr                     => carbonflux_vars%grain_mr_patch                        , & ! Input:  [real(r8) (:)   ]
         annsum_npp                   => carbonflux_vars%annsum_npp_patch                      , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         gpp                          => carbonflux_vars%gpp_before_downreg_patch              , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                       => carbonflux_vars%availc_patch                          , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         xsmrpool_recover             => carbonflux_vars%xsmrpool_recover_patch                , & ! Output: [real(r8) (:)   ]  C flux assigned to recovery of negative cpool (gC/m2/s)
         excess_cflux                 => carbonflux_vars%excess_cflux_patch                    , & ! Output: [real(r8) (:)   ]  C flux not allocated due to downregulation (gC/m2/s)
         plant_calloc                 => carbonflux_vars%plant_calloc_patch                    , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)
         psnsun_to_cpool              => carbonflux_vars%psnsun_to_cpool_patch                 , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool            => carbonflux_vars%psnshade_to_cpool_patch               , & ! Output: [real(r8) (:)   ]

         leaf_curmr                   => carbonflux_vars%leaf_curmr_patch                      , &
         froot_curmr                  => carbonflux_vars%froot_curmr_patch                     , & ! Output: [real(r8) (:)   ]
         livestem_curmr               => carbonflux_vars%livestem_curmr_patch                  , & ! Output: [real(r8) (:)   ]
         livecroot_curmr              => carbonflux_vars%livecroot_curmr_patch                 , & ! Output: [real(r8) (:)   ]
         grain_curmr                  => carbonflux_vars%grain_curmr_patch                     , & ! Output: [real(r8) (:)   ]
         leaf_xsmr                    => carbonflux_vars%leaf_xsmr_patch                       , & ! Output: [real(r8) (:)   ]
         froot_xsmr                   => carbonflux_vars%froot_xsmr_patch                      , & ! Output: [real(r8) (:)   ]
         livestem_xsmr                => carbonflux_vars%livestem_xsmr_patch                   , & ! Output: [real(r8) (:)   ]
         livecroot_xsmr               => carbonflux_vars%livecroot_xsmr_patch                  , & ! Output: [real(r8) (:)   ]
         grain_xsmr                   => carbonflux_vars%grain_xsmr_patch                      , & ! Output: [real(r8) (:)   ]
         cpool_to_xsmrpool            => carbonflux_vars%cpool_to_xsmrpool_patch               , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc               => carbonflux_vars%cpool_to_leafc_patch                  , & ! Output: [real(r8) (:)   ]
         cpool_to_leafc_storage       => carbonflux_vars%cpool_to_leafc_storage_patch          , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc              => carbonflux_vars%cpool_to_frootc_patch                 , & ! Output: [real(r8) (:)   ]
         cpool_to_frootc_storage      => carbonflux_vars%cpool_to_frootc_storage_patch         , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc           => carbonflux_vars%cpool_to_livestemc_patch              , & ! Output: [real(r8) (:)   ]
         cpool_to_livestemc_storage   => carbonflux_vars%cpool_to_livestemc_storage_patch      , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc           => carbonflux_vars%cpool_to_deadstemc_patch              , & ! Output: [real(r8) (:)   ]
         cpool_to_deadstemc_storage   => carbonflux_vars%cpool_to_deadstemc_storage_patch      , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc          => carbonflux_vars%cpool_to_livecrootc_patch             , & ! Output: [real(r8) (:)   ]
         cpool_to_livecrootc_storage  => carbonflux_vars%cpool_to_livecrootc_storage_patch     , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc          => carbonflux_vars%cpool_to_deadcrootc_patch             , & ! Output: [real(r8) (:)   ]
         cpool_to_deadcrootc_storage  => carbonflux_vars%cpool_to_deadcrootc_storage_patch     , & ! Output: [real(r8) (:)   ]
         cpool_to_gresp_storage       => carbonflux_vars%cpool_to_gresp_storage_patch          , & ! Output: [real(r8) (:)   ]  allocation to growth respiration storage (gC/m2/s)
         cpool_to_grainc              => carbonflux_vars%cpool_to_grainc_patch                 , & ! Output: [real(r8) (:)   ]  allocation to grain C (gC/m2/s)
         cpool_to_grainc_storage      => carbonflux_vars%cpool_to_grainc_storage_patch         , & ! Output: [real(r8) (:)   ]  allocation to grain C storage (gC/m2/s)

         sminn_vr                     => nitrogenstate_vars%sminn_vr_col                       , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
         retransn                     => nitrogenstate_vars%retransn_patch                     , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N
         smin_nh4_vr                  => nitrogenstate_vars%smin_nh4_vr_col                    , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NH4
         smin_no3_vr                  => nitrogenstate_vars%smin_no3_vr_col                    , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3

         plant_ndemand                => nitrogenflux_vars%plant_ndemand_patch                 , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_nalloc                 => nitrogenflux_vars%plant_nalloc_patch                  , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)
         avail_retransn               => nitrogenflux_vars%avail_retransn_patch                , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         npool_to_grainn              => nitrogenflux_vars%npool_to_grainn_patch               , & ! Output: [real(r8) (:)   ]  allocation to grain N (gN/m2/s)
         npool_to_grainn_storage      => nitrogenflux_vars%npool_to_grainn_storage_patch       , & ! Output: [real(r8) (:)   ]  allocation to grain N storage (gN/m2/s)
         retransn_to_npool            => nitrogenflux_vars%retransn_to_npool_patch             , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         sminn_to_npool               => nitrogenflux_vars%sminn_to_npool_patch                , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
         npool_to_leafn               => nitrogenflux_vars%npool_to_leafn_patch                , & ! Output: [real(r8) (:)   ]  allocation to leaf N (gN/m2/s)
         npool_to_leafn_storage       => nitrogenflux_vars%npool_to_leafn_storage_patch        , & ! Output: [real(r8) (:)   ]  allocation to leaf N storage (gN/m2/s)
         npool_to_frootn              => nitrogenflux_vars%npool_to_frootn_patch               , & ! Output: [real(r8) (:)   ]  allocation to fine root N (gN/m2/s)
         npool_to_frootn_storage      => nitrogenflux_vars%npool_to_frootn_storage_patch       , & ! Output: [real(r8) (:)   ]  allocation to fine root N storage (gN/m2/s)
         npool_to_livestemn           => nitrogenflux_vars%npool_to_livestemn_patch            , & ! Output: [real(r8) (:)   ]
         npool_to_livestemn_storage   => nitrogenflux_vars%npool_to_livestemn_storage_patch    , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn           => nitrogenflux_vars%npool_to_deadstemn_patch            , & ! Output: [real(r8) (:)   ]
         npool_to_deadstemn_storage   => nitrogenflux_vars%npool_to_deadstemn_storage_patch    , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn          => nitrogenflux_vars%npool_to_livecrootn_patch           , & ! Output: [real(r8) (:)   ]
         npool_to_livecrootn_storage  => nitrogenflux_vars%npool_to_livecrootn_storage_patch   , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn          => nitrogenflux_vars%npool_to_deadcrootn_patch           , & ! Output: [real(r8) (:)   ]
         npool_to_deadcrootn_storage  => nitrogenflux_vars%npool_to_deadcrootn_storage_patch   , & ! Output: [real(r8) (:)   ]
         leafn_to_retransn            => nitrogenflux_vars%leafn_to_retransn_patch             , & ! Output: [real(r8) (:)   ]
         frootn_to_retransn           => nitrogenflux_vars%frootn_to_retransn_patch            , & ! Output: [real(r8) (:)   ]
         livestemn_to_retransn        => nitrogenflux_vars%livestemn_to_retransn_patch         , & ! Output: [real(r8) (:)   ]
         potential_immob              => nitrogenflux_vars%potential_immob_col                 , & ! Output: [real(r8) (:)   ]
         actual_immob                 => nitrogenflux_vars%actual_immob_col                    , & ! Output: [real(r8) (:)   ]
         sminn_to_denit_excess_vr     => nitrogenflux_vars%sminn_to_denit_excess_vr_col        , & ! Output: [real(r8) (:,:) ]
         pot_f_nit_vr                 => nitrogenflux_vars%pot_f_nit_vr_col                    , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => nitrogenflux_vars%pot_f_denit_vr_col                  , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => nitrogenflux_vars%f_nit_vr_col                        , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => nitrogenflux_vars%f_denit_vr_col                      , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         actual_immob_no3_vr          => nitrogenflux_vars%actual_immob_no3_vr_col             , & ! Output: [real(r8) (:,:) ]
         actual_immob_nh4_vr          => nitrogenflux_vars%actual_immob_nh4_vr_col             , & ! Output: [real(r8) (:,:) ]
         n2_n2o_ratio_denit_vr        => nitrogenflux_vars%n2_n2o_ratio_denit_vr_col           , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => nitrogenflux_vars%f_n2o_denit_vr_col                  , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => nitrogenflux_vars%f_n2o_nit_vr_col                    , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr_col          , & ! Output: [real(r8) (:,:) ]
         potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col              , & ! Output: [real(r8) (:,:) ]
         actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col                 , & ! Output: [real(r8) (:,:) ]

         !!! add phosphorus variables  - X. YANG
         sminp_vr                     => phosphorusstate_vars%sminp_vr_col                     , & ! Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral P
         solutionp_vr                 => phosphorusstate_vars%solutionp_vr_col                 , & ! Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral P
         retransp                     => phosphorusstate_vars%retransp_patch                   , & ! Input:  [real(r8) (:)   ]  (gP/m2) plant pool of retranslocated P

         plant_pdemand                => phosphorusflux_vars%plant_pdemand_patch               , & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)
         plant_palloc                 => phosphorusflux_vars%plant_palloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated P flux (gP/m2/s)
         avail_retransp               => phosphorusflux_vars%avail_retransp_patch              , & ! Output: [real(r8) (:)   ]  P flux available from retranslocation pool (gP/m2/s)
         ppool_to_grainp              => phosphorusflux_vars%ppool_to_grainp_patch             , & ! Output: [real(r8) (:)   ]  allocation to grain P (gP/m2/s)
         ppool_to_grainp_storage      => phosphorusflux_vars%ppool_to_grainp_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to grain P storage (gP/m2/s)
         retransp_to_ppool            => phosphorusflux_vars%retransp_to_ppool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated P (gP/m2/s)
         sminp_to_ppool               => phosphorusflux_vars%sminp_to_ppool_patch              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral P uptake (gP/m2/s)
         ppool_to_leafp               => phosphorusflux_vars%ppool_to_leafp_patch              , & ! Output: [real(r8) (:)   ]  allocation to leaf P (gP/m2/s)
         ppool_to_leafp_storage       => phosphorusflux_vars%ppool_to_leafp_storage_patch      , & ! Output: [real(r8) (:)   ]  allocation to leaf P storage (gP/m2/s)
         ppool_to_frootp              => phosphorusflux_vars%ppool_to_frootp_patch             , & ! Output: [real(r8) (:)   ]  allocation to fine root P (gP/m2/s)
         ppool_to_frootp_storage      => phosphorusflux_vars%ppool_to_frootp_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to fine root P storage (gP/m2/s)
         ppool_to_livestemp           => phosphorusflux_vars%ppool_to_livestemp_patch          , & ! Output: [real(r8) (:)   ]
         ppool_to_livestemp_storage   => phosphorusflux_vars%ppool_to_livestemp_storage_patch  , & ! Output: [real(r8) (:)   ]
         ppool_to_deadstemp           => phosphorusflux_vars%ppool_to_deadstemp_patch          , & ! Output: [real(r8) (:)   ]
         ppool_to_deadstemp_storage   => phosphorusflux_vars%ppool_to_deadstemp_storage_patch  , & ! Output: [real(r8) (:)   ]
         ppool_to_livecrootp          => phosphorusflux_vars%ppool_to_livecrootp_patch         , & ! Output: [real(r8) (:)   ]
         ppool_to_livecrootp_storage  => phosphorusflux_vars%ppool_to_livecrootp_storage_patch , & ! Output: [real(r8) (:)   ]
         ppool_to_deadcrootp          => phosphorusflux_vars%ppool_to_deadcrootp_patch         , & ! Output: [real(r8) (:)   ]
         ppool_to_deadcrootp_storage  => phosphorusflux_vars%ppool_to_deadcrootp_storage_patch , & ! Output: [real(r8) (:)   ]
         leafp_to_retransp            => phosphorusflux_vars%leafp_to_retransp_patch           , & ! Output: [real(r8) (:)   ]
         frootp_to_retransp           => phosphorusflux_vars%frootp_to_retransp_patch          , & ! Output: [real(r8) (:)   ]
         livestemp_to_retransp        => phosphorusflux_vars%livestemp_to_retransp_patch       , & ! Output: [real(r8) (:)   ]
         potential_immob_p            => phosphorusflux_vars%potential_immob_p_col             , & ! Output: [real(r8) (:)   ]
         actual_immob_p               => phosphorusflux_vars%actual_immob_p_col                , & ! Output: [real(r8) (:)   ]
         supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr_col        , & ! Output: [real(r8) (:,:) ]
         potential_immob_p_vr         => phosphorusflux_vars%potential_immob_p_vr_col          , & ! Output: [real(r8) (:,:) ]
         actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr_col             , & ! Output: [real(r8) (:,:) ]
         p_allometry                  => cnstate_vars%p_allometry_patch                        , & ! Output: [real(r8) (:)   ]  P allocation index (DIM)
         tempmax_retransp             => cnstate_vars%tempmax_retransp_patch                   , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated P pool (gP/m2)
         annmax_retransp              => cnstate_vars%annmax_retransp_patch                    , & ! Output: [real(r8) (:)   ]  annual max of retranslocated P pool

         c13cf                        => c13_carbonflux_vars                                   , &
         c14cf                        => c14_carbonflux_vars                                   , &

         froot_prof                   => cnstate_vars%froot_prof_patch                         , & ! fine root vertical profile Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530
         fpg_nh4_vr                   => cnstate_vars%fpg_nh4_vr_col                           , &
         fpg_no3_vr                   => cnstate_vars%fpg_no3_vr_col                           , &
         fpg_vr                       => cnstate_vars%fpg_vr_col                               , &
         fpg_p_vr                     => cnstate_vars%fpg_p_vr_col                             , &
         cn_scalar                    => cnstate_vars%cn_scalar                                , &
         cp_scalar                    => cnstate_vars%cp_scalar                                , &
         isoilorder                   => cnstate_vars%isoilorder                               , &

         t_scalar                     => carbonflux_vars%t_scalar_col                          , &
         w_scalar                     => carbonflux_vars%w_scalar_col                          , &
         pnup_pfrootc                 => nitrogenstate_vars%pnup_pfrootc_patch                 , &
         leafn                        => nitrogenstate_vars%leafn_patch                          &
         )
      secondp_vr                   => phosphorusstate_vars%secondp_vr_col
      actual_leafcp                => phosphorusstate_vars%actual_leafcp
      actual_frootcp               => phosphorusstate_vars%actual_frootcp
      actual_livewdcp              => phosphorusstate_vars%actual_livewdcp
      actual_deadwdcp              => phosphorusstate_vars%actual_deadwdcp
      leafp                        => phosphorusstate_vars%leafp_patch
      col_plant_ndemand_vr         => nitrogenflux_vars%col_plant_ndemand_vr
      col_plant_nh4demand_vr       => nitrogenflux_vars%col_plant_nh4demand_vr
      col_plant_no3demand_vr       => nitrogenflux_vars%col_plant_no3demand_vr
      col_plant_pdemand_vr         => nitrogenflux_vars%col_plant_pdemand_vr
      plant_nh4demand_vr_patch     => nitrogenflux_vars%plant_nh4demand_vr_patch
      plant_no3demand_vr_patch     => nitrogenflux_vars%plant_no3demand_vr_patch
      plant_ndemand_vr_patch       => nitrogenflux_vars%plant_ndemand_vr_patch
      plant_pdemand_vr_patch       => phosphorusflux_vars%plant_pdemand_vr_patch
      actual_immob_no3             => nitrogenflux_vars%actual_immob_no3_col
      actual_immob_nh4             => nitrogenflux_vars%actual_immob_nh4_col
      adsorb_to_labilep_vr         => phosphorusflux_vars%adsorb_to_labilep_vr
      desorb_to_solutionp_vr       => phosphorusflux_vars%desorb_to_solutionp_vr
      primp_to_labilep_vr_col      => phosphorusflux_vars%primp_to_labilep_vr_col
      biochem_pmin_vr_col          => phosphorusflux_vars%biochem_pmin_vr_col
      secondp_to_labilep_vr_col    => phosphorusflux_vars%secondp_to_labilep_vr_col
      labilep_to_secondp_vr_col    => phosphorusflux_vars%labilep_to_secondp_vr_col
      labilep_vr                   => phosphorusstate_vars%labilep_vr_col
      benefit_pgpp_pleafc          => nitrogenstate_vars%benefit_pgpp_pleafc_patch

      ! for debug
      plant_n_uptake_flux          => nitrogenflux_vars%plant_n_uptake_flux
      plant_p_uptake_flux          => phosphorusflux_vars%plant_p_uptake_flux

      ! set time steps
      dt = real( get_step_size(), r8 )

      call get_curr_date(yr, mon, day, sec)
      if (spinup_state == 1 .and. yr .gt. nyears_ad_carbon_only) then
        call cnallocate_carbon_only_set(.false.)
      end if

     ! loop over patches to assess the total plant N demand and P demand
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

         ! carbon return of leaf C investment
         benefit_pgpp_pleafc(p) = max(gpp(p)  / max(leafc(p) ,1e-20_r8), 0.0_r8)

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

         cpl = leafcp(ivt(p))
         cpfr = frootcp(ivt(p))
         cplw = livewdcp(ivt(p))
         cpdw = deadwdcp(ivt(p))


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
                     aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -   &
                          (arooti(ivt(p)) - arootf(ivt(p))) *  &
                          min(1._r8, hui(p)/gddmaturity(p))))
                     astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
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

                  if (ivt(p) /= nsoybean .or. astem(p) == astemf(ivt(p)) .or. peaklai(p) == 1._r8) then
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
         ! determine P requirements   -X. YANG

         if (woody(ivt(p)) == 1.0_r8) then
            c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
            n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
            p_allometry(p) = 1._r8/cpl + f1/cpfr + (f3*f4*(1._r8+f2))/cplw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cpdw

         else if (ivt(p) >= npcropmin) then ! skip generic crops
            cng = graincn(ivt(p))
            cpg = graincp(ivt(p))
            c_allometry(p) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
            n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cndw
            p_allometry(p) = 1._r8/cpl + f1/cpfr + f5/cpg + (f3*f4*(1._r8+f2))/cplw + &
                 (f3*(1._r8-f4)*(1._r8+f2))/cpdw

         else
            c_allometry(p) = 1._r8+g1+f1+f1*g1
            n_allometry(p) = 1._r8/cnl + f1/cnfr
            p_allometry(p) = 1._r8/cpl + f1/cpfr
         end if
         plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))
         plant_pdemand(p) = availc(p)*(p_allometry(p)/c_allometry(p))

         ! retranslocated N deployment depends on seasonal cycle of potential GPP
         ! (requires one year run to accumulate demand)

         tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

         ! Adding the following line to carry max retransn info to CN Annual Update
         tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))
         tempmax_retransp(p) = max(tempmax_retransp(p),retransp(p))   !! phosphorus

         ! Beth's code: crops pull from retransn pool only during grain fill;
         !              retransn pool has N from leaves, stems, and roots for
         !              retranslocation

         if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
            avail_retransn(p) = plant_ndemand(p)
            avail_retransp(p) = plant_pdemand(p)
         else if (ivt(p) < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
            avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
            avail_retransp(p) = (annmax_retransp(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
         else
            avail_retransn(p) = 0.0_r8
            avail_retransp(p) = 0.0_r8
         end if

         ! make sure available retrans N doesn't exceed storage
         avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)
         avail_retransp(p) = min(avail_retransp(p), retransp(p)/dt)    !! phosphorus

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

         if (plant_pdemand(p) > avail_retransp(p)) then
            retransp_to_ppool(p) = avail_retransp(p)
         else
            retransp_to_ppool(p) = plant_pdemand(p)
         end if
         plant_pdemand(p) = plant_pdemand(p) - retransp_to_ppool(p)


      end do ! end pft loop

      ! now use the p2c routine to get the column-averaged plant_ndemand
      call p2c(bounds, num_soilc, filter_soilc, &
           plant_ndemand(bounds%begp:bounds%endp), &
           col_plant_ndemand(bounds%begc:bounds%endc))

      !!! add phosphorus
      call p2c(bounds, num_soilc, filter_soilc, &
           plant_pdemand(bounds%begp:bounds%endp), &
           col_plant_pdemand(bounds%begc:bounds%endc))

      !!! Starting resolving N limitation
      !! new subroutines to calculate nuptake_prof & puptake_prof

      !! flux_type%var = local var, used in CNAllocation2
      do fc=1, num_soilc
            c = filter_soilc(fc)
            plant_ndemand_col(c) = col_plant_ndemand(c)
            plant_pdemand_col(c) = col_plant_pdemand(c)
      end do

      ! pflotran will need an input from CN: modified 'sum_ndemand_vr' ('potential_immob' excluded).
      if (use_bgc_interface.and.use_pflotran .and. pf_cmode) then
            do j = 1, nlevdecomp
               do fc=1, num_soilc
                  c = filter_soilc(fc)
                  plant_ndemand_vr_col(c,j) = plant_ndemand_col(c) * nuptake_prof(c,j)
                  plant_pdemand_vr_col(c,j) = plant_pdemand_col(c) * puptake_prof(c,j)
               end do
            end do
      endif

    end associate

 end subroutine CNAllocation1_PlantNPDemand
!------------------------------------------------------------------------------
  subroutine calc_plantN_kineticpar(bounds, num_soilc, filter_soilc         , &
                            num_soilp, filter_soilp                         , &
                            cnstate_vars                                    , &
                            carbonstate_vars                                , &
                            nitrogenstate_vars                              , &
                            phosphorusstate_vars                            , &
                            carbonflux_vars                                 , &
                            PlantMicKinetics_vars                             )
  !
  !DESCRIPTION
  !compute kinetic parameters for nutrient competition
  use clm_varpar       , only:  nlevdecomp !!nlevsoi,
  use pftvarcon        , only:  noveg
  implicit none
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilc
  integer, intent(in) :: filter_soilc(:)
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(cnstate_type), intent(in) :: cnstate_vars
  type(carbonstate_type), intent(in) :: carbonstate_vars
  type(nitrogenstate_type), intent(in) :: nitrogenstate_vars
  type(phosphorusstate_type), intent(in):: phosphorusstate_vars
  type(carbonflux_type), intent(in) :: carbonflux_vars
  type(PlantMicKinetics_type), intent(inout) :: PlantMicKinetics_vars

  real(r8) :: leaf_totc
  real(r8) :: leaf_totn
  real(r8) :: leaf_totp
  integer :: p, c, fc, j

  real(r8), parameter :: cn_stoich_var=0.2    ! variability of CN ratio
  real(r8), parameter :: cp_stoich_var=0.5    ! variability of CP ratio

  associate(                                                                            &
     ivt                          => pft%itype                                        , & ! Input:  [integer  (:) ]  pft vegetation type

     frootcn                      => ecophyscon%frootcn                               , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)
     leafcn                       => ecophyscon%leafcn                                , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
     leafcp                       => ecophyscon%leafcp                                , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)

     vmax_plant_nh4               => ecophyscon%vmax_plant_nh4                        , &
     vmax_plant_no3               => ecophyscon%vmax_plant_no3                        , &
     vmax_plant_p                 => ecophyscon%vmax_plant_p                          , &
     decompmicc_patch_vr          => ecophyscon%decompmicc_patch_vr                   , &
     vmax_minsurf_p_vr            => ecophyscon%vmax_minsurf_p_vr                      , &


!  the following parameter are defined uniformly for all grids
     km_decomp_nh4                => ecophyscon%km_decomp_nh4                         , &
     km_decomp_no3                => ecophyscon%km_decomp_no3                         , &
     km_decomp_p                  => ecophyscon%km_decomp_p                           , &
     km_nit                       => ecophyscon%km_nit                                , &
     km_den                       => ecophyscon%km_den                                , &

     km_plant_nh4                 => ecophyscon%km_plant_nh4                          , &
     km_plant_no3                 => ecophyscon%km_plant_no3                          , &
     km_plant_p                   => ecophyscon%km_plant_p                            , &
! the following parameter is defined based on soil order
     km_minsurf_p_vr              => ecophyscon%km_minsurf_p_vr                       , &

     plant_nh4_vmax_vr_patch      => PlantMicKinetics_vars%plant_nh4_vmax_vr_patch    , & ! Output: [real(r8) (:,:) ] vmax for nh4 uptake
     plant_no3_vmax_vr_patch      => PlantMicKinetics_vars%plant_no3_vmax_vr_patch    , & ! Output: [real(r8) (:,:) ] vmax for nh4 uptake
     plant_p_vmax_vr_patch        => PlantMicKinetics_vars%plant_p_vmax_vr_patch      , & ! Output: [real(r8) (:,:) ] vmax for nh4 uptake
     plant_nh4_km_vr_patch        => PlantMicKinetics_vars%plant_nh4_km_vr_patch      , & !
     plant_no3_km_vr_patch        => PlantMicKinetics_vars%plant_no3_km_vr_patch      , & !
     plant_p_km_vr_patch          => PlantMicKinetics_vars%plant_p_km_vr_patch        , & !
     vmax_minsurf_p_vr_col        => PlantMicKinetics_vars%vmax_minsurf_p_vr_col      , & !
     km_minsurf_p_vr_col          => PlantMicKinetics_vars%km_minsurf_p_vr_col        , &
     km_den_no3_vr_col            => PlantMicKinetics_vars%km_den_no3_vr_col          , &
     km_nit_nh4_vr_col            => PlantMicKinetics_vars%km_nit_nh4_vr_col          , &
     km_decomp_p_vr_col           => PlantMicKinetics_vars%km_decomp_p_vr_col         , &
     km_decomp_nh4_vr_col         => PlantMicKinetics_vars%km_decomp_nh4_vr_col       , &
     km_decomp_no3_vr_col         => PlantMicKinetics_vars%km_decomp_no3_vr_col       , &
     plant_eff_ncompet_b          => PlantMicKinetics_vars%plant_eff_ncompet_b_vr_patch  , &
     plant_eff_pcompet_b          => PlantMicKinetics_vars%plant_eff_pcompet_b_vr_patch  , &
     decomp_eff_ncompet_b         => PlantMicKinetics_vars%decomp_eff_ncompet_b_vr_col, &
     decomp_eff_pcompet_b         => PlantMicKinetics_vars%decomp_eff_pcompet_b_vr_col, &
     minsurf_p_compet             => PlantMicKinetics_vars%minsurf_p_compet_vr_col, &
     isoilorder                   => cnstate_vars%isoilorder                          , &
     cn_scalar                    => cnstate_vars%cn_scalar                           , &
     cp_scalar                    => cnstate_vars%cp_scalar                           , &
     froot_prof                   => cnstate_vars%froot_prof_patch                    , & ! fine root vertical profile Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530
     frootc                       => carbonstate_vars%frootc_patch                    , & ! Input:  [real(r8) (:)   ]
     leafc                        => carbonstate_vars%leafc_patch                     , & ! Input:  [real(r8) (:)   ]
     leafc_storage                => carbonstate_vars%leafc_storage_patch             , &
     leafc_xfer                   => carbonstate_vars%leafc_xfer_patch                , &
     t_scalar                     => carbonflux_vars%t_scalar_col                     , &
     leafn                        => nitrogenstate_vars%leafn_patch                   , &
     leafn_storage                => nitrogenstate_vars%leafn_storage_patch           , &
     leafn_xfer                   => nitrogenstate_vars%leafn_xfer_patch              , &
     leafp                        => phosphorusstate_vars%leafp_patch                 , &
     leafp_storage                => phosphorusstate_vars%leafp_storage_patch         , &
     leafp_xfer                   => phosphorusstate_vars%leafp_xfer_patch              &
  )

  do j = 1, nlevdecomp
    do fc=1,num_soilc
      c = filter_soilc(fc)
      vmax_minsurf_p_vr_col(c,j) = vmax_minsurf_p_vr(isoilorder(c),j)
      km_minsurf_p_vr_col(c,j) = km_minsurf_p_vr(isoilorder(c),j)

      !the following is temporary set using alm, in the future, it will be
      !set through the betr_alm interface
      km_den_no3_vr_col(c,j) = km_den
      km_nit_nh4_vr_col(c,j) = km_nit
      km_decomp_p_vr_col(c,j) = km_decomp_p
      km_decomp_no3_vr_col(c,j) = km_decomp_no3
      km_decomp_nh4_vr_col(c,j) = km_decomp_nh4

      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
        ! scaling factor based on  CN ratio flexibility
          leaf_totc=leafc(p) + leafc_storage(p) + leafc_xfer(p)
          leaf_totn=leafn(p) + leafn_storage(p) + leafn_xfer(p)
          leaf_totp=leafp(p) + leafp_storage(p) + leafp_xfer(p)
          cn_scalar(p) = min(max((leaf_totc/max(leaf_totn, 1e-20_r8) - leafcn(ivt(p))*(1- cn_stoich_var)) / &
                (leafcn(ivt(p)) - leafcn(ivt(p))*(1- cn_stoich_var)),0.0_r8),1.0_r8)

          cp_scalar(p) = min(max((leaf_totc/max(leaf_totp, 1e-20_r8) - leafcp(ivt(p))*(1- cp_stoich_var)) / &
           (leafcp(ivt(p)) - leafcp(ivt(p))*(1- cp_stoich_var)),0.0_r8),1.0_r8)

          plant_nh4_vmax_vr_patch(p,j) = vmax_plant_nh4(ivt(p))* frootc(p) * froot_prof(p,j) * &
                             cn_scalar(p) * t_scalar(c,j)
          plant_no3_vmax_vr_patch(p,j) = vmax_plant_no3(ivt(p)) * frootc(p) * froot_prof(p,j) * &
                             cn_scalar(p) * t_scalar(c,j)
          plant_p_vmax_vr_patch(p,j) = vmax_plant_p(ivt(p)) * frootc(p) * froot_prof(p,j) * &
                             cp_scalar(p) * t_scalar(c,j)

          plant_nh4_km_vr_patch(p,j) = km_plant_nh4(ivt(p))
          plant_no3_km_vr_patch(p,j) = km_plant_no3(ivt(p))
          plant_p_km_vr_patch(p,j) = km_plant_p(ivt(p))

          plant_eff_ncompet_b(p,j) = e_plant_scalar*frootc(p)*froot_prof(p,j)
          plant_eff_pcompet_b(p,j) = e_plant_scalar*frootc(p)*froot_prof(p,j)
        else
          cn_scalar(p) = 1.0_r8
        end if
        !effective n competing decomposers
        decomp_eff_ncompet_b(c,j) = decomp_eff_ncompet_b(c,j) + &
              e_decomp_scalar*decompmicc_patch_vr(ivt(p),j)*pft%wtcol(p)
        !effective p competing decomposers
        decomp_eff_pcompet_b(c,j) = decomp_eff_pcompet_b(c,j) + &
              e_decomp_scalar*decompmicc_patch_vr(ivt(p),j)*pft%wtcol(p)
        !effective p competing mineral surfaces
        minsurf_p_compet(c,j) = 0._r8
      end do
    enddo
  end do
  end associate
 end subroutine calc_plantN_kineticpar

!!-------------------------------------------------------------------------------------------------
  subroutine CNAllocation3_PlantCNPAlloc (bounds            , &
        num_soilc, filter_soilc, num_soilp, filter_soilp    , &
        canopystate_vars                                    , &
        cnstate_vars, carbonstate_vars, carbonflux_vars     , &
        c13_carbonflux_vars, c14_carbonflux_vars            , &
        nitrogenstate_vars, nitrogenflux_vars               , &
        phosphorusstate_vars, phosphorusflux_vars)
    !! PHASE-3 of CNAllocation: start new pft loop to distribute the available N/P between the
    ! competing patches on the basis of relative demand, and allocate C/N/P to new growth and storage

    ! !USES:
    use shr_sys_mod      , only: shr_sys_flush
    use clm_varctl       , only: iulog,cnallocate_carbon_only,cnallocate_carbonnitrogen_only,&
                                 cnallocate_carbonphosphorus_only
!    use pftvarcon        , only: npcropmin, declfact, bfact, aleaff, arootf, astemf
!    use pftvarcon        , only: arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
    use pftvarcon        , only: noveg
    use pftvarcon        , only:  npcropmin, grperc, grpnow
    use clm_varpar       , only:  nlevdecomp !!nlevsoi,
    use clm_varcon       , only: nitrif_n2o_loss_frac, secspday
!    use landunit_varcon  , only: istsoil, istcrop
!    use clm_time_manager , only: get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
!    type(photosyns_type)     , intent(in)    :: photosyns_vars
!    type(crop_type)          , intent(in)    :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
!    !!  add phosphorus  -X.YANG
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    !
    integer :: c,p,j         !!l,pi,                                   !indices
    integer :: fp                                                    !lake filter pft index
    integer :: fc                                                    !lake filter column index
    real(r8):: mr                                                    !maintenance respiration (gC/m2/s)
    real(r8):: f1,f2,f3,f4,f5,g1,g2                                     !allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw                                    !C:N ratios for leaf, fine root, and wood
    real(r8):: fcur                                                  !fraction of current psn displayed as growth
    real(r8):: gresp_storage                                         !temporary variable for growth resp to storage
    real(r8):: nlc                                                   !temporary variable for total new leaf carbon allocation
    real(r8):: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8) cng                                                     !C:N ratio for grain (= cnlw for now; slevis)

    !! Local P variables
    real(r8):: rc_npool, rc, r                                               !Factors for nitrogen pool
    real(r8):: cpl,cpfr,cplw,cpdw,cpg                                    !C:N ratios for leaf, fine root, and wood
    real(r8):: puptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)

    !real(r8) :: allocation_leaf(bounds%begp : bounds%endp)              ! fraction of NPP allocated into leaf
    !real(r8) :: allocation_stem(bounds%begp : bounds%endp)              ! fraction of NPP allocated into stem
    !real(r8) :: allocation_froot(bounds%begp : bounds%endp)              ! fraction of NPP allocated into froot

    real(r8):: N_lim_factor(bounds%begp : bounds%endp)                   ! N stress factor that impact dynamic C allocation
    real(r8):: P_lim_factor(bounds%begp : bounds%endp)                   ! P stress factor that impact dynamic C allocation
    real(r8):: W_lim_factor(bounds%begp : bounds%endp)                  ! water stress factor that impact dynamic C allocation
    real(r8):: nlc_adjust_high  ! adjustment of C allocation to non-structural pools due to CNP imbalance
    real(r8):: cn_stoich_var=0.2    ! variability of CN ratio
    real(r8):: cp_stoich_var=0.5    ! variability of CP ratio
    real(r8):: curmr, curmr_ratio                                    !xsmrpool temporary variables
    real(r8), parameter :: taup = 3600._r8 !turnover of the abstract plant p storage
    real(r8), parameter :: taun = 3600._r8 !turnover of the abstract plant n storage
    !-----------------------------------------------------------------------

    associate(                                                                                 &
         ivt                          => pft%itype                                           , & ! Input:  [integer  (:) ]  pft vegetation type
!
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

         croplive                     => cnstate_vars%croplive_patch                         , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
         aleaf                        => cnstate_vars%aleaf_patch                            , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnstate_vars%astem_patch                            , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
         fpg                          => cnstate_vars%fpg_col                                , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)

         !!! add phosphorus
         leafcp                       => ecophyscon%leafcp                                   , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
         frootcp                      => ecophyscon%frootcp                                  , & ! Input:  [real(r8) (:)   ]  fine root C:P (gC/gP)
         livewdcp                     => ecophyscon%livewdcp                                 , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
         deadwdcp                     => ecophyscon%deadwdcp                                 , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:P (gC/gP)
         graincp                      => ecophyscon%graincp                                  , & ! Input:  [real(r8) (:)   ]  grain C:P (gC/gP)
         fpg_p                        => cnstate_vars%fpg_p_col                              , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)

         c_allometry                  => cnstate_vars%c_allometry_patch                      , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry                  => cnstate_vars%n_allometry_patch                      , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
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

         npool                        => nitrogenstate_vars%npool_patch                      , & ! Input:  [real(r8) (:)   ]  (gN/m3) plant N pool storage
         plant_n_buffer_patch        => nitrogenstate_vars%plant_n_buffer_patch            , & ! Inout:  [real(r8) (:)   ] gN/m2

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

         !!! add phosphorus variables  - X. YANG
         plant_pdemand                => phosphorusflux_vars%plant_pdemand_patch               , & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)
         plant_palloc                 => phosphorusflux_vars%plant_palloc_patch                , & ! Output: [real(r8) (:)   ]  total allocated P flux (gP/m2/s)
         ppool_to_grainp              => phosphorusflux_vars%ppool_to_grainp_patch             , & ! Output: [real(r8) (:)   ]  allocation to grain P (gP/m2/s)
         ppool_to_grainp_storage      => phosphorusflux_vars%ppool_to_grainp_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to grain P storage (gP/m2/s)
         retransp_to_ppool            => phosphorusflux_vars%retransp_to_ppool_patch           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated P (gP/m2/s)
         sminp_to_ppool               => phosphorusflux_vars%sminp_to_ppool_patch              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral P uptake (gP/m2/s)
         ppool_to_leafp               => phosphorusflux_vars%ppool_to_leafp_patch              , & ! Output: [real(r8) (:)   ]  allocation to leaf P (gP/m2/s)
         ppool_to_leafp_storage       => phosphorusflux_vars%ppool_to_leafp_storage_patch      , & ! Output: [real(r8) (:)   ]  allocation to leaf P storage (gP/m2/s)
         ppool_to_frootp              => phosphorusflux_vars%ppool_to_frootp_patch             , & ! Output: [real(r8) (:)   ]  allocation to fine root P (gP/m2/s)
         ppool_to_frootp_storage      => phosphorusflux_vars%ppool_to_frootp_storage_patch     , & ! Output: [real(r8) (:)   ]  allocation to fine root P storage (gP/m2/s)
         ppool_to_livestemp           => phosphorusflux_vars%ppool_to_livestemp_patch          , & ! Output: [real(r8) (:)   ]
         ppool_to_livestemp_storage   => phosphorusflux_vars%ppool_to_livestemp_storage_patch  , & ! Output: [real(r8) (:)   ]
         ppool_to_deadstemp           => phosphorusflux_vars%ppool_to_deadstemp_patch          , & ! Output: [real(r8) (:)   ]
         ppool_to_deadstemp_storage   => phosphorusflux_vars%ppool_to_deadstemp_storage_patch  , & ! Output: [real(r8) (:)   ]
         ppool_to_livecrootp          => phosphorusflux_vars%ppool_to_livecrootp_patch         , & ! Output: [real(r8) (:)   ]
         ppool_to_livecrootp_storage  => phosphorusflux_vars%ppool_to_livecrootp_storage_patch , & ! Output: [real(r8) (:)   ]
         ppool_to_deadcrootp          => phosphorusflux_vars%ppool_to_deadcrootp_patch         , & ! Output: [real(r8) (:)   ]
         ppool_to_deadcrootp_storage  => phosphorusflux_vars%ppool_to_deadcrootp_storage_patch , & ! Output: [real(r8) (:)   ]
         p_allometry                  => cnstate_vars%p_allometry_patch                        , & ! Output: [real(r8) (:)   ]  P allocation index (DIM)

         avail_retransn               => nitrogenflux_vars%avail_retransn_patch                , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         avail_retransp               => phosphorusflux_vars%avail_retransp_patch              , & ! Output: [real(r8) (:)   ]  P flux available from retranslocation pool (gP/m2/s)
         retransn                     => nitrogenstate_vars%retransn_patch                     , &
         retransp                     => phosphorusstate_vars%retransp_patch                   , &
         plant_p_buffer_patch         => phosphorusstate_vars%plant_p_buffer_patch             , & ! Inout:  [real(r8) (:)   ] gN/m2

         laisun                       => canopystate_vars%laisun_patch                         , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                       => canopystate_vars%laisha_patch                         , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index
         leafc                        => carbonstate_vars%leafc_patch                          , &
         leafn                        => nitrogenstate_vars%leafn_patch                        , &
         leafp                        => phosphorusstate_vars%leafp_patch                      , &
         supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr_col          , &
         supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr_col        , &

         ! for debug
!         plant_n_uptake_flux          => nitrogenflux_vars%plant_n_uptake_flux                 , &
!         plant_p_uptake_flux          => phosphorusflux_vars%plant_p_uptake_flux               , &
         leafc_storage                => carbonstate_vars%leafc_storage_patch                  , &
         leafc_xfer                   => carbonstate_vars%leafc_xfer_patch                     , &
         leafn_storage                => nitrogenstate_vars%leafn_storage_patch                , &
         leafn_xfer                   => nitrogenstate_vars%leafn_xfer_patch                   , &
         leafp_storage                => phosphorusstate_vars%leafp_storage_patch              , &
         leafp_xfer                   => phosphorusstate_vars%leafp_xfer_patch                 , &
         annsum_potential_gpp         => cnstate_vars%annsum_potential_gpp_patch               , &
         annmax_retransn              => cnstate_vars%annmax_retransn_patch                    , &
         grain_flag                   => cnstate_vars%grain_flag_patch                         , &
         cn_scalar                    => cnstate_vars%cn_scalar                                , &
         cp_scalar                    => cnstate_vars%cp_scalar                                , &
         annmax_retransp              => cnstate_vars%annmax_retransp_patch                    , &
         cpool_to_xsmrpool            => carbonflux_vars%cpool_to_xsmrpool_patch               , &
         w_scalar                     => carbonflux_vars%w_scalar_col                          , &
         froot_prof                   => cnstate_vars%froot_prof_patch                         , &
         leaf_mr                      => carbonflux_vars%leaf_mr_patch                         , &
         froot_mr                     => carbonflux_vars%froot_mr_patch                        , &
         livestem_mr                  => carbonflux_vars%livestem_mr_patch                     , &
         livecroot_mr                 => carbonflux_vars%livecroot_mr_patch                    , &
         grain_mr                     => carbonflux_vars%grain_mr_patch                        , &
         xsmrpool                     => carbonstate_vars%xsmrpool_patch                       , &
         xsmrpool_recover             => carbonflux_vars%xsmrpool_recover_patch                , &
         leaf_curmr                   => carbonflux_vars%leaf_curmr_patch                      , &
         froot_curmr                  => carbonflux_vars%froot_curmr_patch                     , &
         livestem_curmr               => carbonflux_vars%livestem_curmr_patch                  , &
         livecroot_curmr              => carbonflux_vars%livecroot_curmr_patch                 , &
         grain_curmr                  => carbonflux_vars%grain_curmr_patch                     , &
         leaf_xsmr                    => carbonflux_vars%leaf_xsmr_patch                       , &
         froot_xsmr                   => carbonflux_vars%froot_xsmr_patch                      , &
         livestem_xsmr                => carbonflux_vars%livestem_xsmr_patch                   , &
         livecroot_xsmr               => carbonflux_vars%livecroot_xsmr_patch                  , &
         grain_xsmr                   => carbonflux_vars%grain_xsmr_patch                      , &
         allocation_leaf              => carbonflux_vars%allocation_leaf                       , &
         allocation_stem              => carbonflux_vars%allocation_stem                       , &
         allocation_froot             => carbonflux_vars%allocation_froot                      , &

         c13cf => c13_carbonflux_vars, &
         c14cf => c14_carbonflux_vars  &
         )

!
!    !-------------------------------------------------------------------

      ! start new pft loop to distribute the available N between the
      ! competing patches on the basis of relative demand, and allocate C and N to
      ! new growth and storage

      do fp=1,num_soilp
         p = filter_soilp(fp)
         c = pft%column(p)

             ! 'ECA' or 'MIC' mode
             ! dynamic allocation based on light limitation (more woody growth) vs nutrient limitations (more fine root growth)
             ! set allocation coefficients
             N_lim_factor(p) = cn_scalar(p) ! N stress factor
             P_lim_factor(p) = cp_scalar(p) ! P stress factor

             if (cnallocate_carbon_only()) then
                 N_lim_factor(p) = 0.0_r8
                 P_lim_factor(p) = 0.0_r8
             else if (cnallocate_carbonnitrogen_only()) then
                 P_lim_factor(p) = 0.0_r8
             else if (cnallocate_carbonphosphorus_only()) then
                 N_lim_factor(p) = 0.0_r8
             end if
             W_lim_factor(p) = 0.0_r8!w_scalar(c,1)
             do j = 1 , nlevdecomp
                 W_lim_factor(p) = w_scalar(c,j) * froot_prof(p,j)
             end do
             ! N_lim_factor/P_lim_factor ones: highly limited
             ! N_lim_factor/P_lim_factor zeros: not limited
             ! convert to 1- X, see explanation in dynamic_plant_alloc
             call dynamic_plant_alloc(min(1.0_r8-N_lim_factor(p),1.0_r8-P_lim_factor(p)),W_lim_factor(p), &
                  laisun(p)+laisha(p), allocation_leaf(p), allocation_stem(p), allocation_froot(p), woody(ivt(p)))

             f1 = allocation_froot(p) / allocation_leaf(p)
             f2 = croot_stem(ivt(p))
             f3 = allocation_stem(p) / allocation_leaf(p)

             ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
             ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
             ! There was an error in this formula in previous version, where the coefficient
             ! was 0.004 instead of 0.0025.
             ! This variable allocation is only for trees. Shrubs have a constant
             ! allocation as specified in the pft-physiology file.  The value is also used
             ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
             !if (stem_leaf(ivt(p)) == -1._r8) then
             !    f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
             !else
             !    f3 = stem_leaf(ivt(p))
             !end if

             f4 = flivewd(ivt(p))
             g1 = grperc(ivt(p))
             g2 = grpnow(ivt(p))

             cnl = leafcn(ivt(p))
             cnfr = frootcn(ivt(p))
             cnlw = livewdcn(ivt(p))
             cndw = deadwdcn(ivt(p))

             cpl =  leafcp(ivt(p))
             cpfr = frootcp(ivt(p))
             cplw = livewdcp(ivt(p))
             cpdw = deadwdcp(ivt(p))

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

             sminn_to_npool(p) = plant_n_buffer_patch(p)/taun
             sminp_to_ppool(p) = plant_p_buffer_patch(p)/taup
             plant_n_buffer_patch(p) = plant_n_buffer_patch(p) * (1._r8-dt/taun)
             plant_p_buffer_patch(p) = plant_p_buffer_patch(p) * (1._r8-dt/taun)

             if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
                avail_retransn(p) = retransn(p)/dt
                avail_retransp(p) = retransp(p)/dt
             else if (ivt(p) < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
                avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
                avail_retransp(p) = (annmax_retransp(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
             else
                avail_retransn(p) = 0.0_r8
                avail_retransp(p) = 0.0_r8
             end if

             retransn_to_npool(p) = avail_retransn(p)
             retransp_to_ppool(p) = avail_retransp(p)

             plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)
             plant_palloc(p) = sminp_to_ppool(p) + retransp_to_ppool(p)

             mr = leaf_mr(p) + froot_mr(p)
             if (woody(ivt(p)) == 1.0_r8) then
                mr = mr + livestem_mr(p) + livecroot_mr(p)
             else if (ivt(p) >= npcropmin) then
                if (croplive(p)) mr = mr + livestem_mr(p) + grain_mr(p)
             end if
             ! try to take mr from xsmr storage pool first
             if (xsmrpool(p) > 0) then
                if (mr > 0._r8 .and. (xsmrpool(p) + gpp(p)) <= mr) then
                   curmr = gpp(p)
                   curmr_ratio = curmr / mr
                   availc(p) = 0.0
                else if (mr > 0._r8 .and. (xsmrpool(p) + gpp(p)) > mr .and. xsmrpool(p) <= mr ) then
                   curmr = mr - xsmrpool(p)
                   curmr_ratio = curmr / mr
                   availc(p) = gpp(p) - (mr - xsmrpool(p))
                else if (mr > 0._r8 .and. (xsmrpool(p) + gpp(p)) > mr .and. xsmrpool(p) > mr ) then
                   curmr = 0.0
                   curmr_ratio = curmr / mr
                   availc(p) = gpp(p)
                end if
             else
                if (mr > 0._r8 .and.  gpp(p) <= mr) then
                   curmr = gpp(p)
                   curmr_ratio = curmr / mr
                   availc(p) = 0.0
                else if (mr > 0._r8 .and. gpp(p) > mr ) then
                   curmr = mr
                   curmr_ratio = curmr / mr
                   availc(p) = gpp(p) - mr
                end if
             end if

             ! carbon flux available for allocation
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
             else

                ! bug fix: set to zero otherwise xsmrpool may grow infinitely when:
                ! (1) at one timestep xsmrpool(p) <0, cpool_to_xsmrpool(p) is set a positive value
                ! (2) later on if xsmrpool(p) >0; then cpool_to_xsmrpool(p) will neither be updated by following codes nor re-set to zero
                ! (3) each time step in CNCStateUpdate1 xsmrpool(p) = xsmrpool(p) + cpool_to_xsmrpool(p)*dt
                cpool_to_xsmrpool(p) = 0.0_r8
             end if

             plant_calloc(p) = availc(p)

             ! here no down-regulation on allocatable C here, NP limitation is implemented in leaf-level NP control on GPP
             if (woody(ivt(p)) == 1.0_r8) then
                 c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
                 n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
                     (f3*(1._r8-f4)*(1._r8+f2))/cndw
                 p_allometry(p) = 1._r8/cpl + f1/cpfr + (f3*f4*(1._r8+f2))/cplw + &
                     (f3*(1._r8-f4)*(1._r8+f2))/cpdw

             else if (ivt(p) >= npcropmin) then ! skip generic crops
                 cng = graincn(ivt(p))
                 cpg = graincp(ivt(p))
                 c_allometry(p) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
                 n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
                     (f3*(1._r8-f4)*(1._r8+f2))/cndw
                 p_allometry(p) = 1._r8/cpl + f1/cpfr + f5/cpg + (f3*f4*(1._r8+f2))/cplw + &
                     (f3*(1._r8-f4)*(1._r8+f2))/cpdw

             else
                 c_allometry(p) = 1._r8+g1+f1+f1*g1
                 n_allometry(p) = 1._r8/cnl + f1/cnfr
                 p_allometry(p) = 1._r8/cpl + f1/cpfr
             end if


         ! calculate the amount of new leaf C dictated by these allocation
         ! decisions, and calculate the daily fluxes of C and N to current
         ! growth and storage pools

         ! fcur is the proportion of this day's growth that is displayed now,
         ! the remainder going into storage for display next year through the
         ! transfer pools

         ! recover default coefficient for carbon allocation to leaf,  which is possibly changed due to previous time step allocation adjustment
         nlc = plant_calloc(p) / c_allometry(p)
         ! recover allocation fraction,  which is possibly changed due to previous time step allocation adjustment
         !fcur = fcur2(ivt(p))

            ! under ECA or MIC mode, CNP stoichiometry is flexible
            ! If nutrient is limited, plant will accumulate non-structural carbon hydrate (sink strength limitation)
            ! e.g., in the model if allocatable C is too much, allocate excess C to storage pool, later could be respired
            ! Here, adjust the fraction allocate to structure vs storage pool so that:
            ! CN only mode adjust C allocation to maintain CN ratio within natural variability
            ! CP only mode adjust C allocation to maintain CP ratio within natural variability
            ! CNP mode adjust C allocation to maintain CN and CP ratio within natural variability

            if (cnallocate_carbon_only()) then ! C only mode
               ! nothing to adjust
               nlc_adjust_high = nlc
            else if (cnallocate_carbonnitrogen_only()) then ! CN only mode

               ! maximum amount of C allocated to leaf pool that could be supported by plant N allocated to leaf pool:
               ! plant_nalloc(p) / (n_allometry(p) )/ cnl * (cnl*(1 + cn_stoich_var ) )
               ! maximum amount of C allocated to leaf pool that could be supported by plant P allocated to leaf pool:
               ! plant_palloc(p) / (p_allometry(p) )/ cpl * (cpl* (1 + cp_stoich_var ) )
               ! actual amount of C allocated to leaf pool if no adjustment occur
               ! plant_calloc/c_allometry * x* (x*=1)
               ! adjust fcur* to reduce the C allocated to leaf pool
               ! x* = plant_nalloc(p) / n_allometry(p) * (1 + cn_stoich_var )  /  (plant_calloc/c_allometry)
               ! x* = plant_palloc(p) / p_allometry(p) * (1 + cp_stoich_var )  /  (plant_calloc/c_allometry)


               nlc_adjust_high = plant_nalloc(p) / n_allometry(p) * (1 + cn_stoich_var )  ! upper bound of allocatable C to leaf  to satisfy N allocation
               nlc_adjust_high = nlc_adjust_high + max((leafn(p)+leafn_storage(p) + leafn_xfer(p))* cnl *  (1 + cn_stoich_var ) - &
                  (leafc(p)+leafc_storage(p) + leafc_xfer(p)),0.0_r8)/dt ! upper bound of allocatable C to leaf account for offsetting current leaf N deficit
            else if (cnallocate_carbonphosphorus_only()) then ! CP only mode
               nlc_adjust_high = plant_palloc(p) / p_allometry(p) * (1 + cp_stoich_var )  ! upper bound of allocatable C to leaf  to satisfy P allocation
               nlc_adjust_high = nlc_adjust_high + max((leafp(p)+leafp_storage(p) + leafp_xfer(p))* cpl *  (1 + cp_stoich_var ) - &
                  (leafc(p)+leafc_storage(p) + leafc_xfer(p)),0.0_r8)/dt ! upper bound of allocatable C to leaf account for offsetting current leaf N deficit
            else !  CNP mode
               nlc_adjust_high = min(plant_nalloc(p) / n_allometry(p) * (1 + cn_stoich_var ) + max((leafn(p)+leafn_storage(p) + leafn_xfer(p))* cnl *  (1 + cn_stoich_var ) - &
                  (leafc(p)+leafc_storage(p) + leafc_xfer(p)),0.0_r8)/dt, &
                  plant_palloc(p) / p_allometry(p) * (1 + cp_stoich_var ) + max((leafp(p)+leafp_storage(p) + leafp_xfer(p))* cpl *  (1 + cp_stoich_var ) - &
                  (leafc(p)+leafc_storage(p) + leafc_xfer(p)),0.0_r8)/dt)
            end if

            ! calculate excess carbon
            ! put excess carbon into respiration storage pool (if nlc > nlc_adjust_high)
            nlc = max(nlc  - nlc_adjust_high,0.0_r8)
            cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * fcur * (1 + g1)
            cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * (1._r8 - fcur) * (1 + g1)
            cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f1 * fcur * (1 + g1)
            cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f1 * (1._r8 - fcur) * (1 + g1)
            if (woody(ivt(p)) == 1._r8) then
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * f4 * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * f4 * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * (1._r8 - f4) * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * (1._r8 - f4) * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * f4 * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * f4 * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * (1._r8 - f4) * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur) * (1 + g1)
            end if
            if (ivt(p) >= npcropmin) then ! skip 2 generic crops
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * f4 * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * f4 * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * (1._r8 - f4) * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f3 * (1._r8 - f4) * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * f4 * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * f4 * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * (1._r8 - f4) * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur) * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f5 * fcur * (1 + g1)
               cpool_to_xsmrpool(p)  = cpool_to_xsmrpool(p) + nlc * f5 * (1._r8 -fcur) * (1 + g1)
            end if

            ! updated allocation if necessary
            nlc = min(nlc_adjust_high ,plant_calloc(p) / c_allometry(p) )


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
         ! recover default coefficient for carbon allocation to leaf,  which is possibly changed due to previous time step allocation adjustment
         !nlc = plant_calloc(p) / c_allometry(p)
         ! recover allocation fraction,  which is possibly changed due to previous time step allocation adjustment
         !fcur = fcur2(ivt(p))

            if (cnallocate_carbon_only()) then ! C only mode
               ! nothing to adjust
            else ! CN/ CP/ CNP mode
            !   ! minimum amount of C allocated to structural leaf pool that could be supported by plant N allocated to structural leaf pool:
            !   ! plant_nalloc(p) / (n_allometry(p) )/ cnl * (cnl*(1 - cn_stoich_var ) )*x* (x*=1)
            !   ! minimum amount of C allocated to structural leaf pool that could be supported by plant P allocated to structural leaf pool:
            !   ! plant_palloc(p) / (p_allometry(p) )/ cpl * (cpl* (1 - cp_stoich_var ) )*x* (x*=1)
            !   ! actual amount of C allocated to structural leaf pool if no adjustment occur
            !   ! plant_calloc/c_allometry
            !   ! adjust fcur* to reduce the NP allocated to structural leaf pool
            !   ! x* = (plant_calloc/c_allometry)* fcur /(plant_nalloc(p) / n_allometry(p) * (1 - cn_stoich_var ) )
            !   ! x* = (plant_calloc/c_allometry)* fcur /(plant_palloc(p) / p_allometry(p) * (1 - cp_stoich_var ) )
            !
            !   if (plant_nalloc(p) / n_allometry(p) / cnl * fcur > cpool_to_leafc(p) / (cnl *  (1 - cn_stoich_var ) ) ) then ! excess N
            !      fcur = cpool_to_leafc(p) / (plant_nalloc(p) / n_allometry(p) * (1 - cn_stoich_var ) )
            !   end if
               nlc = plant_nalloc(p) / n_allometry(p)
            end if


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

         ! corresponding P fluxes
         ! recover default coefficient for carbon allocation to leaf,  which is possibly changed due to previous time step allocation adjustment
         !nlc = plant_calloc(p) / c_allometry(p)
         ! recover allocation fraction,  which is possibly changed due to previous time step allocation adjustment
         !fcur = fcur2(ivt(p))

            if (cnallocate_carbon_only()) then ! C only mode
               ! nothing to adjust
            else ! CN/ CP/ CNP mode
            !   if (plant_palloc(p) / p_allometry(p) / cpl * fcur > cpool_to_leafc(p) / (cpl *  (1 - cp_stoich_var ) ) ) then ! excess P
            !      fcur = cpool_to_leafc(p) / (plant_palloc(p) / p_allometry(p) * (1 - cp_stoich_var ) )
            !   end if
               nlc = plant_palloc(p) / p_allometry(p)
            end if


         ppool_to_leafp(p)          = (nlc / cpl) * fcur
         ppool_to_leafp_storage(p)  = (nlc / cpl) * (1._r8 - fcur)
         ppool_to_frootp(p)         = (nlc * f1 / cpfr) * fcur
         ppool_to_frootp_storage(p) = (nlc * f1 / cpfr) * (1._r8 - fcur)
         if (woody(ivt(p)) == 1._r8) then
            ppool_to_livestemp(p)          = (nlc * f3 * f4 / cplw) * fcur
            ppool_to_livestemp_storage(p)  = (nlc * f3 * f4 / cplw) * (1._r8 -fcur)
            ppool_to_deadstemp(p)          = (nlc * f3 * (1._r8 - f4) / cpdw) *fcur
            ppool_to_deadstemp_storage(p)  = (nlc * f3 * (1._r8 - f4) / cpdw) *(1._r8 - fcur)
            ppool_to_livecrootp(p)         = (nlc * f2 * f3 * f4 / cplw) * fcur
            ppool_to_livecrootp_storage(p) = (nlc * f2 * f3 * f4 / cplw) * (1._r8 -fcur)
            ppool_to_deadcrootp(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* fcur
            ppool_to_deadcrootp_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* (1._r8 - fcur)
         end if
         if (ivt(p) >= npcropmin) then ! skip 2 generic crops
            cpg = graincp(ivt(p))
            ppool_to_livestemp(p)          = (nlc * f3 * f4 / cplw) * fcur
            ppool_to_livestemp_storage(p)  = (nlc * f3 * f4 / cplw) * (1._r8 -fcur)
            ppool_to_deadstemp(p)          = (nlc * f3 * (1._r8 - f4) / cpdw) * fcur
            ppool_to_deadstemp_storage(p)  = (nlc * f3 * (1._r8 - f4) / cpdw) *(1._r8 - fcur)
            ppool_to_livecrootp(p)         = (nlc * f2 * f3 * f4 / cplw) * fcur
            ppool_to_livecrootp_storage(p) = (nlc * f2 * f3 * f4 / cplw) * (1._r8 -fcur)
            ppool_to_deadcrootp(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* fcur
            ppool_to_deadcrootp_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* (1._r8 - fcur)
            ppool_to_grainp(p)             = (nlc * f5 / cpg) * fcur
            ppool_to_grainp_storage(p)     = (nlc * f5 / cpg) * (1._r8 -fcur)
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

         ! ECA root NP uptake is based on kinetics, plant CNP stoichiometry can vary even
         ! when certain element is set to not limiting (e.g., P not limiting under CN mode)
         ! additional supplement N/P come from first soil layer
         ! must ensure plant get enough N or P or both to maintain its stoichiometry:
         ! (1) maintain plant PC stoichiometry at optimal ratio under CN mode
         ! (2) maintain plant NC stoichiometry at optimal ratio under CP mode
         ! (3) maintain plant PC/NC stoichiometry at optimal ratios under C mode

             if (cnallocate_carbon_only() .or. cnallocate_carbonphosphorus_only()) then

                 supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) + cpool_to_leafc(p) / cnl - npool_to_leafn(p)
                 supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) + cpool_to_leafc_storage(p) / cnl -  npool_to_leafn_storage(p)
                 supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) + cpool_to_frootc(p) / cnfr - npool_to_frootn(p)
                 supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) + cpool_to_frootc_storage(p) / cnfr- npool_to_frootn_storage(p)

                 npool_to_leafn(p) = cpool_to_leafc(p) / cnl
                 npool_to_leafn_storage(p) =  cpool_to_leafc_storage(p) / cnl
                 npool_to_frootn(p) = cpool_to_frootc(p) / cnfr
                 npool_to_frootn_storage(p) = cpool_to_frootc_storage(p) / cnfr

                 if (woody(ivt(p)) == 1._r8) then
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livestemc(p) / cnlw - npool_to_livestemn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livestemc_storage(p) / cnlw - npool_to_livestemn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadstemc(p) / cndw - npool_to_deadstemn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadstemc_storage(p) / cndw - npool_to_deadstemn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livecrootc(p) / cnlw - npool_to_livecrootn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livecrootc_storage(p) / cnlw - npool_to_livecrootn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadcrootc(p) / cndw - npool_to_deadcrootn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadcrootc_storage(p) / cndw - npool_to_deadcrootn_storage(p)

                     npool_to_livestemn(p)  =  cpool_to_livestemc(p) / cnlw
                     npool_to_livestemn_storage(p) =  cpool_to_livestemc_storage(p) / cnlw
                     npool_to_deadstemn(p) = cpool_to_deadstemc(p) / cndw
                     npool_to_deadstemn_storage(p) = cpool_to_deadstemc_storage(p) / cndw
                     npool_to_livecrootn(p) =  cpool_to_livecrootc(p) / cnlw
                     npool_to_livecrootn_storage(p) = cpool_to_livecrootc_storage(p) / cnlw
                     npool_to_deadcrootn(p) = cpool_to_deadcrootc(p) / cndw
                     npool_to_deadcrootn_storage(p) = cpool_to_deadcrootc_storage(p) / cndw
                 end if
                 if (ivt(p) >= npcropmin) then ! skip 2 generic crops
                     cng = graincn(ivt(p))
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livestemc(p) / cnlw - npool_to_livestemn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livestemc_storage(p) / cnlw - npool_to_livestemn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadstemc(p) / cndw - npool_to_deadstemn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadstemc_storage(p) / cndw - npool_to_deadstemn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livecrootc(p) / cnlw - npool_to_livecrootn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_livecrootc_storage(p) / cnlw - npool_to_livecrootn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadcrootc(p) / cndw - npool_to_deadcrootn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_deadcrootc_storage(p) / cndw - npool_to_deadcrootn_storage(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_grainc(p) / cng - npool_to_grainn(p)
                     supplement_to_sminn_vr(c,1) = supplement_to_sminn_vr(c,1) +  cpool_to_grainc_storage(p) / cng - npool_to_grainn_storage(p)

                     npool_to_livestemn(p) = cpool_to_livestemc(p) / cnlw
                     npool_to_livestemn_storage(p) = cpool_to_livestemc_storage(p) / cnlw
                     npool_to_deadstemn(p) = cpool_to_deadstemc(p) / cndw
                     npool_to_deadstemn_storage(p) = cpool_to_deadstemc_storage(p) / cndw
                     npool_to_livecrootn(p) = cpool_to_livecrootc(p) / cnlw
                     npool_to_livecrootn_storage(p) = cpool_to_livecrootc_storage(p) / cnlw
                     npool_to_deadcrootn(p) = cpool_to_deadcrootc(p) / cndw
                     npool_to_deadcrootn_storage(p) = cpool_to_deadcrootc_storage(p) / cndw
                     npool_to_grainn(p) = cpool_to_grainc(p) / cng
                     npool_to_grainn_storage(p) =  cpool_to_grainc_storage(p) / cng
                 end if

             else if (cnallocate_carbon_only() .or. cnallocate_carbonnitrogen_only()) then

                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_leafc(p) / cpl - ppool_to_leafp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_leafc_storage(p) / cpl - ppool_to_leafp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_frootc(p) / cpfr - ppool_to_frootp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_frootc_storage(p) / cpfr - ppool_to_frootp_storage(p)

                     ppool_to_leafp(p) = cpool_to_leafc(p) / cpl
                     ppool_to_leafp_storage(p) =  cpool_to_leafc_storage(p) / cpl
                     ppool_to_frootp(p) = cpool_to_frootc(p) / cpfr
                     ppool_to_frootp_storage(p) = cpool_to_frootc_storage(p) / cpfr

                 if (woody(ivt(p)) == 1._r8) then
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livestemc(p) / cplw - ppool_to_livestemp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livestemc_storage(p) / cplw - ppool_to_livestemp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadstemc(p) /cpdw - ppool_to_deadstemp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadstemc_storage(p)  / cpdw- ppool_to_deadstemp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livecrootc(p) / cplw - ppool_to_livecrootp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livecrootc_storage(p) / cplw - ppool_to_livecrootp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadcrootc(p) / cpdw - ppool_to_deadcrootp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadcrootc_storage(p) / cpdw - ppool_to_deadcrootp_storage(p)

                     ppool_to_livestemp(p) = cpool_to_livestemc(p) / cplw
                     ppool_to_livestemp_storage(p) = cpool_to_livestemc_storage(p) / cplw
                     ppool_to_deadstemp(p) = cpool_to_deadstemc(p) / cpdw
                     ppool_to_deadstemp_storage(p) = cpool_to_deadstemc_storage(p) / cpdw
                     ppool_to_livecrootp(p) = cpool_to_livecrootc(p) / cplw
                     ppool_to_livecrootp_storage(p) = cpool_to_livecrootc_storage(p) / cplw
                     ppool_to_deadcrootp(p) = cpool_to_deadcrootc(p) / cpdw
                     ppool_to_deadcrootp_storage(p) = cpool_to_deadcrootc_storage(p) / cpdw
                 end if
                 if (ivt(p) >= npcropmin) then ! skip 2 generic crops
                     cpg = graincp(ivt(p))
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livestemc(p) / cplw - ppool_to_livestemp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livestemc_storage(p) / cplw - ppool_to_livestemp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadstemc(p) /cpdw - ppool_to_deadstemp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadstemc_storage(p)  / cpdw- ppool_to_deadstemp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livecrootc(p) / cplw - ppool_to_livecrootp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_livecrootc_storage(p) / cplw - ppool_to_livecrootp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadcrootc(p) / cpdw - ppool_to_deadcrootp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_deadcrootc_storage(p) / cpdw - ppool_to_deadcrootp_storage(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_grainc(p) / cpg - ppool_to_grainp(p)
                     supplement_to_sminp_vr(c,1) = supplement_to_sminp_vr(c,1) +  cpool_to_grainc_storage(p) / cpg - ppool_to_grainp_storage(p)

                     ppool_to_livestemp(p) = cpool_to_livestemc(p) / cplw
                     ppool_to_livestemp_storage(p) = cpool_to_livestemc_storage(p) / cplw
                     ppool_to_deadstemp(p) = cpool_to_deadstemc(p) / cpdw
                     ppool_to_deadstemp_storage(p) = cpool_to_deadstemc_storage(p) / cpdw
                     ppool_to_livecrootp(p) = cpool_to_livecrootc(p) / cplw
                     ppool_to_livecrootp_storage(p) = cpool_to_livecrootc_storage(p) / cplw
                     ppool_to_deadcrootp(p) = cpool_to_deadcrootc(p) / cpdw
                     ppool_to_deadcrootp_storage(p) = cpool_to_deadcrootc_storage(p) / cpdw
                     ppool_to_grainp(p) = cpool_to_grainc(p) / cpg
                     ppool_to_grainp_storage(p) =  cpool_to_grainc_storage(p) / cpg
                 end if

             end if

      end do ! end pft loop

      !----------------------------------------------------------------

    end associate

 end subroutine CNAllocation3_PlantCNPAlloc


!-----------------------------------------------------------------------

    subroutine dynamic_plant_alloc( nutrient_scalar, water_scalar, laindex, alloc_leaf, alloc_stem, alloc_froot, woody)

    ! !DESCRIPTION
    ! Added by Qing Zhu 2015 based on P. Friedlingstein DOI: 10.1046/j.1365-2486.1999.00269.x
    ! allocation coefficients for leaf, stem and root are not fixed
    ! update allocation coefficients based on nutrient and light availability
    ! (1) light limited, allocate more C into stem
    ! (2) nutrient/water limited, allocate more C into root

    ! !USES:

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: nutrient_scalar  ! scalar for nutrient availability
    real(r8), intent(in) :: water_scalar    !  scalar for water availability
    real(r8), intent(in) :: laindex      ! lai
    real(r8), intent(out) :: alloc_leaf
    real(r8), intent(out) :: alloc_stem
    real(r8), intent(out) :: alloc_froot
    real(r8), intent(in) :: woody

    !! variables
    real(r8) :: laindex_max = 8
    real(r8) :: allocmin_leaf = 0.25
    !real(r8) :: allocmax_leaf = 0.5
    real(r8) :: alloc_r0 = 0.25     ! initial allocation to roots for unlimiting conditions
    real(r8) :: alloc_s0 = 0.25     ! initial allocation to stem for unlimiting conditions
    real(r8) :: klight_ex = 0.5    ! light extinction parameter
    real(r8) :: light_scalar       ! scalar for light availability
    real(r8) :: nu_scalar
    real(r8) :: w_scalar

    ! general framework P. Friedlingstein DOI: 10.1046/j.1365-2486.1999.00269.x
    ! allocation to a certain compartment A = sum(X)/(X + Y)
    ! increase resource X availability lead to increase allocation to A
    ! increase resource Y availability lead to decrease allocation to A

    ! for nu_scalar from 0->1, system from high nutrient limited -> non-nutrient limited
    ! nutrient resource availability increase, root allocation decrease
    ! in this case nu_scalar is the availability scalar

    ! light scalar lai high->low, light_scalar 0->1
    ! light availability increase, allocation to wood decrease
    ! define the light availability scalar based on LAI
    light_scalar = exp (-klight_ex * laindex)

    ! adjust scalar for numerical stability purposes
    light_scalar = max( 0.1_r8, min( 1.0_r8, light_scalar ) )
    nu_scalar = max( 0.1_r8, min( 1.0_r8, nutrient_scalar ) )
    w_scalar = max( 0.1_r8, min( 1.0_r8, water_scalar ) )

    ! root allocation
    alloc_froot = alloc_r0 * 3.0_r8 * light_scalar / (light_scalar + 2.0_r8 * min(nu_scalar,w_scalar))
    alloc_froot = min(alloc_froot, 0.4_r8)

    ! stem allocation
    if (woody == 1.0_r8) then
       alloc_stem = alloc_s0 * 3.0_r8 *  min(nu_scalar,w_scalar) / (2.0_r8 * light_scalar + min(nu_scalar,w_scalar))
    else
       alloc_stem = 0.0_r8
    end if
    ! leaf allocation
    alloc_leaf = 1.0_r8 - (alloc_froot + alloc_stem)

    ! adjustment under extreme nutrient/light limitation condition
    if (alloc_leaf < allocmin_leaf) then
       alloc_leaf = allocmin_leaf
       alloc_froot = alloc_froot * (1-allocmin_leaf) / (alloc_froot + alloc_stem)
       alloc_stem = 1.0 - alloc_leaf - alloc_froot
    end if

    ! if lai greater than laimax then no allocation to leaf; leaf allocation goes to stem or fine root
    if (laindex > laindex_max) then
       if (woody == 1.0_r8) then
          alloc_stem = alloc_stem + alloc_leaf - 0.01_r8
       else
          alloc_froot = alloc_froot + alloc_leaf - 0.01_r8
       end if
       alloc_leaf = 0.01_r8
    end if

  end subroutine dynamic_plant_alloc

!-----------------------------------------------------------------------

  subroutine update_plant_stoichiometry(num_soilp, filter_soilp, &
           carbonstate_vars, nitrogenstate_vars, phosphorusstate_vars)

    ! ! DESCRIPTION
    ! Added by Qing Zhu 2015
    ! update leaf, fine root, livewood, deadwood stoichiometry

    ! ! USES:
    use pftvarcon      , only : noveg,npcropmin
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    integer                    , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                    , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(carbonstate_type)     , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type)   , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

    !! variables
    integer :: p,fp                                            !indices

    associate(                                                                          &
        ivt                          => pft%itype                                     , & ! Input:  [integer  (:) ]  pft vegetation type
        woody                        => ecophyscon%woody                              , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)

        leafc                        => carbonstate_vars%leafc_patch                  , &
        frootc                       => carbonstate_vars%frootc_patch                 , &
        livestemc                    => carbonstate_vars%livestemc_patch              , &
        livecrootc                   => carbonstate_vars%livecrootc_patch             , &
        deadstemc                    => carbonstate_vars%deadstemc_patch              , &
        deadcrootc                   => carbonstate_vars%deadcrootc_patch             , &
        grainc                       => carbonstate_vars%grainc_patch                 , &
        leafc_storage                => carbonstate_vars%leafc_storage_patch	      , &
        frootc_storage               => carbonstate_vars%frootc_storage_patch         , &
        livestemc_storage            => carbonstate_vars%livestemc_storage_patch      , &
        deadstemc_storage            => carbonstate_vars%deadstemc_storage_patch      , &
        livecrootc_storage           => carbonstate_vars%livecrootc_storage_patch     , &
        deadcrootc_storage           => carbonstate_vars%deadcrootc_storage_patch     , &
        grainc_storage               => carbonstate_vars%grainc_storage_patch         , &
        leafc_xfer                   => carbonstate_vars%leafc_xfer_patch             , &
        frootc_xfer                  => carbonstate_vars%frootc_xfer_patch            , &
        livestemc_xfer               => carbonstate_vars%livestemc_xfer_patch         , &
        deadstemc_xfer               => carbonstate_vars%deadstemc_xfer_patch         , &
        livecrootc_xfer              => carbonstate_vars%livecrootc_xfer_patch 	      , &
        deadcrootc_xfer              => carbonstate_vars%deadcrootc_xfer_patch 	      , &
        grainc_xfer                  => carbonstate_vars%grainc_xfer_patch            , &

        leafn                        => nitrogenstate_vars%leafn_patch                , &
        frootn                       => nitrogenstate_vars%frootn_patch               , &
        livestemn                    => nitrogenstate_vars%livestemn_patch            , &
        livecrootn                   => nitrogenstate_vars%livecrootn_patch           , &
        deadstemn                    => nitrogenstate_vars%deadstemn_patch            , &
        deadcrootn                   => nitrogenstate_vars%deadcrootn_patch           , &
        grainn                       => nitrogenstate_vars%grainn_patch               , &
        leafn_storage                => nitrogenstate_vars%leafn_storage_patch	      , &
        frootn_storage               => nitrogenstate_vars%frootn_storage_patch       , &
        livestemn_storage            => nitrogenstate_vars%livestemn_storage_patch    , &
        deadstemn_storage            => nitrogenstate_vars%deadstemn_storage_patch    , &
        livecrootn_storage           => nitrogenstate_vars%livecrootn_storage_patch   , &
        deadcrootn_storage           => nitrogenstate_vars%deadcrootn_storage_patch   , &
        grainn_storage               => nitrogenstate_vars%grainn_storage_patch       , &
        leafn_xfer                   => nitrogenstate_vars%leafn_xfer_patch           , &
        frootn_xfer                  => nitrogenstate_vars%frootn_xfer_patch          , &
        livestemn_xfer               => nitrogenstate_vars%livestemn_xfer_patch       , &
        deadstemn_xfer               => nitrogenstate_vars%deadstemn_xfer_patch       , &
        livecrootn_xfer              => nitrogenstate_vars%livecrootn_xfer_patch      , &
        deadcrootn_xfer              => nitrogenstate_vars%deadcrootn_xfer_patch      , &
        grainn_xfer                  => nitrogenstate_vars%grainn_xfer_patch          , &

        leafp                        => phosphorusstate_vars%leafp_patch              , &
        frootp                       => phosphorusstate_vars%frootp_patch             , &
        livestemp                    => phosphorusstate_vars%livestemp_patch          , &
        livecrootp                   => phosphorusstate_vars%livecrootp_patch         , &
        deadstemp                    => phosphorusstate_vars%deadstemp_patch          , &
        deadcrootp                   => phosphorusstate_vars%deadcrootp_patch         , &
        grainp                       => phosphorusstate_vars%grainp_patch             , &
        leafp_storage                => phosphorusstate_vars%leafp_storage_patch      , &
        frootp_storage               => phosphorusstate_vars%frootp_storage_patch     , &
        livestemp_storage            => phosphorusstate_vars%livestemp_storage_patch  , &
        deadstemp_storage            => phosphorusstate_vars%deadstemp_storage_patch  , &
        livecrootp_storage           => phosphorusstate_vars%livecrootp_storage_patch , &
        deadcrootp_storage           => phosphorusstate_vars%deadcrootp_storage_patch , &
        grainp_storage               => phosphorusstate_vars%grainp_storage_patch     , &
        leafp_xfer                   => phosphorusstate_vars%leafp_xfer_patch         , &
        frootp_xfer                  => phosphorusstate_vars%frootp_xfer_patch        , &
        livestemp_xfer               => phosphorusstate_vars%livestemp_xfer_patch     , &
        deadstemp_xfer               => phosphorusstate_vars%deadstemp_xfer_patch     , &
        livecrootp_xfer              => phosphorusstate_vars%livecrootp_xfer_patch    , &
        deadcrootp_xfer              => phosphorusstate_vars%deadcrootp_xfer_patch    , &
        grainp_xfer                  => phosphorusstate_vars%grainp_xfer_patch        , &

        actual_leafcn                => nitrogenstate_vars%actual_leafcn              , &
        actual_frootcn               => nitrogenstate_vars%actual_frootcn             , &
        actual_livewdcn              => nitrogenstate_vars%actual_livewdcn            , &
        actual_deadwdcn              => nitrogenstate_vars%actual_deadwdcn            , &
        actual_graincn               => nitrogenstate_vars%actual_graincn             , &

        actual_leafcp                => phosphorusstate_vars%actual_leafcp            , &
        actual_frootcp               => phosphorusstate_vars%actual_frootcp           , &
        actual_livewdcp              => phosphorusstate_vars%actual_livewdcp          , &
        actual_deadwdcp              => phosphorusstate_vars%actual_deadwdcp          , &
        actual_graincp               => phosphorusstate_vars%actual_graincp             &
        )

    ! loop over patches
    do fp=1,num_soilp
        p = filter_soilp(fp)
        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
            actual_leafcn(p) = (leafc(p) + leafc_storage(p) + leafc_xfer(p)) / max((leafn(p) + leafn_storage(p) + leafn_xfer(p)),1e-20_r8)
            actual_frootcn(p) = (frootc(p) + frootc_storage(p) + frootc_xfer(p)) / max((frootn(p) + frootn_storage(p) + frootn_xfer(p)),1e-20_r8)
            actual_leafcp(p) = (leafc(p) + leafc_storage(p) + leafc_xfer(p)) / max((leafp(p) + leafp_storage(p) + leafp_xfer(p)),1e-20_r8)
            actual_frootcp(p) = (frootc(p) + frootc_storage(p) + frootc_xfer(p)) / max((frootp(p) + frootp_storage(p) + frootp_xfer(p)),1e-20_r8)
            if (woody(ivt(p)) == 1._r8) then
                actual_livewdcn(p) = (livestemc(p) + livecrootc(p) + livestemc_storage(p) + livecrootc_storage(p) + livestemc_xfer(p) + livecrootc_xfer(p)) / &
                                 max((livestemn(p) + livecrootn(p) + livestemn_storage(p) + livecrootn_storage(p) + livestemn_xfer(p) + livecrootn_xfer(p)),1e-20_r8)
                actual_deadwdcn(p) = (deadstemc(p) + deadcrootc(p) + deadstemc_storage(p) + deadcrootc_storage(p) + deadstemc_xfer(p) + deadcrootc_xfer(p)) / &
                                 max((deadstemn(p) + deadcrootn(p) + deadstemn_storage(p) + deadcrootn_storage(p) + deadstemn_xfer(p) + deadcrootn_xfer(p)),1e-20_r8)
                actual_livewdcp(p) = (livestemc(p) + livecrootc(p) + livestemc_storage(p) + livecrootc_storage(p) + livestemc_xfer(p) + livecrootc_xfer(p)) / &
                                 max((livestemp(p) + livecrootp(p) + livestemp_storage(p) + livecrootp_storage(p) + livestemp_xfer(p) + livecrootp_xfer(p)),1e-20_r8)
                actual_deadwdcp(p) = (deadstemc(p) + deadcrootc(p) + deadstemc_storage(p) + deadcrootc_storage(p) + deadstemc_xfer(p) + deadcrootc_xfer(p)) / &
                                 max((deadstemp(p) + deadcrootp(p) + deadstemp_storage(p) + deadcrootp_storage(p) + deadstemp_xfer(p) + deadcrootp_xfer(p)),1e-20_r8)
            end if
            if (ivt(p) >= npcropmin) then
                actual_graincn(p) = (grainc(p) + grainc_storage(p) + grainc_xfer(p)) / max((grainn(p) + grainn_storage(p) + grainn_xfer(p)),1e-20_r8)
                actual_graincp(p) = (grainc(p) + grainc_storage(p) + grainc_xfer(p)) / max((grainp(p) + grainp_storage(p) + grainp_xfer(p)),1e-20_r8)
            end if

        end if
    end do

    end associate

  end subroutine update_plant_stoichiometry

!-----------------------------------------------------------------------

end module CNAllocationBeTRMod
