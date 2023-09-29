module AllocationMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in allocation model for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use elm_varcon          , only : dzsoi_decomp
  use elm_varctl          , only : use_c13, use_c14, spinup_state
  use elm_varctl          , only : nyears_ad_carbon_only
  use elm_varctl          , only : use_fates
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use subgridAveMod       , only : p2c, p2c_1d_filter_parallel
  use CanopyStateType     , only : canopystate_type
  !!! add phosphorus
  use CNStateType               , only : cnstate_type
  use PhotosynthesisType        , only : photosyns_type
  use CropType                  , only : crop_type
  use VegetationPropertiesType  , only : veg_vp
  use LandunitType        , only : lun_pp
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : col_ws
  use ColumnDataType      , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType      , only : col_ns, col_nf, col_ps, col_pf
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_cs, veg_ns, veg_nf, veg_ps, veg_pf
  use VegetationDataType  , only : veg_cf, c13_veg_cf, c14_veg_cf
  ! bgc interface & pflotran module switches
  use elm_varctl          , only: use_elm_interface,use_elm_bgc, use_pflotran, pf_cmode
  use elm_varctl          , only : nu_com
  use SoilStatetype       , only : soilstate_type
  use elm_varctl          , only : NFIX_PTASE_plant
  !!!!!use ELMFatesInterfaceMod  , only : hlm_fates_interface_type
  use elm_varctl      , only: iulog
  use shr_infnan_mod  , only: nan => shr_infnan_nan 
  
  !
  implicit none
  save
  ! pflotran
  private :: calc_nuptake_prof
  private :: calc_puptake_prof
  private :: NAllocationRD
  private :: PAllocationRD
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readCNAllocParams
  public :: AllocationInit         ! Initialization
!  public :: Allocation             ! run method
  !-----------------------------------------------------------------------------------------------------
  ! Allocation is divided into 3 subroutines/phases:
  public :: Allocation1_PlantNPDemand     !Plant N/P Demand;       called in EcosystemDynNoLeaching1
  public :: Allocation2_ResolveNPLimit    !Resolve N/P Limitation; called in SoilLittDecompAlloc
  public :: Allocation3_PlantCNPAlloc     !Plant C/N/P Allocation; called in SoilLittDecompAlloc2
  !-----------------------------------------------------------------------------------------------------
  public :: dynamic_plant_alloc        ! dynamic plant carbon allocation based on different nutrient stress

  type :: AllocParamsType

     real(r8), pointer :: bdnr              => null() ! bulk denitrification rate (1/s)
     real(r8), pointer :: dayscrecover      => null() ! number of days to recover negative cpool
     real(r8), pointer :: compet_plant_no3  => null() ! (unitless) relative compettiveness of plants for NO3
     real(r8), pointer :: compet_plant_nh4  => null() ! (unitless) relative compettiveness of plants for NH4
     real(r8), pointer :: compet_decomp_no3 => null() ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8), pointer :: compet_decomp_nh4 => null() ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8), pointer :: compet_denit      => null() ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8), pointer :: compet_nit        => null() ! (unitless) relative competitiveness of nitrifiers for NH4

  end type AllocParamsType
  !
  ! AllocParamsInst is populated in readCNAllocParams which is called in
  type(AllocParamsType)  ,  public ::  AllocParamsInst
  !$acc declare create(AllocParamsInst)

  !
  ! !PUBLIC DATA MEMBERS:
  character(len=3), parameter, public :: suplnAll='ALL'  ! Supplemental Nitrogen for all PFT's
  character(len=4), parameter, public :: suplnNon='NONE' ! No supplemental Nitrogen
  character(len=15), public :: suplnitro = suplnNon      ! Supplemental Nitrogen mode
  !! add phosphorus  - X. YANG
  character(len=3), parameter, public :: suplpAll='ALL'  ! Supplemental Phosphorus for all PFT's
  character(len=4), parameter, public :: suplpNon='NONE' ! No supplemental Phosphorus
  character(len=15), public :: suplphos = suplpAll    ! Supplemental Phosphorus mode
  !! add competition, - Q. Zhu
  logical,          public :: nu_com_leaf_physiology = .false.
  logical,          public :: nu_com_root_kinetics   = .false.
  logical,          public :: nu_com_phosphatase     = .false.
  logical,          public :: nu_com_nfix            = .false.

  !$acc declare create(nu_com_leaf_physiology)
  !$acc declare create(nu_com_root_kinetics  )
  !$acc declare create(nu_com_phosphatase    )
  !$acc declare create(nu_com_nfix           )
  !
  ! !PRIVATE DATA MEMBERS:
  real(r8)              :: bdnr                 !bulk denitrification rate (1/s)
  real(r8)              :: dayscrecover         !number of days to recover negative cpool
  real(r8), allocatable :: arepr(:)             !reproduction allocation coefficient
  real(r8), allocatable :: aroot(:)             !root allocation coefficient

  !$acc declare create(bdnr                )
  !$acc declare create(dayscrecover        )
  !$acc declare create(arepr(:)            )
  !$acc declare create(aroot(:)            )

  logical :: crop_supln  = .false.    !Prognostic crop receives supplemental Nitrogen
  
  real(r8), allocatable,target :: veg_rootc_bigleaf(:,:)        ! column-level fine-root biomas kgc/m3
  integer,  pointer :: ft_index_bigleaf(:)                      ! array holding the pft index of each competitor

  ! ECA parameters
  ! scaling factor for plant fine root biomass to calculate nutrient carrier enzyme abundance                                         
  real(r8), parameter :: e_plant_scalar  = 0.0000125_r8 
  
  ! scaling factor for plant fine root biomass to calculate nutrient carrier enzyme abundance                                         
  real(r8), parameter :: e_decomp_scalar = 0.05_r8      

  !$acc declare create(e_decomp_scalar)
  !$acc declare create(e_plant_scalar)
  
  !$acc declare copyin(crop_supln)
  !!!!$acc declare create(filter_pcomp(:))
  !$acc declare create(veg_rootc_bigleaf(:,:))
  !$acc declare create(ft_index_bigleaf)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNAllocParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'AllocParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    allocate(AllocParamsInst%bdnr              )
    allocate(AllocParamsInst%dayscrecover      )
    allocate(AllocParamsInst%compet_plant_no3  )
    allocate(AllocParamsInst%compet_plant_nh4  )
    allocate(AllocParamsInst%compet_decomp_no3 )
    allocate(AllocParamsInst%compet_decomp_nh4 )
    allocate(AllocParamsInst%compet_denit      )
    allocate(AllocParamsInst%compet_nit        )
    ! read in parameters
    tString='bdnr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%bdnr=tempr

    tString='dayscrecover'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%dayscrecover=tempr

    tString='compet_plant_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_plant_no3=tempr

    tString='compet_plant_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_plant_nh4=tempr

    tString='compet_decomp_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_decomp_no3=tempr

    tString='compet_decomp_nh4'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_decomp_nh4=tempr

    tString='compet_denit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_denit=tempr

    tString='compet_nit'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    AllocParamsInst%compet_nit=tempr

  end subroutine readCNAllocParams

  !-----------------------------------------------------------------------
  subroutine AllocationInit ( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use elm_varcon      , only: secspday, spval
    use clm_time_manager, only: get_step_size, get_curr_date
    use elm_varpar      , only: crop_prog
    use elm_varctl      , only: iulog
    use elm_varctl      , only : carbon_only
    use elm_varctl      , only : carbonnitrogen_only
    use elm_varctl      , only : carbonphosphorus_only
    use elm_varpar      , only: nlevdecomp
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=),isnan => shr_infnan_isnan
    
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !#fates_py type(hlm_fates_interface_type), intent(in) :: elm_fates  ! This will be needed in soon
                                                             ! to be released features
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'AllocationInit'
    real(r8) :: dt
    integer ::  yr, mon, day, sec
    integer :: f
    integer :: max_comps  ! maximum number of possible plant competitors
    ! elm big-leaf: number of pfts/patches
    ! fates: number of cohorts in the column
    !-----------------------------------------------------------------------

    if ( crop_prog )then
       allocate(arepr(bounds%begp:bounds%endp)); arepr(bounds%begp : bounds%endp) = nan
       allocate(aroot(bounds%begp:bounds%endp)); aroot(bounds%begp : bounds%endp) = nan
    end if

    ! Allocate scratch space for ECA and FATES/ECA

    if (nu_com .eq. 'ECA' .or. nu_com .eq. 'MIC') then
       if (use_fates) then
       else
          max_comps = bounds%endp-bounds%begp+1
          !allocate(filter_pcomp(max_comps)); filter_pcomp(:) = -1
          allocate(ft_index_bigleaf(bounds%begp:bounds%endp)); ft_index_bigleaf(bounds%begp:bounds%endp) = -1
          allocate(veg_rootc_bigleaf(bounds%begp:bounds%endp,1:nlevdecomp)); veg_rootc_bigleaf(bounds%begp:bounds%endp,1:nlevdecomp) = nan
       end if
    end if

    ! Allocate scratch space for ECA and FATES/ECA

    if (nu_com .eq. 'ECA' .or. nu_com .eq. 'MIC') then
       if (.not.use_fates) then
          allocate(ft_index_bigleaf(bounds%begp:bounds%endp)); ft_index_bigleaf(bounds%begp:bounds%endp) = -1
          allocate(veg_rootc_bigleaf(bounds%begp:bounds%endp,1:nlevdecomp)); veg_rootc_bigleaf(bounds%begp:bounds%endp,1:nlevdecomp) = nan
       end if
    end if

    
    ! set time steps
    dt = real( get_step_size(), r8 )

    ! set space-and-time parameters from parameter file
    bdnr         = AllocParamsInst%bdnr * (dt/secspday)
    dayscrecover = AllocParamsInst%dayscrecover

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

    select case(nu_com)
        case('RD') ! relative demand mode, same as CLM-CNP Yang 2014
            nu_com_leaf_physiology = .false.
            nu_com_root_kinetics   = .false.
            nu_com_phosphatase     = .false.
            nu_com_nfix            = .false.
        case('ECA') ! ECA competition version of CLM-CNP
            nu_com_leaf_physiology = .true. ! leaf level physiology must be true if using ECA
            nu_com_root_kinetics   = .true. ! root uptake kinetics must be true if using ECA
            nu_com_phosphatase = .true.     ! new phosphatase activity
            nu_com_nfix = .true.            ! new fixation
        case('MIC') ! MIC outcompete plant version of CLM-CNP
            nu_com_leaf_physiology = .true.
            nu_com_root_kinetics   = .true.
            nu_com_phosphatase = .true.
            nu_com_nfix = .true.
    end select
    ! phosphorus conditions of plants are needed, in order to use new fixation and phosphatase
    ! activity subroutines, under carbon only or carbon nitrogen only mode, fixation and phosphatase
    ! activity are set to false
    if (carbon_only) then
        nu_com_nfix = .false.
        nu_com_phosphatase = .false.
    end if
    if (carbonnitrogen_only) then
        nu_com_phosphatase = .false.
    end if

    call get_curr_date(yr, mon, day, sec)
    if (spinup_state == 1 .and. yr .le. nyears_ad_carbon_only) then
      Carbon_only = .true.
     end if
     write(iulog,*) "CNP variables: "
     write(iulog,*) "carbon only:",carbon_only 
     write(iulog,*) "CP only: ", carbonphosphorus_only
     write(iulog,*) "CN only :", carbonnitrogen_only 
     !$acc update device(carbon_only, carbonnitrogen_only,&
     !$acc carbonphosphorus_only)

  end subroutine AllocationInit

!-------------------------------------------------------------------------------------------------
  subroutine Allocation1_PlantNPDemand (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       photosyns_vars, crop_vars, canopystate_vars, cnstate_vars, dt, yr)
    ! PHASE-1 of Allocation: loop over patches to assess the total plant N demand and P demand
    ! !USES:
    use elm_varctl, only : carbon_only
    use elm_varctl, only : carbonnitrogen_only
    use elm_varctl, only : carbonphosphorus_only
    use elm_varpar, only : nlevdecomp
    use elm_varcon, only : nitrif_n2o_loss_frac, secspday
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer                  , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)  ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    real(r8), intent(in) :: dt
    integer, intent(in) :: yr
    !
    ! !LOCAL VARIABLES:
    !
    integer :: c,p,l,j        !indices
    integer :: fp             !lake filter pft index
    integer :: fc             !lake filter column index
    real(r8):: nuptake_prof(num_soilc, 1:nlevdecomp)
    real(r8):: puptake_prof(num_soilc, 1:nlevdecomp)
    !-----------------------------------------------------------------------
    associate(                &
         nfixation_prof               => cnstate_vars%nfixation_prof_col , &
         smin_no3_vr                  => col_ns%smin_no3_vr              , & ! Input: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3
         smin_nh4_vr                  => col_ns%smin_nh4_vr              , & ! Input: [real(r8) (:,:) ]  (gN/m3) soil mineral NH4
         solutionp_vr                 => col_ps%solutionp_vr       , & ! Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral P
         plant_pdemand_col            => col_pf%plant_pdemand      , & ! Output:  [real(r8) (:,:) ]
         plant_ndemand_col            => col_nf%plant_ndemand      , & ! Output:  [real(r8) (:,:) ]
         plant_ndemand_vr_col         => col_nf%plant_ndemand_vr   , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_vr_col         => col_pf%plant_pdemand_vr   , & ! Output:  [real(r8) (:,:) ]
         plant_ndemand                => veg_nf%plant_ndemand      , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_pdemand                => veg_pf%plant_pdemand        & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)
         )

      ! set time steps
      if (spinup_state == 1 .and. yr .gt. nyears_ad_carbon_only) then
         !$acc serial default(present)
         carbon_only = .false.
         !$acc end serial
      end if

     ! loop over patches to assess the total plant N demand and P demand
      call TotalNPDemand(num_soilp, filter_soilp, photosyns_vars, &
                canopystate_vars, crop_vars, cnstate_vars, dt)
      ! now use the p2c routine to get the column-averaged plant_ndemand
      call p2c_1d_filter_parallel( num_soilc, filter_soilc, &
           plant_ndemand,  plant_ndemand_col)

      !!! add phosphorus
      call p2c_1d_filter_parallel( num_soilc, filter_soilc, &
           plant_pdemand, plant_pdemand_col)

   !!! Starting resolving N limitation
      ! pflotran will need an input from CN: modified 'sum_ndemand_vr' ('potential_immob' excluded).
      if (use_elm_interface.and.use_pflotran .and. pf_cmode) then
         !! new subroutines to calculate nuptake_prof & puptake_prof
         if (nu_com .eq. 'RD') then ! 'RD' : relative demand approach
            call calc_nuptake_prof(num_soilc, filter_soilc, cnstate_vars, nuptake_prof)
            call calc_puptake_prof(num_soilc, filter_soilc, cnstate_vars, puptake_prof)
         end if

            do j = 1, nlevdecomp
               do fc=1, num_soilc
                  c = filter_soilc(fc)
                  plant_ndemand_vr_col(c,j) = plant_ndemand_col(c) * nuptake_prof(fc,j)
                  plant_pdemand_vr_col(c,j) = plant_pdemand_col(c) * puptake_prof(fc,j)
               end do
            end do
      endif

    end associate

 end subroutine Allocation1_PlantNPDemand

 subroutine TotalNPDemand(num_soilp,filter_soilp, photosyns_vars, &
            canopystate_vars, crop_vars, cnstate_vars, dt)
    use pftvarcon   , only: npcropmin, declfact, bfact, aleaff, arootf, astemf, noveg
    use pftvarcon   , only: arooti, fleafi, allconsl, allconss, grperc, grpnow, nsoybean
    use elm_varpar  , only: nlevdecomp
    use elm_varcon  , only: nitrif_n2o_loss_frac, secspday
    !!
    integer  ,intent(in) :: num_soilp
    integer , intent(in) :: filter_soilp(:)
    type(photosyns_type)  , intent(in) :: photosyns_vars
    type(canopystate_type), intent(in) :: canopystate_vars
    type(crop_type)       , intent(in) :: crop_vars
    type(cnstate_type)    , intent(inout) :: cnstate_vars
    real(r8) , intent(in) :: dt
    !! Local Variables
    real(r8):: f5               !grain allocation parameter
    real(r8):: cng              !C:N ratio for grain (= cnlw for now; slevis)
    real(r8):: fleaf            !fraction allocated to leaf
    real(r8):: t1               !temporary variable
    real(r8):: mr                  !maintenance respiration (gC/m2/s)
    real(r8):: f1,f2,f3,f4,g1,g2   !allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw  !C:N ratios for leaf, fine root, and wood

    real(r8):: curmr, curmr_ratio  !xsmrpool temporary variables

    !! Local P variables
    real(r8):: cpl,cpfr,cplw,cpdw,cpg  !C:N ratios for leaf, fine root, and wood

    integer :: ivt, fp, p
    !-----------------------------------------------------------------------
    associate(                                                                                   &
         woody                        => veg_vp%woody               , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
         froot_leaf                   => veg_vp%froot_leaf          , & ! Input:  [real(r8) (:)   ]  allocation parameter: new fine root C per new leaf C (gC/gC)
         croot_stem                   => veg_vp%croot_stem          , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
         stem_leaf                    => veg_vp%stem_leaf           , & ! Input:  [real(r8) (:)   ]  allocation parameter: new stem c per new leaf C (gC/gC)
         flivewd                      => veg_vp%flivewd             , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
         leafcn                       => veg_vp%leafcn              , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
         frootcn                      => veg_vp%frootcn             , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)
         livewdcn                     => veg_vp%livewdcn            , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn                     => veg_vp%deadwdcn            , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
         graincn                      => veg_vp%graincn             , & ! Input:  [real(r8) (:)   ]  grain C:N (gC/gN)
         fleafcn                      => veg_vp%fleafcn             , & ! Input:  [real(r8) (:)   ]  leaf c:n during organ fill
         ffrootcn                     => veg_vp%ffrootcn            , & ! Input:  [real(r8) (:)   ]  froot c:n during organ fill
         fstemcn                      => veg_vp%fstemcn             , & ! Input:  [real(r8) (:)   ]  stem c:n during organ fill

         psnsun                       => photosyns_vars%psnsun_patch        , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         psnsha                       => photosyns_vars%psnsha_patch        , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun                   => photosyns_vars%c13_psnsun_patch    , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c13_psnsha                   => photosyns_vars%c13_psnsha_patch    , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsun                   => photosyns_vars%c14_psnsun_patch    , & ! Input:  [real(r8) (:)   ]  sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
         c14_psnsha                   => photosyns_vars%c14_psnsha_patch    , & ! Input:  [real(r8) (:)   ]  shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

         laisun                       => canopystate_vars%laisun_patch      , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index
         laisha                       => canopystate_vars%laisha_patch      , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index

         hui                          => crop_vars%gddplant_patch           , & ! Input:  [real(r8) (:)   ]  =gdd since planting (gddplant)
         leafout                      => crop_vars%gddtsoi_patch            , & ! Input:  [real(r8) (:)   ]  =gdd from top soil layer temperature

         xsmrpool                     => veg_cs%xsmrpool                    , & ! Input:  [real(r8) (:)   ]  (gC/m2) temporary photosynthate C pool
         leafc                        => veg_cs%leafc                       , & ! Input:  [real(r8) (:)   ]
         frootc                       => veg_cs%frootc                      , & ! Input:  [real(r8) (:)   ]
         livestemc                    => veg_cs%livestemc                   , & ! Input:  [real(r8) (:)   ]
         plant_ndemand_col            => col_nf%plant_ndemand               , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_col            => col_pf%plant_pdemand               , & ! Output:  [real(r8) (:,:) ]
         plant_ndemand_vr_col         => col_nf%plant_ndemand_vr            , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_vr_col         => col_pf%plant_pdemand_vr            , & ! Output:  [real(r8) (:,:) ]

         gddmaturity                  => cnstate_vars%gddmaturity_patch      , & ! Input:  [real(r8) (:)   ]  gdd needed to harvest
         huileaf                      => cnstate_vars%huileaf_patch          , & ! Input:  [real(r8) (:)   ]  heat unit index needed from planting to leaf emergence
         huigrain                     => cnstate_vars%huigrain_patch         , & ! Input:  [real(r8) (:)   ]  same to reach vegetative maturity
         croplive                     => crop_vars%croplive_patch            , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
         peaklai                      => cnstate_vars%peaklai_patch          , & ! Input:  [integer  (:)   ]  1: max allowed lai; 0: not at max
         !lgsf                        => cnstate_vars%lgsf_patch             , & ! Input:  [real(r8) (:)   ]  long growing season factor [0-1]
         aleafi                       => cnstate_vars%aleafi_patch           , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         astemi                       => cnstate_vars%astemi_patch           , & ! Output: [real(r8) (:)   ]  saved allocation coefficient from phase 2
         aleaf                        => cnstate_vars%aleaf_patch            , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
         astem                        => cnstate_vars%astem_patch            , & ! Output: [real(r8) (:)   ]  stem allocation coefficient

         !!! add phosphorus
         leafcp                       => veg_vp%leafcp              , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
         frootcp                      => veg_vp%frootcp             , & ! Input:  [real(r8) (:)   ]  fine root C:P (gC/gP)
         livewdcp                     => veg_vp%livewdcp            , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
         deadwdcp                     => veg_vp%deadwdcp            , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:P (gC/gP)
         graincp                      => veg_vp%graincp             , & ! Input:  [real(r8) (:)   ]  grain C:P (gC/gP)

         grain_flag                   => cnstate_vars%grain_flag_patch              , & ! Output: [real(r8) (:)   ]  1: grain fill stage; 0: not
         c_allometry                  => cnstate_vars%c_allometry_patch             , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
         n_allometry                  => cnstate_vars%n_allometry_patch             , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
         tempsum_potential_gpp        => cnstate_vars%tempsum_potential_gpp_patch   , & ! Output: [real(r8) (:)   ]  temporary annual sum of potential GPP
         tempmax_retransn             => cnstate_vars%tempmax_retransn_patch        , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated N pool (gN/m2)
         annsum_potential_gpp         => cnstate_vars%annsum_potential_gpp_patch    , & ! Output: [real(r8) (:)   ]  annual sum of potential GPP
         annmax_retransn              => cnstate_vars%annmax_retransn_patch         , & ! Output: [real(r8) (:)   ]  annual max of retranslocated N pool

         leaf_mr                      => veg_cf%leaf_mr                         , & ! Input:  [real(r8) (:)   ]
         froot_mr                     => veg_cf%froot_mr                        , & ! Input:  [real(r8) (:)   ]
         livestem_mr                  => veg_cf%livestem_mr                     , & ! Input:  [real(r8) (:)   ]
         livecroot_mr                 => veg_cf%livecroot_mr                    , & ! Input:  [real(r8) (:)   ]
         grain_mr                     => veg_cf%grain_mr                        , & ! Input:  [real(r8) (:)   ]
         annsum_npp                   => veg_cf%annsum_npp                      , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
         gpp                          => veg_cf%gpp_before_downreg              , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
         availc                       => veg_cf%availc                          , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
         xsmrpool_recover             => veg_cf%xsmrpool_recover                , & ! Output: [real(r8) (:)   ]  C flux assigned to recovery of negative cpool (gC/m2/s)
         psnsun_to_cpool              => veg_cf%psnsun_to_cpool                 , & ! Output: [real(r8) (:)   ]
         psnshade_to_cpool            => veg_cf%psnshade_to_cpool               , & ! Output: [real(r8) (:)   ]

         leaf_curmr                   => veg_cf%leaf_curmr                      , &
         froot_curmr                  => veg_cf%froot_curmr                     , & ! Output: [real(r8) (:)   ]
         livestem_curmr               => veg_cf%livestem_curmr                  , & ! Output: [real(r8) (:)   ]
         livecroot_curmr              => veg_cf%livecroot_curmr                 , & ! Output: [real(r8) (:)   ]
         grain_curmr                  => veg_cf%grain_curmr                     , & ! Output: [real(r8) (:)   ]
         leaf_xsmr                    => veg_cf%leaf_xsmr                       , & ! Output: [real(r8) (:)   ]
         froot_xsmr                   => veg_cf%froot_xsmr                      , & ! Output: [real(r8) (:)   ]
         livestem_xsmr                => veg_cf%livestem_xsmr                   , & ! Output: [real(r8) (:)   ]
         livecroot_xsmr               => veg_cf%livecroot_xsmr                  , & ! Output: [real(r8) (:)   ]
         grain_xsmr                   => veg_cf%grain_xsmr                      , & ! Output: [real(r8) (:)   ]
         cpool_to_xsmrpool            => veg_cf%cpool_to_xsmrpool               , & ! Output: [real(r8) (:)   ]

         retransn                     => veg_ns%retransn                     , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant pool of retranslocated N

         plant_ndemand                => veg_nf%plant_ndemand                 , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         avail_retransn               => veg_nf%avail_retransn                , & ! Output: [real(r8) (:)   ]  N flux available from retranslocation pool (gN/m2/s)
         retransn_to_npool            => veg_nf%retransn_to_npool             , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
         leafn_to_retransn            => veg_nf%leafn_to_retransn             , & ! Output: [real(r8) (:)   ]
         frootn_to_retransn           => veg_nf%frootn_to_retransn            , & ! Output: [real(r8) (:)   ]
         livestemn_to_retransn        => veg_nf%livestemn_to_retransn         , & ! Output: [real(r8) (:)   ]

         !!! add phosphorus variables  - X. YANG
         retransp                     => veg_ps%retransp                   , & ! Input:  [real(r8) (:)   ]  (gP/m2) plant pool of retranslocated P

         plant_pdemand                => veg_pf%plant_pdemand               , & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)
         avail_retransp               => veg_pf%avail_retransp              , & ! Output: [real(r8) (:)   ]  P flux available from retranslocation pool (gP/m2/s)
         retransp_to_ppool            => veg_pf%retransp_to_ppool           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated P (gP/m2/s)
         p_allometry                  => cnstate_vars%p_allometry_patch        , & ! Output: [real(r8) (:)   ]  P allocation index (DIM)
         tempmax_retransp             => cnstate_vars%tempmax_retransp_patch   , & ! Output: [real(r8) (:)   ]  temporary annual max of retranslocated P pool (gP/m2)
         annmax_retransp              => cnstate_vars%annmax_retransp_patch    , & ! Output: [real(r8) (:)   ]  annual max of retranslocated P pool
         km_plant_p                   => veg_vp%km_plant_p                     , &
         nfixation_prof               => cnstate_vars%nfixation_prof_col , &
         smin_no3_vr                  => col_ns%smin_no3_vr              , & ! Input: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3
         smin_nh4_vr                  => col_ns%smin_nh4_vr              , & ! Input: [real(r8) (:,:) ]  (gN/m3) soil mineral NH4
         solutionp_vr                 => col_ps%solutionp_vr               , & ! Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral P
         benefit_pgpp_pleafc          => veg_ns%benefit_pgpp_pleafc     &
         )

      !$acc parallel loop independent gang vector default(present) 
      do fp = 1, num_soilp 
         p = filter_soilp(fp) 
         ivt = veg_pp%itype(p)
        ! The input psn (psnsun and psnsha) are expressed per unit LAI
        ! in the sunlit and shaded canopy, respectively. These need to be
        ! scaled by laisun and laisha to get the total gpp for allocation

        ! Note that no associate statement is used for the isotope carbon fluxes below
        ! since they are not always allocated AND nag compiler will complain if you try to
        ! to have an associate statement with unallocated memory

        psnsun_to_cpool(p)   = psnsun(p) * laisun(p) * 12.011e-6_r8
        psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6_r8

        if ( use_c13 ) then
           c13_veg_cf%psnsun_to_cpool(p)   = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
           c13_veg_cf%psnshade_to_cpool(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
        endif

        if ( use_c14 ) then
           c14_veg_cf%psnsun_to_cpool(p)   = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
           c14_veg_cf%psnshade_to_cpool(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
        endif

        gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

        ! carbon return of leaf C investment
        benefit_pgpp_pleafc(p) = max(gpp(p)  / max(leafc(p),1e-20_r8), 0.0_r8)

        ! get the time step total maintenance respiration
        ! These fluxes should already be in gC/m2/s

        mr = leaf_mr(p) + froot_mr(p)
        if (woody(ivt) == 1.0_r8) then
           mr = mr + livestem_mr(p) + livecroot_mr(p)
        else if (ivt >= npcropmin) then
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

        f1 = froot_leaf(ivt)
        f2 = croot_stem(ivt)

        ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
        ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
        ! This variable allocation is only for trees. Shrubs have a constant
        ! allocation as specified in the pft-physiology file.  The value is also used
        ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

        if (stem_leaf(ivt) < 0._r8) then
           if (stem_leaf(ivt) == -1._r8) then
                f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
           else
                f3 = max((-1.0_r8*stem_leaf(ivt)*2.7_r8)/(1.0_r8+exp(-0.004_r8*(annsum_npp(p) - &
                          300.0_r8))) - 0.4_r8, 0.2_r8)
           end if
        else
           f3 = stem_leaf(ivt)
        end if

        f4   = flivewd(ivt)
        g1   = grperc(ivt)
        g2   = grpnow(ivt)
        cnl  = leafcn(ivt)
        cnfr = frootcn(ivt)
        cnlw = livewdcn(ivt)
        cndw = deadwdcn(ivt)

        cpl = leafcp(ivt)
        cpfr = frootcp(ivt)
        cplw = livewdcp(ivt)
        cpdw = deadwdcp(ivt)


        ! calculate f1 to f5 for prog crops following AgroIBIS subr phenocrop

        f5 = 0._r8 ! continued intializations from above

        if (ivt >= npcropmin) then ! skip 2 generic crops

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
                    aroot(p) = max(0._r8, min(1._r8, arooti(ivt) -   &
                        (arooti(ivt) - arootf(ivt)) *  &
                        min(1._r8, hui(p)/gddmaturity(p))))
                    astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
                 else
                    arepr(p) = 0._r8
                    aroot(p) = max(0._r8, min(1._r8, arooti(ivt) -   &
                        (arooti(ivt) - arootf(ivt)) *  &
                        min(1._r8, hui(p)/gddmaturity(p))))
                    fleaf = fleafi(ivt) * (exp(-bfact(ivt)) -         &
                        exp(-bfact(ivt)*hui(p)/huigrain(p))) / &
                        (exp(-bfact(ivt))-1) ! fraction alloc to leaf (from J Norman alloc curve)
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

                 aroot(p) = max(0._r8, min(1._r8, arooti(ivt) - &
                      (arooti(ivt) - arootf(ivt)) * min(1._r8, hui(p)/gddmaturity(p))))
                 if (astemi(p) > astemf(ivt)) then
                    astem(p) = max(0._r8, max(astemf(ivt), astem(p) * &
                        (1._r8 - min((hui(p)-                 &
                        huigrain(p))/((gddmaturity(p)*declfact(ivt))- &
                        huigrain(p)),1._r8)**allconss(ivt) )))
                 end if
                 if (aleafi(p) > aleaff(ivt)) then
                    aleaf(p) = max(1.e-5_r8, max(aleaff(ivt), aleaf(p) * &
                        (1._r8 - min((hui(p)-                    &
                        huigrain(p))/((gddmaturity(p)*declfact(ivt))- &
                        huigrain(p)),1._r8)**allconsl(ivt) )))
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

                 if (ivt /= nsoybean .or. astem(p) == astemf(ivt) .or. peaklai(p) == 1._r8) then
                    if (grain_flag(p) == 0._r8) then
                       t1 = 1 / dt
                       leafn_to_retransn(p) = t1 * ((leafc(p) / leafcn(ivt)) - (leafc(p) / &
                           fleafcn(ivt)))
                       livestemn_to_retransn(p) = t1 * ((livestemc(p) / livewdcn(ivt)) - (livestemc(p) / &
                           fstemcn(ivt)))
                       frootn_to_retransn(p) = 0._r8
                       if (ffrootcn(ivt) > 0._r8) then
                          frootn_to_retransn(p) = t1 * ((frootc(p) / frootcn(ivt)) - (frootc(p) / &
                              ffrootcn(ivt)))
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

        if (woody(ivt) == 1.0_r8) then
           c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
           n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
                (f3*(1._r8-f4)*(1._r8+f2))/cndw
           p_allometry(p) = 1._r8/cpl + f1/cpfr + (f3*f4*(1._r8+f2))/cplw + &
                (f3*(1._r8-f4)*(1._r8+f2))/cpdw

        else if (ivt >= npcropmin) then ! skip generic crops
           cng = graincn(ivt)
           cpg = graincp(ivt)
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

        if (ivt >= npcropmin .and. grain_flag(p) == 1._r8) then
           avail_retransn(p) = plant_ndemand(p)
           avail_retransp(p) = plant_pdemand(p)
        else if (ivt < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
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
      end do 

   end associate

 end subroutine TotalNPDemand

!-------------------------------------------------------------------------------------------------

 subroutine Allocation2_ResolveNPLimit (bounds,num_soilc, filter_soilc, &
                                        cnstate_vars , soilstate_vars, dt )
   ! PHASE-2 of Allocation:  resolving N/P limitation
   ! !USES:
   use elm_varctl      , only : carbon_only          !
   use elm_varctl      , only : carbonnitrogen_only  !
   use elm_varctl      , only : carbonphosphorus_only!
   use pftvarcon       , only : noveg
   use elm_varpar      , only : nlevdecomp, ndecomp_cascade_transitions
   use elm_varcon      , only : nitrif_n2o_loss_frac, secspday
   use elm_varcon      , only : zisoi
   !
   ! !ARGUMENTS:
   type(bounds_type)        , intent(in)    :: bounds
   integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
   integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
   type(cnstate_type)       , intent(inout) :: cnstate_vars
   type(soilstate_type)     , intent(in)    :: soilstate_vars
   real(r8)  ,  intent(in)  :: dt
   !
   ! !LOCAL VARIABLES:
   real(r8) :: fpi_no3_vr(num_soilc,1:nlevdecomp) ! fraction of potential immobilization supplied by no3(no units)
   real(r8) :: fpi_nh4_vr(num_soilc,1:nlevdecomp) ! fraction of potential immobilization supplied by nh4 (no units)

   integer :: c,p,l,j,k ! indices
   integer :: fp        ! lake filter pft index
   integer :: fc        ! lake filter column index

   ! Fractional uptake profiles, that are proportional to root density
   real(r8):: nuptake_prof(num_soilc,1:nlevdecomp)
   real(r8):: puptake_prof(num_soilc,1:nlevdecomp)

   real(r8) :: sum1,sum2,sum3,sum4,sum5,sum6
   real  :: startt, stopt
   integer :: begc, endc 
   !-----------------------------------------------------------------------
   associate(        &
        ivt                          => veg_pp%itype                    , & ! Input:  [integer  (:) ]  pft vegetation type
                                                                            ! new variables due to partition of Allocation to 3 subroutines: BEG
        plant_ndemand_col            => col_nf%plant_ndemand            , & ! Output:  [real(r8) (:,:) ]
        plant_pdemand_col            => col_pf%plant_pdemand            , & ! Output:  [real(r8) (:,:) ]
                                                                            ! new variables due to partition of Allocation to 3 subroutines: END
        fpg                          => cnstate_vars%fpg_col            , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
        fpi                          => cnstate_vars%fpi_col            , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
        fpi_vr                       => cnstate_vars%fpi_vr_col         , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
        fpg_p                        => cnstate_vars%fpg_p_col          , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
        fpi_p                        => cnstate_vars%fpi_p_col          , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
        fpi_p_vr                     => cnstate_vars%fpi_p_vr_col       , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
        fpg_nh4_vr                   => cnstate_vars%fpg_nh4_vr_col     , &
        fpg_no3_vr                   => cnstate_vars%fpg_no3_vr_col     , &
        fpg_vr                       => cnstate_vars%fpg_vr_col         , &
        fpg_p_vr                     => cnstate_vars%fpg_p_vr_col       , &
        nfixation_prof               => cnstate_vars%nfixation_prof_col , &
        sminn_vr                     => col_ns%sminn_vr                 , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
        smin_nh4_vr                  => col_ns%smin_nh4_vr              , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NH4
        smin_no3_vr                  => col_ns%smin_no3_vr              , & ! Output: [real(r8) (:,:) ]  (gN/m3) soil mineral NO3
        potential_immob              => col_nf%potential_immob          , & ! Output: [real(r8) (:)   ]
        actual_immob                 => col_nf%actual_immob             , & ! Output: [real(r8) (:)   ]
        sminn_to_plant               => col_nf%sminn_to_plant           , & ! Output: [real(r8) (:)   ]
        pot_f_nit_vr                 => col_nf%pot_f_nit_vr             , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
        pot_f_denit_vr               => col_nf%pot_f_denit_vr              , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
        f_nit_vr                     => col_nf%f_nit_vr                    , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
        f_denit_vr                   => col_nf%f_denit_vr                  , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
        actual_immob_no3_vr          => col_nf%actual_immob_no3_vr         , & ! Output: [real(r8) (:,:) ]
        actual_immob_nh4_vr          => col_nf%actual_immob_nh4_vr         , & ! Output: [real(r8) (:,:) ]
        smin_no3_to_plant_vr         => col_nf%smin_no3_to_plant_vr        , & ! Output: [real(r8) (:,:) ]
        smin_nh4_to_plant_vr         => col_nf%smin_nh4_to_plant_vr        , & ! Output: [real(r8) (:,:) ]
        n2_n2o_ratio_denit_vr        => col_nf%n2_n2o_ratio_denit_vr       , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
        f_n2o_denit_vr               => col_nf%f_n2o_denit_vr              , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
        f_n2o_nit_vr                 => col_nf%f_n2o_nit_vr                , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
        supplement_to_sminn_vr       => col_nf%supplement_to_sminn_vr      , & ! Output: [real(r8) (:,:) ]
        sminn_to_plant_vr            => col_nf%sminn_to_plant_vr           , & ! Output: [real(r8) (:,:) ]
        potential_immob_vr           => col_nf%potential_immob_vr          , & ! Output: [real(r8) (:,:) ]
        actual_immob_vr              => col_nf%actual_immob_vr             , & ! Output: [real(r8) (:,:) ]
        col_plant_ndemand_vr         => col_nf%col_plant_ndemand_vr        , &
        col_plant_nh4demand_vr       => col_nf%col_plant_nh4demand_vr      , &
        col_plant_no3demand_vr       => col_nf%col_plant_no3demand_vr      , &
        col_plant_pdemand_vr         => col_pf%col_plant_pdemand_vr        , &
        cn_scalar                    => cnstate_vars%cn_scalar             , &
        cp_scalar                    => cnstate_vars%cp_scalar             , &
        t_scalar                     => col_cf%t_scalar                    , &
        plant_nh4demand_vr_patch     => veg_nf%plant_nh4demand_vr          , &
        plant_no3demand_vr_patch     => veg_nf%plant_no3demand_vr          , &
        plant_ndemand_vr_patch       => veg_nf%plant_ndemand_vr            , &
        plant_pdemand_vr_patch       => veg_pf%plant_pdemand_vr            , &
        pnup_pfrootc                 => veg_ns%pnup_pfrootc                , &
        isoilorder                   => cnstate_vars%isoilorder            , &
        sminp_to_plant_patch         => veg_pf%sminp_to_plant              , &
        sminn_to_plant_patch         => veg_nf%sminn_to_plant              , &
        smin_nh4_to_plant_patch      => veg_nf%smin_nh4_to_plant           , &
        smin_no3_to_plant_patch      => veg_nf%smin_no3_to_plant           , &
        actual_immob_no3             => col_nf%actual_immob_no3            , &
        actual_immob_nh4             => col_nf%actual_immob_nh4            , &
        froot_prof                   => cnstate_vars%froot_prof_patch      , & ! fine root vertical profile Zeng, X. 2001. Global vegetation root distribution for land modeling. J. Hydrometeor. 2:525-530
        frootc                       => veg_cs%frootc                      , & ! Input:  [real(r8) (:)   ]
        leafc                        => veg_cs%leafc                       , & ! Input:  [real(r8) (:)   ]
        leafcn                       => veg_vp%leafcn                      , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
        leafcp                       => veg_vp%leafcp                      , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
        vmax_minsurf_p_vr            => veg_vp%vmax_minsurf_p_vr           , &
        leafn                        => veg_ns%leafn                       , &
        vmax_plant_nh4               => veg_vp%vmax_plant_nh4              , &
        vmax_plant_no3               => veg_vp%vmax_plant_no3              , &
        vmax_plant_p                 => veg_vp%vmax_plant_p                , &
        primp_to_labilep_vr_col      => col_pf%primp_to_labilep_vr         , &
        leafp                        => veg_ps%leafp                      , &
        km_decomp_nh4                => veg_vp%km_decomp_nh4              , &
        km_decomp_no3                => veg_vp%km_decomp_no3              , &
        km_decomp_p                  => veg_vp%km_decomp_p                , &
        km_nit                       => veg_vp%km_nit                     , &
        km_den                       => veg_vp%km_den                     , &
        km_plant_nh4                 => veg_vp%km_plant_nh4               , &
        km_plant_no3                 => veg_vp%km_plant_no3               , &
        km_plant_p                   => veg_vp%km_plant_p                 , &
        km_minsurf_p_vr              => veg_vp%km_minsurf_p_vr            , &
        decompmicc_patch_vr          => veg_vp%decompmicc_patch_vr        , &
        adsorb_to_labilep_vr         => col_pf%adsorb_to_labilep_vr       , &
        desorb_to_solutionp_vr       => col_pf%desorb_to_solutionp_vr     , &
        biochem_pmin_vr_col          => col_pf%biochem_pmin_vr            , &
        labilep_vr                   => col_ps%labilep_vr                 , &
        pdep_to_sminp                => col_pf%pdep_to_sminp              , &
        pdep_prof                    => cnstate_vars%pdep_prof_col        , &
        gross_pmin_vr                => col_pf%gross_pmin_vr              , &
                                !! add phosphorus variables  - X. YANG
        solutionp_vr                 => col_ps%solutionp_vr               , & ! Input:  [real(r8) (:,:) ]  (gP/m3) soil mineral P
        potential_immob_p            => col_pf%potential_immob_p           , & ! Output: [real(r8) (:)   ]
        actual_immob_p               => col_pf%actual_immob_p              , & ! Output: [real(r8) (:)   ]
        sminp_to_plant               => col_pf%sminp_to_plant              , & ! Output: [real(r8) (:)   ]
        supplement_to_sminp_vr       => col_pf%supplement_to_sminp_vr      , & ! Output: [real(r8) (:,:) ]
        sminp_to_plant_vr            => col_pf%sminp_to_plant_vr           , & ! Output: [real(r8) (:,:) ]
        potential_immob_p_vr         => col_pf%potential_immob_p_vr        , & ! Output: [real(r8) (:,:) ]
        actual_immob_p_vr            => col_pf%actual_immob_p_vr           , & ! Output: [real(r8) (:,:) ]
        bd                           => soilstate_vars%bd_col              , &
        h2osoi_vol                   => col_ws%h2osoi_vol                  , &
        pmnf_decomp_cascade          => col_nf%pmnf_decomp_cascade         , &
        pmpf_decomp_cascade          => col_pf%pmpf_decomp_cascade         , &
        leafc_storage                => veg_cs%leafc_storage               , &
        leafc_xfer                   => veg_cs%leafc_xfer                  , &
        leafn_storage                => veg_ns%leafn_storage              , &
        leafn_xfer                   => veg_ns%leafn_xfer                 , &
        leafp_storage                => veg_ps%leafp_storage              , &
        leafp_xfer                   => veg_ps%leafp_xfer                 &
        )

     !if (nu_com .eq. 'RD') then ! 'RD' : relative demand approach
      ! Starting resolving N/P limitation
      ! calculate nuptake & puptake profile
      begc = bounds%begc 
      endc = bounds%endc 
      !$acc enter data create(nuptake_prof(1:num_soilc,1:nlevdecomp),puptake_prof(1:num_soilc,1:nlevdecomp),&
      !$acc fpi_no3_vr(:,1:nlevdecomp),fpi_nh4_vr(:,1:nlevdecomp))

      call calc_nuptake_prof(num_soilc, filter_soilc, cnstate_vars, nuptake_prof)
      call calc_puptake_prof(num_soilc, filter_soilc, cnstate_vars, puptake_prof)

   ! ------------------------------------------------------------------------------
   ! PART I.
   ! Determine the boundary conditions for the competitive allocation modules
   ! This is mostly about pointing to either the big-leaf or FATES boundary
   ! conditions.
   ! ------------------------------------------------------------------------------
      call cpu_time(startt)
      !$acc parallel loop independent gang default(present)
      do j = 1, nlevdecomp
         !$acc loop vector independent private(c)
         do fc=1,num_soilc
           c = filter_soilc(fc)
           col_plant_ndemand_vr(c,j) = plant_ndemand_col(c) * nuptake_prof(fc,j)
           col_plant_pdemand_vr(c,j) = plant_pdemand_col(c) * puptake_prof(fc,j)
        end do
      end do
   
   if (nu_com .eq. 'RD') then
   ! Estimate actual allocation rates via Relative Demand
   ! approach (RD)
      ! Starting resolving N limitation !!!
      ! =============================================================
      ! This section is modified, Aug 2015 by Q. Zhu
      ! (1) add nitrogen and phosphorus competition
      ! (2) nitrogen and phosphorus uptake is based on root kinetics
      ! (3) no second pass nutrient uptake for plants
      ! =============================================================
       call NAllocationRD(num_soilc, filter_soilc, &
           col_plant_ndemand_vr(begc:endc,1:nlevdecomp), & ! IN
           potential_immob_vr(begc:endc,1:nlevdecomp),   & ! IN
           AllocParamsInst%compet_plant_nh4,        & ! IN
           AllocParamsInst%compet_decomp_nh4,       & ! IN
           dt,                               & ! IN
           smin_nh4_vr(begc:endc,1:nlevdecomp),     & ! IN
           fpi_nh4_vr(1:num_soilc,1:nlevdecomp),    & ! OUT
           actual_immob_nh4_vr(begc:endc,1:nlevdecomp), & ! OUT
           smin_nh4_to_plant_vr(begc:endc,1:nlevdecomp),& ! OUT
           smin_no3_vr(begc:endc,1:nlevdecomp),     & ! IN
           AllocParamsInst%compet_plant_no3,        & ! IN
           AllocParamsInst%compet_decomp_no3,       & ! IN
           AllocParamsInst%compet_nit,              & ! IN
           AllocParamsInst%compet_denit,            & ! IN
           pot_f_nit_vr(begc:endc,1:nlevdecomp),    & ! IN
           pot_f_denit_vr(begc:endc,1:nlevdecomp),  & ! IN
           fpi_no3_vr(1:num_soilc,1:nlevdecomp),    & ! OUT
           actual_immob_no3_vr(begc:endc,1:nlevdecomp),& ! OUT
           smin_no3_to_plant_vr(begc:endc,1:nlevdecomp),&! OUT
           f_nit_vr(begc:endc,1:nlevdecomp),           & ! OUT
           f_denit_vr(begc:endc,1:nlevdecomp))           ! OUT
   end if
     !$acc parallel loop independent collapse(2) gang vector default(present)
      do j = 1, nlevdecomp
         do fc=1,num_soilc !col_loop
           c = filter_soilc(fc)
           l = col_pp%landunit(c)
           ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
           f_n2o_nit_vr(c,j) = f_nit_vr(c,j) * nitrif_n2o_loss_frac
           f_n2o_denit_vr(c,j) = f_denit_vr(c,j) / (1._r8 + n2_n2o_ratio_denit_vr(c,j))

           ! eliminate any N limitation if we are supplemnting nitrogen
           ! The amount of supplemental N is the deficit for plant-uptake and immobilization,
           ! for both NO3 and Nh4.  However, this quantity, even though it factors in the
           ! deficit for NO3 to meet allocation needs, is only added to the NH4 pool.
           ! Thus, the NH4 fluxes are increased, for itself and as a surrogate to meet the
           ! NO3 flux demands.

           if (carbon_only .or. carbonphosphorus_only) then

              if ( fpi_no3_vr(fc,j) + fpi_nh4_vr(fc,j) < 1._r8 ) then
                 fpi_vr(c,j) = 1._r8
                 fpi_nh4_vr(fc,j) = 1.0_r8 - fpi_no3_vr(fc,j)
                 supplement_to_sminn_vr(c,j) = (potential_immob_vr(c,j) - actual_immob_no3_vr(c,j)) - actual_immob_nh4_vr(c,j)
                 ! update to new values that satisfy demand
                 actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)
              end if

              !if (nu_com .eq. 'RD') then
               if ( smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j) < col_plant_ndemand_vr(c,j)) then
                   supplement_to_sminn_vr(c,j) = supplement_to_sminn_vr(c,j) + &
                        col_plant_ndemand_vr(c,j) - (smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j))

                   ! update to new values that satisfy demand
                   smin_nh4_to_plant_vr(c,j) = col_plant_ndemand_vr(c,j) - smin_no3_to_plant_vr(c,j)
               end if
              !else

           end if

           ! sum up nitrogen limitation to decomposition
           fpi_vr(c,j) = fpi_nh4_vr(fc,j) + fpi_no3_vr(fc,j)

           ! sum up no3 and nh4 fluxes
           sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
           actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
         end do ! col loop
      end do ! nlevdecomp


        ! Starting resolving P limitation !!!
        ! =============================================================
        if (nu_com .eq. 'RD') then
           ! Relative Demand (RD)
           call PAllocationRD(num_soilc,filter_soilc, &
                col_plant_pdemand_vr(begc:endc,1:nlevdecomp), & ! IN
                potential_immob_p_vr(begc:endc,1:nlevdecomp),               & ! IN
                solutionp_vr(begc:endc,1:nlevdecomp),                       & ! IN
                dt,                                      & ! IN
                fpi_p_vr(begc:endc,1:nlevdecomp),                           & ! OUT
                actual_immob_p_vr(begc:endc,1:nlevdecomp),                  & ! OUT
                sminp_to_plant_vr(begc:endc,1:nlevdecomp),                  & ! OUT
                supplement_to_sminp_vr(begc:endc,1:nlevdecomp))               ! OUT
        end if ! end of P competition

         !  resolving N limitation vs. P limitation for decomposition
         !  update (1) actual immobilization for N and P (2) sminn_to_plant and sminp_to_plant
         !  We only resolve co-limitations when are supplementing neither element

     np_bothactive: if ( .not.carbon_only .and.  &
                         .not.carbonphosphorus_only .and. &
                         .not.carbonnitrogen_only ) then

            !$acc enter data create(sum1,sum2,sum3,sum4,sum5,sum6)

           !if (nu_com .eq. 'RD') then
           !$acc parallel loop independent collapse(2) gang worker private(c,sum1,sum2) default(present)
           do fc = 1, num_soilc
              do j = 1, nlevdecomp
                 c= filter_soilc(fc)
                 sum1 = 0._r8
                 sum2 = 0._r8
                 if( fpi_vr(c,j) <= fpi_p_vr(c,j) )then ! more N limited
                    !$acc loop vector reduction(+:sum1)
                    do k = 1, ndecomp_cascade_transitions
                       if (pmnf_decomp_cascade(c,j,k) > 0.0_r8 .and. pmpf_decomp_cascade(c,j,k) > 0.0_r8) then
                           sum1 = sum1 - pmpf_decomp_cascade(c,j,k) *(fpi_p_vr(c,j)-fpi_vr(c,j))
                       end if
                    end do
                    actual_immob_p_vr(c,j) = actual_immob_p_vr(c,j) + sum1
                 else
                    if (fpi_nh4_vr(fc,j) > fpi_p_vr(fc,j)) then ! more P limited
                       !$acc loop vector reduction(+:sum1,sum2)
                       do k = 1, ndecomp_cascade_transitions
                          if (pmnf_decomp_cascade(c,j,k) > 0.0_r8 .and. pmpf_decomp_cascade(c,j,k) > 0.0_r8) then
                              sum1 = sum1 - pmnf_decomp_cascade(c,j,k) * (fpi_nh4_vr(fc,j) - fpi_p_vr(c,j))
                              sum2 = sum2 - pmnf_decomp_cascade(c,j,k) * fpi_no3_vr(fc,j)
                          end if
                       end do
                       actual_immob_nh4_vr(c,j) = actual_immob_nh4_vr(c,j)+ sum1
                       actual_immob_no3_vr(c,j) = actual_immob_no3_vr(c,j)+ sum2
                    else
                       !$acc loop vector reduction(+:sum1)
                       do k = 1, ndecomp_cascade_transitions
                          if (pmnf_decomp_cascade(c,j,k) > 0.0_r8 .and. pmpf_decomp_cascade(c,j,k) > 0.0_r8) then
                             sum1  = actual_immob_no3_vr(c,j) - pmnf_decomp_cascade(c,j,k) * &
                                  (fpi_nh4_vr(fc,j) + fpi_no3_vr(fc,j) - fpi_p_vr(c,j) )
                          end if
                       end do
                       actual_immob_no3_vr(c,j) = actual_immob_no3_vr(c,j) + sum1
                    end if
                 endif
                 ! sum up no3 and nh4 fluxes
                 actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
              end do
           end do ! end col loops
           !else
           !end if
        endif np_bothactive

      if(carbonnitrogen_only)then
         !$acc parallel loop independent collapse(2) gang vector private(c) default(present)
         do fc = 1, num_soilc
           do j = 1, nlevdecomp
             c = filter_soilc(fc)
             actual_immob_p_vr(c,j) = potential_immob_p_vr(c,j) * fpi_vr(c,j)
           end do
        end do
      end if

      if(carbonphosphorus_only)then
         !$acc parallel loop independent collapse(2) gang vector private(c) default(present)
         do fc = 1, num_soilc
           do j = 1, nlevdecomp
             c = filter_soilc(fc)
              actual_immob_vr(c,j) = potential_immob_vr(c,j) * fpi_p_vr(c,j)
           end do
        end do
      end if

      ! sum up plant N/P uptake at column level and patch level
      ! sum up N fluxes to plant after initial competition
      !$acc parallel loop independent gang worker private(c,sum1,sum2) default(present)
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        sum1 = 0._r8
        sum2 = 0._r8
        !$acc loop vector reduction(+:sum1,sum2)
        do j = 1, nlevdecomp
           sum1 = sum1  + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
           sum2 = sum2 + sminp_to_plant_vr(c,j) * dzsoi_decomp(j)
        end do
        sminn_to_plant(c) = sum1
        sminp_to_plant(c) = sum2
     end do !col_loop

     ! sum up N fluxes to immobilization
     !$acc parallel loop independent gang worker private(c,sum1,sum2,sum3,sum4,sum5,sum6) default(present)
     do fc=1,num_soilc
        c = filter_soilc(fc)
        sum1 = 0._r8; sum2 = 0._r8;
        sum3 = 0._r8;sum4 = 0._r8;
        sum5 = 0._r8;sum6 = 0._r8;
        !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5,sum6)
        do j = 1, nlevdecomp
           sum1 = sum1 + actual_immob_vr(c,j) * dzsoi_decomp(j)
           sum2 = sum2 + potential_immob_vr(c,j) * dzsoi_decomp(j)
           sum3 = sum3 + actual_immob_no3_vr(c,j) * dzsoi_decomp(j)
           sum4 = sum4 + actual_immob_nh4_vr(c,j) * dzsoi_decomp(j)
           sum5 = sum5 + actual_immob_p_vr(c,j) * dzsoi_decomp(j)
           sum6 = sum6 + potential_immob_p_vr(c,j) * dzsoi_decomp(j)
        end do
        actual_immob(c)     = sum1
        potential_immob(c)  = sum2
        actual_immob_no3(c) = sum3
        actual_immob_nh4(c) = sum4
        actual_immob_p(c)   = sum5
        potential_immob_p(c)= sum6
     end do

     !$acc parallel loop independent gang vector private(c) default(present)
     do fc=1,num_soilc
        c = filter_soilc(fc)
        ! calculate the fraction of potential growth that can be
        ! acheived with the N available to plants
        if (plant_ndemand_col(c) > 0.0_r8) then
           fpg(c) = sminn_to_plant(c) / plant_ndemand_col(c)
        else
           fpg(c) = 1._r8
        end if

        ! calculate the fraction of immobilization realized (for diagnostic purposes)
        if (potential_immob(c) > 0.0_r8) then
           fpi(c) = actual_immob(c) / potential_immob(c)
        else
           fpi(c) = 1._r8
        end if
     end do ! end of column loops

     !$acc parallel loop independent gang vector private(c) default(present)
     do fc=1,num_soilc
        c = filter_soilc(fc)
        ! calculate the fraction of potential growth that can be
        ! acheived with the P available to plants
        if (plant_pdemand_col(c) > 0.0_r8) then
           fpg_p(c) = sminp_to_plant(c) / plant_pdemand_col(c)
        else
           fpg_p(c) = 1.0_r8
        end if
        ! calculate the fraction of immobilization realized (for diagnostic purposes)
        if (potential_immob_p(c) > 0.0_r8) then
           fpi_p(c) = actual_immob_p(c) / potential_immob_p(c)
        else
           fpi_p(c) = 1.0_r8
        end if
     end do
     !$acc exit data delete(nuptake_prof(:,:),puptake_prof(:,:),fpi_no3_vr(:,:),fpi_nh4_vr(:,:),sum1,sum2,&
     !$acc  sum3, sum4, sum5, sum6)
    end associate
 end subroutine Allocation2_ResolveNPLimit

!-------------------------------------------------------------------------------------------------
  subroutine Allocation3_PlantCNPAlloc ( &
        num_soilc, filter_soilc, num_soilp, filter_soilp    , &
        canopystate_vars                                    , &
        cnstate_vars, crop_vars , dt )
    ! PHASE-3 of Allocation: start new pft loop to distribute the available N/P between the
    ! competing patches on the basis of relative demand, and allocate C/N/P to new growth and storage

    ! !USES:
    use elm_varctl  , only : carbon_only , carbonnitrogen_only ,carbonphosphorus_only!
    use pftvarcon   , only : noveg
    use pftvarcon   , only : npcropmin, grperc, grpnow
    use elm_varpar  , only : nlevdecomp
    use elm_varcon  , only : nitrif_n2o_loss_frac, secspday
    !
    ! !ARGUMENTS:
    integer           , intent(in)    :: num_soilc        ! number of soil columns in filter
    integer           , intent(in)    :: filter_soilc(:)  ! filter for soil columns
    integer           , intent(in)    :: num_soilp        ! number of soil patches in filter
    integer           , intent(in)    :: filter_soilp(:)  ! filter for soil patches

    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars

    type(crop_type)          , intent(inout) :: crop_vars
    real(r8)                  , intent(in)   :: dt
    !
    ! !LOCAL VARIABLES:
    !
    integer :: c,p,j                  !indices
    integer :: fp                     !lake filter pft index
    integer :: fc                     !lake filter column index

    real(r8):: temp_sminn_to_plant(num_soilc)
    real(r8):: temp_sminp_to_plant(num_soilc)
    real(r8) :: sum1,sum2
    !-----------------------------------------------------------------------

    associate(      &
         fpg                          => cnstate_vars%fpg_col         , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpg_p                        => cnstate_vars%fpg_p_col       , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         smin_no3_to_plant_vr         => col_nf%smin_no3_to_plant_vr  , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_vr         => col_nf%smin_nh4_to_plant_vr  , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_patch      => veg_nf%smin_nh4_to_plant     , &
         smin_no3_to_plant_patch      => veg_nf%smin_no3_to_plant     , &
         sminp_to_plant_patch         => veg_pf%sminp_to_plant        , &
         sminn_to_plant_patch         => veg_nf%sminn_to_plant        , &
         sminn_to_npool               => veg_nf%sminn_to_npool        , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)
         sminp_to_ppool               => veg_pf%sminp_to_ppool        , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)

         sminn_to_plant               => col_nf%sminn_to_plant       , & ! Output: [real(r8) (:)   ]
         sminp_to_plant               => col_pf%sminp_to_plant       , & ! Output: [real(r8) (:)   ]

         sminn_to_plant_vr            => col_nf%sminn_to_plant_vr     , & ! Output: [real(r8) (:,:) ]
         sminp_to_plant_vr            => col_pf%sminp_to_plant_vr   , & ! Output: [real(r8) (:,:) ]
         plant_ndemand                => veg_nf%plant_ndemand       , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
         plant_pdemand                => veg_pf%plant_pdemand       , & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)

         plant_n_uptake_flux    => col_nf%plant_n_uptake_flux , &
         plant_p_uptake_flux    => col_pf%plant_p_uptake_flux &
         )

      !-------------------------------------------------------------------
      ! set time steps
      !$acc enter data create(sum1,sum2) 
      !if (nu_com .eq. 'RD') then
      !$acc parallel loop  independent gang worker private(c,sum1,sum2) default(present)
      do fc=1,num_soilc
         sum1=0.0_r8;sum2=0.0_r8;
         c = filter_soilc(fc)
         !$acc loop vector reduction(+:sum1,sum2)
         do p = col_pp%pfti(c), col_pp%pftf(c)
            if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
               sum1 = sum1 + plant_ndemand(p) * fpg(c)*veg_pp%wtcol(p)
                sum2= sum2 + plant_pdemand(p) * fpg_p(c)*veg_pp%wtcol(p)
            end if

         end do
         plant_n_uptake_flux(c) = sum1
         plant_p_uptake_flux(c) = sum2
      end do


      ! start new pft loop to distribute the available N between the
      ! competing patches on the basis of relative demand, and allocate C and N to
      ! new growth and storage

      call DistributeN_RD(num_soilp,filter_soilp,cnstate_vars,crop_vars)

      !----------------------------------------------------------------
      ! now use the p2c routine to update column level soil mineral N and P uptake
      ! based on competition between N and P limitation       - XYANG
      !! Nitrogen
      !if (nu_com .eq. 'RD') then
        !! Phosphorus
        if( .not.carbonphosphorus_only .and. .not.carbonnitrogen_only .and. &
             .not. carbon_only )then

             !$acc enter data create(temp_sminn_to_plant(:), temp_sminp_to_plant(:))
             !$acc parallel loop gang worker vector private(c)
             do fc = 1, num_soilc
               c = filter_soilc(fc)
               temp_sminn_to_plant(fc) = sminn_to_plant(c)
               temp_sminp_to_plant(fc) = sminp_to_plant(c)
            end do

            call p2c_1d_filter_parallel(num_soilc,filter_soilc, &
                sminn_to_npool, sminn_to_plant)

            call p2c_1d_filter_parallel(num_soilc,filter_soilc, &
                sminp_to_ppool, sminp_to_plant )


            !$acc parallel loop gang independent default(present)
            do j = 1, nlevdecomp
               !$acc loop worker vector private(c)
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  if ( temp_sminn_to_plant(fc) > 0._r8) then
                     sminn_to_plant_vr(c,j)    = sminn_to_plant_vr(c,j) * ( sminn_to_plant(c)/temp_sminn_to_plant(fc) )
                     smin_nh4_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) * ( sminn_to_plant(c)/temp_sminn_to_plant(fc) )
                     smin_no3_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) * ( sminn_to_plant(c)/temp_sminn_to_plant(fc) )
                  else
                     sminn_to_plant_vr(c,j)    = 0._r8
                     smin_nh4_to_plant_vr(c,j) = 0._r8
                     smin_no3_to_plant_vr(c,j) = 0._r8
                  endif

                  if ( temp_sminp_to_plant(fc) > 0._r8) then
                     sminp_to_plant_vr(c,j) =  sminp_to_plant_vr(c,j) * ( sminp_to_plant(c)/temp_sminp_to_plant(fc) )
                  else
                     sminp_to_plant_vr(c,j) = 0._r8
                  endif
               end do
            end do
            !$acc exit data delete(temp_sminp_to_plant(:), temp_sminn_to_plant(:))

          end if   ! carbonnitrogenphosphorus

          if(  carbonnitrogen_only  )then

             !$acc enter data create(temp_sminp_to_plant(:))
             !$acc parallel loop independent gang worker vector private(c)
             do fc = 1, num_soilc
                c = filter_soilc(fc)
                temp_sminp_to_plant(fc) = sminp_to_plant(c)
            end do

            call p2c_1d_filter_parallel(num_soilc,filter_soilc, &
                sminp_to_ppool , sminp_to_plant )

            !$acc parallel loop gang independent default(present)
            do j = 1, nlevdecomp
               !$acc loop worker vector private(c)
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  if ( temp_sminp_to_plant(fc) > 0._r8) then
                     sminp_to_plant_vr(c,j) =  sminp_to_plant_vr(c,j) * ( sminp_to_plant(c)/temp_sminp_to_plant(fc) )
                  else
                     sminp_to_plant_vr(c,j) = 0._r8
                  endif
               end do
            end do

            !$acc exit data delete(temp_sminp_to_plant(:))
          end if  ! carbonnitrogen

      !end if ! nu_com .eq. RD
      !----------------------------------------------------------------
      !$acc exit data delete(sum1,sum2) 
    end associate

  end subroutine Allocation3_PlantCNPAlloc

  subroutine DistributeN_RD(num_soilp,filter_soilp,cnstate_vars,crop_vars)
     ! Routine called in Allocation Phase 3
     use pftvarcon   , only : npcropmin, grperc, grpnow
     use elm_varctl  , only : carbon_only , carbonnitrogen_only ,carbonphosphorus_only!
     use pftvarcon   , only : noveg

     !
     integer :: num_soilp 
     integer :: filter_soilp(:)
     type(cnstate_type), intent(inout) :: cnstate_vars
     type(crop_type), intent(in) :: crop_vars
     real(r8):: f1,f2,f3,f4,f5,g1,g2   !allocation parameters
     real(r8):: cnl,cnfr,cnlw,cndw     !C:N ratios for leaf, fine root, and wood
     real(r8):: fcur                   !fraction of current psn displayed as growth
     real(r8):: gresp_storage          !temporary variable for growth resp to storage
     real(r8):: nlc                    !temporary variable for total new leaf carbon allocation
     real(r8):: cng                      !C:N ratio for grain (= cnlw for now; slevis)

     !! Local P variables
     real(r8):: rc, rc_p, r            !Factors for nitrogen pool
     real(r8):: cpl,cpfr,cplw,cpdw,cpg !C:N ratios for leaf, fine root, and wood
     integer :: ivt,fp ,p,c 

     associate(&
        woody                        => veg_vp%woody                         , & ! Input:  [real(r8) (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
        froot_leaf                   => veg_vp%froot_leaf                    , & ! Input:  [real(r8) (:)   ]  allocation parameter: new fine root C per new leaf C (gC/gC)
        croot_stem                   => veg_vp%croot_stem                    , & ! Input:  [real(r8) (:)   ]  allocation parameter: new coarse root C per new stem C (gC/gC)
        stem_leaf                    => veg_vp%stem_leaf                     , & ! Input:  [real(r8) (:)   ]  allocation parameter: new stem c per new leaf C (gC/gC)
        flivewd                      => veg_vp%flivewd                       , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
        leafcn                       => veg_vp%leafcn                        , & ! Input:  [real(r8) (:)   ]  leaf C:N (gC/gN)
        frootcn                      => veg_vp%frootcn                       , & ! Input:  [real(r8) (:)   ]  fine root C:N (gC/gN)
        livewdcn                     => veg_vp%livewdcn                      , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
        deadwdcn                     => veg_vp%deadwdcn                      , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:N (gC/gN)
        fcur2                        => veg_vp%fcur                          , & ! Input:  [real(r8) (:)   ]  allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
        graincn                      => veg_vp%graincn                       , & ! Input:  [real(r8) (:)   ]  grain C:N (gC/gN)
        croplive                     => crop_vars%croplive_patch             , & ! Input:  [logical  (:)   ]  flag, true if planted, not harvested
        aleaf                        => cnstate_vars%aleaf_patch             , & ! Output: [real(r8) (:)   ]  leaf allocation coefficient
        astem                        => cnstate_vars%astem_patch             , & ! Output: [real(r8) (:)   ]  stem allocation coefficient
        fpg                          => cnstate_vars%fpg_col                 , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
        !!! add phosphorus
        leafcp                       => veg_vp%leafcp                        , & ! Input:  [real(r8) (:)   ]  leaf C:P (gC/gP)
        frootcp                      => veg_vp%frootcp                       , & ! Input:  [real(r8) (:)   ]  fine root C:P (gC/gP)
        livewdcp                     => veg_vp%livewdcp                      , & ! Input:  [real(r8) (:)   ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
        deadwdcp                     => veg_vp%deadwdcp                      , & ! Input:  [real(r8) (:)   ]  dead wood (xylem and heartwood) C:P (gC/gP)
        graincp                      => veg_vp%graincp                       , & ! Input:  [real(r8) (:)   ]  grain C:P (gC/gP)
        fpg_p                        => cnstate_vars%fpg_p_col               , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
        c_allometry                  => cnstate_vars%c_allometry_patch       , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
        n_allometry                  => cnstate_vars%n_allometry_patch       , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
        downreg                      => cnstate_vars%downreg_patch           , & ! Output: [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)
        annsum_npp                   => veg_cf%annsum_npp                    , & ! Input:  [real(r8) (:)   ]  annual sum of NPP, for wood allocation
        gpp                          => veg_cf%gpp_before_downreg            , & ! Output: [real(r8) (:)   ]  GPP flux before downregulation (gC/m2/s)
        availc                       => veg_cf%availc                        , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
        excess_cflux                 => veg_cf%excess_cflux                  , & ! Output: [real(r8) (:)   ]  C flux not allocated due to downregulation (gC/m2/s)
        plant_calloc                 => veg_cf%plant_calloc                  , & ! Output: [real(r8) (:)   ]  total allocated C flux (gC/m2/s)
        psnsun_to_cpool              => veg_cf%psnsun_to_cpool               , & ! Output: [real(r8) (:)   ]
        psnshade_to_cpool            => veg_cf%psnshade_to_cpool             , & ! Output: [real(r8) (:)   ]
        cpool_to_leafc               => veg_cf%cpool_to_leafc                , & ! Output: [real(r8) (:)   ]
        cpool_to_leafc_storage       => veg_cf%cpool_to_leafc_storage        , & ! Output: [real(r8) (:)   ]
        cpool_to_frootc              => veg_cf%cpool_to_frootc               , & ! Output: [real(r8) (:)   ]
        cpool_to_frootc_storage      => veg_cf%cpool_to_frootc_storage       , & ! Output: [real(r8) (:)   ]
        cpool_to_livestemc           => veg_cf%cpool_to_livestemc            , & ! Output: [real(r8) (:)   ]
        cpool_to_livestemc_storage   => veg_cf%cpool_to_livestemc_storage    , & ! Output: [real(r8) (:)   ]
        cpool_to_deadstemc           => veg_cf%cpool_to_deadstemc            , & ! Output: [real(r8) (:)   ]
        cpool_to_deadstemc_storage   => veg_cf%cpool_to_deadstemc_storage    , & ! Output: [real(r8) (:)   ]
        cpool_to_livecrootc          => veg_cf%cpool_to_livecrootc           , & ! Output: [real(r8) (:)   ]
        cpool_to_livecrootc_storage  => veg_cf%cpool_to_livecrootc_storage   , & ! Output: [real(r8) (:)   ]
        cpool_to_deadcrootc          => veg_cf%cpool_to_deadcrootc           , & ! Output: [real(r8) (:)   ]
        cpool_to_deadcrootc_storage  => veg_cf%cpool_to_deadcrootc_storage   , & ! Output: [real(r8) (:)   ]
        cpool_to_gresp_storage       => veg_cf%cpool_to_gresp_storage        , & ! Output: [real(r8) (:)   ]  allocation to growth respiration storage (gC/m2/s)
        cpool_to_grainc              => veg_cf%cpool_to_grainc               , & ! Output: [real(r8) (:)   ]  allocation to grain C (gC/m2/s)
        cpool_to_grainc_storage      => veg_cf%cpool_to_grainc_storage       , & ! Output: [real(r8) (:)   ]  allocation to grain C storage (gC/m2/s)
        npool                        => veg_ns%npool                        , & ! Input:  [real(r8) (:)   ]  (gN/m2) plant N pool storage
        plant_ndemand                => veg_nf%plant_ndemand               , & ! Output: [real(r8) (:)   ]  N flux required to support initial GPP (gN/m2/s)
        plant_nalloc                 => veg_nf%plant_nalloc                , & ! Output: [real(r8) (:)   ]  total allocated N flux (gN/m2/s)
        npool_to_grainn              => veg_nf%npool_to_grainn             , & ! Output: [real(r8) (:)   ]  allocation to grain N (gN/m2/s)
        npool_to_grainn_storage      => veg_nf%npool_to_grainn_storage     , & ! Output: [real(r8) (:)   ]  allocation to grain N storage (gN/m2/s)
        retransn_to_npool            => veg_nf%retransn_to_npool           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
        sminn_to_npool               => veg_nf%sminn_to_npool              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral N uptake (gN/m2/s)

        npool_to_leafn               => veg_nf%npool_to_leafn              , & ! Output: [real(r8) (:)   ]  allocation to leaf N (gN/m2/s)
        npool_to_leafn_storage       => veg_nf%npool_to_leafn_storage      , & ! Output: [real(r8) (:)   ]  allocation to leaf N storage (gN/m2/s)
        npool_to_frootn              => veg_nf%npool_to_frootn             , & ! Output: [real(r8) (:)   ]  allocation to fine root N (gN/m2/s)
        npool_to_frootn_storage      => veg_nf%npool_to_frootn_storage     , & ! Output: [real(r8) (:)   ]  allocation to fine root N storage (gN/m2/s)
        npool_to_livestemn           => veg_nf%npool_to_livestemn          , & ! Output: [real(r8) (:)   ]
        npool_to_livestemn_storage   => veg_nf%npool_to_livestemn_storage  , & ! Output: [real(r8) (:)   ]
        npool_to_deadstemn           => veg_nf%npool_to_deadstemn          , & ! Output: [real(r8) (:)   ]
        npool_to_deadstemn_storage   => veg_nf%npool_to_deadstemn_storage  , & ! Output: [real(r8) (:)   ]
        npool_to_livecrootn          => veg_nf%npool_to_livecrootn         , & ! Output: [real(r8) (:)   ]
        npool_to_livecrootn_storage  => veg_nf%npool_to_livecrootn_storage , & ! Output: [real(r8) (:)   ]
        npool_to_deadcrootn          => veg_nf%npool_to_deadcrootn         , & ! Output: [real(r8) (:)   ]
        npool_to_deadcrootn_storage  => veg_nf%npool_to_deadcrootn_storage , & ! Output: [real(r8) (:)   ]
        !!! add phosphorus variables  - X. YANG
        ppool                        => veg_ps%ppool                      , & ! Input: [real(r8)       ] Plant non-structural P storage (gP/m2)
        plant_pdemand                => veg_pf%plant_pdemand               , & ! Output: [real(r8) (:)   ]  P flux required to support initial GPP (gP/m2/s)
        plant_palloc                 => veg_pf%plant_palloc                , & ! Output: [real(r8) (:)   ]  total allocated P flux (gP/m2/s)
        ppool_to_grainp              => veg_pf%ppool_to_grainp             , & ! Output: [real(r8) (:)   ]  allocation to grain P (gP/m2/s)
        ppool_to_grainp_storage      => veg_pf%ppool_to_grainp_storage     , & ! Output: [real(r8) (:)   ]  allocation to grain P storage (gP/m2/s)
        retransp_to_ppool            => veg_pf%retransp_to_ppool           , & ! Output: [real(r8) (:)   ]  deployment of retranslocated P (gP/m2/s)
        sminp_to_ppool               => veg_pf%sminp_to_ppool              , & ! Output: [real(r8) (:)   ]  deployment of soil mineral P uptake (gP/m2/s)
        ppool_to_leafp               => veg_pf%ppool_to_leafp              , & ! Output: [real(r8) (:)   ]  allocation to leaf P (gP/m2/s)
        ppool_to_leafp_storage       => veg_pf%ppool_to_leafp_storage      , & ! Output: [real(r8) (:)   ]  allocation to leaf P storage (gP/m2/s)
        ppool_to_frootp              => veg_pf%ppool_to_frootp             , & ! Output: [real(r8) (:)   ]  allocation to fine root P (gP/m2/s)
        ppool_to_frootp_storage      => veg_pf%ppool_to_frootp_storage     , & ! Output: [real(r8) (:)   ]  allocation to fine root P storage (gP/m2/s)
        ppool_to_livestemp           => veg_pf%ppool_to_livestemp          , & ! Output: [real(r8) (:)   ]
        ppool_to_livestemp_storage   => veg_pf%ppool_to_livestemp_storage  , & ! Output: [real(r8) (:)   ]
        ppool_to_deadstemp           => veg_pf%ppool_to_deadstemp          , & ! Output: [real(r8) (:)   ]
        ppool_to_deadstemp_storage   => veg_pf%ppool_to_deadstemp_storage  , & ! Output: [real(r8) (:)   ]
        ppool_to_livecrootp          => veg_pf%ppool_to_livecrootp         , & ! Output: [real(r8) (:)   ]
        ppool_to_livecrootp_storage  => veg_pf%ppool_to_livecrootp_storage , & ! Output: [real(r8) (:)   ]
        ppool_to_deadcrootp          => veg_pf%ppool_to_deadcrootp         , & ! Output: [real(r8) (:)   ]
        ppool_to_deadcrootp_storage  => veg_pf%ppool_to_deadcrootp_storage , & ! Output: [real(r8) (:)   ]
        p_allometry                  => cnstate_vars%p_allometry_patch      & ! Output: [real(r8) (:)   ]  P allocation index (DIM)
        )
      
     print *, "num_soilp:", num_soilp 
      !$acc parallel loop independent gang vector default(present)
      do fp = 1, num_soilp 
         p= filter_soilp(fp) 
         c= veg_pp%column(p)
         ivt = veg_pp%itype(p)
            ! set some local allocation variables
            f1 = froot_leaf(ivt)
            f2 = croot_stem(ivt)

            ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
            ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
            ! There was an error in this formula in previous version, where the coefficient
            ! was 0.004 instead of 0.0025.
            ! This variable allocation is only for trees. Shrubs have a constant
            ! allocation as specified in the pft-physiology file.  The value is also used
            ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

            if (stem_leaf(ivt) < 0._r8) then
                if (stem_leaf(ivt) == -1._r8) then
                    f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
                else
                    f3 = max((-1.0_r8*stem_leaf(ivt)*2.7_r8)/(1.0_r8+exp(-0.004_r8*(annsum_npp(p) - &
                              300.0_r8))) - 0.4_r8, 0.2_r8)
                end if
            else
                f3 = stem_leaf(ivt)
            end if

            f4 = flivewd(ivt)
            g1 = grperc(ivt)
            g2 = grpnow(ivt)
            cnl = leafcn(ivt)
            cnfr = frootcn(ivt)
            cnlw = livewdcn(ivt)
            cndw = deadwdcn(ivt)

            cpl = leafcp(ivt)
            cpfr = frootcp(ivt)
            cplw = livewdcp(ivt)
            cpdw = deadwdcp(ivt)

            fcur = fcur2(ivt)

            if (ivt >= npcropmin) then ! skip 2 generic crops
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

            if (veg_vp%nstor(ivt) > 1e-6_r8) then
              !N pool modification
              sminn_to_npool(p) = plant_ndemand(p) * min(fpg(c), fpg_p(c))
              sminp_to_ppool(p) = plant_pdemand(p) * min(fpg(c), fpg_p(c))

              rc   = veg_vp%nstor(ivt) * max(annsum_npp(p) * n_allometry(p) / c_allometry(p), 0.01_r8)
              rc_p = veg_vp%nstor(ivt) * max(annsum_npp(p) * p_allometry(p) / c_allometry(p), 0.01_r8)

              if ( carbon_only  .or.  carbonphosphorus_only ) then
                r = 1.0_r8
              else
                r  = max(1._r8,rc/max(npool(p), 1e-15_r8))
              end if
              plant_nalloc(p) = (plant_ndemand(p) + retransn_to_npool(p)) / r

              if ( carbon_only  .or.  carbonnitrogen_only ) then
                r = 1.0_r8
              else
                r  = max(1._r8,rc_p/max(ppool(p), 1e-15_r8))
              end if
              plant_palloc(p) = (plant_pdemand(p) + retransp_to_ppool(p)) / r

            else
              sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
              sminp_to_ppool(p) = plant_pdemand(p) * fpg_p(c)

              plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)
              plant_palloc(p) = sminp_to_ppool(p) + retransp_to_ppool(p)
            end if

            ! calculate the associated carbon allocation, and the excess
            ! carbon flux that must be accounted for through downregulation

            if( .not.carbonphosphorus_only .and. .not.carbonnitrogen_only .and. .not.carbon_only )then
                if( plant_nalloc(p) * (c_allometry(p)/n_allometry(p)) < &
                    plant_palloc(p) * (c_allometry(p)/p_allometry(p)) )then

                    plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
                    plant_palloc(p) = plant_nalloc(p) * (p_allometry(p)/n_allometry(p))
                    !in case of strong N limitation, and plant_palloc(p) < retransp_to_ppool(p)
                    if (veg_vp%nstor(ivt) < 1e-6_r8) then
                        sminp_to_ppool(p) = max(plant_palloc(p) - retransp_to_ppool(p),0.0_r8)
                        retransp_to_ppool(p) = min(plant_palloc(p), retransp_to_ppool(p))
                    end if
                else
                    plant_calloc(p) = plant_palloc(p) * (c_allometry(p)/p_allometry(p))
                    plant_nalloc(p) = plant_palloc(p) * (n_allometry(p)/p_allometry(p))
                    ! in case of strong P limitation, and plant_nalloc(p) < retransn_to_npool(p)
                    if (veg_vp%nstor(ivt) < 1e-6_r8) then
                        sminn_to_npool(p) = max(plant_nalloc(p) - retransn_to_npool(p), 0.0_r8)
                        retransn_to_npool(p) = min(plant_nalloc(p) , retransn_to_npool(p))
                    end if
                endif
            endif

            if(carbonphosphorus_only .or. carbon_only )then
                plant_calloc(p) = plant_palloc(p) * (c_allometry(p)/p_allometry(p))
            endif

            if(carbonnitrogen_only .or. carbon_only )then
                plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
                plant_palloc(p) = plant_calloc(p) * (p_allometry(p)/c_allometry(p))
                if (veg_vp%nstor(ivt) < 1e-6_r8) then
                    sminp_to_ppool(p) = max(plant_palloc(p) - retransp_to_ppool(p), 0.0_r8)
                end if
            endif

            excess_cflux(p) = availc(p) - plant_calloc(p)

            ! reduce gpp fluxes due to N limitation
            if (gpp(p) > 0.0_r8) then
                downreg(p) = excess_cflux(p)/gpp(p)

                if (veg_vp%br_xr(veg_pp%itype(p)) > 1e-9_r8) then
                    !Excess carbon goes to temporary NSC pool instead of
                    !instantaneous downregulation
                    psnsun_to_cpool(p) = psnsun_to_cpool(p)
                    psnshade_to_cpool(p) = psnshade_to_cpool(p)
                    if ( use_c13 ) then
                        c13_veg_cf%psnsun_to_cpool(p) = c13_veg_cf%psnsun_to_cpool(p)
                        c13_veg_cf%psnshade_to_cpool(p) = c13_veg_cf%psnshade_to_cpool(p)
                    endif

                    if ( use_c14 ) then
                        c14_veg_cf%psnsun_to_cpool(p) = c14_veg_cf%psnsun_to_cpool(p)
                        c14_veg_cf%psnshade_to_cpool(p) = c14_veg_cf%psnshade_to_cpool(p)
                    endif

                else
                    psnsun_to_cpool(p)   = psnsun_to_cpool(p)  *(1._r8 - downreg(p))
                    psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1._r8 - downreg(p))
                    if ( use_c13 ) then
                        c13_veg_cf%psnsun_to_cpool(p)   = c13_veg_cf%psnsun_to_cpool(p)  *(1._r8 - downreg(p))
                        c13_veg_cf%psnshade_to_cpool(p) = c13_veg_cf%psnshade_to_cpool(p)*(1._r8 - downreg(p))
                    endif

                    if ( use_c14 ) then
                        c14_veg_cf%psnsun_to_cpool(p)   = c14_veg_cf%psnsun_to_cpool(p)  *(1._r8 - downreg(p))
                        c14_veg_cf%psnshade_to_cpool(p) = c14_veg_cf%psnshade_to_cpool(p)*(1._r8 - downreg(p))
                    endif
                endif
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
       !fcur = fcur2(ivt)

       cpool_to_leafc(p)          = nlc * fcur
       cpool_to_leafc_storage(p)  = nlc * (1._r8 - fcur)
       cpool_to_frootc(p)         = nlc * f1 * fcur
       cpool_to_frootc_storage(p) = nlc * f1 * (1._r8 - fcur)
       if (woody(ivt) == 1._r8) then
          cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
          cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
          cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
          cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
          cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
          cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
          cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
          cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops
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
       !fcur = fcur2(ivt)

       npool_to_leafn(p)          = (nlc / cnl) * fcur
       npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur)
       npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
       npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
       if (woody(ivt) == 1._r8) then
          npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
          npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
          npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
          npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
          npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
          npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
          npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
          npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops
          cng = graincn(ivt)
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
       !fcur = fcur2(ivt)
       ppool_to_leafp(p)          = (nlc / cpl) * fcur
       ppool_to_leafp_storage(p)  = (nlc / cpl) * (1._r8 - fcur)
       ppool_to_frootp(p)         = (nlc * f1 / cpfr) * fcur
       ppool_to_frootp_storage(p) = (nlc * f1 / cpfr) * (1._r8 - fcur)
       if (woody(ivt) == 1._r8) then
          ppool_to_livestemp(p)          = (nlc * f3 * f4 / cplw) * fcur
          ppool_to_livestemp_storage(p)  = (nlc * f3 * f4 / cplw) * (1._r8 -fcur)
          ppool_to_deadstemp(p)          = (nlc * f3 * (1._r8 - f4) / cpdw) *fcur
          ppool_to_deadstemp_storage(p)  = (nlc * f3 * (1._r8 - f4) / cpdw) *(1._r8 - fcur)
          ppool_to_livecrootp(p)         = (nlc * f2 * f3 * f4 / cplw) * fcur
          ppool_to_livecrootp_storage(p) = (nlc * f2 * f3 * f4 / cplw) * (1._r8 -fcur)
          ppool_to_deadcrootp(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* fcur
          ppool_to_deadcrootp_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cpdw)* (1._r8 - fcur)
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops
          cpg = graincp(ivt)
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
       if (woody(ivt) == 1._r8) then
          gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
          gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)
          gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
          gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
       end if
       if (ivt >= npcropmin) then ! skip 2 generic crops
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
      end do 

     end associate
  end subroutine DistributeN_RD

  ! ======================================================================================

  subroutine NAllocationRD( num_soilc, filter_soilc, &
       col_plant_ndemand_vr,   &! IN (j)
       potential_immob_vr,  &    ! IN (j)
       compet_plants_nh4,   &    ! IN
       compet_decomp_nh4,   &    ! IN
       dt,                  &    ! IN
       smin_nh4_vr,         &    ! IN (j)
       fpi_nh4_vr,          &    ! OUT (:)
       actual_immob_nh4_vr, &    ! OUT (:)
       smin_nh4_to_plant_vr, &   ! OUT (:)
       smin_no3_vr,          &   ! IN (j)
       compet_plants_no3,    &   ! IN
       compet_decomp_no3,    &   ! IN
       compet_nit,           &   ! IN
       compet_denit,         &   ! IN
       pot_f_nit_vr,         &   ! IN (j)
       pot_f_denit_vr,       &   ! IN (j)
       fpi_no3_vr,           &   ! OUT (j)
       actual_immob_no3_vr,  &   ! OUT (j)
       smin_no3_to_plant_vr, &   ! OUT (j)
       f_nit_vr,             &   ! OUT (j)
       f_denit_vr)               ! OUT (j)

    use elm_varpar, only : nlevdecomp
    ! Arguments
    integer , intent(in)  :: num_soilc 
    integer , intent(in)  :: filter_soilc(:) 
    real(r8), intent(in)  :: col_plant_ndemand_vr(:,:)    ! How much N all plants demand as group [g/m3]
    real(r8), intent(in)  :: potential_immob_vr(:,:)      ! potential N immobilization [g/m3/s]
    real(r8), intent(in)  :: compet_plants_nh4       ! relative competability of plants (unitless)
    real(r8), intent(in)  :: compet_decomp_nh4       ! relative competability of decomposers (unitless)
    real(r8), intent(in)  :: dt                      ! timestep [seconds]
    real(r8), intent(in)  :: smin_nh4_vr(:,:)             ! mineralized nh4 [g/m3]
    real(r8), intent(inout) :: fpi_nh4_vr (:,:)             ! fraction of potential immobilization supplied by nh4 (no units)
    real(r8), intent(inout) :: actual_immob_nh4_vr(:,:)     ! actual nh4 immobilization [g/m3/s]
    real(r8), intent(inout) :: smin_nh4_to_plant_vr(:,:)    ! nh4 flux to plant competitors [g/m3/s]

    ! Optional (for NO3)
    real(r8), intent(in)  :: smin_no3_vr(:,:)             ! mineralized no3 [g/m3]
    real(r8), intent(in)  :: compet_plants_no3       ! relative competability of plants (unitless)
    real(r8), intent(in)  :: compet_decomp_no3       ! relative competability of decomposers (unitless)
    real(r8), intent(in)  :: compet_nit              ! relative competitiveness of nitrifiers for NH4
    real(r8), intent(in)  :: compet_denit            ! relative competitiveness of denitrifiers for NO3
    real(r8), intent(in)  :: pot_f_nit_vr(:,:)            ! potential soil nitrification flux [g/m3/s]
    real(r8), intent(in)  :: pot_f_denit_vr(:,:)          ! potential soil denitrification flux [g/m3/s]
    real(r8), intent(inout) :: fpi_no3_vr(:,:)            ! fraction of potential immobilization supplied by NO3
    real(r8), intent(inout) :: actual_immob_no3_vr(:,:)   ! actual no3 immobilization [g/m3/s]
    real(r8), intent(inout) :: smin_no3_to_plant_vr(:,:)  ! no3 flux to plant competitors [g/m3/s]
    real(r8), intent(inout) :: f_nit_vr(:,:)              ! soil nitrification flux [g/m3/s]
    real(r8), intent(inout) :: f_denit_vr(:,:)            ! soil denitrification flux [g/m3/s]

    ! Locals
    real(r8) :: sum_nh4_demand        ! Total nh4 demand over all competitors
    real(r8) :: sum_nh4_demand_scaled ! Total nh4 demand, but scaled by competitivness
    real(r8) :: sum_no3_demand        ! "" no3
    real(r8) :: sum_no3_demand_scaled ! "" no3
    integer  :: j,fc,c                ! soil decomp layer loop

    !$acc parallel loop independent gang vector default(present) collapse(2) 
    do j = 1, nlevdecomp
       do fc=1,num_soilc !col_loop
         c = filter_soilc(fc)

         sum_nh4_demand        = col_plant_ndemand_vr(c,j) + potential_immob_vr(c,j) + pot_f_nit_vr(c,j)
         sum_nh4_demand_scaled = col_plant_ndemand_vr(c,j) * compet_plants_nh4 + &
                                     potential_immob_vr(c,j)*compet_decomp_nh4 + pot_f_nit_vr(c,j)*compet_nit

         if (sum_nh4_demand*dt < smin_nh4_vr(c,j)) then
             ! NH4 availability is not limiting immobilization or plant
             ! uptake, and all can proceed at their potential rates
             fpi_nh4_vr(fc,j) = 1.0_r8
             actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
             smin_nh4_to_plant_vr(c,j) = col_plant_ndemand_vr(c,j)
             f_nit_vr(c,j) = pot_f_nit_vr(c,j)

         else

            ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
            ! plant growth demands, so these three demands compete for available
            ! soil mineral NH4 resource.
            if (sum_nh4_demand > 0.0_r8 .and. smin_nh4_vr(c,j) > 0.0_r8 &
                 .and. sum_nh4_demand_scaled > 0.0_r8) then
               actual_immob_nh4_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(potential_immob_vr(c,j)* &
                    compet_decomp_nh4 / sum_nh4_demand_scaled), potential_immob_vr(c,j))
               smin_nh4_to_plant_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*&
                    (col_plant_ndemand_vr(c,j)*compet_plants_nh4 / sum_nh4_demand_scaled), &
                    col_plant_ndemand_vr(c,j))
               f_nit_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nit_vr(c,j)*compet_nit / &
                    sum_nh4_demand_scaled), pot_f_nit_vr(c,j))
            else
               actual_immob_nh4_vr(c,j) = 0.0_r8
               smin_nh4_to_plant_vr(c,j) = 0.0_r8
               f_nit_vr(c,j) = 0.0_r8
            end if

            if (potential_immob_vr(c,j) > 0.0_r8) then
               fpi_nh4_vr(fc,j) = actual_immob_nh4_vr(c,j) / potential_immob_vr(c,j)
            else
               fpi_nh4_vr(fc,j)= 0.0_r8
            end if

         end if    ! if (sum_nh4_demand*dt < smin_nh4_vr) then
       !
       ! If we passed in parameters and mineralized no3, then
       ! we are free to calculate competitive allocation rates on it
       ! ------------------------------------------------------------------------
       
       ! next compete for no3
       sum_no3_demand = (col_plant_ndemand_vr(c,j)-smin_nh4_to_plant_vr(c,j)) + &
            (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)

       sum_no3_demand_scaled = (col_plant_ndemand_vr(c,j)-smin_nh4_to_plant_vr(c,j)) &
            * compet_plants_no3 + (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 &
            + pot_f_denit_vr(c,j)*compet_denit

       if (sum_no3_demand*dt < smin_no3_vr(c,j)) then

          ! NO3 availability is not limiting immobilization or plant
          ! uptake, and all can proceed at their potential rates
          fpi_no3_vr(fc,j) = 1.0_r8 -  fpi_nh4_vr(fc,j)
          actual_immob_no3_vr(c,j) = (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
          smin_no3_to_plant_vr(c,j) = (col_plant_ndemand_vr(c,j)-smin_nh4_to_plant_vr(c,j))
          f_denit_vr(c,j) = pot_f_denit_vr(c,j)

       else

          ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
          ! plant growth demands, so these three demands compete for available
          ! soil mineral NO3 resource.
          if (sum_no3_demand > 0.0_r8 .and. smin_no3_vr(c,j) > 0.0_r8 &
               .and. sum_no3_demand_scaled > 0.0_r8) then
             actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)- &
                  actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled), &
                  potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
             smin_no3_to_plant_vr(c,j) = min((smin_no3_vr(c,j)/dt) * &
                  ((col_plant_ndemand_vr(c,j)-smin_nh4_to_plant_vr(c,j)) * &
                  compet_plants_no3 / sum_no3_demand_scaled), &
                  col_plant_ndemand_vr(c,j)-smin_nh4_to_plant_vr(c,j))
             f_denit_vr(c,j) =  min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / &
                  sum_no3_demand_scaled), pot_f_denit_vr(c,j))
          else
             actual_immob_no3_vr(c,j) = 0.0_r8
             smin_no3_to_plant_vr(c,j) = 0.0_r8
             f_denit_vr (c,j)= 0.0_r8
          end if

          if (potential_immob_vr(c,j) > 0.0_r8) then
             fpi_no3_vr(fc,j) = actual_immob_no3_vr(c,j) / potential_immob_vr(c,j)
          else
             fpi_no3_vr(fc,j) = 0.0_r8
          end if
       end if ! if sum_no3_demand*dt < smin_no3_vr
     end do    ! j = 1,nlevdecomp
    end do 
  end subroutine NAllocationRD

  ! ======================================================================================

  subroutine PAllocationRD( num_soilc, filter_soilc, &
       col_plant_pdemand_vr, &    ! IN
       potential_immob_p_vr, &    ! IN 
       solutionp_vr,         &    ! IN 
       dt,                   &    ! IN
       fpi_p_vr,             &    ! OUT 
       actual_immob_p_vr,    &    ! OUT 
       sminp_to_plant_vr,    &    ! OUT 
       supplement_to_sminp_vr)    ! OUT 

    use elm_varctl       , only:  carbon_only, carbonnitrogen_only
    use elm_varpar, only : nlevdecomp

    ! Arguments
    integer , intent(in) :: num_soilc
    integer , intent(in) :: filter_soilc(:)
    real(r8), intent(in) :: col_plant_pdemand_vr(:,:) ! demand on phos, all plant grouped [g/m3]
    real(r8), intent(in) :: potential_immob_p_vr(:,:)  ! potential P immobilization [g/m3/s]
    real(r8), intent(in) :: solutionp_vr(:,:)          ! soil mineral P   [g/m3]
    real(r8), intent(in) :: dt      ! timestep in seconds
    real(r8), intent(inout) :: fpi_p_vr(:,:)             ! fraction of potential immobilization supplied by p
    real(r8), intent(inout) :: actual_immob_p_vr(:,:)    ! actual P immobilization [g/m3/s]
    real(r8), intent(inout) :: sminp_to_plant_vr(:,:)    ! P flux to plant competitors [g/m3/s]
    real(r8), intent(inout) :: supplement_to_sminp_vr(:,:)

    ! Locals
    real(r8) :: sum_pdemand          ! Total phos demand over all competitors
    integer :: j, fc ,c 
   !$acc parallel loop independent collapse(2) gang vector default(present)
    do j = 1, nlevdecomp
      do fc = 1, num_soilc
         c = filter_soilc(fc)
          
         sum_pdemand = col_plant_pdemand_vr(c,j) + potential_immob_p_vr(c,j)

        if (sum_pdemand*dt < solutionp_vr(c,j)) then

            ! P availability is not limiting immobilization or plant
            ! uptake, and both can proceed at their potential rates
            fpi_p_vr(c,j) = 1.0_r8
            actual_immob_p_vr(c,j) = potential_immob_p_vr(c,j)
            sminp_to_plant_vr(c,j) = col_plant_pdemand_vr(c,j)

         elseif(carbon_only .or. carbonnitrogen_only    ) then

            fpi_p_vr(c,j) = 1.0_r8
            actual_immob_p_vr(c,j) = potential_immob_p_vr(c,j)
            sminp_to_plant_vr(c,j) =  col_plant_pdemand_vr(c,j)
            supplement_to_sminp_vr(c,j) = sum_pdemand - (solutionp_vr(c,j)/dt)

         else
            ! P availability can not satisfy the sum of immobilization and
            ! plant growth demands, so these two demands compete for
            ! available soil mineral solution P resource.
            if (sum_pdemand > 0.0_r8 .and. solutionp_vr(c,j) >0._r8) then
               actual_immob_p_vr(c,j) = (solutionp_vr(c,j)/dt)*(potential_immob_p_vr(c,j) / sum_pdemand)
            else
               actual_immob_p_vr(c,j) = 0.0_r8
            end if

            if (potential_immob_p_vr(c,j) > 0.0_r8) then
               fpi_p_vr(c,j) = actual_immob_p_vr(c,j) / potential_immob_p_vr(c,j)
            else
               fpi_p_vr(c,j) = 0.0_r8
            end if

            sminp_to_plant_vr(c,j) = max( 0._r8,(solutionp_vr(c,j)/dt) - actual_immob_p_vr(c,j) )
         end if
      end do 
   end do 
  end subroutine PAllocationRD

  !-------------------------------------------------------------------------------------------------
  subroutine calc_nuptake_prof(num_soilc, filter_soilc, cnstate_vars, nuptake_prof)
     ! bgc interface & pflotran:
     ! nuptake_prof is used in Allocation1, 2, 3
     ! !USES:
     use elm_varpar       , only: nlevdecomp
     ! !ARGUMENTS:
     integer              , intent(in)    :: num_soilc        ! number of soil columns in filter
     integer              , intent(in)    :: filter_soilc(:)  ! filter for soil columns
     type(cnstate_type)   , intent(in)    :: cnstate_vars
     real(r8)             , intent(inout) :: nuptake_prof(num_soilc, 1:nlevdecomp)

     integer :: c,j,fc                                            !indices
     real(r8):: sminn_tot(num_soilc)
     real(r8):: sminn_vr_loc, sum1
     !-----------------------------------------------------------------------

     associate( &
          nfixation_prof               => cnstate_vars%nfixation_prof_col      , & ! Input: [real(r8) (:,:) ]
          smin_no3_vr                  => col_ns%smin_no3_vr                  , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
          smin_nh4_vr                  => col_ns%smin_nh4_vr                    & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
          )

          ! column loops to resolve plant/heterotroph competition for mineral N

          !NOTE:  This data creation can be pipelined if need be
          !$acc enter data create(sminn_tot(1:num_soilc))
         if(not (use_pflotran .and. pf_cmode) ) then

          !$acc parallel loop independent gang worker default(present) private(sum1,c)
          do fc=1,num_soilc
             sum1 = 0._r8
             c = filter_soilc(fc)
             !$acc loop vector reduction(+:sum1) private(sminn_vr_loc)
             do j = 1, nlevdecomp
                sminn_vr_loc = smin_no3_vr(c,j) + smin_nh4_vr(c,j)
                sum1 = sum1 + sminn_vr_loc * dzsoi_decomp(j)
             end do
             sminn_tot(fc) = sum1
          end do

          !$acc parallel loop independent gang default(present)
          do j = 1, nlevdecomp
             !$acc loop worker vector independent private(c,sminn_vr_loc)
             do fc=1,num_soilc
                c = filter_soilc(fc)
                sminn_vr_loc = smin_no3_vr(c,j) + smin_nh4_vr(c,j)
                if (sminn_tot(fc)  >  0._r8) then
                    !original:  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
                   nuptake_prof(fc,j) = sminn_vr_loc / sminn_tot(fc)
                else
                   nuptake_prof(fc,j) = nfixation_prof(c,j)
                end if
             end do
          end do

          end if
          if(use_pflotran .and. pf_cmode) then
             !$acc parallel loop independent gang worker default(present) private(sum1,c)
            do fc=1,num_soilc
               sum1 = 0._r8
               c = filter_soilc(fc)
               !$acc loop vector reduction(+:sum1)
               do j = 1, nlevdecomp
                  sminn_vr_loc = smin_no3_vr(c,j) + smin_nh4_vr(c,j)
                  sum1 = sum1 + sminn_vr_loc * dzsoi_decomp(j) &
                           *(nfixation_prof(c,j)*dzsoi_decomp(j))         ! weighted by froot fractions in annual max. active layers
             end do
             sminn_tot(fc) = sum1
            end do

            do j = 1, nlevdecomp
              do fc=1,num_soilc
                 c = filter_soilc(fc)
                 sminn_vr_loc = smin_no3_vr(c,j) + smin_nh4_vr(c,j)

                 if (sminn_tot(fc)  >  0._r8) then
                    nuptake_prof(fc,j) = sminn_vr_loc/ sminn_tot(fc) &
                          *(nfixation_prof(c,j)*dzsoi_decomp(j))         ! weighted by froot fractions in annual max. active layers
                 else
                    nuptake_prof(fc,j) = nfixation_prof(c,j)
                 end if
              end do
           end do
       end if
       !$acc exit data delete(sminn_tot(1:num_soilc))

     end associate

  end subroutine calc_nuptake_prof

!-------------------------------------------------------------------------------------------------
subroutine calc_puptake_prof(num_soilc, filter_soilc, cnstate_vars, puptake_prof)
  ! bgc interface & pflotran:
  ! puptake_prof is used in Allocation1, 2, & 3
  ! !USES:
  use elm_varpar       , only: nlevdecomp
  ! !ARGUMENTS:
  integer                  , intent(in)    :: num_soilc        ! number of soil columns in filter
  integer                  , intent(in)    :: filter_soilc(:)  ! filter for soil columns
  type(cnstate_type)       , intent(in)    :: cnstate_vars
  real(r8)                 , intent(inout) :: puptake_prof(num_soilc, 1:nlevdecomp)

  integer :: c,j,fc                                            !indices
  real(r8):: solutionp_tot(num_soilc), sum1

  !-----------------------------------------------------------------------
  associate( &
       nfixation_prof               => cnstate_vars%nfixation_prof_col , & ! Output: [real(r8) (:,:) ]
       solutionp_vr                 => col_ps%solutionp_vr  & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N
       )
       ! column loops to resolve plant/heterotroph competition for mineral N
       ! init sminn_tot
       ! do fc=1,num_soilc
       !    solutionp_tot(fc) = 0.
       ! end do
       !$acc enter data create(solutionp_tot(1:num_soilc))

       !$acc parallel loop independent gang worker default(present) private(c,sum1)
       do fc=1,num_soilc
          c = filter_soilc(fc)
          sum1 = 0._r8
          !$acc loop vector reduction(+:sum1)
          do j = 1, nlevdecomp
             sum1 = sum1 + solutionp_vr(c,j) * dzsoi_decomp(j)
          end do
          solutionp_tot(fc) = sum1
       end do

       !$acc parallel loop independent gang default(present)
       do j = 1, nlevdecomp
          !$acc loop vector independent private(c)
          do fc=1,num_soilc
             c = filter_soilc(fc)
             !!! add P demand calculation
             if (solutionp_tot(fc)  >  0.) then
                puptake_prof(fc,j) = solutionp_vr(c,j) / solutionp_tot(fc)
             else
                puptake_prof(fc,j) = nfixation_prof(c,j)      ! need modifications
             endif

          end do
       end do

       !$acc exit data delete(solutionp_tot(:))

  end associate

end subroutine calc_puptake_prof

!-----------------------------------------------------------------------
    subroutine dynamic_plant_alloc( nutrient_scalar, water_scalar, laindex, alloc_leaf, alloc_stem, alloc_froot, woody)

    ! !DESCRIPTION
    ! Added by Qing Zhu 2015 based on P. Friedlingstein DOI: 10.1046/j.1365-2486.1999.00269.x
    ! allocation coefficients for leaf, stem and root are not fixed
    ! update allocation coefficients based on nutrient and light availability
    ! (1) light limited, allocate more C into stem
    ! (2) nutrient/water limited, allocate more C into root

    ! !USES:
      !$acc routine seq
    use pftvarcon      , only : laimax

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
    if (laindex > laimax) then
       if (woody == 1.0_r8) then
          alloc_stem = alloc_stem + alloc_leaf/2._r8 - 0.005_r8
          alloc_froot = alloc_froot + alloc_leaf/2._r8 - 0.005_r8
       else
          alloc_froot = alloc_froot + alloc_leaf - 0.01_r8
       end if
       alloc_leaf = 0.01_r8
    end if

  end subroutine dynamic_plant_alloc

!-----------------------------------------------------------------------

end module AllocationMod
