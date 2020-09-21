module DecompCascadeCNMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coeffiecients used in the decomposition cascade submodel.  
  ! This uses the CN parameters as in CLMCN 4.0
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_TKFRZ
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : nlevsoi, nlevgrnd, nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar             , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, spinup_state, anoxia, use_lch4, use_vertsoilc, use_fates
  use elm_varcon             , only : zsoi
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use SharedParamsMod      , only : ParamsShareInst, anoxia_wtsat, nlev_soildecomp_standard 
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use SoilStateType          , only : soilstate_type
  use CanopyStateType        , only : canopystate_type
  use TemperatureType        , only : temperature_type 
  use ch4Mod                 , only : ch4_type
  use ColumnType             , only : col_pp   
  use ColumnDataType         , only : col_es, col_cf  
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readDecompCNParams
  public :: init_decompcascade_cn
  public :: decomp_rate_constants_cn

  type, private :: DecompCNParamsType
     real(r8):: cn_s1_cn        !C:N for SOM 1
     real(r8):: cn_s2_cn        !C:N for SOM 2
     real(r8):: cn_s3_cn        !C:N for SOM 3
     real(r8):: cn_s4_cn        !C:N for SOM 4

     real(r8):: np_s1_new_cn        !C:P for SOM 1
     real(r8):: np_s2_new_cn        !C:P for SOM 2
     real(r8):: np_s3_new_cn        !C:P for SOM 3
     real(r8):: np_s4_new_cn        !C:P for SOM 4

     real(r8):: cp_s1_new_cn        !C:P for SOM 1
     real(r8):: cp_s2_new_cn        !C:P for SOM 2
     real(r8):: cp_s3_new_cn        !C:P for SOM 3
     real(r8):: cp_s4_new_cn        !C:P for SOM 4

     real(r8):: rf_l1s1_cn      !respiration fraction litter 1 -> SOM 1
     real(r8):: rf_l2s2_cn      !respiration fraction litter 2 -> SOM 2
     real(r8):: rf_l3s3_cn      !respiration fraction litter 3 -> SOM 3
     real(r8):: rf_s1s2_cn      !respiration fraction SOM 1 -> SOM 2
     real(r8):: rf_s2s3_cn      !respiration fraction SOM 2 -> SOM 3
     real(r8):: rf_s3s4_cn      !respiration fraction SOM 3 -> SOM 4

     real(r8) :: cwd_fcel_cn    !cellulose fraction for CWD
     real(r8) :: cwd_flig_cn    !

     real(r8) :: k_l1_cn        !decomposition rate for litter 1
     real(r8) :: k_l2_cn        !decomposition rate for litter 2
     real(r8) :: k_l3_cn        !decomposition rate for litter 3
     real(r8) :: k_s1_cn        !decomposition rate for SOM 1
     real(r8) :: k_s2_cn        !decomposition rate for SOM 2
     real(r8) :: k_s3_cn        !decomposition rate for SOM 3
     real(r8) :: k_s4_cn        !decomposition rate for SOM 4

     real(r8) :: k_frag_cn      !fragmentation rate for CWD
     real(r8) :: minpsi_cn      !minimum soil water potential for heterotrophic resp

     integer  :: nsompools = 4 
     integer  :: nlitpools = 3
     integer  :: ncwdpools = 1  
     real(r8), allocatable :: spinup_vector(:) ! multipliers for soil decomp during accelerated spinup

     
  end type DecompCNParamsType

  type(DecompCNParamsType),private ::  DecompCNParamsInst

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readDecompCNParams ( ncid )
    !
    ! !USES:
    use ncdio_pio , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !CALLED FROM:   readParamsMod.F90::CNParamsReadFile
    !
    ! !REVISION HISTORY:
    !  Dec 3 2012 : Created by S. Muszala
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'DecompCNParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading

    !EOP
    !-----------------------------------------------------------------------

    ! These are not read off of netcdf file
    allocate(DecompCNParamsInst%spinup_vector(DecompCNParamsInst%nsompools+DecompCNParamsInst%nlitpools+ &
        DecompCNParamsInst%ncwdpools))
    !DecompCNParamsInst%spinup_vector(:) = (/ 1.0_r8, 1.0_r8, 5.0_r8, 70.0_r8 /)
    !These will be set below

    ! Read off of netcdf file
    tString='cn_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cn_s1_cn=tempr

    tString='cn_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cn_s2_cn=tempr

    tString='cn_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cn_s3_cn=tempr

    tString='cn_s4'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cn_s4_cn=tempr

!!! read in phosphorus variables - X. YANG
    tString='np_s1_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%np_s1_new_cn=tempr

    tString='np_s2_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%np_s2_new_cn=tempr

    tString='np_s3_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%np_s3_new_cn=tempr

    tString='np_s4_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%np_s4_new_cn=tempr

    tString='rf_l1s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_l1s1_cn=tempr

    tString='rf_l2s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_l2s2_cn=tempr

    tString='rf_l3s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_l3s3_cn=tempr

    tString='rf_s1s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_s1s2_cn=tempr

    tString='rf_s2s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_s2s3_cn=tempr

    tString='rf_s3s4'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%rf_s3s4_cn=tempr

    tString='cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cwd_fcel_cn=tempr

    tString='k_l1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_l1_cn=tempr

    tString='k_l2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_l2_cn=tempr

    tString='k_l3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_l3_cn=tempr

    tString='k_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_s1_cn=tempr

    tString='k_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_s2_cn=tempr

    tString='k_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_s3_cn=tempr

    tString='k_s4'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_s4_cn=tempr

    tString='k_frag'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%k_frag_cn=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%minpsi_cn=tempr 

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    DecompCNParamsInst%cwd_flig_cn=tempr

    DecompCNParamsInst%spinup_vector(1) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_l1_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(2) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_l2_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(3) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_l3_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(4) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_s1_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(5) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_s2_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(6) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_s3_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(7) = max(1.0_r8, 1.0_r8 / (DecompCNParamsInst%k_s4_cn * 365.0_r8))
    DecompCNParamsInst%spinup_vector(8) = max(1.0_r8, 0.5_r8 * 1.0_r8 / (DecompCNParamsInst%k_frag_cn * 365.0_r8))

  end subroutine readDecompCNParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_cn(bounds, cnstate_vars)
    !
    ! !DESCRIPTION:
    !  initialize rate constants and decomposition pathways for the BGC model originally implemented in CLM-CN
    !  written by C. Koven based on original CLM4 decomposition cascade by P. Thornton
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds  
    type(cnstate_type), intent(inout) :: cnstate_vars
    !
    !-- properties of each pathway along decomposition cascade 
    !-- properties of each decomposing pool
    real(r8) :: rf_l1s1      !respiration fraction litter 1 -> SOM 1
    real(r8) :: rf_l2s2      !respiration fraction litter 2 -> SOM 2
    real(r8) :: rf_l3s3      !respiration fraction litter 3 -> SOM 3
    real(r8) :: rf_s1s2      !respiration fraction SOM 1 -> SOM 2
    real(r8) :: rf_s2s3      !respiration fraction SOM 2 -> SOM 3
    real(r8) :: rf_s3s4      !respiration fraction SOM 3 -> SOM 4
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    real(r8) :: cn_s4
    real(r8) :: np_s1_new
    real(r8) :: np_s2_new
    real(r8) :: np_s3_new
    real(r8) :: np_s4_new

    integer :: i_litr1
    integer :: i_litr2
    integer :: i_litr3
    integer :: i_soil1
    integer :: i_soil2
    integer :: i_soil3
    integer :: i_soil4
    integer :: i_atm
    integer :: i_l1s1
    integer :: i_l2s2
    integer :: i_l3s3
    integer :: i_s1s2
    integer :: i_s2s3
    integer :: i_s3s4
    integer :: i_s4atm
    integer :: i_cwdl2
    integer :: i_cwdl3
    !-----------------------------------------------------------------------

    associate(&
         rf_decomp_cascade              =>    cnstate_vars%rf_decomp_cascade_col                , & ! Output:  [real(r8)           (:,:,:) ]  respired fraction in decomposition step (frac)       
         pathfrac_decomp_cascade        =>    cnstate_vars%pathfrac_decomp_cascade_col          , & ! Output:  [real(r8)           (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)
         
         cascade_step_name              =>    decomp_cascade_con%cascade_step_name              , & ! Output:   [character(len=8)  (:)     ]  name of transition                               
         cascade_donor_pool             =>    decomp_cascade_con%cascade_donor_pool             , & ! Output:   [integer           (:)     ]  which pool is C taken from for a given decomposition step 
         cascade_receiver_pool          =>    decomp_cascade_con%cascade_receiver_pool          , & ! Output:   [integer           (:)     ]  which pool is C added to for a given decomposition step   
         floating_cn_ratio_decomp_pools =>    decomp_cascade_con%floating_cn_ratio_decomp_pools , & ! Output:   [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
         floating_cp_ratio_decomp_pools =>    decomp_cascade_con%floating_cp_ratio_decomp_pools , & ! Output:   [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
         decomp_pool_name_restart       =>    decomp_cascade_con%decomp_pool_name_restart       , & ! Output:   [character(len=8)  (:)     ]  name of pool for restart files                   
         decomp_pool_name_history       =>    decomp_cascade_con%decomp_pool_name_history       , & ! Output:   [character(len=8)  (:)     ]  name of pool for history files                   
         decomp_pool_name_long          =>    decomp_cascade_con%decomp_pool_name_long          , & ! Output:   [character(len=20) (:)     ]  name of pool for netcdf long names              
         decomp_pool_name_short         =>    decomp_cascade_con%decomp_pool_name_short         , & ! Output:   [character(len=8)  (:)     ]  name of pool for netcdf short names              
         is_litter                      =>    decomp_cascade_con%is_litter                      , & ! Output:   [logical           (:)     ]  TRUE => pool is a litter pool                             
         is_soil                        =>    decomp_cascade_con%is_soil                        , & ! Output:   [logical           (:)     ]  TRUE => pool is a soil pool                               
         is_cwd                         =>    decomp_cascade_con%is_cwd                         , & ! Output:   [logical           (:)     ]  TRUE => pool is a cwd pool                                
         initial_cn_ratio               =>    decomp_cascade_con%initial_cn_ratio               , & ! Output:   [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
         initial_cp_ratio               =>    decomp_cascade_con%initial_cp_ratio               , & ! Output:   [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
         initial_stock                  =>    decomp_cascade_con%initial_stock                  , & ! Output:   [real(r8)          (:)     ]  initial concentration for seeding at spinup              
         is_metabolic                   =>    decomp_cascade_con%is_metabolic                   , & ! Output:   [logical           (:)     ]  TRUE => pool is metabolic material                        
         is_cellulose                   =>    decomp_cascade_con%is_cellulose                   , & ! Output:   [logical           (:)     ]  TRUE => pool is cellulose                                 
         is_lignin                      =>    decomp_cascade_con%is_lignin                      , & ! Output:   [logical           (:)     ]  TRUE => pool is lignin                                    
         spinup_factor                  =>    decomp_cascade_con%spinup_factor                    & ! Output:   [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )
      
      !------- time-constant coefficients ---------- !
      ! set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0)
      cn_s1=DecompCNParamsInst%cn_s1_cn
      cn_s2=DecompCNParamsInst%cn_s2_cn
      cn_s3=DecompCNParamsInst%cn_s3_cn
      cn_s4=DecompCNParamsInst%cn_s4_cn
       
      ! set soil organic matter C:P ratios -X. YANG 
      np_s1_new=DecompCNParamsInst%np_s1_new_cn
      np_s2_new=DecompCNParamsInst%np_s2_new_cn
      np_s3_new=DecompCNParamsInst%np_s3_new_cn
      np_s4_new=DecompCNParamsInst%np_s4_new_cn

      ! set respiration fractions for fluxes between compartments
      ! (from Biome-BGC v4.2.0)
      rf_l1s1=DecompCNParamsInst%rf_l1s1_cn
      rf_l2s2=DecompCNParamsInst%rf_l2s2_cn
      rf_l3s3=DecompCNParamsInst%rf_l3s3_cn
      rf_s1s2=DecompCNParamsInst%rf_s1s2_cn
      rf_s2s3=DecompCNParamsInst%rf_s2s3_cn
      rf_s3s4=DecompCNParamsInst%rf_s3s4_cn

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel=DecompCNParamsInst%cwd_fcel_cn
      cwd_flig=DecompCNParamsInst%cwd_flig_cn

      !-------------------  list of pools and their attributes  ------------

      i_litr1 = i_met_lit
      floating_cn_ratio_decomp_pools (i_litr1) = .true.
      floating_cp_ratio_decomp_pools (i_litr1) = .true.
      decomp_pool_name_restart       (i_litr1) = 'litr1'
      decomp_pool_name_history       (i_litr1) = 'LITR1'
      decomp_pool_name_long          (i_litr1) = 'litter 1'
      decomp_pool_name_short         (i_litr1) = 'L1'
      is_litter                      (i_litr1) = .true.
      is_soil                        (i_litr1) = .false.
      is_cwd                         (i_litr1) = .false.
      initial_cn_ratio               (i_litr1) = 90._r8
      initial_cp_ratio               (i_litr1) = 900._r8
      initial_stock                  (i_litr1) = 0._r8
      is_metabolic                   (i_litr1) = .true.
      is_cellulose                   (i_litr1) = .false.
      is_lignin                      (i_litr1) = .false.

      i_litr2 = i_cel_lit
      floating_cn_ratio_decomp_pools (i_litr2) = .true.
      floating_cp_ratio_decomp_pools (i_litr2) = .true.
      decomp_pool_name_restart       (i_litr2) = 'litr2'
      decomp_pool_name_history       (i_litr2) = 'LITR2'
      decomp_pool_name_long          (i_litr2) = 'litter 2'
      decomp_pool_name_short         (i_litr2) = 'L2'
      is_litter                      (i_litr2) = .true.
      is_soil                        (i_litr2) = .false.
      is_cwd                         (i_litr2) = .false.
      initial_cn_ratio               (i_litr2) = 90._r8
      initial_cp_ratio               (i_litr2) = 900._r8
      initial_stock                  (i_litr2) = 0._r8
      is_metabolic                   (i_litr2) = .false.
      is_cellulose                   (i_litr2) = .true.
      is_lignin                      (i_litr2) = .false.

      i_litr3 = i_lig_lit
      floating_cn_ratio_decomp_pools (i_litr3) = .true.
      floating_cp_ratio_decomp_pools (i_litr3) = .true.
      decomp_pool_name_restart       (i_litr3) = 'litr3'
      decomp_pool_name_history       (i_litr3) = 'LITR3'
      decomp_pool_name_long          (i_litr3) = 'litter 3'
      decomp_pool_name_short         (i_litr3) = 'L3'
      is_litter                      (i_litr3) = .true.
      is_soil                        (i_litr3) = .false.
      is_cwd                         (i_litr3) = .false.
      initial_cn_ratio               (i_litr3) = 90._r8
      initial_cp_ratio               (i_litr3) = 900._r8
      initial_stock                  (i_litr3) = 0._r8
      is_metabolic                   (i_litr3) = .false.
      is_cellulose                   (i_litr3) = .false.
      is_lignin                      (i_litr3) = .true.

      if (.not. use_fates) then
         floating_cn_ratio_decomp_pools (i_cwd) = .true.
         floating_cp_ratio_decomp_pools (i_cwd) = .true.
         decomp_pool_name_restart       (i_cwd) = 'cwd'
         decomp_pool_name_history       (i_cwd) = 'CWD'
         decomp_pool_name_long          (i_cwd) = 'coarse woody debris'
         decomp_pool_name_short         (i_cwd) = 'CWD'
         is_litter                      (i_cwd) = .false.
         is_soil                        (i_cwd) = .false.
         is_cwd                         (i_cwd) = .true.
         initial_cn_ratio               (i_cwd) = 500._r8
         initial_cp_ratio               (i_cwd) = 5000._r8
         initial_stock                  (i_cwd) = 0._r8
         is_metabolic                   (i_cwd) = .false.
         is_cellulose                   (i_cwd) = .false.
         is_lignin                      (i_cwd) = .false.
      end if

      if ( .not. use_fates ) then
         i_soil1 = 5
      else
         i_soil1 = 4
      endif
      floating_cn_ratio_decomp_pools (i_soil1) = .false.
      floating_cp_ratio_decomp_pools (i_soil1) = .true.
      decomp_pool_name_restart       (i_soil1) = 'soil1'
      decomp_pool_name_history       (i_soil1) = 'SOIL1'
      decomp_pool_name_long          (i_soil1) = 'soil 1'
      decomp_pool_name_short         (i_soil1) = 'S1'
      is_litter                      (i_soil1) = .false.
      is_soil                        (i_soil1) = .true.
      is_cwd                         (i_soil1) = .false.
      initial_cn_ratio               (i_soil1) = cn_s1
      initial_cp_ratio               (i_soil1) = cn_s1*np_s1_new 
      initial_stock                  (i_soil1) = 0._r8
      is_metabolic                   (i_soil1) = .false.
      is_cellulose                   (i_soil1) = .false.
      is_lignin                      (i_soil1) = .false.

      if ( .not. use_fates ) then
         i_soil2 = 6
      else
         i_soil2 = 5
      endif
      floating_cn_ratio_decomp_pools (i_soil2) = .false.
      floating_cp_ratio_decomp_pools (i_soil2) = .true.
      decomp_pool_name_restart       (i_soil2) = 'soil2'
      decomp_pool_name_history       (i_soil2) = 'SOIL2'
      decomp_pool_name_long          (i_soil2) = 'soil 2'
      decomp_pool_name_short         (i_soil2) = 'S2'
      is_litter                      (i_soil2) = .false.
      is_soil                        (i_soil2) = .true.
      is_cwd                         (i_soil2) = .false.
      initial_cn_ratio               (i_soil2) = cn_s2
      initial_cp_ratio               (i_soil2) = cn_s2*np_s2_new
      initial_stock                  (i_soil2) = 0._r8
      is_metabolic                   (i_soil2) = .false.
      is_cellulose                   (i_soil2) = .false.
      is_lignin                      (i_soil2) = .false.

      if ( .not. use_fates ) then
         i_soil3 = 7
      else
         i_soil3 = 6
      endif
      floating_cn_ratio_decomp_pools (i_soil3) = .false.
      floating_cp_ratio_decomp_pools (i_soil3) = .true.
      decomp_pool_name_restart       (i_soil3) = 'soil3'
      decomp_pool_name_history       (i_soil3) = 'SOIL3'
      decomp_pool_name_long          (i_soil3) = 'soil 3'
      decomp_pool_name_short         (i_soil3) = 'S3'
      is_litter                      (i_soil3) = .false.
      is_soil                        (i_soil3) = .true.
      is_cwd                         (i_soil3) = .false.
      initial_cn_ratio               (i_soil3) = cn_s3
      initial_cp_ratio               (i_soil3) = cn_s3*np_s3_new
      initial_stock                  (i_soil3) = 0._r8
      is_metabolic                   (i_soil3) = .false.
      is_cellulose                   (i_soil3) = .false.
      is_lignin                      (i_soil3) = .false.

      if ( .not. use_fates ) then
         i_soil4 = 8
      else
         i_soil4 = 7
      endif
      floating_cn_ratio_decomp_pools   (i_soil4) = .false.
      floating_cp_ratio_decomp_pools   (i_soil4) = .true.
      decomp_pool_name_restart         (i_soil4) = 'soil4'
      decomp_pool_name_history         (i_soil4) = 'SOIL4'
      decomp_pool_name_long            (i_soil4) = 'soil 4'
      decomp_pool_name_short           (i_soil4) = 'S4'
      is_litter                        (i_soil4) = .false.
      is_soil                          (i_soil4) = .true.
      is_cwd                           (i_soil4) = .false.
      initial_cn_ratio                 (i_soil4) = cn_s4
      initial_cp_ratio                 (i_soil4) = cn_s4*np_s4_new
      initial_stock                    (i_soil4) = 10._r8
      is_metabolic                     (i_soil4) = .false.
      is_cellulose                     (i_soil4) = .false.
      is_lignin                        (i_soil4) = .false.

      i_atm = 0  !! for terminal pools (i.e. 100% respiration)
      floating_cn_ratio_decomp_pools   (i_atm) = .false.
      floating_cp_ratio_decomp_pools   (i_atm) = .false.
      decomp_pool_name_restart         (i_atm) = 'atmosphere'
      decomp_pool_name_history         (i_atm) = 'atmosphere'
      decomp_pool_name_long            (i_atm) = 'atmosphere'
      decomp_pool_name_short           (i_atm) = ''
      is_litter                        (i_atm) = .true.
      is_soil                          (i_atm) = .false.
      is_cwd                           (i_atm) = .false.
      initial_cn_ratio                 (i_atm) = 0._r8
      initial_cp_ratio                 (i_atm) = 0._r8
      initial_stock                    (i_atm) = 0._r8
      is_metabolic                     (i_atm) = .false.
      is_cellulose                     (i_atm) = .false.
      is_lignin                        (i_atm) = .false.

      spinup_factor(i_litr1) = DecompCNParamsInst%spinup_vector(1)
      spinup_factor(i_litr2) = DecompCNParamsInst%spinup_vector(2)
      spinup_factor(i_litr3) = DecompCNParamsInst%spinup_vector(3)
      if (.not. use_fates) then
         spinup_factor(i_cwd) =   DecompCNParamsInst%spinup_vector(8)
      end if
      spinup_factor(i_soil1) = DecompCNParamsInst%spinup_vector(4)
      spinup_factor(i_soil2) = DecompCNParamsInst%spinup_vector(5)
      spinup_factor(i_soil3) = DecompCNParamsInst%spinup_vector(6)
      spinup_factor(i_soil4) = DecompCNParamsInst%spinup_vector(7)


      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      cascade_step_name(i_l1s1) = 'L1S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = rf_l1s1
      cascade_donor_pool(i_l1s1) = i_litr1
      cascade_receiver_pool(i_l1s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = 1.0_r8

      i_l2s2 = 2
      cascade_step_name(i_l2s2) = 'L2S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s2) = rf_l2s2
      cascade_donor_pool(i_l2s2) = i_litr2
      cascade_receiver_pool(i_l2s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s2) = 1.0_r8

      i_l3s3 = 3
      cascade_step_name(i_l3s3) = 'L3S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s3) = rf_l3s3
      cascade_donor_pool(i_l3s3) = i_litr3
      cascade_receiver_pool(i_l3s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s3) = 1.0_r8

      i_s1s2 = 4
      cascade_step_name(i_s1s2) = 'S1S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = rf_s1s2
      cascade_donor_pool(i_s1s2) = i_soil1
      cascade_receiver_pool(i_s1s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = 1.0_r8

      i_s2s3 = 5
      cascade_step_name(i_s2s3) = 'S2S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      cascade_donor_pool(i_s2s3) = i_soil2
      cascade_receiver_pool(i_s2s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = 1.0_r8

      i_s3s4 = 6 
      cascade_step_name(i_s3s4) = 'S3S4'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s4) = rf_s3s4
      cascade_donor_pool(i_s3s4) = i_soil3
      cascade_receiver_pool(i_s3s4) = i_soil4
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s4) = 1.0_r8

      i_s4atm = 7
      cascade_step_name(i_s4atm) = 'S4'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s4atm) = 1.
      cascade_donor_pool(i_s4atm) = i_soil4
      cascade_receiver_pool(i_s4atm) = i_atm
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s4atm) = 1.0_r8

      if (.not. use_fates) then
         i_cwdl2 = 8
         cascade_step_name(i_cwdl2) = 'CWDL2'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = 0._r8
         cascade_donor_pool(i_cwdl2) = i_cwd
         cascade_receiver_pool(i_cwdl2) = i_litr2
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fcel
         
         i_cwdl3 = 9
         cascade_step_name(i_cwdl3) = 'CWDL3'
         rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = 0._r8
         cascade_donor_pool(i_cwdl3) = i_cwd
         cascade_receiver_pool(i_cwdl3) = i_litr3
         pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = cwd_flig
      end if

    end associate

   end subroutine init_decompcascade_cn

   !-----------------------------------------------------------------------
   subroutine decomp_rate_constants_cn(bounds, &
        num_soilc, filter_soilc, &
        canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
     !
     ! !DESCRIPTION:
     ! calculate rate constants and decomposition pathways for the BGC model 
     ! originally implemented in CLM-CN  
     ! written by C. Koven based on original CLM4 decomposition cascade by P. Thornton
     !
     ! !USES:
     use clm_time_manager, only : get_step_size, get_nstep, get_curr_date
     use elm_varcon      , only : secspday
     use clm_varpar      , only : i_cwd
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds          
     integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
     integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
     type(canopystate_type) , intent(in)    :: canopystate_vars
     type(soilstate_type)   , intent(in)    :: soilstate_vars
     type(temperature_type) , intent(in)    :: temperature_vars 
     type(ch4_type)         , intent(in)    :: ch4_vars
     type(carbonflux_type)  , intent(inout) :: carbonflux_vars
     type(cnstate_type)     , intent(inout) :: cnstate_vars

     !
     ! !LOCAL VARIABLES:
     real(r8):: dt                           ! decomp timestep (seconds)   
     real(r8):: dtd                          ! decomp timestep (days)
     real(r8):: frw(bounds%begc:bounds%endc) ! rooting fraction weight
     real(r8), allocatable:: fr(:,:)         ! column-level rooting fraction by soil depth
     real(r8):: minpsi, maxpsi               ! limits for soil water scalar for decomp
     real(r8):: psi                          ! temporary soilpsi for water scalar
     real(r8):: rate_scalar                  ! combined rate scalar for decomp
     real(r8):: k_l1                         ! decomposition rate constant litter 1
     real(r8):: k_l2                         ! decomposition rate constant litter 2
     real(r8):: k_l3                         ! decomposition rate constant litter 3
     real(r8):: k_s1                         ! decomposition rate constant SOM 1
     real(r8):: k_s2                         ! decomposition rate constant SOM 2
     real(r8):: k_s3                         ! decomposition rate constant SOM 3
     real(r8):: k_s4                         ! decomposition rate constant SOM 4
     real(r8):: k_frag                       ! fragmentation rate constant CWD
     real(r8):: ck_l1                        ! corrected decomposition rate constant litter 1
     real(r8):: ck_l2                        ! corrected decomposition rate constant litter 2
     real(r8):: ck_l3                        ! corrected decomposition rate constant litter 3
     real(r8):: ck_s1                        ! corrected decomposition rate constant SOM 1
     real(r8):: ck_s2                        ! corrected decomposition rate constant SOM 2
     real(r8):: ck_s3                        ! corrected decomposition rate constant SOM 3
     real(r8):: ck_s4                        ! corrected decomposition rate constant SOM 4
     real(r8):: ck_frag                      ! corrected fragmentation rate constant CWD
     real(r8):: cwdc_loss                    ! fragmentation rate for CWD carbon (gC/m2/s)
     real(r8):: cwdn_loss                    ! fragmentation rate for CWD nitrogen (gN/m2/s)
     integer :: i_litr1
     integer :: i_litr2
     integer :: i_litr3
     integer :: i_soil1
     integer :: i_soil2
     integer :: i_soil3
     integer :: i_soil4
     integer :: c, fc, j, k, l
     real(r8):: Q10                          ! temperature dependence
     real(r8):: froz_q10                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
     real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
     real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
     real(r8) :: mino2lim                    ! minimum anaerobic decomposition rate as a
     integer :: year, mon, day, sec          ! fraction of potential aerobic rate
     !-----------------------------------------------------------------------

     associate(                                             &
          dz             => col_pp%dz                        , & ! Input:  [real(r8) (:,:)   ]  soil layer thickness (m)                               

          sucsat         => soilstate_vars%sucsat_col     , & ! Input:  [real(r8) (:,:)   ]  minimum soil suction (mm)                              
          soilpsi        => soilstate_vars%soilpsi_col    , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          

          alt_indx       => canopystate_vars%alt_indx_col , & ! Input:  [integer  (:)     ]  current depth of thaw                                     

          t_soisno       => col_es%t_soisno , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       

          o2stress_sat   => ch4_vars%o2stress_sat_col     , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
          o2stress_unsat => ch4_vars%o2stress_unsat_col   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
          finundated     => ch4_vars%finundated_col       , & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
          
          t_scalar       => col_cf%t_scalar  , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
          w_scalar       => col_cf%w_scalar  , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
          o_scalar       => col_cf%o_scalar  , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
          decomp_k       => col_cf%decomp_k  , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec) 
          decomp_k_pools => decomp_cascade_con%decomp_k_pools  & !(0: ndecomp_pools)    ! pflotran (0 for atm. co2)
          )

       mino2lim = ParamsShareInst%mino2lim

       ! set time steps
       dt = real( get_step_size(), r8 )
       dtd = dt/secspday

       ! set initial base rates for decomposition mass loss (1/day)
       ! (from Biome-BGC v4.2.0, using three SOM pools)
       ! Value inside log function is the discrete-time values for a
       ! daily time step model, and the result of the log function is
       ! the corresponding continuous-time decay rate (1/day), following
       ! Olson, 1963.
       k_l1=DecompCNParamsInst%k_l1_cn
       k_l2=DecompCNParamsInst%k_l2_cn
       k_l3=DecompCNParamsInst%k_l3_cn

       k_s1=DecompCNParamsInst%k_s1_cn
       k_s2=DecompCNParamsInst%k_s2_cn
       k_s3=DecompCNParamsInst%k_s3_cn
       k_s4=DecompCNParamsInst%k_s4_cn

       k_frag=DecompCNParamsInst%k_frag_cn

       ! calculate the new discrete-time decay rate for model timestep
       k_l1 = 1.0_r8-exp(-k_l1*dtd)
       k_l2 = 1.0_r8-exp(-k_l2*dtd)
       k_l3 = 1.0_r8-exp(-k_l3*dtd)

       k_s1 = 1.0_r8-exp(-k_s1*dtd)
       k_s2 = 1.0_r8-exp(-k_s2*dtd)
       k_s3 = 1.0_r8-exp(-k_s3*dtd)
       k_s4 = 1.0_r8-exp(-k_s4*dtd)

       k_frag = 1.0_r8-exp(-k_frag*dtd)

       minpsi=DecompCNParamsInst%minpsi_cn

       Q10 = ParamsShareInst%Q10_hr

       ! set "froz_q10" parameter
       froz_q10  = ParamsShareInst%froz_q10

       if (use_vertsoilc) then
          ! Set "decomp_depth_efolding" parameter
          decomp_depth_efolding = ParamsShareInst%decomp_depth_efolding
       end if

       i_litr1 = 1
       i_litr2 = 2
       i_litr3 = 3
       if (.not.use_fates) then
          i_soil1 = 5
          i_soil2 = 6
          i_soil3 = 7
          i_soil4 = 8
       else
          i_soil1 = 4
          i_soil2 = 5
          i_soil3 = 6
          i_soil4 = 7
       end if

       ! pflotran:beg---saving kd (NOT ad scaled) for passing to pflotran bgc decomposition sandboxes
       decomp_k_pools(i_litr1) = k_l1 / dt
       decomp_k_pools(i_litr2) = k_l2 / dt
       decomp_k_pools(i_litr3) = k_l3 / dt
       if (.not.use_fates) decomp_k_pools(i_cwd) = k_frag / dt
       decomp_k_pools(i_soil1) = k_s1 / dt
       decomp_k_pools(i_soil2) = k_s2 / dt
       decomp_k_pools(i_soil3) = k_s3 / dt
       decomp_k_pools(i_soil4) = k_s4 / dt
       ! pflotran:end

       ! The following code implements the acceleration part of the AD spinup
       ! algorithm, by multiplying all of the SOM decomposition base rates by 10.0.

       if ( spinup_state .eq. 1 ) then
          k_l1 = k_l1 * DecompCNParamsInst%spinup_vector(1)
          k_l2 = k_l2 * DecompCNParamsInst%spinup_vector(2)
          k_l3 = k_l3 * DecompCNParamsInst%spinup_vector(3)
          k_s1 = k_s1 * DecompCNParamsInst%spinup_vector(4)
          k_s2 = k_s2 * DecompCNParamsInst%spinup_vector(5)
          k_s3 = k_s3 * DecompCNParamsInst%spinup_vector(6)
          k_s4 = k_s4 * DecompCNParamsInst%spinup_vector(7)
          if (.not.use_fates) k_frag = k_frag * DecompCNParamsInst%spinup_vector(8)
       endif

       !--- time dependent coefficients-----!
       if ( nlevdecomp .eq. 1 ) then

          ! calculate function to weight the temperature and water potential scalars
          ! for decomposition control.  


          ! the following normalizes values in fr so that they
          ! sum to 1.0 across top nlevdecomp levels on a column
          frw(bounds%begc:bounds%endc) = 0._r8
          nlev_soildecomp_standard=5
          allocate(fr(bounds%begc:bounds%endc,nlev_soildecomp_standard))
          do j=1,nlev_soildecomp_standard
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                frw(c) = frw(c) + dz(c,j)
             end do
          end do
          do j = 1,nlev_soildecomp_standard
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                if (frw(c) /= 0._r8) then
                   fr(c,j) = dz(c,j) / frw(c)
                else
                   fr(c,j) = 0._r8
                end if
             end do
          end do

          ! calculate rate constant scalar for soil temperature
          ! assuming that the base rate constants are assigned for non-moisture
          ! limiting conditions at 25 C. 
          ! Peter Thornton: 3/13/09
          ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
          ! as part of the modifications made to improve the seasonal cycle of 
          ! atmospheric CO2 concentration in global simulations. This does not impact
          ! the base rates at 25 C, which are calibrated from microcosm studies.
          do j = 1,nlev_soildecomp_standard
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                if (j==1) t_scalar(c,:) = 0._r8
                !! use separate (possibly equal) t funcs above and below freezing point
                !! t_scalar(c,1)=t_scalar(c,1) + (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
                if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                   t_scalar(c,1)=t_scalar(c,1) + &
                        (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
                else
                   t_scalar(c,1)=t_scalar(c,1) + &
                        (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))*fr(c,j)
                endif
             end do
          end do

          ! calculate the rate constant scalar for soil water content.
          ! Uses the log relationship with water potential given in
          ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
          ! a comparison of models. Ecology, 68(5):1190-1200.
          ! and supported by data in
          ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
          ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

          do j = 1,nlev_soildecomp_standard
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                if (j==1) w_scalar(c,:) = 0._r8
                maxpsi = sucsat(c,j) * (-9.8e-6_r8)
                psi = min(soilpsi(c,j),maxpsi)
                ! decomp only if soilpsi is higher than minpsi
                if (psi > minpsi) then
                   w_scalar(c,1) = w_scalar(c,1) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
                end if
             end do
          end do

          if (use_lch4) then
             if (anoxia_wtsat) then ! Adjust for saturated fraction if unfrozen.
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   if (alt_indx(c) >= nlev_soildecomp_standard .and. t_soisno(c,1) > SHR_CONST_TKFRZ) then
                      w_scalar(c,1) = w_scalar(c,1)*(1._r8 - finundated(c)) + finundated(c)
                   end if
                end do
             end if
          end if

          if (use_lch4) then
             ! Calculate ANOXIA
             if (anoxia) then
                ! Check for anoxia w/o LCH4 now done in controlMod.

                do j = 1,nlev_soildecomp_standard
                   do fc = 1,num_soilc
                      c = filter_soilc(fc)

                      if (j==1) o_scalar(c,:) = 0._r8

                      if (.not. anoxia_wtsat) then
                         o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
                      else
                         o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * &
                              (max(o2stress_unsat(c,j), mino2lim)*(1._r8 - finundated(c)) + &
                              max(o2stress_sat(c,j), mino2lim)*finundated(c) )
                      end if
                   end do
                end do
             else
                o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
             end if
          else
             o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
          end if

          deallocate(fr)

       else
          
          ! calculate rate constant scalar for soil temperature
          ! assuming that the base rate constants are assigned for non-moisture
          ! limiting conditions at 25 C. 
          ! Peter Thornton: 3/13/09
          ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
          ! as part of the modifications made to improve the seasonal cycle of 
          ! atmospheric CO2 concentration in global simulations. This does not impact
          ! the base rates at 25 C, which are calibrated from microcosm studies.

          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                !! use separate (possibly equal) t funcs above and below freezing point
                !! t_scalar(c,j)= (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
                if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                   t_scalar(c,j)= (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
                else
                   t_scalar(c,j)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
                endif
             end do
          end do


          ! calculate the rate constant scalar for soil water content.
          ! Uses the log relationship with water potential given in
          ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
          ! a comparison of models. Ecology, 68(5):1190-1200.
          ! and supported by data in
          ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
          ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

          do j = 1,nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                maxpsi = sucsat(c,j) * (-9.8e-6_r8)
                psi = min(soilpsi(c,j),maxpsi)
                ! decomp only if soilpsi is higher than minpsi
                if (psi > minpsi) then
                   w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
                else
                   w_scalar(c,j) = 0._r8
                end if
                if (use_lch4) then
                   if (anoxia_wtsat .and. t_soisno(c,j) > SHR_CONST_TKFRZ) then ! wet area will have w_scalar of 1 if unfrozen
                      w_scalar(c,j) = w_scalar(c,j)*(1._r8 - finundated(c)) + finundated(c)
                   end if
                end if
             end do
          end do

       end if

       if (use_lch4) then
          ! Calculate ANOXIA
          ! Check for anoxia w/o LCH4 now done in controlMod.

          if (anoxia) then
             do j = 1,nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)

                   if (.not. anoxia_wtsat) then
                      o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
                   else
                      o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim) * (1._r8 - finundated(c)) + &
                           max(o2stress_sat(c,j), mino2lim) * finundated(c)
                   end if
                end do
             end do
          else
             o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
          end if
       else
          o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
       end if

       if (use_vertsoilc) then
          ! add a term to reduce decomposition rate at depth
          ! for now used a fixed e-folding depth
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding)
             end do
          end do
       end if

       call get_curr_date(year, mon, day, sec)
       !Calcluate location and depth-specific acceleration factors
       do fc=1,num_soilc
           c = filter_soilc(fc)
           if (year < 20 .and. spinup_state == 1) then 
               cnstate_vars%scalaravg_col(c,:) = 0._r8
           else if (year < 40 .and. spinup_state == 1) then  
               cnstate_vars%scalaravg_col(c,:) = cnstate_vars%scalaravg_col(c,:) + &
                     (t_scalar(c,4) * w_scalar(c,4) * o_scalar(c,4) * depth_scalar(c,4) ) &
                     * dt / (86400._r8 * 365._r8 * 20._r8)
           else
               if (cnstate_vars%scalaravg_col(c,4) < 1.0e-3) then 
                    cnstate_vars%scalaravg_col(c,:) = 1.0_r8
               end if
           end if
       end do
  
       if (use_vertsoilc) then
          do j = 1,nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
             end do
          end do
          if (.not.use_fates) then
             do j = 1,nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
                end do
             end do
          end if
       else
          do j = 1,nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
             end do
          end do
          if (.not.use_fates) then
             do j = 1,nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
                end do
             end do
          end if
       end if
       

       if (spinup_state == 1 .and. year >= 40) then 
         !adjust decomposition factors based on scalar factors from first 20 years of simulation
         do j=1,nlevdecomp
           do fc = 1, num_soilc
             c = filter_soilc(fc)
             if ( decomp_cascade_con%spinup_factor(i_litr1) > 1._r8) decomp_k(c,j,i_litr1) = decomp_k(c,j,i_litr1)  &
               / cnstate_vars%scalaravg_col(c,j) 
             if ( decomp_cascade_con%spinup_factor(i_litr2) > 1._r8) decomp_k(c,j,i_litr2) = decomp_k(c,j,i_litr2)  &
               / cnstate_vars%scalaravg_col(c,j) 
             if ( decomp_cascade_con%spinup_factor(i_litr3) > 1._r8) decomp_k(c,j,i_litr3) = decomp_k(c,j,i_litr3)  &
                  / cnstate_vars%scalaravg_col(c,j) 
             if ( .not. use_fates ) then
                if ( decomp_cascade_con%spinup_factor(i_cwd)   > 1._r8) decomp_k(c,j,i_cwd)   = decomp_k(c,j,i_cwd)    &
                     / cnstate_vars%scalaravg_col(c,j) 
             endif
             if ( decomp_cascade_con%spinup_factor(i_soil1) > 1._r8) decomp_k(c,j,i_soil1) = decomp_k(c,j,i_soil1)  &
               / cnstate_vars%scalaravg_col(c,j) 
             if ( decomp_cascade_con%spinup_factor(i_soil2) > 1._r8) decomp_k(c,j,i_soil2) = decomp_k(c,j,i_soil2)  &
               / cnstate_vars%scalaravg_col(c,j) 
             if ( decomp_cascade_con%spinup_factor(i_soil3) > 1._r8) decomp_k(c,j,i_soil3) = decomp_k(c,j,i_soil3)  &
               / cnstate_vars%scalaravg_col(c,j) 
             if ( decomp_cascade_con%spinup_factor(i_soil4) > 1._r8) decomp_k(c,j,i_soil4) = decomp_k(c,j,i_soil4)  &
               / cnstate_vars%scalaravg_col(c,j) 
           end do
         end do
       end if    
     end associate
   end subroutine decomp_rate_constants_cn

end module DecompCascadeCNMod
