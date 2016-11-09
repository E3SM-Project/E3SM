module CNDecompCascadeBGCMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Sets the coeffiecients used in the decomposition cascade submodel.  
  ! This uses the CENTURY/BGC parameters
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_TKFRZ
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : nlevsoi, nlevgrnd, nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar             , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl             , only : iulog, spinup_state, anoxia, use_lch4, use_vertsoilc
  use clm_varcon             , only : zsoi
  use decompMod              , only : bounds_type
  use abortutils             , only : endrun
  use CNSharedParamsMod      , only : CNParamsShareInst, anoxia_wtsat, nlev_soildecomp_standard 
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use SoilStateType          , only : soilstate_type
  use CanopyStateType        , only : canopystate_type
  use TemperatureType        , only : temperature_type 
  use ch4Mod                 , only : ch4_type
  use ColumnType             , only : col                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: init_decompcascade_bgc
  public :: readCNDecompBgcParams
  public :: decomp_rate_constants_bgc
  !
  ! !PUBLIC DATA MEMBERS 
  logical , public :: normalize_q10_to_century_tfunc = .true.! do we normalize the century decomp. rates so that they match the CLM Q10 at a given tep?
  logical , public :: use_century_tfunc = .false.
  real(r8), public :: normalization_tref = 15._r8            ! reference temperature for normalizaion (degrees C)
  !
  ! !PRIVATE DATA MEMBERS 
  type, private :: CNDecompBgcParamsType
     real(r8):: cn_s1_bgc     !C:N for SOM 1
     real(r8):: cn_s2_bgc     !C:N for SOM 2
     real(r8):: cn_s3_bgc     !C:N for SOM 3
     
     real(r8):: np_s1_new_bgc  !C:P for SOM 1
     real(r8):: np_s2_new_bgc  !C:P for SOM 2
     real(r8):: np_s3_new_bgc  !C:P for SOM 3

     real(r8):: cp_s1_new_bgc        !C:P for SOM 1
     real(r8):: cp_s2_new_bgc        !C:P for SOM 2
     real(r8):: cp_s3_new_bgc        !C:P for SOM 3



     real(r8):: rf_l1s1_bgc   !respiration fraction litter 1 -> SOM 1
     real(r8):: rf_l2s1_bgc
     real(r8):: rf_l3s2_bgc

     real(r8):: rf_s2s1_bgc    
     real(r8):: rf_s2s3_bgc    
     real(r8):: rf_s3s1_bgc    

     real(r8):: rf_cwdl2_bgc 
     real(r8):: rf_cwdl3_bgc

     real(r8):: tau_l1_bgc    ! turnover time of  litter 1 (yr)
     real(r8):: tau_l2_l3_bgc ! turnover time of  litter 2 and litter 3 (yr)
     real(r8):: tau_s1_bgc    ! turnover time of  SOM 1 (yr)
     real(r8):: tau_s2_bgc    ! turnover time of  SOM 2 (yr)
     real(r8):: tau_s3_bgc    ! turnover time of  SOM 3 (yr)
     real(r8):: tau_cwd_bgc   ! corrected fragmentation rate constant CWD

     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD
     real(r8) :: cwd_flig_bgc !

     real(r8) :: k_frag_bgc   !fragmentation rate for CWD
     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     
     integer  :: nsompools = 3
     real(r8),allocatable :: spinup_vector(:) ! multipliers for soil decomp during accelerated spinup

  end type CNDecompBgcParamsType

  type(CNDecompBgcParamsType),private ::  CNDecompBgcParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNDecompBgcParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompBgcParamsType'
    character(len=100) :: errCode = 'Error reading in CN const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! These are not read off of netcdf file
    allocate(CNDecompBgcParamsInst%spinup_vector(CNDecompBgcParamsInst%nsompools))
    CNDecompBgcParamsInst%spinup_vector(:) = (/ 1.0_r8, 15.0_r8, 675.0_r8 /)

    ! Read off of netcdf file
    tString='tau_l1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_l1_bgc=tempr

    tString='tau_l2_l3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_l2_l3_bgc=tempr

    tString='tau_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_s1_bgc=tempr

    tString='tau_s2'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_s2_bgc=tempr

    tString='tau_s3'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_s3_bgc=tempr

    tString='tau_cwd'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%tau_cwd_bgc=tempr

    tString='cn_s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%cn_s1_bgc=tempr

    tString='cn_s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%cn_s2_bgc=tempr

    tString='cn_s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%cn_s3_bgc=tempr

!!! read in phosphorus variables - note that these NP ratio parameters for BGC  will have
!!! to be added in the parameter file
 
    tString='np_s1_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%np_s1_new_bgc=tempr

    tString='np_s2_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%np_s2_new_bgc=tempr

    tString='np_s3_new'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%np_s3_new_bgc=tempr

    tString='rf_l1s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_l1s1_bgc=tempr

    tString='rf_l2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_l2s1_bgc=tempr

    tString='rf_l3s2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_l3s2_bgc=tempr   

    tString='rf_s2s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_s2s1_bgc=tempr

    tString='rf_s2s3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_s2s3_bgc=tempr

    tString='rf_s3s1_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_s3s1_bgc=tempr

    tString='rf_cwdl2_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_cwdl2_bgc=tempr

    tString='rf_cwdl3_bgc'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%rf_cwdl3_bgc=tempr

    tString='cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%cwd_fcel_bgc=tempr

    tString='k_frag'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%k_frag_bgc=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%minpsi_bgc=tempr 

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompBgcParamsInst%cwd_flig_bgc=tempr 

  end subroutine readCNDecompBgcParams

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_bgc(bounds, cnstate_vars, soilstate_vars)
    !
    ! !DESCRIPTION:
    !  initialize rate constants and decomposition pathways following the decomposition cascade of the BGC model.
    !  written by C. Koven 
    !
    ! !USES:
    use clm_time_manager , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds  
    type(cnstate_type)   , intent(inout) :: cnstate_vars
    type(soilstate_type) , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES
    !-- properties of each decomposing pool
    real(r8) :: rf_l1s1
    real(r8) :: rf_l2s1
    real(r8) :: rf_l3s2
    !real(r8) :: rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: rf_s1s2(:,:)
    real(r8), allocatable :: rf_s1s3(:,:)
    real(r8) :: rf_s2s1
    real(r8) :: rf_s2s3
    real(r8) :: rf_s3s1
    real(r8) :: rf_cwdl2
    real(r8) :: rf_cwdl3
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    real(r8) :: np_s1_new
    real(r8) :: np_s2_new
    real(r8) :: np_s3_new
    !real(r8) :: f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
    !real(r8) :: f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8), allocatable :: f_s1s2(:,:)
    real(r8), allocatable :: f_s1s3(:,:)
    real(r8) :: f_s2s1
    real(r8) :: f_s2s3

    integer :: i_litr1
    integer :: i_litr2
    integer :: i_litr3
    integer :: i_soil1
    integer :: i_soil2
    integer :: i_soil3
    integer :: i_l1s1
    integer :: i_l2s1
    integer :: i_l3s2
    integer :: i_s1s2
    integer :: i_s1s3
    integer :: i_s2s1
    integer :: i_s2s3
    integer :: i_s3s1
    integer :: i_cwdl2
    integer :: i_cwdl3

    integer  :: c, j    ! indices
    real(r8) :: t       ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                 &
         rf_decomp_cascade              => cnstate_vars%rf_decomp_cascade_col                , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)       
         pathfrac_decomp_cascade        => cnstate_vars%pathfrac_decomp_cascade_col          , & ! Input:  [real(r8)          (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)
         cellsand                       => soilstate_vars%cellsand_col                       , & ! Input:  [real(r8)          (:,:)   ]  column 3D sand                                         
         
         cascade_step_name              => decomp_cascade_con%cascade_step_name              , & ! Output: [character(len=8)  (:)     ]  name of transition                               
         cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool             , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step 
         cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool          , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step   
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
         floating_cp_ratio_decomp_pools => decomp_cascade_con%floating_cp_ratio_decomp_pools , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:P ratio                          
         decomp_pool_name_restart       => decomp_cascade_con%decomp_pool_name_restart       , & ! Output: [character(len=8)  (:)     ]  name of pool for restart files                   
         decomp_pool_name_history       => decomp_cascade_con%decomp_pool_name_history       , & ! Output: [character(len=8)  (:)     ]  name of pool for history files                   
         decomp_pool_name_long          => decomp_cascade_con%decomp_pool_name_long          , & ! Output: [character(len=20) (:)     ]  name of pool for netcdf long names              
         decomp_pool_name_short         => decomp_cascade_con%decomp_pool_name_short         , & ! Output: [character(len=8)  (:)     ]  name of pool for netcdf short names              
         is_litter                      => decomp_cascade_con%is_litter                      , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool                             
         is_soil                        => decomp_cascade_con%is_soil                        , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool                               
         is_cwd                         => decomp_cascade_con%is_cwd                         , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool                                
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio               , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
         initial_cp_ratio               => decomp_cascade_con%initial_cp_ratio               , & ! Output: [real(r8)          (:)     ]  c:p ratio for initialization of pools                    
         initial_stock                  => decomp_cascade_con%initial_stock                  , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup              
         is_metabolic                   => decomp_cascade_con%is_metabolic                   , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material                        
         is_cellulose                   => decomp_cascade_con%is_cellulose                   , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose                                 
         is_lignin                      => decomp_cascade_con%is_lignin                      , & ! Output: [logical           (:)     ]  TRUE => pool is lignin                                    
         spinup_factor                  => decomp_cascade_con%spinup_factor                    & ! Output: [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
         )

      allocate(rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp))
      allocate(f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp))

      !------- time-constant coefficients ---------- !
      ! set soil organic matter compartment C:N ratios
      cn_s1 = CNDecompBgcParamsInst%cn_s1_bgc
      cn_s2 = CNDecompBgcParamsInst%cn_s2_bgc
      cn_s3 = CNDecompBgcParamsInst%cn_s3_bgc

      ! set soil organic matter C:P ratios -X. YANG 
      np_s1_new=CNDecompBgcParamsInst%np_s1_new_bgc
      np_s2_new=CNDecompBgcParamsInst%np_s2_new_bgc
      np_s3_new=CNDecompBgcParamsInst%np_s3_new_bgc

      ! set respiration fractions for fluxes between compartments
      rf_l1s1 = CNDecompBgcParamsInst%rf_l1s1_bgc
      rf_l2s1 = CNDecompBgcParamsInst%rf_l2s1_bgc
      rf_l3s2 = CNDecompBgcParamsInst%rf_l3s2_bgc
      rf_s2s1 = CNDecompBgcParamsInst%rf_s2s1_bgc
      rf_s2s3 = CNDecompBgcParamsInst%rf_s2s3_bgc
      rf_s3s1 = CNDecompBgcParamsInst%rf_s3s1_bgc

      rf_cwdl2 = CNDecompBgcParamsInst%rf_cwdl2_bgc
      rf_cwdl3 = CNDecompBgcParamsInst%rf_cwdl3_bgc

      ! set the cellulose and lignin fractions for coarse woody debris
      cwd_fcel = CNDecompBgcParamsInst%cwd_fcel_bgc
      cwd_flig = CNDecompBgcParamsInst%cwd_flig_bgc

        ! set path fractions
      f_s2s1 = 0.42_r8/(0.45_r8)
      f_s2s3 = 0.03_r8/(0.45_r8)

      ! some of these are dependent on the soil texture properties
      do c = bounds%begc, bounds%endc
         do j = 1, nlevdecomp
            t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
            f_s1s2(c,j) = 1._r8 - .004_r8 / (1._r8 - t)
            f_s1s3(c,j) = .004_r8 / (1._r8 - t)
            rf_s1s2(c,j) = t
            rf_s1s3(c,j) = t
         end do
      end do

      !-------------------  list of pools and their attributes  ------------
      i_litr1                                  = i_met_lit
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

      i_litr2                                  = i_cel_lit
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

      i_litr3                                  = i_lig_lit
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

      ! CWD
      floating_cn_ratio_decomp_pools (i_cwd)   = .true.
      floating_cp_ratio_decomp_pools (i_cwd)   = .true.
      decomp_pool_name_restart       (i_cwd)   = 'cwd'
      decomp_pool_name_history       (i_cwd)   = 'CWD'
      decomp_pool_name_long          (i_cwd)   = 'coarse woody debris'
      decomp_pool_name_short         (i_cwd)   = 'CWD'
      is_litter                      (i_cwd)   = .false.
      is_soil                        (i_cwd)   = .false.
      is_cwd                         (i_cwd)   = .true.
      initial_cn_ratio               (i_cwd)   = 90._r8
      initial_cp_ratio               (i_cwd)   = 900._r8
      initial_stock                  (i_cwd)   = 0._r8
      is_metabolic                   (i_cwd)   = .false.
      is_cellulose                   (i_cwd)   = .false.
      is_lignin                      (i_cwd)   = .false.

      i_soil1                                  = 5
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
      initial_stock                  (i_soil1) = 20._r8
      is_metabolic                   (i_soil1) = .false.
      is_cellulose                   (i_soil1) = .false.
      is_lignin                      (i_soil1) = .false.

      i_soil2                                  = 6
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
      initial_stock                  (i_soil2) = 20._r8
      is_metabolic                   (i_soil2) = .false.
      is_cellulose                   (i_soil2) = .false.
      is_lignin                      (i_soil2) = .false.

      i_soil3                                  = 7
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
      initial_stock                  (i_soil3) = 20._r8
      is_metabolic                   (i_soil3) = .false.
      is_cellulose                   (i_soil3) = .false.
      is_lignin                      (i_soil3) = .false.

      spinup_factor(i_litr1) = 1._r8
      spinup_factor(i_litr2) = 1._r8
      spinup_factor(i_litr3) = 1._r8
      spinup_factor(i_cwd) = 1._r8
      spinup_factor(i_soil1) = CNDecompBgcParamsInst%spinup_vector(1)
      spinup_factor(i_soil2) = CNDecompBgcParamsInst%spinup_vector(2)
      spinup_factor(i_soil3) = CNDecompBgcParamsInst%spinup_vector(3)


      !----------------  list of transitions and their time-independent coefficients  ---------------!
      i_l1s1 = 1
      cascade_step_name(i_l1s1) = 'L1S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = rf_l1s1
      cascade_donor_pool(i_l1s1) = i_litr1
      cascade_receiver_pool(i_l1s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1s1) = 1.0_r8

      i_l2s1 = 2
      cascade_step_name(i_l2s1) = 'L2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1) = rf_l2s1
      cascade_donor_pool(i_l2s1) = i_litr2
      cascade_receiver_pool(i_l2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2s1)= 1.0_r8

      i_l3s2 = 3
      cascade_step_name(i_l3s2) = 'L3S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = rf_l3s2
      cascade_donor_pool(i_l3s2) = i_litr3
      cascade_receiver_pool(i_l3s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l3s2) = 1.0_r8

      i_s1s2 = 4
      cascade_step_name(i_s1s2) = 'S1S2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = rf_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s2) = i_soil1
      cascade_receiver_pool(i_s1s2) = i_soil2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s2) = f_s1s2(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s1s3 = 5
      cascade_step_name(i_s1s3) = 'S1S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = rf_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)
      cascade_donor_pool(i_s1s3) = i_soil1
      cascade_receiver_pool(i_s1s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1s3) = f_s1s3(bounds%begc:bounds%endc,1:nlevdecomp)

      i_s2s1 = 6
      cascade_step_name(i_s2s1) = 'S2S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = rf_s2s1
      cascade_donor_pool(i_s2s1) = i_soil2
      cascade_receiver_pool(i_s2s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = f_s2s1

      i_s2s3 = 7 
      cascade_step_name(i_s2s3) = 'S2S3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = rf_s2s3
      cascade_donor_pool(i_s2s3) = i_soil2
      cascade_receiver_pool(i_s2s3) = i_soil3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s3) = f_s2s3

      i_s3s1 = 8
      cascade_step_name(i_s3s1) = 'S3S1'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = rf_s3s1
      cascade_donor_pool(i_s3s1) = i_soil3
      cascade_receiver_pool(i_s3s1) = i_soil1
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8

      i_cwdl2 = 9
      cascade_step_name(i_cwdl2) = 'CWDL2'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
      cascade_donor_pool(i_cwdl2) = i_cwd
      cascade_receiver_pool(i_cwdl2) = i_litr2
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = cwd_fcel

      i_cwdl3 = 10
      cascade_step_name(i_cwdl3) = 'CWDL3'
      rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
      cascade_donor_pool(i_cwdl3) = i_cwd
      cascade_receiver_pool(i_cwdl3) = i_litr3
      pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl3) = cwd_flig
      
      deallocate(rf_s1s2)
      deallocate(rf_s1s3)
      deallocate(f_s1s2)
      deallocate(f_s1s3)

    end associate

  end subroutine init_decompcascade_bgc

  !-----------------------------------------------------------------------
  subroutine decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
       canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    !  calculate rate constants and decomposition pathways for teh CENTURY decomposition cascade model
    !  written by C. Koven based on original CLM4 decomposition cascade
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
    use shr_const_mod    , only : SHR_CONST_PI
    use clm_varcon       , only : secspday
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
    !
    ! !LOCAL VARIABLES:
    real(r8):: frw(bounds%begc:bounds%endc) ! rooting fraction weight
    real(r8), allocatable:: fr(:,:)         ! column-level rooting fraction by soil depth
    real(r8):: minpsi, maxpsi               ! limits for soil water scalar for decomp
    real(r8):: psi                          ! temporary soilpsi for water scalar
    real(r8):: rate_scalar                  ! combined rate scalar for decomp
    real(r8):: k_l1                         ! decomposition rate constant litter 1 (1/sec)
    real(r8):: k_l2_l3                      ! decomposition rate constant litter 2 and litter 3 (1/sec)
    real(r8):: k_s1                         ! decomposition rate constant SOM 1 (1/sec)
    real(r8):: k_s2                         ! decomposition rate constant SOM 2 (1/sec)
    real(r8):: k_s3                         ! decomposition rate constant SOM 3 (1/sec)
    real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
    real(r8):: tau_l1                       ! turnover time of  litter 1 (yr)
    real(r8):: tau_l2_l3                    ! turnover time of  litter 2 and litter 3 (yr)
    real(r8):: tau_l3                       ! turnover time of  litter 3 (yr)
    real(r8):: tau_s1                       ! turnover time of  SOM 1 (yr)
    real(r8):: tau_s2                       ! turnover time of  SOM 2 (yr)
    real(r8):: tau_s3                       ! turnover time of  SOM 3 (yr)
    real(r8):: tau_cwd                      ! corrected fragmentation rate constant CWD
    real(r8):: cwdc_loss                    ! fragmentation rate for CWD carbon (gC/m2/s)
    real(r8):: cwdn_loss                    ! fragmentation rate for CWD nitrogen (gN/m2/s)
    real(r8):: Q10                          ! temperature dependence
    real(r8):: froz_q10                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
    integer :: i_litr1
    integer :: i_litr2
    integer :: i_litr3
    integer :: i_soil1
    integer :: i_soil2
    integer :: i_soil3
    integer :: c, fc, j, k, l
    real(r8):: catanf                       ! hyperbolic temperature function from CENTURY
    real(r8):: catanf_30                    ! reference rate at 30C
    real(r8):: t1                           ! temperature argument
    real(r8):: normalization_factor         ! factor by which to offset the decomposition rates frm century to a q10 formulation
    real(r8):: days_per_year                ! days per year
    real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
    real(r8):: mino2lim                     !minimum anaerobic decomposition rate
    !-----------------------------------------------------------------------

    !----- CENTURY T response function
    catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

    associate(                                             &
         sucsat         => soilstate_vars%sucsat_col     , & ! Input:  [real(r8) (:,:)   ]  minimum soil suction (mm)                              
         soilpsi        => soilstate_vars%soilpsi_col    , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          

         alt_indx       => canopystate_vars%alt_indx_col , & ! Input:  [integer  (:)     ]  current depth of thaw                                     

         t_soisno       => temperature_vars%t_soisno_col , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       

         o2stress_sat   => ch4_vars%o2stress_sat_col     , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         o2stress_unsat => ch4_vars%o2stress_unsat_col   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
         finundated     => ch4_vars%finundated_col       , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
         
         t_scalar       => carbonflux_vars%t_scalar_col  , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
         w_scalar       => carbonflux_vars%w_scalar_col  , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
         o_scalar       => carbonflux_vars%o_scalar_col  , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
         decomp_k       => carbonflux_vars%decomp_k_col    & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)             
         )

      mino2lim = CNParamsShareInst%mino2lim

      if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
         call endrun(msg='ERROR: cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true'//&
              errMsg(__FILE__, __LINE__))
      endif

      days_per_year = get_days_per_year()

      ! the belowground parameters from century
      tau_l1 = 1./18.5
      tau_l2_l3 = 1./4.9
      tau_s1 = 1./7.3
      tau_s2 = 1./0.2
      tau_s3 = 1./.0045

      ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
      tau_cwd  = 1./0.3

      ! Todo:  FIX(SPM,032414) - the explicit divide gives different results than when that
      ! value is placed in the parameters netcdf file.  To get bfb, keep the 
      ! divide in source.

      !tau_l1 = CNDecompBgcParamsInst%tau_l1_bgc
      !tau_l2_l3 = CNDecompBgcParamsInst%tau_l2_l3_bgc
      !tau_s1 = CNDecompBgcParamsInst%tau_s1_bgc
      !tau_s2 = CNDecompBgcParamsInst%tau_s2_bgc
      !tau_s3 = CNDecompBgcParamsInst%tau_s3_bgc

      !set turnover rate of coarse woody debris
      !tau_cwd = CNDecompBgcParamsInst%tau_cwd_bgc

      ! set "Q10" parameter
      Q10 = CNParamsShareInst%Q10

      ! set "froz_q10" parameter
      froz_q10  = CNParamsShareInst%froz_q10 

      ! Set "decomp_depth_efolding" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      ! translate to per-second time constant
      k_l1 = 1._r8    / (secspday * days_per_year * tau_l1)
      k_l2_l3 = 1._r8 / (secspday * days_per_year * tau_l2_l3)
      k_s1 = 1._r8    / (secspday * days_per_year * tau_s1)
      k_s2 = 1._r8    / (secspday * days_per_year * tau_s2)
      k_s3 = 1._r8    / (secspday * days_per_year * tau_s3)
      k_frag = 1._r8  / (secspday * days_per_year * tau_cwd)

      ! calc ref rate
      catanf_30 = catanf(30._r8)
      ! The following code implements the acceleration part of the AD spinup algorithm

      if ( spinup_state .eq. 1 ) then
         k_s1 = k_s1 * CNDecompBgcParamsInst%spinup_vector(1)
         k_s2 = k_s2 * CNDecompBgcParamsInst%spinup_vector(2)
         k_s3 = k_s3 * CNDecompBgcParamsInst%spinup_vector(3)
      endif

      i_litr1 = 1
      i_litr2 = 2
      i_litr3 = 3
      i_soil1 = 5
      i_soil2 = 6
      i_soil3 = 7


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
               frw(c) = frw(c) + col%dz(c,j)
            end do
         end do
         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (frw(c) /= 0._r8) then
                  fr(c,j) = col%dz(c,j) / frw(c)
               else
                  fr(c,j) = 0._r8
               end if
            end do
         end do

         if ( .not. use_century_tfunc ) then
            ! calculate rate constant scalar for soil temperature
            ! assuming that the base rate constants are assigned for non-moisture
            ! limiting conditions at 25 C. 

            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (j==1) t_scalar(c,:) = 0._r8
                  if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c,1)=t_scalar(c,1) + &
                          (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
                  else
                     t_scalar(c,1)=t_scalar(c,1) + &
                          (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))*fr(c,j)
                  endif
               end do
            end do

         else
            ! original century uses an arctangent function to calculate the temperature dependence of decomposition
            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (j==1) t_scalar(c,:) = 0._r8

                  t_scalar(c,1)=t_scalar(c,1) +max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30*fr(c,j),0.01_r8)
               end do
            end do

         endif

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

         minpsi = -10.0_r8;

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
            if (anoxia_wtsat) then ! Adjust for saturated fraction if unfrozen
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

         if ( .not. use_century_tfunc ) then
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
                  if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
                     t_scalar(c,j)= (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
                  else
                     t_scalar(c,j)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
                  endif
               end do
            end do

         else

            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  t_scalar(c,j)= max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30, 0.01_r8)
               end do
            end do

         endif

         ! calculate the rate constant scalar for soil water content.
         ! Uses the log relationship with water potential given in
         ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
         ! a comparison of models. Ecology, 68(5):1190-1200.
         ! and supported by data in
         ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
         ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

         minpsi = -10.0_r8;
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

      end if

      if ( normalize_q10_to_century_tfunc ) then
         ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
         normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               t_scalar(c,j) = t_scalar(c,j) * normalization_factor
            end do
         end do
      endif

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

      if (use_vertsoilc) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               decomp_k(c,j,i_litr1) = k_l1    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_cwd)   = k_frag  * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil1) = k_s1    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil2) = k_s2    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil3) = k_s3    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
            end do
         end do
      else
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               decomp_k(c,j,i_litr1) = k_l1    * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_cwd)   = k_frag  * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil1) = k_s1    * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil2) = k_s2    * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
               decomp_k(c,j,i_soil3) = k_s3    * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
            end do
         end do
      end if

    end associate

 end subroutine decomp_rate_constants_bgc

end module CNDecompCascadeBGCMod
