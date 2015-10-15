module BGCCenturyParMod
  !
  ! !DESCRIPTION:
  !  parameterization module for century bgc
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  implicit none

  public :: readCentDecompBgcParams
  public :: readCentNitrifDenitrifParams
  public :: readCentCNAllocParams

  type, private :: CNNitrifDenitrifParamsType
     real(r8) :: k_nitr_max               !  maximum nitrification rate constant (1/s)
     real(r8) :: surface_tension_water    !  surface tension of water(J/m^2), Arah an and Vinten 1995
     real(r8) :: rij_kro_a                !  Arah and Vinten 1995)
     real(r8) :: rij_kro_alpha            !  parameter to calculate anoxic fraction of soil  (Arah and Vinten 1995)
     real(r8) :: rij_kro_beta             !  (Arah and Vinten 1995)
     real(r8) :: rij_kro_gamma            !  (Arah and Vinten 1995)
     real(r8) :: rij_kro_delta            !  (Arah and Vinten 1995)
  end type CNNitrifDenitrifParamsType

  type(CNNitrifDenitrifParamsType), protected ::  CNNitrifDenitrifParamsInst


  type :: NutrientCompetitionParamsType

     real(r8) :: dayscrecover      ! number of days to recover negative cpool
     real(r8) :: compet_plant_no3  ! (unitless) relative compettiveness of plants for NO3
     real(r8) :: compet_plant_nh4  ! (unitless) relative compettiveness of plants for NH4
     real(r8) :: compet_decomp_no3 ! (unitless) relative competitiveness of immobilizers for NO3
     real(r8) :: compet_decomp_nh4 ! (unitless) relative competitiveness of immobilizers for NH4
     real(r8) :: compet_denit      ! (unitless) relative competitiveness of denitrifiers for NO3
     real(r8) :: compet_nit        ! (unitless) relative competitiveness of nitrifiers for NH4
  end type NutrientCompetitionParamsType

  ! NutrientCompetitionParamsInst is populated in readCNAllocParams which is called in
  type(NutrientCompetitionParamsType),protected ::  NutrientCompetitionParamsInst


  type, private :: CNDecompBgcParamsType
     real(r8)             :: cn_s1_bgc        !C:N for SOM 1
     real(r8)             :: cn_s2_bgc        !C:N for SOM 2
     real(r8)             :: cn_s3_bgc        !C:N for SOM 3

     real(r8)             :: rf_l1s1_bgc      !respiration fraction litter 1 -> SOM 1
     real(r8)             :: rf_l2s1_bgc
     real(r8)             :: rf_l3s2_bgc

     real(r8)             :: rf_s2s1_bgc
     real(r8)             :: rf_s2s3_bgc
     real(r8)             :: rf_s3s1_bgc

     real(r8)             :: rf_cwdl2_bgc
     real(r8)             :: rf_cwdl3_bgc

     real(r8)             :: tau_l1_bgc       ! turnover time of  litter 1 (yr)
     real(r8)             :: tau_l2_l3_bgc    ! turnover time of  litter 2 and litter 3 (yr)
     real(r8)             :: tau_s1_bgc       ! turnover time of  SOM 1 (yr)
     real(r8)             :: tau_s2_bgc       ! turnover time of  SOM 2 (yr)
     real(r8)             :: tau_s3_bgc       ! turnover time of  SOM 3 (yr)
     real(r8)             :: tau_cwd_bgc      ! corrected fragmentation rate constant CWD

     real(r8)             :: k_decay_lit1
     real(r8)             :: k_decay_lit2
     real(r8)             :: k_decay_lit3
     real(r8)             :: k_decay_som1
     real(r8)             :: k_decay_som2
     real(r8)             :: k_decay_som3
     real(r8)             :: k_decay_cwd

     real(r8)             :: cwd_fcel_bgc     !cellulose fraction for CWD
     real(r8)             :: cwd_flig_bgc     !

     real(r8)             :: k_frag_bgc       !fragmentation rate for CWD
     real(r8)             :: minpsi_bgc       !minimum soil water potential for heterotrophic resp

     integer              :: nsompools = 3

     real(r8),allocatable :: spinup_vector(:) ! multipliers for soil decomp during accelerated spinup

  end type CNDecompBgcParamsType

  type(CNDecompBgcParamsType),protected ::  CNDecompBgcParamsInst


contains

  !-------------------------------------------------------------------------------
  subroutine readCentNitrifDenitrifParams ( ncid )
    !
    ! !DESCRIPTION:
    ! read in nitrification denitrification parameters:

    ! !USES:
    use ncdio_pio        , only : file_desc_t,ncd_io
    use clm_varcon       , only : secspday
    use clm_time_manager , only : get_days_per_year
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNNitrifDenitrifParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in constants
    !
    tString='k_nitr_max'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%k_nitr_max=tempr

    tString='surface_tension_water'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%surface_tension_water=tempr

    tString='rij_kro_a'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_a=tempr

    tString='rij_kro_alpha'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_alpha=tempr

    tString='rij_kro_beta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_beta=tempr

    tString='rij_kro_gamma'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_gamma=tempr

    tString='rij_kro_delta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_delta=tempr

  end subroutine readCentNitrifDenitrifParams

  !-----------------------------------------------------------------------
  subroutine readCentDecompBgcParams ( ncid, nelms, betrtracer_vars )
    !
    ! !DESCRIPTION:
    ! read in decomposition parameters for century bgc
    !
    ! !USES:
    use ncdio_pio              , only: file_desc_t,ncd_io
    use clm_varcon             , only : secspday
    use clm_varctl             , only : spinup_state
    use clm_varpar             , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
    use clm_time_manager       , only : get_days_per_year
    use BeTRTracerType         , only : BeTRTracer_Type
    use CNDecompCascadeConType , only : decomp_cascade_con
    !
    ! !ARGUMENTS:
    type(file_desc_t)       , intent(inout) :: ncid   ! pio netCDF file id
    type(BeTRTracer_Type)   , intent(inout) :: betrtracer_vars
    integer                 , intent(in) :: nelms
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompBgcParamsType'
    character(len=100) :: errCode = 'Error reading in CN const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading

    real(r8) :: tau_l1
    real(r8) :: tau_l2_l3
    real(r8) :: tau_s1
    real(r8) :: tau_s2
    real(r8) :: tau_s3
    real(r8) :: days_per_year
    real(r8) :: tau_cwd
    real(r8) :: cn_s1
    real(r8) :: cn_s2
    real(r8) :: cn_s3
    integer  :: i_litr1
    integer  :: i_litr2
    integer  :: i_litr3
    integer  :: i_soil1
    integer  :: i_soil2
    integer  :: i_soil3
    integer  :: ii, jj, kk
    !-----------------------------------------------------------------------
    associate(                                                                                 & !
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio
         decomp_pool_name_restart       => decomp_cascade_con%decomp_pool_name_restart       , & ! Output: [character(len=8)  (:)     ]  name of pool for restart files
         decomp_pool_name_history       => decomp_cascade_con%decomp_pool_name_history       , & ! Output: [character(len=8)  (:)     ]  name of pool for history files
         decomp_pool_name_long          => decomp_cascade_con%decomp_pool_name_long          , & ! Output: [character(len=20) (:)     ]  name of pool for netcdf long names
         decomp_pool_name_short         => decomp_cascade_con%decomp_pool_name_short         , & ! Output: [character(len=8)  (:)     ]  name of pool for netcdf short names
         is_litter                      => decomp_cascade_con%is_litter                      , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool
         is_soil                        => decomp_cascade_con%is_soil                        , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool
         is_cwd                         => decomp_cascade_con%is_cwd                         , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio               , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools
         initial_stock                  => decomp_cascade_con%initial_stock                  , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup
         is_metabolic                   => decomp_cascade_con%is_metabolic                   , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material
         is_cellulose                   => decomp_cascade_con%is_cellulose                   , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose
         is_lignin                      => decomp_cascade_con%is_lignin                      , & ! Output: [logical           (:)     ]  TRUE => pool is lignin
         spinup_factor                  => decomp_cascade_con%spinup_factor                    & ! Output: [real(r8)          (:)
         )
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

    !------- time-constant coefficients ---------- !
    ! set soil organic matter compartment C:N ratios
    cn_s1 = CNDecompBgcParamsInst%cn_s1_bgc
    cn_s2 = CNDecompBgcParamsInst%cn_s2_bgc
    cn_s3 = CNDecompBgcParamsInst%cn_s3_bgc

    !-------------------  list of pools and their attributes  ------------
    i_litr1                                 = i_met_lit
    floating_cn_ratio_decomp_pools(i_litr1) = .true.
    decomp_pool_name_restart(i_litr1)       = 'litr1'
    decomp_pool_name_history(i_litr1)       = 'LITR1'
    decomp_pool_name_long(i_litr1)          = 'litter 1'
    decomp_pool_name_short(i_litr1)         = 'L1'
    is_litter(i_litr1)                      = .true.
    is_soil(i_litr1)                        = .false.
    is_cwd(i_litr1)                         = .false.
    initial_cn_ratio(i_litr1)               = 90._r8
    initial_stock(i_litr1)                  = 0._r8
    is_metabolic(i_litr1)                   = .true.
    is_cellulose(i_litr1)                   = .false.
    is_lignin(i_litr1)                      = .false.

    i_litr2                                 = i_cel_lit
    floating_cn_ratio_decomp_pools(i_litr2) = .true.
    decomp_pool_name_restart(i_litr2)       = 'litr2'
    decomp_pool_name_history(i_litr2)       = 'LITR2'
    decomp_pool_name_long(i_litr2)          = 'litter 2'
    decomp_pool_name_short(i_litr2)         = 'L2'
    is_litter(i_litr2)                      = .true.
    is_soil(i_litr2)                        = .false.
    is_cwd(i_litr2)                         = .false.
    initial_cn_ratio(i_litr2)               = 90._r8
    initial_stock(i_litr2)                  = 0._r8
    is_metabolic(i_litr2)                   = .false.
    is_cellulose(i_litr2)                   = .true.
    is_lignin(i_litr2)                      = .false.

    i_litr3                                 = i_lig_lit
    floating_cn_ratio_decomp_pools(i_litr3) = .true.
    decomp_pool_name_restart(i_litr3)       = 'litr3'
    decomp_pool_name_history(i_litr3)       = 'LITR3'
    decomp_pool_name_long(i_litr3)          = 'litter 3'
    decomp_pool_name_short(i_litr3)         = 'L3'
    is_litter(i_litr3)                      = .true.
    is_soil(i_litr3)                        = .false.
    is_cwd(i_litr3)                         = .false.
    initial_cn_ratio(i_litr3)               = 90._r8
    initial_stock(i_litr3)                  = 0._r8
    is_metabolic(i_litr3)                   = .false.
    is_cellulose(i_litr3)                   = .false.
    is_lignin(i_litr3)                      = .true.

    ! CWD
    floating_cn_ratio_decomp_pools(i_cwd)   = .true.
    decomp_pool_name_restart(i_cwd)         = 'cwd'
    decomp_pool_name_history(i_cwd)         = 'CWD'
    decomp_pool_name_long(i_cwd)            = 'coarse woody debris'
    decomp_pool_name_short(i_cwd)           = 'CWD'
    is_litter(i_cwd)                        = .false.
    is_soil(i_cwd)                          = .false.
    is_cwd(i_cwd)                           = .true.
    initial_cn_ratio(i_cwd)                 = 90._r8
    initial_stock(i_cwd)                    = 0._r8
    is_metabolic(i_cwd)                     = .false.
    is_cellulose(i_cwd)                     = .false.
    is_lignin(i_cwd)                        = .false.

    i_soil1                                 = 5
    floating_cn_ratio_decomp_pools(i_soil1) = .false.
    decomp_pool_name_restart(i_soil1)       = 'soil1'
    decomp_pool_name_history(i_soil1)       = 'SOIL1'
    decomp_pool_name_long(i_soil1)          = 'soil 1'
    decomp_pool_name_short(i_soil1)         = 'S1'
    is_litter(i_soil1)                      = .false.
    is_soil(i_soil1)                        = .true.
    is_cwd(i_soil1)                         = .false.
    initial_cn_ratio(i_soil1)               = cn_s1
    initial_stock(i_soil1)                  = 20._r8
    is_metabolic(i_soil1)                   = .false.
    is_cellulose(i_soil1)                   = .false.
    is_lignin(i_soil1)                      = .false.

    i_soil2                                 = 6
    floating_cn_ratio_decomp_pools(i_soil2) = .false.
    decomp_pool_name_restart(i_soil2)       = 'soil2'
    decomp_pool_name_history(i_soil2)       = 'SOIL2'
    decomp_pool_name_long(i_soil2)          = 'soil 2'
    decomp_pool_name_short(i_soil2)         = 'S2'
    is_litter(i_soil2)                      = .false.
    is_soil(i_soil2)                        = .true.
    is_cwd(i_soil2)                         = .false.
    initial_cn_ratio(i_soil2)               = cn_s2
    initial_stock(i_soil2)                  = 20._r8
    is_metabolic(i_soil2)                   = .false.
    is_cellulose(i_soil2)                   = .false.
    is_lignin(i_soil2)                      = .false.

    i_soil3                                 = 7
    floating_cn_ratio_decomp_pools(i_soil3) = .false.
    decomp_pool_name_restart(i_soil3)       = 'soil3'
    decomp_pool_name_history(i_soil3)       = 'SOIL3'
    decomp_pool_name_long(i_soil3)          = 'soil 3'
    decomp_pool_name_short(i_soil3)         = 'S3'
    is_litter(i_soil3)                      = .false.
    is_soil(i_soil3)                        = .true.
    is_cwd(i_soil3)                         = .false.
    initial_cn_ratio(i_soil3)               = cn_s3
    initial_stock(i_soil3)                  = 20._r8
    is_metabolic(i_soil3)                   = .false.
    is_cellulose(i_soil3)                   = .false.
    is_lignin(i_soil3)                      = .false.

    spinup_factor(i_litr1)                  = 1._r8
    spinup_factor(i_litr2)                  = 1._r8
    spinup_factor(i_litr3)                  = 1._r8
    spinup_factor(i_cwd)                    = 1._r8
    spinup_factor(i_soil1)                  = CNDecompBgcParamsInst%spinup_vector(1)
    spinup_factor(i_soil2)                  = CNDecompBgcParamsInst%spinup_vector(2)
    spinup_factor(i_soil3)                  = CNDecompBgcParamsInst%spinup_vector(3)

    tau_l1                                  = 1./18.5
    tau_l2_l3                               = 1./4.9
    tau_s1                                  = 1./7.3
    tau_s2                                  = 1./0.2
    tau_s3                                  = 1./.0045

    ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
    tau_cwd                                 = 1./0.3
    days_per_year                           = get_days_per_year()

    CNDecompBgcParamsInst%k_decay_lit1=1._r8/(secspday * days_per_year * tau_l1)        ![1/s]
    CNDecompBgcParamsInst%k_decay_lit2=1._r8/(secspday * days_per_year * tau_l2_l3)
    CNDecompBgcParamsInst%k_decay_lit3=1._r8/(secspday * days_per_year * tau_l2_l3)
    CNDecompBgcParamsInst%k_decay_som1=1._r8/(secspday * days_per_year * tau_s1)
    CNDecompBgcParamsInst%k_decay_som2=1._r8/(secspday * days_per_year * tau_s2)
    CNDecompBgcParamsInst%k_decay_som3=1._r8/(secspday * days_per_year * tau_s3)
    CNDecompBgcParamsInst%k_decay_cwd =1._r8/(secspday * days_per_year * tau_cwd)


    kk = 1
    betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) = 1._r8
    betrtracer_vars%tracer_solid_passive_diffus_thc_group(kk) = 1.e-30_r8

    if ( spinup_state .eq. 1 ) then
       CNDecompBgcParamsInst%k_decay_som1 = CNDecompBgcParamsInst%k_decay_som1 * CNDecompBgcParamsInst%spinup_vector(1)
       CNDecompBgcParamsInst%k_decay_som2 = CNDecompBgcParamsInst%k_decay_som2 * CNDecompBgcParamsInst%spinup_vector(2)
       CNDecompBgcParamsInst%k_decay_som3 = CNDecompBgcParamsInst%k_decay_som3 * CNDecompBgcParamsInst%spinup_vector(3)

       ii=i_soil1
       kk = 2
       do jj = 1, nelms
          betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) = &
               betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) * spinup_factor(ii)
       enddo

       ii=i_soil2
       kk = 3
       do jj = 1, nelms
          betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) = &
               betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) * spinup_factor(ii)
       enddo

       ii=i_soil3
       kk = 4
       do jj = 1, nelms
          betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) = &
               betrtracer_vars%tracer_solid_passive_diffus_scal_group(kk) * spinup_factor(ii)
       enddo

    endif
  end associate

end subroutine readCentDecompBgcParams

!-----------------------------------------------------------------------
subroutine readCentCNAllocParams ( ncid )
  !
  ! !DESCRIPTION:
  ! read in allocation parameters.
  !
  ! !USES:
  use ncdio_pio , only : file_desc_t,ncd_io

  ! !ARGUMENTS:
  implicit none
  type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
  !
  ! !LOCAL VARIABLES:
  character(len=32)  :: subname = 'readCentCNAllocParams'
  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv   ! has variable been read in or not
  real(r8)           :: tempr   ! temporary to read in parameter
  character(len=100) :: tString ! temp. var for reading
  !-----------------------------------------------------------------------

  ! read in parameters

  tString='compet_plant_no3'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_plant_no3=tempr

  tString='compet_plant_nh4'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_plant_nh4=tempr

  tString='compet_decomp_no3'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_decomp_no3=tempr

  tString='compet_decomp_nh4'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_decomp_nh4=tempr

  tString='compet_denit'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_denit=tempr

  tString='compet_nit'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
  NutrientCompetitionParamsInst%compet_nit=tempr


end subroutine readCentCNAllocParams


end module BGCCenturyParMod
