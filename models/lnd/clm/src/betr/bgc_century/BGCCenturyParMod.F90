module BGCCenturyParMod
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use abortutils                         , only : endrun
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
implicit none

  public :: readCentDecompBgcParams
  public :: readCentNitrifDenitrifParams


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


  type, private :: CNDecompBgcParamsType
     real(r8) :: cn_s1_bgc     !C:N for SOM 1
     real(r8) :: cn_s2_bgc     !C:N for SOM 2
     real(r8) :: cn_s3_bgc     !C:N for SOM 3

     real(r8) :: rf_l1s1_bgc   !respiration fraction litter 1 -> SOM 1
     real(r8) :: rf_l2s1_bgc
     real(r8) :: rf_l3s2_bgc

     real(r8) :: rf_s2s1_bgc    
     real(r8) :: rf_s2s3_bgc    
     real(r8) :: rf_s3s1_bgc    

     real(r8) :: rf_cwdl2_bgc 
     real(r8) :: rf_cwdl3_bgc

     real(r8) :: tau_l1_bgc    ! turnover time of  litter 1 (yr)
     real(r8) :: tau_l2_l3_bgc ! turnover time of  litter 2 and litter 3 (yr)
     real(r8) :: tau_s1_bgc    ! turnover time of  SOM 1 (yr)
     real(r8) :: tau_s2_bgc    ! turnover time of  SOM 2 (yr)
     real(r8) :: tau_s3_bgc    ! turnover time of  SOM 3 (yr)
     real(r8) :: tau_cwd_bgc   ! corrected fragmentation rate constant CWD

     real(r8) :: k_decay_lit1
     real(r8) :: k_decay_lit2
     real(r8) :: k_decay_lit3
     real(r8) :: k_decay_som1
     real(r8) :: k_decay_som2
     real(r8) :: k_decay_som3
     real(r8) :: k_decay_cwd
     
     real(r8) :: cwd_fcel_bgc !cellulose fraction for CWD
     real(r8) :: cwd_flig_bgc !

     real(r8) :: k_frag_bgc   !fragmentation rate for CWD
     real(r8) :: minpsi_bgc   !minimum soil water potential for heterotrophic resp
     
     integer  :: nsompools = 3
     real(r8),allocatable :: spinup_vector(:) ! multipliers for soil decomp during accelerated spinup

  end type CNDecompBgcParamsType

  type(CNDecompBgcParamsType),protected ::  CNDecompBgcParamsInst

  
  contains


!-------------------------------------------------------------------------------
  subroutine readCentNitrifDenitrifParams ( ncid )
    !
    use ncdio_pio    , only: file_desc_t,ncd_io
   
    use clm_varcon   , only : secspday
    use clm_time_manager , only : get_days_per_year
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNNitrifDenitrifParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
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
  subroutine readCentDecompBgcParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    use clm_varcon   , only : secspday
    use clm_time_manager , only : get_days_per_year
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
    
    real(r8) :: tau_l1  
    real(r8) :: tau_l2_l3 
    real(r8) :: tau_s1  
    real(r8) :: tau_s2  
    real(r8) :: tau_s3     
    real(r8) :: days_per_year
    real(r8) :: tau_cwd
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

    tau_l1 = 1./18.5
    tau_l2_l3 = 1./4.9
    tau_s1 = 1./7.3
    tau_s2 = 1./0.2
    tau_s3 = 1./.0045
    ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
    tau_cwd  = 1./0.3    
    days_per_year = get_days_per_year()    

    CNDecompBgcParamsInst%k_decay_lit1=1._r8/(secspday * days_per_year * tau_l1)        ![1/s]
    CNDecompBgcParamsInst%k_decay_lit2=1._r8/(secspday * days_per_year * tau_l2_l3)
    CNDecompBgcParamsInst%k_decay_lit3=1._r8/(secspday * days_per_year * tau_l2_l3)
    CNDecompBgcParamsInst%k_decay_som1=1._r8/(secspday * days_per_year * tau_s1)
    CNDecompBgcParamsInst%k_decay_som2=1._r8/(secspday * days_per_year * tau_s2)
    CNDecompBgcParamsInst%k_decay_som3=1._r8/(secspday * days_per_year * tau_s3)
    CNDecompBgcParamsInst%k_decay_cwd =1._r8/(secspday * days_per_year * tau_cwd)    
  end subroutine readCentDecompBgcParams
end module BGCCenturyParMod
