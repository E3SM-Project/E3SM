module CNSharedParamsMod

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save

  ! CNParamsShareInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public  :: CNParamsShareType
      real(r8) :: Q10         ! temperature dependence
      real(r8) :: minpsi      ! minimum soil water potential for heterotrophic resp	  
      real(r8) :: cwd_fcel    ! cellulose fraction of coarse woody debris
      real(r8) :: cwd_flig    ! lignin fraction of coarse woody debris
      real(r8) :: froz_q10    ! separate q10 for frozen soil respiration rates
      real(r8) :: decomp_depth_efolding ! e-folding depth for reduction in decomposition (m) 
      real(r8) :: mino2lim    ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
      real(r8) :: organic_max ! organic matter content (kg/m3) where soil is assumed to act like peat
  end type CNParamsShareType

  type(CNParamsShareType),protected :: CNParamsShareInst

  logical, public :: anoxia_wtsat = .false.
  integer, public :: nlev_soildecomp_standard = 5

  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine CNParamsReadShared(ncid)
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'CNParamsReadShared'
    character(len=100) :: errCode = '-Error reading in CN and BGC shared params file. Var:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! netcdf read here
    !
    tString='q10_mr'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%Q10=tempr

    tString='minpsi_hr'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%minpsi=tempr 

    tString='cwd_fcel'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%cwd_fcel=tempr

    tString='cwd_flig'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%cwd_flig=tempr 

    tString='froz_q10'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%froz_q10=tempr   

    tString='decomp_depth_efolding'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%decomp_depth_efolding=tempr  

    tString='mino2lim'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%mino2lim=tempr 
    !CNParamsShareInst%mino2lim=0.2_r8 

    tString='organic_max'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNParamsShareInst%organic_max=tempr

  end subroutine CNParamsReadShared

end module CNSharedParamsMod
