module EDSharedParamsMod

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none

  ! EDParamsShareInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public  :: EDParamsShareType
      real(r8) :: Q10      ! temperature dependence
      real(r8) :: froz_q10 ! separate q10 for frozen soil respiration rates
  end type EDParamsShareType

  type(EDParamsShareType), protected :: EDParamsShareInst

  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine EDParamsReadShared(ncid)
    !
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    character(len=32)  :: subname = 'EDParamsReadShared'
    character(len=100) :: errCode = '-Error reading in ED shared params file. Var:'
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
    EDParamsShareInst%Q10=tempr

    tString='froz_q10'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    EDParamsShareInst%froz_q10=tempr   

  end subroutine EDParamsReadShared
  
end module EDSharedParamsMod
