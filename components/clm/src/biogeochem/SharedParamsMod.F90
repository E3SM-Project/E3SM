module SharedParamsMod

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save

  ! ParamsShareInst.  PGI wants the type decl. public but the instance
  ! is indeed protected.  A generic private statement at the start of the module
  ! overrides the protected functionality with PGI

  type, public :: ParamsShareType

      real(r8), pointer :: Q10_mr                => null() ! temperature dependence for maintenance respiraton
      real(r8), pointer :: Q10_hr                => null() ! temperature dependence for heterotrophic respiration
      real(r8), pointer :: minpsi                => null() ! minimum soil water potential for heterotrophic resp
      real(r8), pointer :: cwd_fcel              => null() ! cellulose fraction of coarse woody debris
      real(r8), pointer :: cwd_flig              => null() ! lignin fraction of coarse woody debris
      real(r8), pointer :: froz_q10              => null() ! separate q10 for frozen soil respiration rates
      real(r8), pointer :: decomp_depth_efolding => null() ! e-folding depth for reduction in decomposition (m)
      real(r8), pointer :: mino2lim              => null() ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
      real(r8), pointer :: organic_max           => null() ! organic matter content (kg/m3) where soil is assumed to act like peat

  end type ParamsShareType

  type(ParamsShareType),protected :: ParamsShareInst

  !$acc declare create(ParamsShareInst)
  logical, public :: anoxia_wtsat = .false.
  integer, public :: nlev_soildecomp_standard = 5
  !$acc declare create(anoxia_wtsat)
  !$acc declare create(nlev_soildecomp_standard)

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
   subroutine ParamsReadShared(ncid)
     !
     use ncdio_pio   , only : file_desc_t,ncd_io
     use abortutils  , only : endrun
     use shr_log_mod , only : errMsg => shr_log_errMsg
     !
     implicit none
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     !
     character(len=32)  :: subname = 'ParamsReadShared'
     character(len=100) :: errCode = '-Error reading in CN and BGC shared params file. Var:'
     logical            :: readv ! has variable been read in or not
     real(r8)           :: tempr ! temporary to read in parameter
     character(len=100) :: tString ! temp. var for reading
     !-----------------------------------------------------------------------
     !
     ! netcdf read here
     !
     allocate(ParamsShareInst%Q10_mr               )
     allocate(ParamsShareInst%Q10_hr               )
     allocate(ParamsShareInst%minpsi               )
     allocate(ParamsShareInst%cwd_fcel             )
     allocate(ParamsShareInst%cwd_flig             )
     allocate(ParamsShareInst%froz_q10             )
     allocate(ParamsShareInst%decomp_depth_efolding)
     allocate(ParamsShareInst%mino2lim             )
     allocate(ParamsShareInst%organic_max          )
     tString='q10_mr'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%Q10_mr=tempr

     tString='q10_hr'
     call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%Q10_hr=tempr


     tString='minpsi_hr'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%minpsi=tempr

     tString='cwd_fcel'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%cwd_fcel=tempr

     tString='cwd_flig'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%cwd_flig=tempr

     tString='froz_q10'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%froz_q10=tempr

     tString='decomp_depth_efolding'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%decomp_depth_efolding=tempr

     tString='mino2lim'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%mino2lim=tempr
     !ParamsShareInst%mino2lim=0.2_r8

     tString='organic_max'
     call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
     if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
     ParamsShareInst%organic_max=tempr

   end subroutine ParamsReadShared

end module SharedParamsMod
