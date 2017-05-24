module paramUtilMod
   !
   ! module that deals with reading parameter files
   !
   use shr_kind_mod , only: r8 => shr_kind_r8
   implicit none
   save
   private

   interface readNcdio
      module procedure readNcdioScalar
      module procedure readNcdioArray1d
      module procedure readNcdioArray2d
   end interface

   public :: readNcdioScalar
   public :: readNcdioArray1d
   public :: readNcdioArray2d

   public :: readNcdio

contains
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine readNcdioScalar(ncid, varName, callingName, retVal)
     !
     ! read the netcdf file...generic, could be used for any parameter read
     !
     use abortutils   , only : endrun
     use ncdio_pio    , only : file_desc_t,ncd_io

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     character(len=*), intent(in)    :: varName ! variable we are reading
     character(len=*), intent(in)    :: callingName ! calling routine
     real(r8), intent(inout) :: retVal

     ! local vars
     character(len=32)  :: subname = 'readNcdio::'
     character(len=100) :: errCode = ' - Error reading.  Var: '
     logical            :: readv     ! has variable been read in or not

     !
     ! netcdf read here
     !

     call ncd_io(varname=trim(varName),data=retVal, flag='read', ncid=ncid, readvar=readv)

     if ( .not. readv ) then
        call endrun(trim(callingName)//trim(subname)//trim(errCode)//trim(varName))
     endif

  end subroutine readNcdioScalar
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine readNcdioArray1d(ncid, varName, callingName, retVal) 
     !
     ! read the netcdf file...generic, could be used for any parameter read
     !
     use abortutils   , only : endrun
     use ncdio_pio    , only : file_desc_t,ncd_io

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     character(len=*), intent(in)    :: varName ! variable we are reading
     character(len=*), intent(in)    :: callingName ! calling routine
     real(r8), intent(inout) :: retVal( 1: )

     ! local vars
     character(len=32)  :: subname = 'readNcdio::'
     character(len=100) :: errCode = ' - Error reading.  Var: '
     logical            :: readv     ! has variable been read in or not

     !
     ! netcdf read here
     !

     call ncd_io(varname=trim(varName),data=retVal, flag='read', ncid=ncid, readvar=readv)

     if ( .not. readv ) then
        call endrun(trim(callingName)//trim(subname)//trim(errCode)//trim(varName))
     endif

  end subroutine readNcdioArray1d
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  subroutine readNcdioArray2d(ncid, varName, callingName, retVal) 
     !
     ! read the netcdf file...generic, could be used for any parameter read
     !
     use abortutils   , only : endrun
     use ncdio_pio    , only : file_desc_t,ncd_io

     implicit none

     ! arguments
     type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
     character(len=*), intent(in)    :: varName ! variable we are reading
     character(len=*), intent(in)    :: callingName ! calling routine
     real(r8), intent(inout) :: retVal( 1: , :)

     ! local vars
     character(len=32)  :: subname = 'readNcdio::'
     character(len=100) :: errCode = ' - Error reading.  Var: '
     logical            :: readv     ! has variable been read in or not

     !
     ! netcdf read here
     !

     call ncd_io(varname=trim(varName),data=retVal, flag='read', ncid=ncid, readvar=readv)

     if ( .not. readv ) then
        call endrun(trim(callingName)//trim(subname)//trim(errCode)//trim(varName))
     endif

  end subroutine readNcdioArray2d
  !-----------------------------------------------------------------------

end module paramUtilMod
