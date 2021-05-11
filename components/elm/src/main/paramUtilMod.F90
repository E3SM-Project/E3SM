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
      module procedure readNcdioScalarCheckDimensions
      module procedure readNcdioArray1dCheckDimensions
      module procedure readNcdioArray2dCheckDimensions
   end interface readNcdio

   public :: readNcdioScalar
   public :: readNcdioArray1d
   public :: readNcdioArray2d
   public :: readNcdioScalarCheckDimensions
   public :: readNcdioArray1dCheckDimensions
   public :: readNcdioArray2dCheckDimensions

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
   !-----------------------------------------------------------------------
   !
   !-----------------------------------------------------------------------
   subroutine readNcdioScalarCheckDimensions(ncid, varName, expected_numDims, expected_dimNames, &
         callingName, retVal) 
      !
      ! read the netcdf file...generic, could be used for any parameter read
      !
      use abortutils   , only : endrun
      use ncdio_pio    , only : file_desc_t

      implicit none

      ! arguments
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      character(len=*), intent(in)    :: varName ! variable we are reading
      integer, intent(in) :: expected_numDims
      character(len=*), intent(in)    :: expected_dimNames(:) ! expected dimension name
      character(len=*), intent(in)    :: callingName ! calling routine
      real(r8), intent(inout) :: retVal

      ! local vars
      character(len=32)  :: subname = 'readNcdio::'
      character(len=100) :: errCode = ' - Error reading.  Var: '

      !
      ! netcdf read here
      !
      call checkDimensions(ncid, varName, expected_numDims, expected_dimNames, subname)
      call readNcdio(ncid, varName, callingName, retVal)

   end subroutine readNcdioScalarCheckDimensions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !
   !-----------------------------------------------------------------------
   subroutine readNcdioArray1dCheckDimensions(ncid, varName, expected_numDims, expected_dimNames, &
         callingName, retVal) 
      !
      ! read the netcdf file...generic, could be used for any parameter read
      !
      use abortutils   , only : endrun
      use ncdio_pio    , only : file_desc_t

      implicit none

      ! arguments
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      character(len=*), intent(in)    :: varName ! variable we are reading
      integer, intent(in) :: expected_numDims
      character(len=*), intent(in)    :: expected_dimNames(:) ! expected dimension name
      character(len=*), intent(in)    :: callingName ! calling routine
      real(r8), intent(inout) :: retVal( 1: )

      ! local vars
      character(len=32)  :: subname = 'readNcdio::'
      character(len=100) :: errCode = ' - Error reading.  Var: '
      !
      ! netcdf read here
      !
      call checkDimensions(ncid, varName, expected_numDims, expected_dimNames, subname)
      call readNcdio(ncid, varName, callingName, retVal)

   end subroutine readNcdioArray1dCheckDimensions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !
   !-----------------------------------------------------------------------
   subroutine readNcdioArray2dCheckDimensions(ncid, varName, expected_numDims, expected_dimNames, &
         callingName, retVal) 
      !
      ! read the netcdf file...generic, could be used for any parameter read
      !
      use abortutils   , only : endrun
      use ncdio_pio    , only : file_desc_t

      implicit none

      ! arguments
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      character(len=*), intent(in)    :: varName ! variable we are reading
      integer, intent(in) :: expected_numDims
      character(len=*), intent(in)    :: expected_dimNames(:) ! expected dimension name
      character(len=*), intent(in)    :: callingName ! calling routine
      real(r8), intent(inout) :: retVal(1:, : )

      ! local vars
      character(len=32)  :: subname = 'readNcdio::'
      character(len=100) :: errCode = ' - Error reading.  Var: '
      !
      ! netcdf read here
      !
      call checkDimensions(ncid, varName, expected_numDims, expected_dimNames, subname)
      call readNcdio(ncid, varName, callingName, retVal)

   end subroutine readNcdioArray2dCheckDimensions
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   !
   !-----------------------------------------------------------------------
   subroutine checkDimensions(ncid, varName, expected_numDims, expected_dimNames, callingName) 
      !
      ! Assert that the expected number of dimensions and dimension
      ! names for a variable match the actual names on the file.
      !
      use abortutils   , only : endrun
      use ncdio_pio    , only : file_desc_t, var_desc_t, check_var, ncd_inqvdname, ncd_inqvdims

      implicit none

      ! arguments
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      character(len=*), intent(in)    :: varName ! variable we are reading
      integer,          intent(in)    :: expected_numDims ! number of expected dimensions on the variable
      character(len=*), intent(in)    :: expected_dimNames(:) ! expected dimension names
      character(len=*), intent(in)    :: callingName ! calling routine
      integer                         :: error_num

      ! local vars
      character(len=32)  :: subname = 'checkDimensions::'
      type(Var_desc_t)   :: var_desc        ! variable descriptor
      logical            :: readvar        ! whether the variable was found
      character(len=100) :: received_dimName
      integer            :: d, num_dims
      character(len=256) :: msg

      call check_var(ncid, varName, var_desc, readvar)
      if (readvar) then
         call ncd_inqvdims(ncid, num_dims, var_desc)
         if (num_dims /= expected_numDims) then
            write(msg, *) trim(callingName)//trim(subname)//trim(varname)//":: expected number of dimensions = ", &
                  expected_numDims, "   num dimensions received from file = ", num_dims
            call endrun(msg)
         end if
         do d = 1, num_dims
            received_dimName = ''
            call ncd_inqvdname(ncid, varname=trim(varName), dimnum=d, dname=received_dimName, err_code=error_num)
            if (trim(expected_dimNames(d)) /= trim(received_dimName)) then
               write(msg, *) trim(callingName)//trim(subname)//trim(varname)//":: dimension ", d, &
                     " expected dimension name '"//trim(expected_dimNames(d))//&
                     "' dimension name received from file '"//trim(received_dimName)//"'." 
               call endrun(msg)
            end if
         end do
      end if

   end subroutine checkDimensions
   !-----------------------------------------------------------------------
end module paramUtilMod
