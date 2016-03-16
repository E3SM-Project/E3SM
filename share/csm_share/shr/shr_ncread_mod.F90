!===============================================================================
! SVN $Id: shr_ncread_mod.F90 29597 2011-08-04 02:24:04Z erik $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_150116/shr/shr_ncread_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_ncread_mod -- semi-generic netCDF file reader
!
! !DESCRIPTION:
! Reads netcdf stuff off a file
! \newline
! General Usage:
!    check = shr_ncread_varExists('myfile','sst')
!    call shr_ncread_varDimNum('myfile','sst',ndims)
!    call shr_ncread_varDimSize('myfile','sst','lon',nsize)
!    call shr_ncread_varDimSize('myfile','sst',  2  ,nsize)
!    call shr_ncread_varDimSizes('myfile','sst',ns1,ns2,ns3)
!    call shr_ncread_dimSize('myfile','lon',nsize)
!    call shr_ncread_domain('myfile','xc',lon,'yc',lat,'mask',imask,'area',area)
!    call shr_ncread_tField('myfile',6,'sst',a2d)
!    call shr_ncread_tField('myfile',6,'sst',a2d,'xc','yc','time')
!    call shr_ncread_tField('myfile',1,'zlev',a1d)
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d)
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d,dim1='lon',dim2='lat',dim3='time',dim3i=21)
!    call shr_ncread_setAbort(.true.)
!    call shr_ncread_setDebug(1)
! \newline
!
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!     2005-Apr-21 - B. Kauffman, J. Schramm, M. Vertenstein  - first design
!
! !INTERFACE: ------------------------------------------------------------------

module shr_ncread_mod

! !USES:

   use shr_string_mod  ! string methods
   use shr_kind_mod    ! kinds
   use shr_sys_mod     ! shared system calls
   use shr_file_mod    ! file methods
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use netcdf

   implicit none

   private              ! everything is default private

! !PUBLIC TYPES:

   ! no public data types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_ncread_varExists
   public :: shr_ncread_varDimNum
   public :: shr_ncread_varDimSize
   public :: shr_ncread_varDimSizes
   public :: shr_ncread_dimSize
   public :: shr_ncread_domain
   public :: shr_ncread_tField
   public :: shr_ncread_Field4dG
   public :: shr_ncread_handleErr
   public :: shr_ncread_setAbort
   public :: shr_ncread_setDebug
   public :: shr_ncread_open
   public :: shr_ncread_close

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

   interface shr_ncread_varDimSize  ; module procedure &
      shr_ncread_varDimSizeName, &
      shr_ncread_varDimSizeID
   end interface

   interface shr_ncread_dimSize  ; module procedure &
      shr_ncread_dimSizeName
   end interface

   interface shr_ncread_tField  ; module procedure &
      shr_ncread_tField2dR8, &
      shr_ncread_tField1dR8, &
      shr_ncread_tField2dIN, &
      shr_ncread_tField1dIN
   end interface

   logical             ,save      :: doabort = .true.
   integer(SHR_KIND_IN),save      :: debug = 0

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_varExists -- return logical for existance of var
!
! !DESCRIPTION:
! Return logical if variable name exists on file
! \newline
! General Usage:
!   check = shr_ncread_varExists('myfile','sst')
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_ncread_varExists(fileName, varName)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)   :: fileName ! nc file name
   character(*),intent(in)   :: varName  ! name of variable

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: fid
   integer(SHR_KIND_IN) :: vid
   integer(SHR_KIND_IN) :: debug0
   integer(SHR_KIND_IN) :: rCode

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_varExists)"
   character(*),parameter :: F00     = "('(shr_ncread_varExists) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_varExists) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   !--- turn off debug writing ---
   debug0 = debug
   call shr_ncread_setDebug(0)

   shr_ncread_varExists = .false.
   call shr_ncread_open(fileName,fid,rCode)
   rcode = nf90_inq_varid(fid,trim(varName),vid)
   if (rcode == nf90_noerr) shr_ncread_varExists = .true.
   call shr_ncread_close(fid,rCode)

   !--- reset debug code ---
   call shr_ncread_setDebug(debug0)

end function shr_ncread_varExists
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_varDimNum -- return num of dimensions of a variable
!
! !DESCRIPTION:
! Returns the number of dimensions of a named variable
! \newline
! General Usage:
!    call shr_ncread_varDimNum('myfile','sst',ndims)
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_varDimNum(fileName, varName, ns, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fileName ! nc file name
   character(*)        ,intent(in)           :: varName  ! name of variable
   integer(SHR_KIND_IN),intent(out)          :: ns       ! number of dims of var
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: fid
   integer(SHR_KIND_IN) :: vid
   integer(SHR_KIND_IN) :: rCode

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_varDimNum)"
   character(*),parameter :: F00     = "('(shr_ncread_varDimNum) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_varDimNum) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   call shr_ncread_open(fileName,fid,rCode)

   !--- read variable info ---
   rcode = nf90_inq_varid(fid,trim(varName),vid)
   call shr_ncread_handleErr(rCode, subName//" ERROR inq varid")
   rcode = nf90_inquire_variable(fid,vid,ndims=ns)
   call shr_ncread_handleErr(rCode, subName//" ERROR inq var")
   if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) trim(varName)//' has dims = ',ns

   call shr_ncread_close(fid,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_varDimNum
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_varDimSizeName -- return var dim size by dim name
!
! !DESCRIPTION:
! Returns the size of a dimension of a variable, both dimension and 
! variable are named.
! \newline
! General Usage:
!    call shr_ncread_varDimSize('myfile','sst','lon',nsize)
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_varDimSizeName(fileName, varName, dimName,  ns, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fileName ! nc file name
   character(*)        ,intent(in)           :: varName  ! name of variable
   character(*)        ,intent(in)           :: dimName  ! name of dimension
   integer(SHR_KIND_IN),intent(out)          :: ns       ! number of dims of var
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: rCode

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_varDimSizeName)"
   character(*),parameter :: F00     = "('(shr_ncread_varDimSizeName) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_varDimSizeName) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   call shr_ncread_dimSizeName(fileName,dimName,ns,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_varDimSizeName
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_varDimSizeID -- return var dim size by dim number
!
! !DESCRIPTION:
! Returns the size of a dimension of a variable where the variable is
! named and the dimension is numbered.
! \newline
! General Usage:
!    call shr_ncread_varDimSize('myfile','sst',2,nsize)
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_varDimSizeID(fileName, varName, dnum,  ns, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fileName ! nc file name
   character(*)        ,intent(in)           :: varName  ! name of variable
   integer(SHR_KIND_IN),intent(in)           :: dnum     ! dim number in var
   integer(SHR_KIND_IN),intent(out)          :: ns       ! size of dim in var
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: fid                   ! file id
   integer(SHR_KIND_IN)   :: vid                   ! var id
   integer(SHR_KIND_IN)   :: ndims                 ! number of dims
   character(SHR_KIND_CS) :: dimName               ! dim name
   integer(SHR_KIND_IN),allocatable :: dids(:)     ! dim ids
   integer(SHR_KIND_IN)   :: rCode                 ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_varDimSizeID)"
   character(*),parameter :: F00     = "('(shr_ncread_varDimSizeID) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_varDimSizeID) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   call shr_ncread_open(fileName,fid,rCode)

   rCode = nf90_inq_varid(fid,trim(varName),vid)
   call shr_ncread_handleErr(rCode,subName//' ERROR inq varid vid')
   rCode = nf90_inquire_variable(fid,vid,ndims=ndims)
   call shr_ncread_handleErr(rCode,subName//' ERROR inquire variable ndims')
   allocate(dids(ndims))
   rCode = nf90_inquire_variable(fid,vid,dimids=dids)
   call shr_ncread_handleErr(rCode,subName//' ERROR inquire variable dimids')
   rcode = nf90_inquire_dimension(fid,dids(dnum),name=dimName,len=ns)
   call shr_ncread_handleErr(rCode, subName//" ERROR inquire dimension")
   if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) trim(dimName)//' dimension has size = ',ns

   deallocate(dids)
   call shr_ncread_close(fid,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_varDimSizeID
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_varDimSizes -- return var dim sizes 
!
! !DESCRIPTION:
! Returns the dimension sizes of a named variable using optional arguments.
! Each optional argument represents a numbered dimension.
! /newline
! General Usage:
!    call shr_ncread_varDimSizes('myfile','sst',ns1,ns2,ns3)
! /newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_varDimSizes(fileName, varName, n1, n2, n3, n4, n5, n6, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fileName ! nc file name
   character(*)        ,intent(in)           :: varName  ! name of variable
   integer(SHR_KIND_IN),intent(out),optional :: n1       ! size of dim1 in var
   integer(SHR_KIND_IN),intent(out),optional :: n2       ! size of dim2 in var
   integer(SHR_KIND_IN),intent(out),optional :: n3       ! size of dim3 in var
   integer(SHR_KIND_IN),intent(out),optional :: n4       ! size of dim4 in var
   integer(SHR_KIND_IN),intent(out),optional :: n5       ! size of dim5 in var
   integer(SHR_KIND_IN),intent(out),optional :: n6       ! size of dim6 in var
   Integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN),parameter :: maxn = 6     ! max number of dims available
   integer(SHR_KIND_IN) :: n                      ! counter
   integer(SHR_KIND_IN) :: fid                    ! file id
   integer(SHR_KIND_IN) :: vid                    ! variable id
   integer(SHR_KIND_IN) :: ndims                  ! number of dims
   integer(SHR_KIND_IN),allocatable :: dids(:)    ! dimids
   integer(SHR_KIND_IN),allocatable :: ns(:)      ! size of dims
   integer(SHR_KIND_IN) :: rCode                  ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_varDimSizes)"
   character(*),parameter :: F00     = "('(shr_ncread_varDimSizes) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_varDimSizes) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   call shr_ncread_open(fileName,fid,rCode)

   rCode = nf90_inq_varid(fid,trim(varName),vid)
   call shr_ncread_handleErr(rCode,subName//' ERROR inq varid vid')
   rCode = nf90_inquire_variable(fid,vid,ndims=ndims)
   call shr_ncread_handleErr(rCode,subName//' ERROR inquire variable ndims')
   allocate(dids(ndims))
   allocate(ns(maxn))
   rCode = nf90_inquire_variable(fid,vid,dimids=dids)
   call shr_ncread_handleErr(rCode,subName//' ERROR inquire variable dimids')

   !--- get dim sizes for all dims or to maxn, default result is 1 ---
   ns = 1
   do n=1,min(ndims,maxn)
     rcode = nf90_inquire_dimension(fid,dids(n),len=ns(n))
     call shr_ncread_handleErr(rCode, subName//" ERROR inquire dimension")
   enddo

   call shr_ncread_close(fid,rCode)

   !--- copy to output optional arguments ---
   if (present(n1)) then
     n1 = ns(1)
   endif
   if (present(n2)) then
     n2 = ns(2)
   endif
   if (present(n3)) then
     n3 = ns(3)
   endif
   if (present(n4)) then
     n4 = ns(4)
   endif
   if (present(n5)) then
     n5 = ns(5)
   endif
   if (present(n6)) then
     n6 = ns(6)
   endif

   deallocate(dids)
   deallocate(ns)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_varDimSizes
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_dimSizeName -- return size of dimension
!
! !DESCRIPTION:
! Returns the size of a named dimension
! \newline
! General Usage:
!    call shr_ncread_dimSize('myfile','lon',nsize)
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_dimSizeName(fileName, dimName, ns, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fileName ! nc file name
   character(*)        ,intent(in)           :: dimName  ! name of dimension
   integer(SHR_KIND_IN),intent(out)          :: ns       ! size of dimension
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: fid       ! file id
   integer(SHR_KIND_IN) :: did       ! dim id
   integer(SHR_KIND_IN) :: rCode     ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_dimSizeName)"
   character(*),parameter :: F00     = "('(shr_ncread_dimSizeName) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_dimSizeName) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   call shr_ncread_open(fileName,fid,rCode)

   !--- read coordinate dimensions ---
   rcode = nf90_inq_dimid        (fid, trim(dimName), did)  ! size of dimension
   call shr_ncread_handleErr(rCode, subName//" ERROR inq dimid")
   rcode = nf90_inquire_dimension(fid,did,len=ns)
   call shr_ncread_handleErr(rCode, subName//" ERROR inquire dimension")
   if (debug > 1 .and. s_loglev > 0) write(s_logunit,F01) trim(dimName)//' dimension has size = ',ns

   call shr_ncread_close(fid,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_dimSizeName
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_domain -- read in domain info from a file
!
! !DESCRIPTION:
!     Read in domain information from a file.  The subroutine is designed
!     specially for certain criteria in CCSM.  longitude, latitude, mask,
!     a area arrays will be read in from the cdf file for each variable
!     name input.  Across the subroutine interface, all arrays are 2d,
!     and each is real*8 except mask which is an integer array.  Within
!     the netcdf file, other scenarios are possible.
! note:
! o always returns 2d lat/lon arrays even if data is 1d in netCDF file
! o works if lat & lon are dimensions or variables
! o assumes area and mask are variables without a time dimension
! o mask is an integer array, all others are real*8
! o mask is read as real*8 array then copied via nint
! o assumes arrays are already allocated by the caller
!
! \newline
! General Usage:
!    call shr_ncread_domain('myfile','xc',lon,'yc',lat,'mask',imask,'area',area)
! \newline
! !REVISION HISTORY:
!     2005-Apr-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_domain(fn, lonName,  lon,  latName,  lat, &
                             &  maskName, mask, areaName, area, &
                             &  fracName, frac, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)          :: fn        ! nc file name
   character(*)        ,intent(in)          :: lonName   ! name of longitude var
   real(SHR_KIND_R8)   ,intent(out)         :: lon(:,:)  ! longitudes
   character(*)        ,intent(in)          :: latName   ! name of latitude var
   real(SHR_KIND_R8)   ,intent(out)         :: lat(:,:)  ! latitudes
   character(*)        ,intent(in) ,optional:: maskName  ! name of mask var
   integer(SHR_KIND_IN),intent(out),optional:: mask(:,:) ! domain mask
   character(*)        ,intent(in) ,optional:: areaName  ! name of area var
   real(SHR_KIND_R8)   ,intent(out),optional:: area(:,:) ! cell area
   character(*)        ,intent(in) ,optional:: fracName  ! name of frac var
   real(SHR_KIND_R8)   ,intent(out),optional:: frac(:,:) ! cell frac
   integer(SHR_KIND_IN),intent(out),optional:: rc        ! return code

!EOP
   !----- local -----
   real(SHR_KIND_R8),allocatable :: A4d(:,:,:,:)  ! local 4d array
   real(SHR_KIND_R8),allocatable :: P2d(:,:)      ! pointer to 2d arrays
   character(SHR_KIND_CS)        :: varName       ! var name
   integer(SHR_KIND_IN)          :: nflds         ! number of flds to read
   integer(SHR_KIND_IN)          :: n,i,j         ! counters
   integer(SHR_KIND_IN)          :: ndim,nd1,nd2  ! dims and size of 2 dims for cdf field
   integer(SHR_KIND_IN)          ::      pd1,pd2  ! size of 2 dims for P2d
   integer(SHR_KIND_IN)          :: rCode         ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_domain)"
   character(*),parameter :: F00     = "('(shr_ncread_domain) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_domain) ',2a,3i6,2x,a)"
   character(*),parameter :: F02     = "('(shr_ncread_domain) ',a,i6)"
   character(*),parameter :: F03     = "('(shr_ncread_domain) ',a,2i6)"
   character(*),parameter :: F04     = "('(shr_ncread_domain) ',a,2g17.8)"

   logical :: readmask
   logical :: readarea 
   logical :: readfrac 
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   rCode = 0

   nflds = 2 ! manditory fields

   if (present(maskName).and.present(mask)) then
      nflds = nflds + 1
   else if ((     present(maskName) .and. .not.present(mask)) .or. &
        (.not.present(maskName) .and.      present(mask))) then
      write(s_logunit,F00) ' ERROR: maskName and mask must both be present or not '
      call shr_ncread_abort(subName//' ERROR subroutine arguments, mask')
   end if

   if (present(areaName).and.present(area)) then
      nflds = nflds + 1
   else if ((     present(areaName) .and. .not.present(area)) .or. &
        (.not.present(areaName) .and.      present(area))) then
      write(s_logunit,F00) ' ERROR: areaName and area must both be present or not '
      call shr_ncread_abort(subName//' ERROR subroutine arguments, area')
   end if

   if (present(fracName).and.present(frac)) then
      nflds = nflds + 1
   else if ((     present(fracName) .and. .not.present(frac)) .or. &
        (.not.present(fracName) .and.      present(frac))) then
      write(s_logunit,F00) ' ERROR: fracName and frac must both be present or not '
      call shr_ncread_abort(subName//' ERROR subroutine arguments, frac')
   end if

   ! --- two fields hardwired ---
   readmask = .true.
   readarea = .true.
   readfrac = .true.
   do n=1,nflds
     if (n == 1) then
        varName = trim(lonName)
        allocate(P2d(size(lon,1),size(lon,2)))
     elseif (n == 2) then
        varName = trim(latName)
        allocate(P2d(size(lat,1),size(lat,2)))
     elseif (n > 2) then
        if (present(maskName) .and. readmask) then 
           varName = trim(maskName)
           !--- since mask in an integer, allocate P2d and copy back later ---
           allocate(P2d(size(mask,1),size(mask,2)))
           readmask = .false.
        else if (present(areaName) .and. readarea) then
           varName = trim(areaName)
           allocate(P2d(size(area,1),size(area,2)))
           readarea = .false.
        else if (present(fracName) .and. readfrac) then
           varName = trim(fracName)
           allocate(P2d(size(frac,1),size(frac,2)))
           readfrac = .false.
        endif
     end if

     if (.not.shr_ncread_varExists(fn,varName)) &
       call shr_ncread_abort(subName//' ERROR var does not exist '//trim(varName))

     !--- get size of input array ---
     pd1 = size(P2d,1)
     pd2 = size(P2d,2)

     !--- get var dims and check ---
     call shr_ncread_varDimNum(fn,varName,ndim)
     if (n > 2 .and. ndim /= 2) then
       write(s_logunit,F02) 'ERROR '//trim(varName)//' ndim = ',ndim
       call shr_ncread_abort(subName//' ERROR ndim must be 2 for '//trim(varName))
     elseif (ndim < 1 .or. ndim > 2) then
       write(s_logunit,F02) 'ERROR '//trim(varName)//' ndim = ',ndim
       call shr_ncread_abort(subName//' ERROR ndim must be 1 or 2 for '//trim(varName))
     endif
     nd1 = 1
     nd2 = 1
     if (ndim > 0)  call shr_ncread_varDimSize(fn,varName,1,nd1)
     if (ndim > 1)  call shr_ncread_varDimSize(fn,varName,2,nd2)

     !--- error check dimensions, special case for 1d lat  ---
     if (n == 2 .and. ndim == 1) then
       if ( nd1 /= pd2) then
         write(s_logunit,F03) ' nd1 pd2 error ',nd1,pd2
         call shr_ncread_abort(subName//' ERROR nd1 pd2 error')
       endif
     elseif (ndim > 0 .and. nd1 /= pd1) then
       write(s_logunit,F03) ' nd1 pd1 error ',nd1,pd1
       call shr_ncread_abort(subName//' ERROR nd1 pd1 error')
     endif
     if (ndim > 1 .and. nd2 /= pd2) then
       write(s_logunit,F03) ' nd2 pd2 error ',nd2,pd2
       call shr_ncread_abort(subName//' ERROR nd2 pd2 error')
     endif

     !--- allocate A4d and read ---
     allocate(A4d(nd1,nd2,1,1))
     A4d = 0.0_SHR_KIND_R8
     call shr_ncread_field4dG(fn,varName,rfld=A4d)

     !--- copy into P2d as appropriate ---
     do j = 1,pd2
     do i = 1,pd1
       if (n == 2 .and. ndim == 1) then
         P2d(i,j) = A4d(j,1,1,1)
       elseif (ndim == 1) then
         P2d(i,j) = A4d(i,1,1,1)
       else
         P2d(i,j) = A4d(i,j,1,1)
       endif
     enddo
     enddo

     !--- copy into mask R8 to IN ---
     if (n == 1) then
       lon(:,:) = P2d(:,:)
     elseif (n == 2) then
       lat(:,:) = P2d(:,:)
     elseif (n == 3) then
       mask(:,:) = nint(P2d(:,:))
     elseif (n == 4) then
       area(:,:) = P2d(:,:)
     elseif (n == 5) then
       frac(:,:) = P2d(:,:)
     endif

     !--- clean up ---
     deallocate(A4d)
     deallocate(P2d,stat=rCode)
!     nullify(P2d)
     
   enddo

   if (debug > 1 .and. s_loglev > 0) then
     write(s_logunit,F04) 'min/max lon  ',minval(lon),maxval(lon)
     write(s_logunit,F04) 'min/max lat  ',minval(lat),maxval(lat)
     write(s_logunit,F04) 'min/max mask ',minval(mask),maxval(mask)
     write(s_logunit,F04) 'min/max area ',minval(area),maxval(area)
     if (nflds >= 5 .and. s_loglev > 0) write(s_logunit,F04) 'min/max frac ',minval(frac),maxval(frac)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_ncread_domain
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tField2dR8 -- read in field data from a file
!
! !DESCRIPTION:
!     Read in field data from a netcdf file.  This is a special routine
!     built specificallly for CCSM.  The idea is to read a snapshot of
!     (possibly) time-varying data from a netcdf file.  The array is a
!     2d real*8 field in this case.  Inputs are filename, timeslice 
!     (integer), and variable name.  Optional inputs include the
!     time dimension name and the 2 dimension names for the array.
!     If dim1 is sent as an optional argument, dim2 must also be sent.
!     Otherwise, the time dimension is assumed to be the third
!     dimension and the first 2 dimensions are associated with the
!     2d array.
!
! \newline
! General Usage:
!    call shr_ncread_tField('myfile',6,'sst',a2d)
! \newline
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tField2dR8(fn, tIndex, fldName, fld, dim1, dim2, tName, fidi, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   real(SHR_KIND_R8)   ,intent(out)          :: fld(:,:) ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: dim2     ! name of dim2 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
   integer(SHR_KIND_IN),intent(in) ,optional :: fidi     ! file id
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   real(SHR_KIND_R8),allocatable :: lfld(:,:,:,:)   ! local 4d array
   integer(SHR_KIND_IN) :: rCode                    ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tField2dR8)"
   character(*),parameter :: F00     = "('(shr_ncread_tField2dR8) ',4a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
   allocate(lfld(size(fld,1),size(fld,2),1,1))

   if (present(dim1).and.present(dim2).and.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(dim1).and.present(dim2)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3=tName,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3=tName,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (.not.present(dim1).and..not.present(dim2).and..not.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   else
     call shr_ncread_abort(subName//' ERROR argument combination not supported')
   endif

   fld(:,:) = lfld(:,:,1,1)
   deallocate(lfld)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tField2dR8
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tField1dR8 -- read in field data from a file
!
! !DESCRIPTION:
!     Read in field data from a netcdf file.  This is a special routine
!     built specificallly for CCSM.  The idea is to read a snapshot of
!     (possibly) time-varying data from a netcdf file.  The array is a
!     1d real*8 field in this case.  Inputs are filename, timeslice 
!     (integer), and variable name.  Optional inputs include the
!     time dimension name and the dimension name for the array.
!     Otherwise, the time dimension is assumed to be the second
!     dimension and the first dimension is associated with the
!     1d array.
! \newline
! General Usage:
!    call shr_ncread_tField('myfile',1,'zlev',a1d)
! \newline
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tField1dR8(fn, tIndex, fldName, fld, dim1, tName, fidi, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   real(SHR_KIND_R8)   ,intent(out)          :: fld(:)   ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
   integer(SHR_KIND_IN),intent(in) ,optional :: fidi     ! file id
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   real(SHR_KIND_R8),allocatable :: lfld(:,:,:,:)    ! local 4d array
   integer(SHR_KIND_IN) :: rCode                     ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tField1dR8)"
   character(*),parameter :: F00     = "('(shr_ncread_tField1dR8) ',4a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   allocate(lfld(size(fld,1),1,1,1))

   if (present(dim1).and.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(dim1)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2=tName,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2=tName,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (.not.present(dim1).and..not.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   else
     call shr_ncread_abort(subName//' ERROR argument combination not supported')
   endif

   fld(:) = lfld(:,1,1,1)
   deallocate(lfld)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tField1dR8
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tField2dIN -- read in field data from a file
!
! !DESCRIPTION:
!     Read in field data from a netcdf file.  This is a special routine
!     built specificallly for CCSM.  The idea is to read a snapshot of
!     (possibly) time-varying data from a netcdf file.  The array is a
!     2d integer field in this case.  Inputs are filename, timeslice 
!     (integer), and variable name.  Optional inputs include the
!     time dimension name and the 2 dimension names for the array.
!     If dim1 is sent as an optional argument, dim2 must also be sent.
!     Otherwise, the time dimension is assumed to be the third
!     dimension and the first 2 dimensions are associated with the
!     2d array.
!
! \newline
! General Usage:
!    call shr_ncread_tField('myfile',1,'index',i2d)
! \newline
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tField2dIN(fn, tIndex, fldName, fld, dim1, dim2, tName, fidi, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   integer(SHR_KIND_IN),intent(out)          :: fld(:,:) ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: dim2     ! name of dim2 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
   integer(SHR_KIND_IN),intent(in) ,optional :: fidi     ! file id
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN),allocatable :: lfld(:,:,:,:)  ! local 4d array
   integer(SHR_KIND_IN) :: rCode                      ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tField2dIN)"
   character(*),parameter :: F00     = "('(shr_ncread_tField2dIN) ',4a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   allocate(lfld(size(fld,1),size(fld,2),1,1))

   if (present(dim1).and.present(dim2).and.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(dim1).and.present(dim2)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3=tName,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3=tName,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (.not.present(dim1).and..not.present(dim2).and..not.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3i=tIndex,fidi=fidi,rc=rCode)
     endif
   else
     call shr_ncread_abort(subName//' ERROR argument combination not supported')
   endif

   fld(:,:) = lfld(:,:,1,1)
   deallocate(lfld)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tField2dIN
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tField1dIN -- read in field data from a file
!
! !DESCRIPTION:
!     Read in field data from a netcdf file.  This is a special routine
!     built specificallly for CCSM.  The idea is to read a snapshot of
!     (possibly) time-varying data from a netcdf file.  The array is a
!     1d integer field in this case.  Inputs are filename, timeslice 
!     (integer), and variable name.  Optional inputs include the
!     time dimension name and the dimension name for the array.
!     Otherwise, the time dimension is assumed to be the second
!     dimension and the first dimension is associated with the
!     1d array.
! \newline
! General Usage:
!    call shr_ncread_tField('myfile',1,'klev',a1d)
! \newline
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tField1dIN(fn, tIndex, fldName, fld, dim1, tName, fidi, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   integer(SHR_KIND_IN),intent(out)          :: fld(:)   ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
   integer(SHR_KIND_IN),intent(in) ,optional :: fidi     ! file id
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN),allocatable :: lfld(:,:,:,:)  ! local 4d array
   integer(SHR_KIND_IN) :: rCode                      ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tField1dIN)"
   character(*),parameter :: F00     = "('(shr_ncread_tField1dIN) ',4a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   allocate(lfld(size(fld,1),1,1,1))

   if (present(dim1).and.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(dim1)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2=tName,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2=tName,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   elseif (.not.present(dim1).and..not.present(tName)) then
     if (.not.present(fidi)) then
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2i=tIndex,rc=rCode)
     else
       call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2i=tIndex,fidi=fidi,rc=rCode)
     endif
   else
     call shr_ncread_abort(subName//' ERROR argument combination not supported')
   endif

   fld(:) = lfld(:,1,1,1)
   deallocate(lfld)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tField1dIN
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_field4dG -- read in field data from a file
!
! !DESCRIPTION:
!     Read in field data from a cdf file, fld is 4d in this case
!     This subroutine supports the reading of 1d, 2d, 3d, or 4d
!        data through the interface as long as the calling argument
!        is explicitly 4d.  The user may need to invoke a temporary
!        4d pointer or array to use this subroutine.
!     Can read in a subset of data from a netcdf file that's up
!        to 6 dimensions large.
!     Supports real*8 and integer arrays, must specify either rfld 
!        or ifld in optional arguments
!     dimN are the dimension names associated with the 4d input array,
!        if N>4, this represents dimensions outside a 4d array which can
!        be optionally set to a specific index using dimNi
!     dimNi set the index to be used for the dimn dimension name
!
! \newline
! General Usage:
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d)
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d,dim1='lon',dim2='lat',dim3='time',dim3i=21)
!    call shr_ncread_field4dG('myfile','tracer',rfld=a4d,dim5='tracer_n',dim5i=3)
! \newline
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_field4dG(fn, fldName, rfld, ifld, &
     dim1, dim1i, dim2, dim2i, dim3, dim3i, dim4, dim4i, &
     dim5, dim5i, dim6, dim6i, fidi, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   character(*)        ,intent(in)           :: fldName  ! name of field
   real(SHR_KIND_R8)   ,intent(out),optional :: rfld(:,:,:,:) ! field array
   integer(SHR_KIND_IN),intent(out),optional :: ifld(:,:,:,:) ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim1i    ! dim1 index
   character(*)        ,intent(in) ,optional :: dim2     ! name of dim2 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim2i    ! dim2 index
   character(*)        ,intent(in) ,optional :: dim3     ! name of dim3 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim3i    ! dim3 index
   character(*)        ,intent(in) ,optional :: dim4     ! name of dim4 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim4i    ! dim4 index
   character(*)        ,intent(in) ,optional :: dim5     ! name of dim5 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim5i    ! dim5 index
   character(*)        ,intent(in) ,optional :: dim6     ! name of dim6 in fld
   integer(SHR_KIND_IN),intent(in) ,optional :: dim6i    ! dim6 index
   integer(SHR_KIND_IN),intent(in) ,optional :: fidi     ! file id
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP
   !----- local -----
   integer(SHR_KIND_IN),parameter :: maxd = 4         ! max num of dims of array
   integer(SHR_KIND_IN) :: fid                        ! file id
   integer(SHR_KIND_IN) :: vid                        ! var id
   integer(SHR_KIND_IN) :: xtype                      ! var type
   integer(SHR_KIND_IN) :: ndims                      ! number of dims
   integer(SHR_KIND_IN) :: n,n1,n2,n3,n4,k            ! counters
   integer(SHR_KIND_IN)  ,allocatable :: dimid(:)     ! dimension ids for array
   integer(SHR_KIND_IN)  ,allocatable :: dids(:)      ! dimension ids for cdf
   integer(SHR_KIND_IN)  ,allocatable :: start(:)     ! cdf start array
   integer(SHR_KIND_IN)  ,allocatable :: count(:)     ! cdf count array
   integer(SHR_KIND_IN)  ,allocatable :: len(:)       ! size of dim
   character(SHR_KIND_CS),allocatable :: name(:)      ! name of dim
   real(SHR_KIND_R8)     ,allocatable :: rin(:,:)     ! local 2d array
   integer(SHR_KIND_IN)  ,allocatable :: iin(:,:)     ! local 2d array
   integer(SHR_KIND_IN)  ,allocatable :: start2d(:)   ! start for 2d local array
   integer(SHR_KIND_IN)  ,allocatable :: count2d(:)   ! count for 2d local array
   logical :: found                                   ! search logical
   integer(SHR_KIND_IN) :: rCode                      ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_field4dG) "
   character(*),parameter :: F00   = "('(shr_ncread_field4dG) ',4a)"
   character(*),parameter :: F01   = "('(shr_ncread_field4dG) ',2a,3i6,2x,a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   !--- check that rfld or ifld is present ---
   if (present(rfld).and.present(ifld)) then
     call shr_ncread_abort(subName//'both rfld and ifld should not be sent')
   endif
   if (.not.present(rfld).and..not.present(ifld)) then
     call shr_ncread_abort(subName//'either rfld or ifld must be sent')
   endif

   if (.not.present(fidi)) then
     call shr_ncread_open(fn,fid,rCode)
   else
     fid = fidi
   endif

   !--- get variable id and ndims for vid
   rCode = nf90_inq_varid(fid,trim(fldName),vid)
   call shr_ncread_handleErr(rCode,subName//'inq varid vid: '//trim(fldName))
   rCode = nf90_inquire_variable(fid,vid,xtype=xtype,ndims=ndims)
   call shr_ncread_handleErr(rCode,subName//'inquire variable ndims: '//trim(fldName))

   !--- allocate locals
   n4 = max(ndims,maxd)
   allocate(dimid  (n4)) ; dimid   = 0
   allocate(dids   (n4)) ; dids    = 0
   allocate(name   (n4)) ; name    = ' '
   allocate(len    (n4)) ; len     = 1
   allocate(start  (n4)) ; start   = 1
   allocate(count  (n4)) ; count   = 1
   allocate(start2d(n4)) ; start2d = 1
   allocate(count2d(n4)) ; count2d = 1

   !--- get dimension info for vid
   rCode = nf90_inquire_variable(fid,vid,dimids=dids)
   call shr_ncread_handleErr(rCode,subName//'inquire variable dids: '//trim(fldName))
   do n=1,ndims
     rCode = nf90_inquire_dimension(fid,dids(n),name=name(n),len=len(n))
     call shr_ncread_handleErr(rCode,subName//'inquire dimension len: '//trim(fldName))
   enddo

   !--- set dimid from dim
   if (present(dim1)) then
     do n=1,ndims
       if (trim(dim1) == trim(name(n))) dimid(1) = n
     enddo
   endif
   if (present(dim2)) then
     do n=1,ndims
       if (trim(dim2) == trim(name(n))) dimid(2) = n
     enddo
   endif
   if (present(dim3)) then
     do n=1,ndims
       if (trim(dim3) == trim(name(n))) dimid(3) = n
     enddo
   endif
   if (present(dim4)) then
     do n=1,ndims
       if (trim(dim4) == trim(name(n))) dimid(4) = n
     enddo
   endif
   if (present(dim5)) then
     do n=1,ndims
       if (trim(dim5) == trim(name(n))) dimid(5) = n
     enddo
   endif
   if (present(dim6)) then
     do n=1,ndims
       if (trim(dim6) == trim(name(n))) dimid(6) = n
     enddo
   endif

   !--- set dimid for non user set dimension based on what's left
   do n1=1,max(maxd,ndims)
     k = 1
     do while (dimid(n1) == 0)
       found = .false.
       do n2 = 1,maxd
         if (dimid(n2) == k) found = .true.
       enddo
       if (found) then
         k = k + 1
       else
         dimid(n1) = k
       endif
     enddo
   enddo

   !--- set count to len if n exists in variable, otherwise set to 1
   do n1=1,maxd
     if (dimid(n1) <= ndims) then
       count(dimid(n1)) = len(dimid(n1))
     else
       count(dimid(n1)) = 1
     endif
   enddo

   !--- modify start and count from user inputs
   if (present(dim1i)) then
     if (dim1i < 1 .or. dim1i > len(dimid(1))) &
       call shr_ncread_abort(subName//'dim1i setting: '//trim(fldName))
     start(dimid(1)) = dim1i
     count(dimid(1)) = 1
   endif
   if (present(dim2i)) then
     if (dim2i < 1 .or. dim2i > len(dimid(2))) &
       call shr_ncread_abort(subName//'dim2i setting: '//trim(fldName))
     start(dimid(2)) = dim2i
     count(dimid(2)) = 1
   endif
   if (present(dim3i)) then
     if (dim3i < 1 .or. dim3i > len(dimid(3))) &
       call shr_ncread_abort(subName//'dim3i setting: '//trim(fldName))
     start(dimid(3)) = dim3i
     count(dimid(3)) = 1
   endif
   if (present(dim4i)) then
     if (dim4i < 1 .or. dim4i > len(dimid(4))) &
       call shr_ncread_abort(subName//'dim4i setting: '//trim(fldName))
     start(dimid(4)) = dim4i
     count(dimid(4)) = 1
   endif
   if (present(dim5i)) then
     if (dim5i < 1 .or. dim5i > len(dimid(5))) &
       call shr_ncread_abort(subName//'dim5i setting: '//trim(fldName))
     start(dimid(5)) = dim5i
     count(dimid(5)) = 1
   endif
   if (present(dim6i)) then
     if (dim6i < 1 .or. dim6i > len(dimid(6))) &
       call shr_ncread_abort(subName//'dim6i setting: '//trim(fldName))
     start(dimid(6)) = dim6i
     count(dimid(6)) = 1
   endif

   !--- error check, fld size must match variable size
   do n=1,maxd
     if (present(rfld)) then
       if (size(rfld,n) /= count(dimid(n))) then
         call shr_ncread_abort(subName//'fld size does not agree with count: '//trim(fldName))
       endif
     endif
     if (present(ifld)) then
       if (size(ifld,n) /= count(dimid(n))) then
         call shr_ncread_abort(subName//'fld size does not agree with count: '//trim(fldName))
       endif
     endif
   enddo

   !--- fill fld, prepare both int and real arrays, just in case
   !--- use rin/iin and transpose if needed
   if (dimid(1) > dimid(2)) then
     allocate(rin(count(dimid(2)),count(dimid(1))))
     allocate(iin(count(dimid(2)),count(dimid(1))))
   else
     allocate(rin(count(dimid(1)),count(dimid(2))))
     allocate(iin(count(dimid(1)),count(dimid(2))))
   endif
   start2d = start
   count2d = count
   count2d(dimid(3)) = 1
   count2d(dimid(4)) = 1
   do n4 = 1,count(dimid(4))
   do n3 = 1,count(dimid(3))
     start2d(dimid(3)) = n3 + start(dimid(3)) - 1
     start2d(dimid(4)) = n4 + start(dimid(4)) - 1
     if (present(rfld)) then
       rCode = nf90_get_var(fid,vid,rin,start=start2d,count=count2d)
     elseif (present(ifld)) then
       rCode = nf90_get_var(fid,vid,iin,start=start2d,count=count2d)
     endif
     call shr_ncread_handleErr(rCode,subName//'get var: '//trim(fldName))

!     if (debug > 1 .and. s_loglev > 0) then
!       write(s_logunit,*) subName,' size rfld',size(rfld,1),size(rfld,2), &
!                                       size(rfld,3),size(rfld,4)
!       write(s_logunit,*) subName,' size ifld',size(ifld,1),size(ifld,2), &
!                                       size(ifld,3),size(ifld,4)
!       write(s_logunit,*) subName,' size rin',size(rin,1),size(rin,2)
!       write(s_logunit,*) subName,' size iin',size(iin,1),size(iin,2)
!       write(s_logunit,*) subName,' dimid ',dimid
!       write(s_logunit,*) subName,' start ',start
!       write(s_logunit,*) subName,' count ',count
!       write(s_logunit,*) subName,' start2d ',start2d
!       write(s_logunit,*) subName,' count2d ',count2d
!       write(s_logunit,*) subName,' min/max rin ',minval(rin),maxval(rin)
!       write(s_logunit,*) subName,' min/max iin ',minval(iin),maxval(iin)
!     endif
     do n2 = 1,count(dimid(2))
       do n1 = 1,count(dimid(1))
         if (dimid(1) > dimid(2)) then
             if (present(rfld)) then
               rfld(n1,n2,n3,n4) = rin(n2,n1)
             elseif (present(ifld)) then
               ifld(n1,n2,n3,n4) = iin(n2,n1)
             endif
         else
             if (present(rfld)) then
               rfld(n1,n2,n3,n4) = rin(n1,n2)
             elseif (present(ifld)) then
               ifld(n1,n2,n3,n4) = iin(n1,n2)
             endif
         endif
       enddo
     enddo
   enddo
   enddo
   deallocate(rin)
   deallocate(iin)

   deallocate(dimid)
   deallocate(dids)
   deallocate(start)
   deallocate(count)
   deallocate(name)
   deallocate(len)
   deallocate(start2d)
   deallocate(count2d)
   if (.not.present(fidi)) then
     call shr_ncread_close(fid,rCode)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_ncread_field4dG
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_open -- Open netcdf file
!
! !DESCRIPTION:
!   Open netcdf file
!
! \newline
! General Usage:
!    call shr_ncread_open('myfile',fid)
! \newline
! !REVISION HISTORY:
!     2005-May-01 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_ncread_open(fileName,fid,rCode)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)          :: fileName
   integer(SHR_KIND_IN),intent(out) :: fid
   integer(SHR_KIND_IN),intent(out) :: rCode

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n
   logical                :: exists

   !----- formats -----
   character(*),parameter :: subName =   '(shr_ncread_open) '
   character(*),parameter :: F00     = "('(shr_ncread_open) ',4a)"

!----------------------------------------------------------------------------
! Notes: simply opens the file, does not acquire from anywhere (eg. mss:)
!----------------------------------------------------------------------------

   !--- verify the file exists ---
   inquire(file=trim(fileName),exist=exists)
   if (.not.exists) then
      if (s_loglev > 0) write(s_logunit,F00) "ERROR: file does not exist: ",trim(fileName)
      call shr_ncread_handleErr(rCode,subName//"ERROR: file does not exist: "//trim(fileName))
   end if

   !--- open the data file ---
   if (debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'open netCDF data file: ',trim(fileName)
   rCode = nf90_open(fileName,nf90_nowrite,fid)
   call shr_ncread_handleErr(rCode, subName//"ERROR opening input data file")

end subroutine shr_ncread_open
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_close -- Close netcdf file
!
! !DESCRIPTION:
!   Close  netcdf file
!
! \newline
! General Usage:
!    call shr_ncread_close(fid)
! \newline
! !REVISION HISTORY:
!     2005-May-01 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_ncread_close(fid,rCode)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: fid
   integer(SHR_KIND_IN),intent(out) :: rCode

!EOP

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_close)"
   character(*),parameter :: F00     = "('(shr_ncread_close) ',4a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   !--- close the data file ---
   if (debug > 1 .and. s_loglev > 0) write(s_logunit,F00) 'close netCDF input data file '
   rCode = nf90_close(fid)
   call shr_ncread_handleErr(rCode, subName//" ERROR closing input data file")

end subroutine shr_ncread_close
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_handleErr -- Print netCDF error message
!
! !DESCRIPTION:
!   Print the error message corresponding to the netCDF error status
!
! \newline
! General Usage:
!   call shr_ncread_handleErr(rCode,' check in xx call in subroutine yy ')
! \newline
! !REVISION HISTORY:
!     2005-Jan-31 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_handleErr(rCode, str)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent (in) :: rCode
   character(*)        ,intent (in) :: str

!EOP

   !----- formats -----
   character(*),parameter :: F00     = "('(shr_ncread_handleErr) ',4a)" 
   
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (rCode /= nf90_noerr) then
      write(s_logunit,F00) "netCDF error: ",trim(nf90_strerror(rCode))
      call shr_ncread_abort(str)
   end if

end subroutine shr_ncread_handleErr
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_setAbort -- Set local shr_ncread abort flag
!
! !DESCRIPTION:
!     Set local shr_ncread abort flag, true = abort, false = print and continue
! \newline
! General Usage:
!    call shr\_ncread\_setAbort(.false.)
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_ncread_setAbort') "
  character(*),parameter :: F00     = "('(shr_ncread_setAbort) ',a) "

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  doabort = flag

end subroutine shr_ncread_setAbort
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_setDebug -- Set local shr_ncread debug level
!
! !DESCRIPTION:
!     Set local shr_ncread debug level, 0 = production
! \newline
! General Usage:
!    call shr\_ncread\_setDebug(2)
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('shr_ncread_setDebug') "
  character(*),parameter :: F00     = "('(shr_ncread_setDebug) ',a) "

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  debug = iflag

end subroutine shr_ncread_setDebug
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_abort -- local abort call
!
! !DESCRIPTION:
!     local abort call
! \newline
! General Usage:
!    call shr\_ncread\_abort(' ERROR in subroutine xyz ')
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(IN) :: string

!EOP

  !--- local ---
  character(SHR_KIND_CL) :: lstring
  character(*),parameter :: subName =   "(shr_ncread_abort)"
  character(*),parameter :: F00     = "('(shr_ncread_abort) ',a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  lstring = ''
  if (present(string)) lstring = string

  if (doabort) then
    call shr_sys_abort(lstring)
  else
    write(s_logunit,F00) ' no abort:'//trim(lstring)
  endif

end subroutine shr_ncread_abort
!===============================================================================
!===============================================================================
end module shr_ncread_mod

