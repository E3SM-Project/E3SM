!===============================================================================
! SVN $Id: shr_ncread_mod.F90 1052 2006-05-25 22:36:57Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_ncread_mod.F90 $
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
!    call shr_ncread_attribute('myfile','sst','units',attribute)
!    call shr_ncread_tCoord('myfile','time',tvar,units,calendar)
!    call shr_ncread_tCoord('myfile','time',dates,secs)
!    call shr_ncread_domain('myfile','xc',lon,'yc',lat,'mask',imask,'area',area)
!    call shr_ncread_tField('myfile',6,'sst',a2d)
!    call shr_ncread_tField('myfile',6,'sst',a2d,'xc','yc','time')
!    call shr_ncread_tField('myfile',1,'zlev',a1d)
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d)
!    call shr_ncread_field4dG('myfile','sst',rfld=a4d,dim1='lon',dim2='lat',dim3='time',dim3i=21)
!    call shr_ncread_print('myfile')
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
   use shr_cal_mod     ! calendar
   use netcdf

   implicit none

   private              ! everything is default private

! !PUBLIC TYPES:

   ! no public data types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_ncread_varExists
   public :: shr_ncread_attExists
   public :: shr_ncread_varDimNum
   public :: shr_ncread_varDimSize
   public :: shr_ncread_varDimSizes
   public :: shr_ncread_dimSize
   public :: shr_ncread_attribute
   public :: shr_ncread_tCoord
   public :: shr_ncread_domain
   public :: shr_ncread_tField
   public :: shr_ncread_Field4dG
   public :: shr_ncread_print
   public :: shr_ncread_setAbort
   public :: shr_ncread_setDebug

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

   interface shr_ncread_tCoord  ; module procedure &
      shr_ncread_tCoordRC, &
      shr_ncread_tCoordII
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
! !IROUTINE: shr_ncread_attExists -- return logical for existance of att
!
! !DESCRIPTION:
! Returns logical if attribute exists on netcdf file
! Returns True if attribute exists, does not return the attribute
! \newline
! General Usage:
!   check = shr_ncread_attExists('myfile','sst','units')
! \newline
! !REVISION HISTORY:
!     2005-May-21 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function shr_ncread_attExists(fileName, varName, attName)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)   :: fileName ! nc file name
   character(*),intent(in)   :: varName  ! name of variable
   character(*),intent(in)   :: attName  ! name of attribute

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: fid
   integer(SHR_KIND_IN) :: vid
   integer(SHR_KIND_IN) :: debug0
   integer(SHR_KIND_IN) :: rCode

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_attExists)"
   character(*),parameter :: F00     = "('(shr_ncread_attExists) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_attExists) ',a,i6)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   !--- turn off debug writing ---
   debug0 = debug
   call shr_ncread_setDebug(0)

   shr_ncread_attExists = .false.
   call shr_ncread_open(fileName,fid,rCode)
   rCode = nf90_inq_varid(fid,trim(varName),vid)
   if (rCode == nf90_noerr) then
     rCode = nf90_inquire_attribute(fid, vid, trim(attName))
     if (rCode == nf90_noerr) shr_ncread_attExists = .true.
   endif

   call shr_ncread_close(fid,rCode)

   !--- reset debug code ---
   call shr_ncread_setDebug(debug0)

end function shr_ncread_attExists
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
   if (debug > 1) write(6,F01) trim(varName)//' has dims = ',ns

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
   if (debug > 1) write(6,F01) trim(dimName)//' dimension has size = ',ns

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
   if (debug > 1) write(6,F01) trim(dimName)//' dimension has size = ',ns

   call shr_ncread_close(fid,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_dimSizeName
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_attribute -- Get attribute
!
! !DESCRIPTION:
! Returns the character string associated with a particular attribute.
! The attribute is specified for a file, a varable name (or GLOBAL) and
! an attribute name.
!
! \newline
! General Usage:
!    call shr_ncread_attribute('myfile','time','units',attrib)
! \newline
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_ncread_attribute(fn,vName,aName,attrib,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   character(*)        ,intent(in)           :: vName    ! name of variable
   character(*)        ,intent(in)           :: aName    ! name of attribute
   character(*)        ,intent(out)          :: attrib   ! attribute
   integer(SHR_KIND_IN),intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   integer(SHR_KIND_IN)   :: n              ! loop index
   integer(SHR_KIND_IN)   :: xtype          ! cdf type
   integer(SHR_KIND_IN)   :: len            ! datatype size
   integer(SHR_KIND_IN)   :: fid            ! file id
   integer(SHR_KIND_IN)   :: vid            ! variable id
   integer(SHR_KIND_IN)   :: rCode          ! rc

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_attribute)"
   character(*),parameter :: F00     = "('(shr_ncread_attribute) ',8a)"
   character(*),parameter :: eos     = "[end-of-string]"

!-------------------------------------------------------------------------------
! Note:
! o there is a problem with netCDF attribute strings -- depending on how the
!   file was created, attribute strings may or may not be terminated with an
!   ASCII "NUL" character (decimal value = 0, hex value = 0).  This is the
!   c-language string terminator, but F90 doesn't treat this NUL char as white 
!   space or a string terminator, and cannot parse strings with a NUL char.
!   (This has been my experience on bluesky (NCAR IBM SP4, circa 2005).  My 
!   solution then is to replace any NUL char string terminator with a blank space.
!   This seems like a reasonable and safe thing to do.  - B. Kauffman, Jun 2005
!-------------------------------------------------------------------------------

   attrib = ' '
   rCode = 0
   call shr_ncread_open(fn,fid,rCode)

   rCode = nf90_inq_varid(fid, trim(vName), vid)
   call shr_ncread_handleErr(rCode,subName//' nf90_inq_var')

   rCode = nf90_inquire_attribute(fid, vid, trim(aName), xtype, len)
   call shr_ncread_handleErr(rCode,subName//' nf90_inq_att')
   if (xtype == NF90_CHAR) then
      rCode = nf90_get_att(fid, vid, trim(aName), attrib)
      n = len_trim(attrib)
      if (ichar(attrib(n:n)) == 0 ) then
         if (debug>0) then
            write(6,F00) 'removed null char from end of attribute...'
            write(6,F00) 'orig: ',trim(vName),':',trim(aName),' = ',trim(attrib),eos
         end if
         attrib(n:n) = ' '
         if (debug>0) then
            write(6,F00) 'new : ',trim(vName),':',trim(aName),' = ',trim(attrib),eos
         end if
      else   
         if (debug>0) then
            write(6,F00) 'read: ',trim(vName),':',trim(aName),' = ',trim(attrib),eos
         end if
      end if 
      call shr_ncread_handleErr(rCode,subName//' nf90_get_att attrib')
   else
      write(6,F00) 'attribute: '//trim(vName)//' '//trim(aName)//' not char'
      call shr_ncread_abort(subName//' attribute '//trim(aName)//' not char')
   endif

   call shr_ncread_close(fid,rCode)

   if (present(rc)) rc = rCode

end subroutine shr_ncread_attribute
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tCoordRC -- read in tCoord variable from a file
!
! !DESCRIPTION:
!     Read in tCoord variable from a file.  Given a filename and
!     a time variable name, will return the array of time and
!     the units and calendar attributes in character string.
! \newline
! General Usage:
!    call shr_ncread_tCoord('myfile','time',tvar,units,calendar)
! \newline
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tCoordRC(fn, tName,  tvar, units, calendar, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)          :: fn        ! nc file name
   character(*)        ,intent(in)          :: tName     ! name of time var
   real(SHR_KIND_R8)   ,intent(out)         :: tvar(:)   ! time array
   character(*)        ,intent(out)         :: units     ! units attribute
   character(*)        ,intent(out)         :: calendar  ! calendar attribute
   integer(SHR_KIND_IN),intent(out),optional:: rc        ! return code

!EOP
   !----- local -----
   real(SHR_KIND_R8),allocatable :: A1d(:)        ! local 1d array
   character(SHR_KIND_CL)        :: string        ! local string var
   integer(SHR_KIND_IN)          :: i             ! counters
   integer(SHR_KIND_IN)          :: ndim,nd1,pd1  ! dims and dim sizes
   integer(SHR_KIND_IN)          :: rCode         ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tCoordRC)"
   character(*),parameter :: F00     = "('(shr_ncread_tCoordRC) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_tCoordRC) ',2a,3i6,2x,a)"
   character(*),parameter :: F02     = "('(shr_ncread_tCoordRC) ',a,i6)"
   character(*),parameter :: F03     = "('(shr_ncread_tCoordRC) ',a,2i6)"
   character(*),parameter :: F04     = "('(shr_ncread_tCoordRC) ',a,2g17.8)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0

   if (.not.shr_ncread_varExists(fn,tName)) &
     call shr_ncread_abort(subName//' ERROR var does not exist '//trim(tName))

   !--- get size of input array ---
   pd1 = size(tvar,1)

   !--- get time dims and check ---
   call shr_ncread_varDimNum(fn,tName,ndim)
   if (ndim /= 1) then
     write(6,F02) 'ERROR '//trim(tName)//' ndim = ',ndim
     call shr_ncread_abort(subName//' ERROR ndim must be 1 for '//trim(tName))
   endif
   call shr_ncread_varDimSize(fn,tName,1,nd1)

   !--- error check dimensions ---
   if ( nd1 > pd1) then
     write(6,F03) ' nd1 pd1 error ',nd1,pd1
     call shr_ncread_abort(subName//' ERROR nd1 pd1 error')
   endif

   call shr_ncread_attribute(fn,tName,'units',units,rc=rCode)
   if (shr_ncread_attExists(fn,tName,'calendar')) then
     call shr_ncread_attribute(fn,tName,'calendar',calendar,rc=rCode)
   else
     calendar = "gregorian"    ! CF-1.0 default value
   endif

   call shr_string_leftAlign(units)
   call shr_string_leftAlign(calendar)

   allocate(A1d(nd1))
   call shr_ncread_tfield(fn,1,tName,A1d,rc=rCode)
   do i=1,nd1
     tvar(i) = A1d(i)
   enddo

   deallocate(A1d)

   if (debug > 1) then
     write(6,F04) 'min/max tvar  ',minval(tvar),maxval(tvar)
   endif

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tCoordRC

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_tCoordII -- read in tCoord variable from a file
!
! !DESCRIPTION:
!     Read in tCoord variable from a file.  Given a filename and time
!     variable name, returns a date and seconds array.  This array is
!     generated by a calendar function based on the time array,
!     units attribute, and calendar attribute on the cdf file.
! \newline
! General Usage:
!    call shr_ncread_tCoord('myfile','time',dates,secs)
! \newline
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_tCoordII(fn, tName,  dates, secs, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)          :: fn        ! nc file name
   character(*)        ,intent(in)          :: tName     ! name of time var
   integer(SHR_KIND_IN),intent(out)         :: dates(:)  ! date array
   integer(SHR_KIND_IN),intent(out)         :: secs(:)   ! seconds array
   integer(SHR_KIND_IN),intent(out),optional:: rc        ! return code

!EOP

   !----- local -----
   real(SHR_KIND_R8),allocatable :: tvar(:)       ! time variable
   character(SHR_KIND_CL)        :: cfUnits       ! CF-1.0 units attribute
   character(SHR_KIND_CL)        :: cfCalendar    ! CF-1.0 calendar attribute
   integer(SHR_KIND_IN)          :: n             ! counters
   integer(SHR_KIND_IN)          :: nd,ns,nmax    ! dim sizes
   character(SHR_KIND_CS)        :: units         ! time units (days,secs,...)
   integer(SHR_KIND_IN)          :: bdate         ! base date: calendar date
   real(SHR_KIND_R8)             :: bsec          ! base date: elapsed secs
   integer(SHR_KIND_IN)          :: ndate         ! calendar date of time value
   real(SHR_KIND_R8)             :: nsec          ! elapsed secs on calendar date
   integer(SHR_KIND_IN)          :: rCode         ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_tCoordII)"
   character(*),parameter :: F00     = "('(shr_ncread_tCoordII) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_tCoordII) ',a,i8)"
   character(*),parameter :: F02     = "('(shr_ncread_tCoordII) ',a,g17.7)"
   character(*),parameter :: eos     = "[end-of-string]"

!-------------------------------------------------------------------------------
! time coordinate values in the data file must be converted to actual calendar
! dates & elpased seconds by parsing the associated CF-1.0 time unit string, eg. 
! time = 15 and units = "days since 2002-01-01 00:00:00" 
! imply the actual calendar date is 2002-01-16 00:00:00  
!-------------------------------------------------------------------------------

   rCode = 0

   nd   = size(dates)
   ns   = size(secs)
   nmax = min(nd,ns)
   allocate(tVar(nmax))

   call shr_ncread_tCoordRC(fn,tName,tVar,cfUnits,cfCalendar,rCode)
   call shr_string_parseCFtunit(cfUnits,units,bdate,bsec)

   if (debug > 0) then
      write(6,F00) ' units      = ',trim(units)      ,eos
      write(6,F01) ' bdate      = ',bdate
      write(6,F02) ' bsec       = ',bsec
      write(6,F00) ' cfUnits    = ',trim(cfUnits)    ,eos
      write(6,F00) ' cfCalendar = ',trim(cfCalendar) ,eos
   endif

   do n = 1,nMax
      call shr_cal_advDate(tVar(n),units,bDate,bSec,nDate,nSec,cfCalendar)
      dates(n) = nDate
      secs (n) = nSec
   enddo

   deallocate(tVar) ! F90 may not dealloc this local array

   if (present(rc)) rc = rCode

end subroutine shr_ncread_tCoordII

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
   character(*)        ,intent(in)          :: maskName  ! name of mask var
   integer(SHR_KIND_IN),intent(out)         :: mask(:,:) ! domain mask
   character(*)        ,intent(in)          :: areaName  ! name of area var
   real(SHR_KIND_R8)   ,intent(out)         :: area(:,:) ! cell area
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

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   rCode = 0

   if (.not.present(fracName) .and. .not.present(frac)) then
      nflds = 4
   elseif (present(fracName).and.present(frac)) then
      nflds = 5
   else
      write(6,F00) ' ERROR: fracName and frac must both be present or not '
      call shr_ncread_abort(subName//' ERROR subroutine arguments, frac')
   endif

   ! --- four fields hardwired ---
   do n=1,nflds
     if (n == 1) then
       varName = trim(lonName)
       allocate(P2d(size(lon,1),size(lon,2)))
     elseif (n == 2) then
       varName = trim(latName)
       allocate(P2d(size(lat,1),size(lat,2)))
     elseif (n == 3) then
       varName = trim(maskName)
       !--- since mask in an integer, allocate P2d and copy back later ---
       allocate(P2d(size(mask,1),size(mask,2)))
     elseif (n == 4) then
       varName = trim(areaName)
       allocate(P2d(size(area,1),size(area,2)))
     elseif (n == 5) then
       varName = trim(fracName)
       allocate(P2d(size(frac,1),size(frac,2)))
     endif

     if (.not.shr_ncread_varExists(fn,varName)) &
       call shr_ncread_abort(subName//' ERROR var does not exist '//trim(varName))

     !--- get size of input array ---
     pd1 = size(P2d,1)
     pd2 = size(P2d,2)

     !--- get var dims and check ---
     call shr_ncread_varDimNum(fn,varName,ndim)
     if (n > 2 .and. ndim /= 2) then
       write(6,F02) 'ERROR '//trim(varName)//' ndim = ',ndim
       call shr_ncread_abort(subName//' ERROR ndim must be 2 for '//trim(varName))
     elseif (ndim < 1 .or. ndim > 2) then
       write(6,F02) 'ERROR '//trim(varName)//' ndim = ',ndim
       call shr_ncread_abort(subName//' ERROR ndim must be 1 or 2 for '//trim(varName))
     endif
     nd1 = 1
     nd2 = 1
     if (ndim > 0)  call shr_ncread_varDimSize(fn,varName,1,nd1)
     if (ndim > 1)  call shr_ncread_varDimSize(fn,varName,2,nd2)

     !--- error check dimensions, special case for 1d lat  ---
     if (n == 2 .and. ndim == 1) then
       if ( nd1 /= pd2) then
         write(6,F03) ' nd1 pd2 error ',nd1,pd2
         call shr_ncread_abort(subName//' ERROR nd1 pd2 error')
       endif
     elseif (ndim > 0 .and. nd1 /= pd1) then
       write(6,F03) ' nd1 pd1 error ',nd1,pd1
       call shr_ncread_abort(subName//' ERROR nd1 pd1 error')
     endif
     if (ndim > 1 .and. nd2 /= pd2) then
       write(6,F03) ' nd2 pd2 error ',nd2,pd2
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

   if (debug > 1) then
     write(6,F04) 'min/max lon  ',minval(lon),maxval(lon)
     write(6,F04) 'min/max lat  ',minval(lat),maxval(lat)
     write(6,F04) 'min/max mask ',minval(mask),maxval(mask)
     write(6,F04) 'min/max area ',minval(area),maxval(area)
     if (nflds >= 5) write(6,F04) 'min/max frac ',minval(frac),maxval(frac)
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

subroutine shr_ncread_tField2dR8(fn, tIndex, fldName, fld, dim1, dim2, tName, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   real(SHR_KIND_R8)   ,intent(out)          :: fld(:,:) ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: dim2     ! name of dim2 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
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
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,rc=rCode)
   elseif (present(dim1).and.present(dim2)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,rc=rCode)
   elseif (present(tName)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3=tName,dim3i=tIndex,rc=rCode)
   elseif (.not.present(dim1).and..not.present(dim2).and..not.present(tName)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim3i=tIndex,rc=rCode)
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

subroutine shr_ncread_tField1dR8(fn, tIndex, fldName, fld, dim1, tName, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   real(SHR_KIND_R8)   ,intent(out)          :: fld(:)   ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
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
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,rc=rCode)
   elseif (present(dim1)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim1=dim1,dim2i=tIndex,rc=rCode)
   elseif (present(tName)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2=tName,dim2i=tIndex,rc=rCode)
   elseif (.not.present(dim1).and..not.present(tName)) then
     call shr_ncread_field4dG(fn,fldName,rfld=lfld,dim2i=tIndex,rc=rCode)
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

subroutine shr_ncread_tField2dIN(fn, tIndex, fldName, fld, dim1, dim2, tName, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   integer(SHR_KIND_IN),intent(out)          :: fld(:,:) ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: dim2     ! name of dim2 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
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
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3=tName,dim3i=tIndex,rc=rCode)
   elseif (present(dim1).and.present(dim2)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=dim2,dim3i=tIndex,rc=rCode)
   elseif (present(tName)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3=tName,dim3i=tIndex,rc=rCode)
   elseif (.not.present(dim1).and..not.present(dim2).and..not.present(tName)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim3i=tIndex,rc=rCode)
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

subroutine shr_ncread_tField1dIN(fn, tIndex, fldName, fld, dim1, tName, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)           :: fn       ! nc file name
   integer(SHR_KIND_IN),intent(in)           :: tIndex   ! time-coord index
   character(*)        ,intent(in)           :: fldName  ! name of field
   integer(SHR_KIND_IN),intent(out)          :: fld(:)   ! field array
   character(*)        ,intent(in) ,optional :: dim1     ! name of dim1 in fld
   character(*)        ,intent(in) ,optional :: tName    ! name of tIndex dim
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
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2=tName,dim2i=tIndex,rc=rCode)
   elseif (present(dim1)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim1=dim1,dim2i=tIndex,rc=rCode)
   elseif (present(tName)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2=tName,dim2i=tIndex,rc=rCode)
   elseif (.not.present(dim1).and..not.present(tName)) then
     call shr_ncread_field4dG(fn,fldName,ifld=lfld,dim2i=tIndex,rc=rCode)
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
     dim5, dim5i, dim6, dim6i, rc)

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

   call shr_ncread_open(fn,fid,rCode)

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

!     if (debug > 1) then
!       write(6,*) subName,' size rfld',size(rfld,1),size(rfld,2), &
!                                       size(rfld,3),size(rfld,4)
!       write(6,*) subName,' size ifld',size(ifld,1),size(ifld,2), &
!                                       size(ifld,3),size(ifld,4)
!       write(6,*) subName,' size rin',size(rin,1),size(rin,2)
!       write(6,*) subName,' size iin',size(iin,1),size(iin,2)
!       write(6,*) subName,' dimid ',dimid
!       write(6,*) subName,' start ',start
!       write(6,*) subName,' count ',count
!       write(6,*) subName,' start2d ',start2d
!       write(6,*) subName,' count2d ',count2d
!       write(6,*) subName,' min/max rin ',minval(rin),maxval(rin)
!       write(6,*) subName,' min/max iin ',minval(iin),maxval(iin)
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
   call shr_ncread_close(fid,rCode)

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
   character(SHR_KIND_CL) :: localFN
   logical                :: exists

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_open)"
   character(*),parameter :: F00     = "('(shr_ncread_open) ',4a)"

!----------------------------------------------------------------------------
! Notes:
! o as per shr_file_get(), the file name format is expected to be
!   fileName = [location:][directory path]localFn
!   eg. "foobar.nc"  "/home/user/foobar.nc"  "mss:/USER/foobar.nc"
!----------------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! get the input data file into the cwd
   !-------------------------------------------------------------------------
   n = shr_string_lastIndex(fileName,"/")
   if (n==0) n = index(fileName,":")
   if (n==0) then
      localFn = fileName
   else
      localFn = fileName(n+1: len_trim(fileName) )
   end if
   inquire(file=trim(localFn),exist=exists)
   if (.not.exists) call shr_file_get(rCode,localFn,fileName)

   !--- open the data file ---
   if (debug > 1) write(6,F00) 'open netCDF input data file: ',trim(localFn)
   rCode = nf90_open(localFn,nf90_nowrite,fid)
   call shr_ncread_handleErr(rCode, subName//" ERROR opening input data file")

end subroutine shr_ncread_open
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_attrbs -- Write Attribute info about netcdf file
!
! !DESCRIPTION:
!   Write attribute info about netcdf file, for a variable or global
!
! \newline
! General Usage:
!    call shr_ncread_attrbs(fid,nf90_global)
!    call shr_ncread_attrbs(fid,'sst')
! \newline
! !REVISION HISTORY:
!     2005-May-01 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_ncread_attrbs(fid,vn)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: fid
   integer(SHR_KIND_IN),intent(in)  :: vn

!EOP

   !----- local -----
   integer(SHR_KIND_IN),parameter :: strlen = SHR_KIND_CL
   character(len=strlen)  :: cvalue         ! long string
   integer(SHR_KIND_IN)   :: natt
   integer(SHR_KIND_IN)   :: xtype          ! cdf type
   integer(SHR_KIND_IN)   :: an             ! counter
   integer(SHR_KIND_IN)   :: len            ! datatype size
   integer(SHR_KIND_IN)   :: rCode          ! rc
   character(SHR_KIND_CS) :: name           ! name

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_attrbs)"
   character(*),parameter :: F00     = "('(shr_ncread_attrbs) ',4a)"
   character(*),parameter :: F04     = "('(shr_ncread_attrbs) ',4x,a,i4,2a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (vn == nf90_global) then
     rCode = nf90_inquire(fid,nAttributes=natt)
   else
     rCode = nf90_inquire_variable(fid,vn,nAtts=natt)
   endif

   do an=1,natt
     rCode = nf90_inq_attname(fid,VN,an,name)
     call shr_ncread_handleErr(rCode,subName//' nf90_inq_attname')
     rCode = nf90_inquire_attribute(fid, VN, trim(name), xtype, len)
     call shr_ncread_handleErr(rCode,subName//' nf90_inq_att')
     if (xtype == NF90_CHAR) then
       if (len < strlen) then
         cvalue = ' '
         rCode = nf90_get_att(fid,VN,trim(name),cvalue)
         call shr_ncread_handleErr(rCode,subName//' nf90_get_att cvalue')
         write(6,F04) 'attribute: ',an,' '//trim(name)//':',trim(cvalue)
       else
         write(6,F04) 'attribute: ',an,' '//trim(name)//':','*** too long ***'
       endif
     else
       write(6,F04) 'attribute: ',an,' '//trim(name)//':',' *** not char ***'
     endif
   enddo

end subroutine shr_ncread_attrbs
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
   if (debug > 1) write(6,F00) 'close netCDF input data file '
   rCode = nf90_close(fid)
   call shr_ncread_handleErr(rCode, subName//" ERROR closing input data file")

end subroutine shr_ncread_close
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_ncread_print -- Print info about netcdf file
!
! !DESCRIPTION:
!   Print info about netcdf file
!
! \newline
! General Usage:
!    call shr_ncread_print('myfile')
! \newline
! !REVISION HISTORY:
!     2005-May-01 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_ncread_print(fileName, rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent (in) :: fileName
   integer(SHR_KIND_IN),optional,intent (in) :: rc

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: fid,vid,did     ! file, var, and dim id
   integer(SHR_KIND_IN) :: nvar,ndim,natt  ! numver of vars, dims, atts
   integer(SHR_KIND_IN) :: vn,dn,an        ! counters for var, dim, att
   integer(SHR_KIND_IN) :: xtype           ! field type
   integer(SHR_KIND_IN) :: len             ! size
   character(SHR_KIND_CS) :: name          ! name
   integer(SHR_KIND_IN) :: debug0          ! debug holder
   integer(SHR_KIND_IN) :: rCode           ! error code

   !----- formats -----
   character(*),parameter :: subName = "(shr_ncread_print)"
   character(*),parameter :: F00     = "('(shr_ncread_print) ',4a)"
   character(*),parameter :: F01     = "('(shr_ncread_print) ',2x,a,i4,a,i4,a)"
   character(*),parameter :: F02     = "('(shr_ncread_print) ',2x,2a,i4,2i12,i8)"
   character(*),parameter :: F03     = "('(shr_ncread_print) ',2x,2a,i4,2g16.7,i8)"
   character(*),parameter :: F04     = "('(shr_ncread_print) ',2x,a,i4,2a)"
   character(*),parameter :: F05     = "('(shr_ncread_print) ',2x,a,3i6)"
   character(*),parameter :: F06     = "('(shr_ncread_print) ',2x,a,i4,a,i4,2a,i4)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
   
   debug0 = debug
   call shr_ncread_setDebug(2)

   write(6,F00) fileName

   call shr_ncread_open(fileName,fid,rCode)

   rCode = nf90_inquire(fid,ndim,nvar,natt)
   call shr_ncread_handleErr(rCode,subName//' nf90_inquire')

   write(6,F05) 'ndim,nvar,natt: ',ndim,nvar,natt

   call shr_ncread_attrbs(fid,nf90_global)

   do dn=1,ndim
     rcode = nf90_inquire_dimension(fid,dn,name,len)
     call shr_ncread_handleErr(rCode,subName//' nf90_inquire_dim')
     write(6,F01) 'dimension: ',dn,' '//trim(name)//'(',len,')'
   enddo

   do vn=1,nvar
     rcode = nf90_inquire_variable(fid,vn,name,xtype,len,nAtts=natt)
     write(6,F06) 'variable: ',vn,' '//trim(name),len,' dims',' xtype=',xtype
     call shr_ncread_attrbs(fid,vn)
   enddo

   call shr_ncread_close(fid,rCode)

   call shr_ncread_setDebug(debug0)

end subroutine shr_ncread_print
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
      write(6,F00) "netCDF error: ",trim(nf90_strerror(rCode))
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
    write(6,F00) ' no abort:'//trim(lstring)
  endif

end subroutine shr_ncread_abort
!===============================================================================
!===============================================================================
end module shr_ncread_mod

