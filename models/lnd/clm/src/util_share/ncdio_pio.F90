module ncdio_pio

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdio_pioMod
!
! !DESCRIPTION:
! Generic interfaces to write fields to netcdf files for CLM
!
! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use spmdMod        , only : masterproc, mpicom, iam, npes,  &
                              MPI_REAL8, MPI_INTEGER, MPI_LOGICAL
  use clmtype        , only : grlnd, nameg, namel, namec, namep
  use clm_varcon     , only : spval,ispval
  use clm_varctl     , only : single_column, iulog
  use shr_sys_mod    , only : shr_sys_flush
  use shr_file_mod   , only : shr_file_getunit, shr_file_freeunit
  use shr_string_mod , only : shr_string_toUpper
  use abortutils     , only : endrun
  use decompMod      , only : get_clmlevel_gsize,get_clmlevel_gsmap
  use perf_mod       , only : t_startf, t_stopf
  use fileutils      , only : getavu, relavu
  use nanMod         , only : nan
  use mct_mod
  use pio
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: check_var          ! determine if variable is on netcdf file
  public :: check_dim          ! validity check on dimension
  public :: ncd_pio_openfile   ! open a file
  public :: ncd_pio_createfile ! create a new file
  public :: ncd_pio_closefile  ! close a file
  public :: ncd_pio_init       ! called from clm_comp
  public :: ncd_enddef         ! end define mode
  public :: ncd_putatt         ! put attribute
  public :: ncd_defdim         ! define dimension
  public :: ncd_inqdid         ! inquire dimension id
  public :: ncd_inqdname       ! inquire dimension name
  public :: ncd_inqdlen        ! inquire dimension length
  public :: ncd_inqfdims       ! inquire file dimnesions 
  public :: ncd_defvar         ! define variables
  public :: ncd_inqvid         ! inquire variable id
  public :: ncd_inqvname       ! inquire variable name
  public :: ncd_inqvdims       ! inquire variable ndims
  public :: ncd_inqvdids       ! inquire variable dimids
  public :: ncd_io             ! write local data

  integer,parameter,public :: ncd_int       = pio_int
  integer,parameter,public :: ncd_log       =-pio_int
  integer,parameter,public :: ncd_float     = pio_real
  integer,parameter,public :: ncd_double    = pio_double
  integer,parameter,public :: ncd_char      = pio_char
  integer,parameter,public :: ncd_global    = pio_global
  integer,parameter,public :: ncd_write     = pio_write
  integer,parameter,public :: ncd_nowrite   = pio_nowrite
  integer,parameter,public :: ncd_clobber   = pio_clobber
  integer,parameter,public :: ncd_noclobber = pio_noclobber
  integer,parameter,public :: ncd_nofill    = pio_nofill
  integer,parameter,public :: ncd_unlimited = pio_unlimited

  ! PIO types needed for ncdio_pio interface calls
  public file_desc_t
  public var_desc_t

!
! !REVISION HISTORY:
!
!
! !PRIVATE MEMBER FUNCTIONS:
!

  interface ncd_putatt
     module procedure ncd_putatt_int
     module procedure ncd_putatt_real
     module procedure ncd_putatt_char
  end interface

  interface ncd_defvar
     module procedure ncd_defvar_bynf
     module procedure ncd_defvar_bygrid
  end interface

  interface ncd_io 
     module procedure ncd_io_log_var0_nf
     module procedure ncd_io_int_var0_nf
     module procedure ncd_io_real_var0_nf
     module procedure ncd_io_int_var1_nf
     module procedure ncd_io_log_var1_nf
     module procedure ncd_io_real_var1_nf
     module procedure ncd_io_char_var1_nf
     module procedure ncd_io_int_var2_nf
     module procedure ncd_io_real_var2_nf
     module procedure ncd_io_char_var2_nf
     module procedure ncd_io_char_var3_nf
     module procedure ncd_io_char_varn_strt_nf
     module procedure ncd_io_int_var1
     module procedure ncd_io_log_var1
     module procedure ncd_io_real_var1
     module procedure ncd_io_int_var2
     module procedure ncd_io_real_var2
     module procedure ncd_io_int_var3
     module procedure ncd_io_real_var3
     module procedure ncd_io_real_var3_nf
  end interface

  private :: ncd_getiodesc      ! obtain iodesc
  private :: scam_field_offsets ! get offset to proper lat/lon gridcell for SCAM

  integer,parameter,private :: debug = 0             ! local debug level

  integer , parameter  , public  :: max_string_len = 256     ! length of strings
  real(r8), parameter  , public  :: fillvalue = 1.e36_r8     ! fill value for netcdf fields

  integer, public :: io_type

  type(iosystem_desc_t), pointer, public  :: pio_subsystem

  type iodesc_plus_type
     character(len=64) :: name
     type(IO_desc_t)   :: iodesc
     integer           :: type
     integer           :: ndims
     integer           :: dims(4)
     integer           :: dimids(4) 
  end type iodesc_plus_type
  integer,parameter      ,private :: max_iodesc = 100
  integer                ,private :: num_iodesc = 0
  type(iodesc_plus_type) ,private, target :: iodesc_list(max_iodesc)

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_pio_init
!
! !INTERFACE:
  subroutine ncd_pio_init()
!
! !DESCRIPTION:
! Initial PIO
!
! !USES:
    use shr_pio_mod     , only : shr_pio_getiosys, shr_pio_getiotype
    use clm_varctl      , only : inst_name
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
    character(len=*),parameter :: subname='ncd_pio_init' ! subroutine name
!-----------------------------------------------------------------------
!EOP
    PIO_subsystem => shr_pio_getiosys(inst_name)
    io_type       =  shr_pio_getiotype(inst_name)

  end subroutine ncd_pio_init

  !=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_pio_openfile
!
! !INTERFACE:
  subroutine ncd_pio_openfile(file, fname, mode)
!
! !DESCRIPTION:
! Open a NetCDF PIO file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: file   ! Output PIO file handle
    character(len=*) , intent(in)    :: fname  ! Input filename to open
    integer          , intent(in)    :: mode   ! file mode
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
    integer :: ierr
    character(len=*),parameter :: subname='ncd_pio_openfile' ! subroutine name
!-----------------------------------------------------------------------
!EOP
    ierr = pio_openfile(pio_subsystem, file, io_type, fname, mode)

    if(ierr/= PIO_NOERR) then
       call endrun(subname//'ERROR: Failed to open file')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    end if

  end subroutine ncd_pio_openfile

  !=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_pio_closefile
!
! !INTERFACE:
  subroutine ncd_pio_closefile(file)
!
! !DESCRIPTION:
! Close a NetCDF PIO file
!
! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: file   ! PIO file handle to close
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!
! !LOCAL VARIABLES:
!-----------------------------------------------------------------------
!EOP

    call pio_closefile(file)

  end subroutine ncd_pio_closefile

  !=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_pio_createfile
!
! !INTERFACE:
  subroutine ncd_pio_createfile(file, fname)
!
! !DESCRIPTION:
! Create a new NetCDF file with PIO
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: file    ! PIO file descriptor
    character(len=*),  intent(in)    :: fname   ! File name to create
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
    integer :: ierr
    character(len=*),parameter :: subname='ncd_pio_createfile' ! subroutine name
!-----------------------------------------------------------------------
!EOP
    ierr = pio_createfile(pio_subsystem, file, io_type, fname, ior(PIO_CLOBBER,PIO_64BIT_OFFSET))

    if(ierr/= PIO_NOERR) then
       call endrun( subname//' ERROR: Failed to open file to write: '//trim(fname))
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened file ', trim(fname),  ' to write', file%fh
    end if

  end subroutine ncd_pio_createfile

  !=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var
!
! !INTERFACE:
  subroutine check_var(ncid, varname, vardesc, readvar, print_err )
!
! !DESCRIPTION:
! Check if variable is on netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)  :: ncid      ! PIO file descriptor
    character(len=*) , intent(in)     :: varname   ! Varible name to check
    type(Var_desc_t) , intent(out)    :: vardesc   ! Output variable descriptor
    logical          , intent(out)    :: readvar   ! If variable exists or not
    logical, optional, intent(in)     :: print_err ! If should print about error
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ret     ! return value
    logical :: log_err ! if should log error
    character(len=*),parameter :: subname='check_var' ! subroutine name
!-----------------------------------------------------------------------


    if ( present(print_err) )then
       log_err = print_err
    else
       log_err = .true.
    end if
    readvar = .true.
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ret = PIO_inq_varid (ncid, varname, vardesc)
    if (ret /= PIO_noerr) then
       readvar = .false.
       if (masterproc .and. log_err) &
           write(iulog,*) subname//': variable ',trim(varname),' is not on dataset'
    end if
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    
  end subroutine check_var

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dim
!
! !INTERFACE:
  subroutine check_dim(ncid, dimname, value)
!
! !DESCRIPTION:
! Validity check on dimension
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(in) :: ncid      ! PIO file handle
    character(len=*), intent(in) :: dimname   ! Dimension name
    integer, intent(in)          :: value     ! Expected dimension size
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: dimid, dimlen    ! temporaries
    integer :: status           ! error code      
    character(len=*),parameter :: subname='check_dim' ! subroutine name
!-----------------------------------------------------------------------

    status = pio_inq_dimid (ncid, trim(dimname), dimid)
    status = pio_inq_dimlen (ncid, dimid, dimlen)
    if (dimlen /= value) then
       write(iulog,*) subname//' ERROR: mismatch of input dimension ',dimlen, &
            ' with expected value ',value,' for variable ',trim(dimname)
       call endrun()
    end if

  end subroutine check_dim
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_enddef
!
! !INTERFACE:
  subroutine ncd_enddef(ncid)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status   ! error status
    character(len=*),parameter :: subname='ncd_enddef' ! subroutine name
!-----------------------------------------------------------------------

    status = PIO_enddef(ncid)

  end subroutine ncd_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdid
!
! !INTERFACE:
  subroutine ncd_inqdid(ncid,name,dimid,dimexist)
!
! !DESCRIPTION:
! inquire on a dimension id
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! netcdf file id
    character(len=*), intent(in) :: name      ! dimension name
    integer         , intent(out):: dimid     ! dimension id
    logical,optional, intent(out):: dimexist  ! if this dimension exists or not
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    if ( present(dimexist) )then
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    end if
    status = PIO_inq_dimid(ncid,name,dimid)
    if ( present(dimexist) )then
       if ( status == PIO_NOERR)then
          dimexist = .true.
       else
          dimexist = .false.
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    end if

  end subroutine ncd_inqdid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdlen
!
! !INTERFACE:
  subroutine ncd_inqdlen(ncid,dimid,len,name)
!
! !DESCRIPTION:
! enddef netcdf file
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid       ! netcdf file id
    integer          , intent(inout) :: dimid      ! dimension id
    integer          , intent(out)   :: len        ! dimension len
    character(len=*), optional, intent(in) :: name ! dimension name
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    if ( present(name) )then
       call ncd_inqdid(ncid,name,dimid)
    end if
    len = -1
    status = PIO_inq_dimlen(ncid,dimid,len)

  end subroutine ncd_inqdlen

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqdname
!
! !INTERFACE:
  subroutine ncd_inqdname(ncid,dimid,dname)
!
! !DESCRIPTION:
! inquire dim name
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    integer          , intent(in) :: dimid     ! dimension id
    character(len=*) , intent(out):: dname     ! dimension name
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_inq_dimname(ncid,dimid,dname)

  end subroutine ncd_inqdname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqfdims
!
! !INTERFACE:
  subroutine ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

!
! !ARGUMENTS:
    type(file_desc_t), intent(inout):: ncid
    logical          , intent(out)  :: isgrid2d
    integer          , intent(out)  :: ni
    integer          , intent(out)  :: nj
    integer          , intent(out)  :: ns
! !LOCAL VARIABLES:
!EOP
    integer  :: dimid                                ! netCDF id
    integer  :: ier                                  ! error status 
    character(len=32) :: subname = 'surfrd_filedims' ! subroutine name
!-----------------------------------------------------------------------

    if (single_column) then
       ni = 1
       nj = 1
       ns = 1
       isgrid2d = .true.
       RETURN
    end if

    ni = 0
    nj = 0

    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ier = pio_inq_dimid (ncid, 'lon', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
    ier = pio_inq_dimid (ncid, 'lat', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)

    ier = pio_inq_dimid (ncid, 'lsmlon', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
    ier = pio_inq_dimid (ncid, 'lsmlat', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)

    ier = pio_inq_dimid (ncid, 'ni', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, ni)
    ier = pio_inq_dimid (ncid, 'nj', dimid)
    if (ier == PIO_NOERR) ier = pio_inq_dimlen(ncid, dimid, nj)

    ier = pio_inq_dimid (ncid, 'gridcell', dimid)
    if (ier == PIO_NOERR) then
       ier = pio_inq_dimlen(ncid, dimid, ni)
       if (ier == PIO_NOERR) nj = 1
    end if

    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

    if (ni == 0 .or. nj == 0) then
       write(iulog,*) trim(subname),' ERROR: ni,nj = ',ni,nj,' cannot be zero '
       call endrun()
    end if

    if (nj == 1) then
       isgrid2d = .false.
    else
       isgrid2d = .true.
    end if

    ns = ni*nj

  end subroutine ncd_inqfdims

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvid
!
! !INTERFACE:
  subroutine ncd_inqvid(ncid,name,varid,vardesc,readvar)
!
! !DESCRIPTION:
! Inquire on a variable ID
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: name      ! variable name
    integer          , intent(out)   :: varid     ! variable id
    type(Var_desc_t) , intent(out)   :: vardesc   ! variable descriptor
    logical, optional, intent(out)   :: readvar   ! does variable exist
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ret               ! return code
    character(len=*),parameter :: subname='ncd_inqvid' ! subroutine name
!-----------------------------------------------------------------------

    if (present(readvar)) then
       readvar = .false.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ret = PIO_inq_varid(ncid,name,vardesc)
       if (ret /= PIO_noerr) then
          if (masterproc) write(iulog,*) subname//': variable ',trim(name),' is not on dataset'
          readvar = .false.
       else
          readvar = .true.
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    else
       ret = PIO_inq_varid(ncid,name,vardesc)
    endif
    varid = vardesc%varid
 
  end subroutine ncd_inqvid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdims
!
! !INTERFACE:
  subroutine ncd_inqvdims(ncid,ndims,vardesc)
!
! !DESCRIPTION:
! inquire variable dimensions
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
    integer          , intent(out)  :: ndims     ! variable ndims
    type(Var_desc_t) , intent(inout):: vardesc   ! variable descriptor
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    ndims = -1
    status = PIO_inq_varndims(ncid,vardesc,ndims)

  end subroutine ncd_inqvdims

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvname
!
! !INTERFACE:
  subroutine ncd_inqvname(ncid,varid,vname,vardesc)
!
! !DESCRIPTION:
! inquire variable name
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in)   :: ncid      ! netcdf file id
    integer          , intent(in)   :: varid     ! variable id
    character(len=*) , intent(out)  :: vname     ! variable vname
    type(Var_desc_t) , intent(inout):: vardesc   ! variable descriptor
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    vname = ''
    status = PIO_inq_varname(ncid,vardesc,vname)

  end subroutine ncd_inqvname

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_inqvdids
!
! !INTERFACE:
  subroutine ncd_inqvdids(ncid,dids,vardesc)
!
! !DESCRIPTION:
! inquire variable dimension ids
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(in)  :: ncid      ! netcdf file id
    integer         ,intent(out)  :: dids(:)   ! variable dids
    type(Var_desc_t),intent(inout):: vardesc   ! variable descriptor
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    dids = -1
    status = PIO_inq_vardimid(ncid,vardesc,dids)

  end subroutine ncd_inqvdids

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_int
!
! !INTERFACE:
  subroutine ncd_putatt_int(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! put integer attributes
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    integer          ,intent(in)    :: value     ! netcdf attrib value
    integer,optional ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_int

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_char
!
! !INTERFACE:
  subroutine ncd_putatt_char(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! put character attributes
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    character(len=*) ,intent(in)    :: value     ! netcdf attrib value
    integer,optional ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = PIO_put_att(ncid,varid,trim(attrib),value)

  end subroutine ncd_putatt_char

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_putatt_real
!
! !INTERFACE:
  subroutine ncd_putatt_real(ncid,varid,attrib,value,xtype)
!
! !DESCRIPTION:
! put real attributes
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid      ! netcdf file id
    integer          ,intent(in)    :: varid     ! netcdf var id
    character(len=*) ,intent(in)    :: attrib    ! netcdf attrib
    real(r8)         ,intent(in)    :: value     ! netcdf attrib value
    integer          ,intent(in)    :: xtype     ! netcdf data type
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
    real*4  :: value4
!-----------------------------------------------------------------------

    value4 = value

    if (xtype == pio_double) then
       status = PIO_put_att(ncid,varid,trim(attrib),value)
    else
       status = PIO_put_att(ncid,varid,trim(attrib),value4)
    endif

  end subroutine ncd_putatt_real

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defdim
!
! !INTERFACE:
  subroutine ncd_defdim(ncid,attrib,value,dimid)
!
! !DESCRIPTION:
! define dimension
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(in) :: ncid      ! netcdf file id
    character(len=*) , intent(in) :: attrib    ! netcdf attrib
    integer          , intent(in) :: value     ! netcdf attrib value
    integer          , intent(out):: dimid     ! netcdf dimension id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status
!-----------------------------------------------------------------------

    status = pio_def_dim(ncid,attrib,value,dimid)

  end subroutine ncd_defdim

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar_bynf
!
! !INTERFACE:
  subroutine ncd_defvar_bynf(ncid, varname, xtype, ndims, dimid, varid, &
                             long_name, units, cell_method, missing_value, fill_value, &
                             imissing_value, ifill_value, comment, flag_meanings, &
                             flag_values, nvalid_range )
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                  ! netcdf file id
    character(len=*) , intent(in)  :: varname                 ! variable name
    integer          , intent(in)  :: xtype                   ! external type
    integer          , intent(in)  :: ndims                   ! number of dims
    integer          , intent(inout) :: varid                 ! returned var id
    integer          , intent(in), optional :: dimid(:)       ! dimids
    character(len=*) , intent(in), optional :: long_name      ! attribute
    character(len=*) , intent(in), optional :: units          ! attribute
    character(len=*) , intent(in), optional :: cell_method    ! attribute
    character(len=*) , intent(in), optional :: comment        ! attribute
    character(len=*) , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)         , intent(in), optional :: missing_value  ! attribute for real
    real(r8)         , intent(in), optional :: fill_value     ! attribute for real
    integer          , intent(in), optional :: imissing_value ! attribute for int
    integer          , intent(in), optional :: ifill_value    ! attribute for int
    integer          , intent(in), optional :: flag_values(:)  ! attribute for int
    integer          , intent(in), optional :: nvalid_range(2)  ! attribute for int
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n                   ! indices
    integer :: ldimid(4)           ! local dimid
    integer :: dimid0(1)           ! local dimid
    integer :: status              ! error status 
    integer :: lxtype              ! local external type (in case logical variable)
    type(var_desc_t)   :: vardesc  ! local vardesc
    character(len=128) :: dimname  ! temporary
    character(len=256) :: str      ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bynf' ! subroutine name
!-----------------------------------------------------------------------

    varid = -1

    dimid0 = 0
    ldimid = 0
    if (present(dimid)) then
       ldimid(1:ndims) = dimid(1:ndims)
    else   ! ndims must be zero if dimid not present
       if (ndims /= 0) then
          write(iulog,*) subname//' ERROR: dimid not supplied and ndims ne 0 ',trim(varname),ndims
          call endrun()
       endif
    endif
       
    if ( xtype == ncd_log )then
       lxtype = ncd_int
    else
       lxtype = xtype
    end if
    if (masterproc .and. debug > 1) then
       write(iulog,*) 'Error in defining variable = ', trim(varname)
       write(iulog,*) subname//' ',trim(varname),lxtype,ndims,ldimid(1:ndims)
    endif
    
    if (ndims >  0) then 
       status = pio_inq_dimname(ncid,ldimid(ndims),dimname)
    end if

    ! Define variable
    if (present(dimid)) then
       status = PIO_def_var(ncid,trim(varname),lxtype,dimid(1:ndims),vardesc)
    else
       status = PIO_def_var(ncid,trim(varname),lxtype,dimid0        ,vardesc)
    endif
    varid = vardesc%varid

    !
    ! Add attributes
    !
    if (present(long_name)) then
       call ncd_putatt(ncid, varid, 'long_name', trim(long_name))
    end if
    if (present(flag_values)) then
       status = PIO_put_att(ncid,varid,'flag_values',flag_values)
       if ( .not. present(flag_meanings)) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call endrun( subname//" ERROR:: flag_values set -- but not flag_meanings" )
       end if
    end if
    if (present(flag_meanings)) then
       if ( .not. present(flag_values)) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call endrun( subname//" ERROR:: flag_meanings set -- but not flag_values" )
       end if
       if ( size(flag_values) /= size(flag_meanings) ) then
          write(iulog,*) 'Error in defining variable = ', trim(varname)
          call endrun( subname//" ERROR:: flag_meanings and flag_values dimension different")
       end if
       str = flag_meanings(1)
       do n = 1, size(flag_meanings)
          if ( index(flag_meanings(n), ' ') /= 0 )then
             write(iulog,*) 'Error in defining variable = ', trim(varname)
             call endrun( subname//" ERROR:: flag_meanings has an invalid space in it" )
          end if
          if ( n > 1 ) str = trim(str)//" "//flag_meanings(n)
       end do
       status = PIO_put_att(ncid,varid,'flag_meanings', trim(str) )
    end if
    if (present(comment)) then
       call ncd_putatt(ncid, varid, 'comment', trim(comment))
    end if
    if (present(units)) then
       call ncd_putatt(ncid, varid, 'units', trim(units))
    end if
    if (present(cell_method)) then
       str = 'time: ' // trim(cell_method)
       call ncd_putatt(ncid, varid, 'cell_methods', trim(str))
    end if
    if (present(fill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', fill_value, lxtype)
    end if
    if (present(missing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', missing_value, lxtype)
    end if
    if (present(ifill_value)) then
       call ncd_putatt(ncid, varid, '_FillValue', ifill_value, lxtype)
    end if
    if (present(imissing_value)) then
       call ncd_putatt(ncid, varid, 'missing_value', imissing_value, lxtype)
    end if
    if (present(nvalid_range)) then
       status = PIO_put_att(ncid,varid,'valid_range', nvalid_range )
    end if
    if ( xtype == ncd_log )then
       status = PIO_put_att(ncid,varid,'flag_values',     (/0, 1/) )
       status = PIO_put_att(ncid,varid,'flag_meanings',  "FALSE TRUE" )
       status = PIO_put_att(ncid,varid,'valid_range',    (/0, 1/) )
    end if

  end subroutine ncd_defvar_bynf

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_defvar_bygrid
!
! !INTERFACE:
  subroutine ncd_defvar_bygrid(ncid, varname, xtype, &
                               dim1name, dim2name, dim3name, dim4name, dim5name, &
                               long_name, units, cell_method, missing_value, fill_value, &
                               imissing_value, ifill_value, switchdim, comment, &
                               flag_meanings, flag_values, nvalid_range )
!
! !DESCRIPTION:
!  Define a netcdf variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                 ! netcdf file id
    character(len=*), intent(in)  :: varname                 ! variable name
    integer         , intent(in)  :: xtype                   ! external type
    character(len=*), intent(in), optional :: dim1name       ! dimension name
    character(len=*), intent(in), optional :: dim2name       ! dimension name
    character(len=*), intent(in), optional :: dim3name       ! dimension name
    character(len=*), intent(in), optional :: dim4name       ! dimension name
    character(len=*), intent(in), optional :: dim5name       ! dimension name
    character(len=*), intent(in), optional :: long_name      ! attribute
    character(len=*), intent(in), optional :: units          ! attribute
    character(len=*), intent(in), optional :: cell_method    ! attribute
    character(len=*), intent(in), optional :: comment        ! attribute
    character(len=*), intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)        , intent(in), optional :: missing_value  ! attribute for real
    real(r8)        , intent(in), optional :: fill_value     ! attribute for real
    integer         , intent(in), optional :: imissing_value ! attribute for int
    integer         , intent(in), optional :: ifill_value    ! attribute for int
    logical         , intent(in), optional :: switchdim      ! true=> permute dim1 and dim2 for output
    integer         , intent(in), optional :: flag_values(:)  ! attribute for int
    integer         , intent(in), optional :: nvalid_range(2)  ! attribute for int
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: n              ! indices
    integer :: ndims          ! dimension counter
    integer :: dimid(5)       ! dimension ids
    integer :: varid          ! variable id
    integer :: itmp           ! temporary
    character(len=256) :: str ! temporary
    character(len=*),parameter :: subname='ncd_defvar_bygrid' ! subroutine name
!-----------------------------------------------------------------------

    dimid(:) = 0

    ! Determine dimension ids for variable

    if (present(dim1name)) call ncd_inqdid(ncid, dim1name, dimid(1))
    if (present(dim2name)) call ncd_inqdid(ncid, dim2name, dimid(2))
    if (present(dim3name)) call ncd_inqdid(ncid, dim3name, dimid(3))
    if (present(dim4name)) call ncd_inqdid(ncid, dim4name, dimid(4))
    if (present(dim5name)) call ncd_inqdid(ncid, dim5name, dimid(5))

    ! Permute dim1 and dim2 if necessary

    if (present(switchdim)) then
       itmp = dimid(2)
       dimid(2) = dimid(1)
       dimid(1) = itmp
    end if

    ! Define variable

    ndims = 0
    if (present(dim1name)) then
       do n = 1, size(dimid)
          if (dimid(n) /= 0) ndims = ndims + 1
       end do
    end if

    call ncd_defvar_bynf(ncid,varname,xtype,ndims,dimid,varid, &
                         long_name=long_name, units=units, cell_method=cell_method, &
                         missing_value=missing_value, fill_value=fill_value, &
                         imissing_value=imissing_value, ifill_value=ifill_value, &
                         comment=comment, flag_meanings=flag_meanings, &
                         flag_values=flag_values, nvalid_range=nvalid_range )

  end subroutine ncd_defvar_bygrid

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_log_var0_nf
!
! !INTERFACE:
  subroutine ncd_io_log_var0_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global integer variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    logical          , intent(inout) :: data      ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(1), count(1)   ! output bounds
    integer :: status               ! error code
    integer :: idata                ! raw integer data
    logical :: varpresent           ! if true, variable is on tape
    integer :: temp(1)              ! temporary
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_log_var0_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: scalar var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
          status = pio_get_var(ncid, varid, idata)
          if (      idata == 0 )then
             data = .false.
          else if ( idata == 1 )then
             data = .true.
          else
             call endrun( subname// &
             ' ERROR: bad integer value for logical data' )
          end if
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = nt
          count(1) = 1
       else
          start(1) = 1
          count(1) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       if ( data )then
          temp(1) = 1
       else
          temp(1) = 0
       end if
       status = pio_put_var(ncid, varid, start, count, temp)

    endif   ! flag

  end subroutine ncd_io_log_var0_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var0_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var0_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global integer variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    integer          , intent(inout) :: data      ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(1), count(1)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    integer :: temp(1)              ! temporary
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_int_var0_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, varid, data)
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: scalar var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = nt
          count(1) = 1
       else
          start(1) = 1
          count(1) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       temp(1) = data
       status = pio_put_var(ncid, varid, start, count, temp)

    endif   ! flag

  end subroutine ncd_io_int_var0_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var0_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var0_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global real variable
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    real(r8)         , intent(inout) :: data      ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(1), count(1)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    real(r8):: temp(1)              ! temporary                
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_real_var0_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, vardesc, data)
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: scalar var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
       endif
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       if (present(nt))      then
          start(1) = nt
          count(1) = 1
       else
          start(1) = 1
          count(1) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       temp(1) = data
       status = pio_put_var(ncid, varid, start, count, temp)

    endif   ! flag

  end subroutine ncd_io_real_var0_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var1_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global integer array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    integer          , intent(inout) :: data(:)   ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(2), count(2)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_int_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, varid, data)
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: 1D var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          count(1) = size(data)
          start(2) = nt
          count(2) = 1
       else
          start(1) = 1
          count(1) = size(data)
          start(2) = 1
          count(2) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   ! flag

  end subroutine ncd_io_int_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_log_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_log_var1_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global integer array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid      ! netcdf file id
    character(len=*) , intent(in)    :: flag      ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname   ! variable name
    logical          , intent(inout) :: data(:)   ! raw data
    logical, optional, intent(out)   :: readvar   ! was var read?
    integer, optional, intent(in)    :: nt        ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(2), count(2)   ! output bounds
    integer :: status               ! error code
    integer, pointer :: idata(:)    ! Temporary integer data to send to file
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_log_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          allocate( idata(size(data)) ) 
          status = pio_get_var(ncid, varid, idata)
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: 1D var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
          data = (idata == 1)
          if ( any(idata /= 0 .and. idata /= 1) )then
             call endrun( subname// &
             ' ERROR: read in bad integer value(s) for logical data' )
          end if
          deallocate( idata )
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          count(1) = size(data)
          start(2) = nt
          count(2) = 1
       else
          start(1) = 1
          count(1) = size(data)
          start(2) = 1
          count(2) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       allocate( idata(size(data)) ) 
       where( data )
          idata = 1
       elsewhere
          idata = 0
       end where
       status = pio_put_var(ncid, varid, start, count, idata)
       deallocate( idata )

    endif   ! flag

  end subroutine ncd_io_log_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var1_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global real array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                ! netcdf file id
    character(len=*) , intent(in)    :: flag                ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname             ! variable name
    real(r8)         , intent(inout) :: data(:)             ! raw data
    logical          , optional, intent(out):: readvar      ! was var read?
    integer          , optional, intent(in) :: nt           ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(2), count(2)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_real_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, varid, data)
          if (single_column .and. present(posNOTonfile) ) then
             if ( .not. posNOTonfile )then
                call endrun( subname// &
                ' ERROR: 1D var is NOT compatable with posNOTonfile = .false.' )
             end if
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          start(2) = nt
          count(1) = size(data)
          count(2) = 1
       else
          start(1) = 1
          start(2) = 1
          count(1) = size(data)
          count(2) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   ! flag

  end subroutine ncd_io_real_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_char_var1_nf
!
! !INTERFACE:
  subroutine ncd_io_char_var1_nf(varname, data, flag, ncid, readvar, nt )
!
! !DESCRIPTION:
! netcdf I/O of global char array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    character(len=*) , intent(inout) :: data             ! raw data
    logical          , optional, intent(out):: readvar   ! was var read?
    integer          , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                   ! netCDF variable id
    integer :: m                       ! indices
    integer :: start(2), count(2)      ! output bounds
    integer :: status                  ! error code
    logical :: varpresent              ! if true, variable is on tape
    character(len=32) :: vname         ! variable error checking
    character(len=1)  :: tmpString(128)! temp for manipulating output string
    type(var_desc_t)  :: vardesc       ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_char_var1_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          status = pio_get_var(ncid, varid, data)
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid  (ncid, varname, varid, vardesc)

       if (present(nt))      then
          count(1) = len(data)
          count(2) = 1
          do m = 1,count(1)
             tmpString(m:m) = data(m:m)
          end do
          if ( count(1) > size(tmpString) )then
             write(iulog,*) subname//' ERROR: input string size is too large:'//trim(data)
          end if
          start(1) = 1
          start(2) = nt
          status = pio_put_var(ncid, varid, start=start, count=count, &
                               ival=tmpString(1:count(1)) )
       else
          status = pio_put_var(ncid, varid, data )
       end if

    endif   ! flag

  end subroutine ncd_io_char_var1_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var2_nf
!
! !INTERFACE:
  subroutine ncd_io_int_var2_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global integer 2D array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid                 ! netcdf file id
    character(len=*) , intent(in)    :: flag                 ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname              ! variable name
    integer          , intent(inout) :: data(:,:)            ! raw data
    logical          , optional, intent(out):: readvar       ! was var read?
    integer          , optional, intent(in) :: nt            ! time sample index
    logical          , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    logical :: found                ! if true, found lat/lon dims on file
    character(len=*),parameter :: subname='ncd_io_int_var2_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',vardesc,&
                  start,count,found=found,posNOTonfile=posNOTonfile)
             if ( found )then
                status = pio_get_var(ncid, varid, start, count, data)
             else
                status = pio_get_var(ncid, varid, data)
             end if
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          start(2) = 1
          start(3) = nt
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       else
          start(1) = 1
          start(2) = 1
          start(3) = 1
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       end if
       call ncd_inqvid(ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_int_var2_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var2_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var2_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global real 2D  array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: flag                ! 'read' or 'write'
    character(len=*), intent(in)    :: varname             ! variable name
    real(r8)        , intent(inout) :: data(:,:)           ! raw data
    logical         , optional, intent(out):: readvar      ! was var read?
    integer         , optional, intent(in) :: nt           ! time sample index
    logical         , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    logical :: found                ! if true, found lat/lon dims on file
    character(len=*),parameter :: subname='ncd_io_real_var2_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',vardesc,start,count,&
                                     found=found,posNOTonfile=posNOTonfile)
             if ( found )then
                status = pio_get_var(ncid, varid, start, count, data)
             else
                status = pio_get_var(ncid, varid, data)
             end if
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          start(2) = 1
          start(3) = nt
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       else
          start(1) = 1
          start(2) = 1
          start(3) = 1
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_real_var2_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_char_var2_nf
!
! !INTERFACE:
  subroutine ncd_io_char_var2_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global character array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: flag                ! 'read' or 'write'
    character(len=*), intent(in)    :: varname             ! variable name
    character(len=*), intent(inout) :: data(:)             ! raw data
    logical         , optional, intent(out):: readvar      ! was var read?
    integer         , optional, intent(in) :: nt           ! time sample index
    logical         , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(3), count(3)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    logical :: found                ! if true, found lat/lon dims on file
    character(len=*),parameter :: subname='ncd_io_char_var2_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          data   = ' '
          status = pio_get_var(ncid, varid, data)
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid  (ncid, varname, varid, vardesc)
       if (present(nt))      then
          start(1) = 1
          start(2) = 1
          start(3) = nt
          count(1) = len(data)
          count(2) = size(data)
          count(3) = 1
          status = pio_put_var(ncid, varid, start, count, data)
       else
          status = pio_put_var(ncid, varid, data)
       end if

    endif   

  end subroutine ncd_io_char_var2_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_char_var3_nf
!
! !INTERFACE:
  subroutine ncd_io_char_var3_nf(varname, data, flag, ncid, readvar, nt )
!
! !DESCRIPTION:
! netcdf I/O of global character array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    character(len=*), intent(in)    :: varname          ! variable name
    character(len=*), intent(inout) :: data(:,:)        ! raw data
    logical         , optional, intent(out):: readvar   ! was var read?
    integer         , optional, intent(in) :: nt        ! time sample index
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(4), count(4)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    character(len=*),parameter :: subname='ncd_io_char_var3_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          data   = ' '
          status = pio_get_var(ncid, varid, data)
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid  (ncid, varname, varid, vardesc)
       if (present(nt))      then
          start(:) = 1
          start(4) = nt
          count(1) = len(data)
          count(2) = size(data,dim=1)
          count(3) = size(data,dim=2)
          count(4) = 1
          status = pio_put_var(ncid, varid, start, count, data)
       else
          status = pio_put_var(ncid, varid, data)
       end if

    endif   

  end subroutine ncd_io_char_var3_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_char_varn_strt_nf
!
! !INTERFACE:
  subroutine ncd_io_char_varn_strt_nf(vardesc, data, flag, ncid, &
                                      start )
!
! !DESCRIPTION:
! netcdf I/O of global character array with start indices input
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid             ! netcdf file id
    character(len=*), intent(in)    :: flag             ! 'read' or 'write'
    type(var_desc_t), intent(in)    :: vardesc          ! local vardesc pointer
    character(len=*), intent(inout) :: data             ! raw data for this index
    integer         , intent(in)    :: start(:)         ! output bounds
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status               ! error code
    character(len=*),parameter :: subname='ncd_io_char_varn_strt_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       status = pio_get_var(ncid, vardesc, start, data )

    elseif (flag == 'write') then

       status = pio_put_var(ncid, vardesc, start, data )

    endif   

  end subroutine ncd_io_char_varn_strt_nf

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var1
!
! !INTERFACE:
  subroutine ncd_io_int_var1(varname, data, dim1name, &
                             flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! I/O for 1d integer field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    integer          , pointer       :: data(:)          ! local decomposition data
    character(len=*) , intent(in)    :: dim1name         ! dimension name
    integer          , optional, intent(in) :: nt        ! time sample index
    logical          , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary
    integer           :: n          ! index      
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    logical           :: varpresent ! if true, variable is on tape
    integer           :: xtype      ! netcdf data type
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var1' ! subroutine name
!-----------------------------------------------------------------------

    clmlevel = dim1name

    if (masterproc .and. debug > 1) then
       write(iulog,*) subname//' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid,clmlevel,vardesc,start,count)
             if (trim(clmlevel) == grlnd) then
                if (present(nt)) then
                   start(3) = nt
                   count(3) = 1
                end if
             else
                if (present(nt)) then
                   start(2) = nt
                   count(2) = 1
                end if
             end if
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_dimname(ncid,dids(ndims),dimname)
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                  PIO_INT, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            PIO_INT, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=0)

    else

       if (masterproc) then
          write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_int_var1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_log_var1
!
! !INTERFACE:
  subroutine ncd_io_log_var1(varname, data, dim1name, &
                             flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! I/O for 1d integer field
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid             ! netcdf file id
    character(len=*) , intent(in)    :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)    :: varname          ! variable name
    logical          , pointer       :: data(:)          ! local decomposition data
    character(len=*) , intent(in)    :: dim1name         ! dimension name
    integer          , optional, intent(in) :: nt        ! time sample index
    logical          , optional, intent(out):: readvar   ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
!EOP
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary
    integer           :: n          ! index      
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    integer, pointer :: idata(:)    ! Temporary integer data to send to file
    logical           :: varpresent ! if true, variable is on tape
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_log_var1' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc .and. debug > 1) then
       write(iulog,*) subname//' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    clmlevel = dim1name

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          allocate( idata(size(data)) ) 
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid,clmlevel,vardesc,start,count)
             if (trim(clmlevel) == grlnd) then
                if (present(nt)) then
                   start(3) = nt
                   count(3) = 1
                end if
             else
                if (present(nt)) then
                   start(2) = nt
                   count(2) = 1
                end if
             end if
             status = pio_get_var(ncid, varid, start, count, idata)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_dimname(ncid,dids(ndims),dimname)
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                  PIO_INT, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, idata, status)
          end if
          data = (idata == 1)
          if ( any(idata /= 0 .and. idata /= 1) )then
             call endrun( subname// &
             ' ERROR: read in bad integer value(s) for logical data' )
          end if
          deallocate( idata )
       end if
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            PIO_INT, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       allocate( idata(size(data)) ) 
       where( data )
          idata = 1
       elsewhere
          idata = 0
       end where
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, idata, status, fillval=0)
       deallocate( idata )

    else

       if (masterproc) then
          write(iulog,*) subname//' ERROR: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_log_var1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var1
!
! !INTERFACE:
  subroutine ncd_io_real_var1(varname, data, dim1name, &
                              flag, ncid, nt, readvar, cnvrtnan2fill)
!
! !DESCRIPTION:
! I/O for 1d real field
!
! !USES:
    use shr_infnan_mod  , only : shr_infnan_isnan
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                 ! netcdf file id
    character(len=*), intent(in)  :: flag                   ! 'read' or 'write'
    character(len=*), intent(in)  :: varname                ! variable name
    real(r8)        , pointer     :: data(:)                ! local decomposition data
    character(len=*), intent(in)  :: dim1name               ! dimension name
    integer         , optional, intent(in) :: nt            ! time sample index
    logical         , optional, intent(in) :: cnvrtnan2fill ! true => convert any NaN's to _FillValue (spval)
    logical         , optional, intent(out):: readvar       ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary
    integer           :: iodnum     ! iodesc num in list
    integer           :: varid      ! varid
    integer           :: ndims      ! ndims for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: start(3)   ! netcdf start index
    integer           :: count(3)   ! netcdf count index
    integer           :: status     ! error code  
    logical           :: varpresent ! if true, variable is on tape
    integer                , pointer  :: compDOF(:)
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var1' ! subroutine name
!-----------------------------------------------------------------------

    clmlevel = dim1name

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    endif

    if ( present(cnvrtnan2fill) )then
       if (.not. cnvrtnan2fill) then
          call endrun( subname//' ERROR: cnvrtnan2fill present but NOT set to true -- MUST set it to TRUE if used' )
       endif
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid,clmlevel,vardesc,start,count)
             if (trim(clmlevel) == grlnd) then
                if (present(nt)) then
                   start(3) = nt
                   count(3) = 1
                end if
             else
                if (present(nt)) then
                   start(2) = nt
                   count(2) = 1
                end if
             end if
             status = pio_get_var(ncid, varid, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_dimname(ncid,dids(ndims),dimname)
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                  PIO_DOUBLE, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
          if ( present(cnvrtnan2fill) )then
             do n = lbound(data,dim=1), ubound(data,dim=1)
                if ( data(n) == spval )then
                   data(n) = nan
                end if
             end do
          end if
       end if
       if (present(readvar)) readvar = varpresent
       
    elseif (flag == 'write') then
       
       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc, dids)
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            PIO_DOUBLE, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       if ( present(cnvrtnan2fill) )then
          do n = lbound(data,dim=1), ubound(data,dim=1)
             if ( shr_infnan_isnan( data(n) ) )then
                data(n) = spval
             end if
          end do
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)
       
    else
       
       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif
       
    endif

  end subroutine ncd_io_real_var1

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var2
!
! !INTERFACE:
  subroutine ncd_io_int_var2(varname, data, dim1name, &
                         lowerb2, upperb2, flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer          , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer           :: ndim1,ndim2       
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var2' ! subroutine name
!-----------------------------------------------------------------------

    clmlevel = dim1name      

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    if (flag == 'read') then
       
       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, vardesc, start, count)
             status = pio_get_var(ncid, vardesc, start, count, data)
             if (trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                if (present(nt)) then
                   start(4) = nt
                   count(4) = 1
                end if
             else
                count(2) = size(data,dim=2)
                if (present(nt)) then
                   start(3) = nt
                   count(3) = 1
                end if
             end if
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                     PIO_INT, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
               PIO_INT, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=0)

    else
       
       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif
       
    endif

  end subroutine ncd_io_int_var2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var2
!
! !INTERFACE:
  subroutine ncd_io_real_var2(varname, data, dim1name, &
                               lowerb2, upperb2, flag, ncid, nt, readvar, switchdim, cnvrtnan2fill)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
    use shr_infnan_mod  , only : shr_infnan_isnan
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)         , pointer     :: data(:,:)       ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    integer, optional, intent(in)  :: lowerb2,upperb2 ! lower and upper bounds of second dimension
    logical, optional, intent(in)  :: cnvrtnan2fill   ! true => convert any NaN's to _FillValue (spval)
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
    logical, optional, intent(in)  :: switchdim       ! true=> permute dim1 and dim2 for output
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), pointer :: temp(:,:)
    integer           :: ndim1,ndim2       
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: dimid      ! dimension id
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(4)   ! netcdf start index
    integer           :: count(4)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    integer           :: i,j      
    integer           :: lb1,lb2
    integer           :: ub1,ub2
    integer           :: itemp      ! integer temporary if switchdim is on
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var2' ! subroutine name
    logical           :: found = .false.
!-----------------------------------------------------------------------

    clmlevel = dim1name      

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    if ( present(cnvrtnan2fill) )then
       if (.not. cnvrtnan2fill) then
          call endrun( subname//' ERROR: cnvrtnan2fill present but NOT set to true -- MUST set it to TRUE if used' )
       endif
    end if

    lb1 = lbound(data, dim=1)
    ub1 = ubound(data, dim=1)
    lb2 = lbound(data, dim=2)
    ub2 = ubound(data, dim=2)
    if (present(switchdim)) then
       if (present(lowerb2)) lb2 = lowerb2
       if (present(upperb2)) ub2 = upperb2
       allocate(temp(lb2:ub2,lb1:ub1))
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, vardesc, start, count)
                if (present(nt)) then 
                status = pio_inq_varndims(ncid, vardesc, ndims)
                status = pio_inq_vardimid(ncid,vardesc, dids)
                do i=1,ndims
                   status = pio_inq_dimname(ncid, dids(i), dimname)
                   if ('time' == trim(dimname)) then
                      start(ndims) = nt
                      count(ndims) = 1
                      found=.true.
                end if
                end do
                if (.not.found) then 
                   write(iulog,*) trim(subname),' ERROR: time dimension not present for var', varname
                   call endrun()
                end if
             end if
             if (present(switchdim)) then
                status = pio_get_var(ncid, vardesc, start, count, temp)
                do j = lb2,ub2
                   do i = lb1,ub1
                      data(i,j) = temp(j,i) 
                   end do
                end do
             else
             status = pio_get_var(ncid, vardesc, start, count, data)
             endif
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             if (present(switchdim)) then
                call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                     PIO_DOUBLE, iodnum, switchdim=.true.)
             else
                call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                     PIO_DOUBLE, iodnum)
             end if
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             if (present(switchdim)) then
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status)
                do j = lb2,ub2
                do i = lb1,ub1
                   data(i,j) = temp(j,i) 
                end do
                end do
             else
                call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
             end if
          end if
          if ( present(cnvrtnan2fill) )then
             do j = lb2,ub2
             do i = lb1,ub1
                if ( data(i,j) == spval )then
                   data(i,j) = nan
                end if
             end do
             end do
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       if (present(switchdim)) then
          call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
               PIO_DOUBLE, iodnum, switchdim=.true.)
       else
          call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
               PIO_DOUBLE, iodnum)
       end if
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       if (present(switchdim)) then
          do j = lb2,ub2
          do i = lb1,ub1
             temp(j,i) = data(i,j)
          end do
          end do
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, temp, status, fillval=spval)
       else
          call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)
       end if
       if ( present(cnvrtnan2fill) )then
          do j = lb2,ub2
          do i = lb1,ub1
             if ( shr_infnan_isnan(data(i,j)) )then
                data(i,j) = spval
             end if
          end do
          end do
       end if

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

    if (present(switchdim)) then
       deallocate(temp)
    end if

  end subroutine ncd_io_real_var2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_int_var3
!
! !INTERFACE:
  subroutine ncd_io_int_var3(varname, data, dim1name, &
                               flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    integer          , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer           :: ndim1,ndim2       
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(5)   ! netcdf start index
    integer           :: count(5)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_int_var3' ! subroutine name
!-----------------------------------------------------------------------

    clmlevel = dim1name      

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, vardesc, start, count)
             if (trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                count(4) = size(data,dim=3)
                if (present(nt)) then
                   start(5) = nt
                   count(5) = 1
                end if
             else
                count(2) = size(data,dim=2)
                count(3) = size(data,dim=3)
                if (present(nt)) then
                   start(4) = nt
                   count(4) = 1
                end if
             end if
             status = pio_get_var(ncid, vardesc, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid,vardesc, dids)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                  PIO_INT, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            PIO_INT, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_int_var3

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var3
!
! !INTERFACE:
  subroutine ncd_io_real_var3(varname, data, dim1name, &
                               flag, ncid, nt, readvar)
!
! !DESCRIPTION:
! Netcdf i/o of 2d initial integer field out to netCDF file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid          ! netcdf file id
    character(len=*) , intent(in)  :: flag            ! 'read' or 'write'
    character(len=*) , intent(in)  :: varname         ! variable name
    real(r8)         , pointer     :: data(:,:,:)     ! local decomposition input data
    character(len=*) , intent(in)  :: dim1name        ! dimension 1 name
    integer, optional, intent(in)  :: nt              ! time sample index
    logical, optional, intent(out) :: readvar         ! true => variable is on initial dataset (read only)
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8), allocatable :: temp(:,:)
    integer               :: ndim1,ndim2       
    character(len=8)  :: clmlevel   ! clmlevel
    character(len=32) :: dimname    ! temporary      
    integer           :: status     ! error status
    integer           :: ndims      ! ndims total for var
    integer           :: ndims_iod  ! ndims iodesc for var
    integer           :: varid      ! varid
    integer           :: n          ! index      
    integer           :: dims(4)    ! dim sizes       
    integer           :: dids(4)    ! dim ids
    integer           :: iodnum     ! iodesc num in list
    integer           :: start(5)   ! netcdf start index
    integer           :: count(5)   ! netcdf count index
    logical           :: varpresent ! if true, variable is on tape
    type(iodesc_plus_type) , pointer  :: iodesc_plus
    type(var_desc_t)                  :: vardesc
    character(len=*),parameter :: subname='ncd_io_real_var3' ! subroutine name
!-----------------------------------------------------------------------

    clmlevel = dim1name      

    if (masterproc .and. debug > 1) then
       write(iulog,*) trim(subname),' ',trim(flag),' ',trim(varname),' ',trim(clmlevel)
    end if

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             start(:) = 1
             count(:) = 1
             call scam_field_offsets(ncid, clmlevel, vardesc, start, count)
             if (trim(clmlevel) == grlnd) then
                count(3) = size(data,dim=2)
                count(4) = size(data,dim=3)
                if (present(nt)) then 
                   start(5) = nt
                   count(5) = 1
                end if
             else
                count(2) = size(data,dim=2)
                count(3) = size(data,dim=3)
                if (present(nt)) then
                   start(4) = nt
                   count(4) = 1
                end if
             end if
             status = pio_get_var(ncid, vardesc, start, count, data)
          else
             status = pio_inq_varndims(ncid, vardesc, ndims)
             status = pio_inq_vardimid(ncid, vardesc, dids)
             status = pio_inq_dimname(ncid, dids(ndims), dimname)
             if (ndims == 0) then
                write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
                call endrun()
             end if
             if ('time' == trim(dimname)) then
                ndims_iod = ndims - 1
             else
                ndims_iod = ndims
             end if
             do n = 1,ndims_iod
                status = pio_inq_dimlen(ncid,dids(n),dims(n))
             enddo
             call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
                  PIO_DOUBLE, iodnum)
             iodesc_plus => iodesc_list(iodnum)
             if (present(nt)) then
                call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
             end if
             call pio_read_darray(ncid, vardesc, iodesc_plus%iodesc, data, status)
          end if
       end if
       if (present(readvar)) readvar = varpresent

    else if (flag == 'write') then

       call ncd_inqvid(ncid, varname ,varid, vardesc)
       status = pio_inq_varndims(ncid, vardesc, ndims)
       status = pio_inq_vardimid(ncid, vardesc , dids)
       if (ndims == 0) then
          write(iulog,*) trim(subname),' ERROR: ndims must be greater than 0'
          call endrun()
       end if
       status = pio_inq_dimname(ncid,dids(ndims),dimname)
       if ('time' == trim(dimname)) then
          ndims_iod = ndims - 1
       else
          ndims_iod = ndims
       end if
       do n = 1,ndims_iod
          status = pio_inq_dimlen(ncid,dids(n),dims(n))
       enddo
       call ncd_getiodesc(ncid, clmlevel, ndims_iod, dims(1:ndims_iod), dids(1:ndims_iod), &
            PIO_DOUBLE, iodnum)
       iodesc_plus => iodesc_list(iodnum)
       if (present(nt)) then
          call pio_setframe(vardesc, int(nt,kind=PIO_Offset))
       end if
       call pio_write_darray(ncid, vardesc, iodesc_plus%iodesc, data, status, fillval=spval)

    else

       if (masterproc) then
          write(iulog,*) subname,' error: unsupported flag ',trim(flag)
          call endrun()
       endif

    endif

  end subroutine ncd_io_real_var3


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ncd_io_real_var3_nf
!
! !INTERFACE:
  subroutine ncd_io_real_var3_nf(varname, data, flag, ncid, readvar, nt, posNOTonfile)
!
! !DESCRIPTION:
! netcdf I/O of global real 3D  array
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: flag                ! 'read' or 'write'
    character(len=*), intent(in)    :: varname             ! variable name
    real(r8)        , intent(inout) :: data(:,:,:)           ! raw data
    logical         , optional, intent(out):: readvar      ! was var read?
    integer         , optional, intent(in) :: nt           ! time sample index
    logical         , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: varid                ! netCDF variable id
    integer :: start(4), count(4)   ! output bounds
    integer :: status               ! error code
    logical :: varpresent           ! if true, variable is on tape
    character(len=32) :: vname      ! variable error checking
    type(var_desc_t)  :: vardesc    ! local vardesc pointer
    logical :: found                ! if true, found lat/lon dims on file
    character(len=*),parameter :: subname='ncd_io_real_var3_nf'
!-----------------------------------------------------------------------

    if (flag == 'read') then

       call ncd_inqvid(ncid, varname, varid, vardesc, readvar=varpresent)
       if (varpresent) then
          if (single_column) then
             call scam_field_offsets(ncid,'undefined',vardesc, start,count, &
                                     found=found,posNOTonfile=posNOTonfile)
             if ( found )then
                status = pio_get_var(ncid, varid, start, count, data)
             else
                status = pio_get_var(ncid, varid, data)
             end if
          else
             status = pio_get_var(ncid, varid, data)
          endif
       endif
       if (present(readvar)) readvar = varpresent

    elseif (flag == 'write') then

       if (present(nt))      then
          start(1) = 1
          start(2) = 1
          start(3) = nt
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       else
          start(1) = 1
          start(2) = 1
          start(3) = 1
          count(1) = size(data, dim=1)
          count(2) = size(data, dim=2)
          count(3) = 1
       end if
       call ncd_inqvid  (ncid, varname, varid, vardesc)
       status = pio_put_var(ncid, varid, start, count, data)

    endif   

  end subroutine ncd_io_real_var3_nf

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine scam_field_offsets
!
! !INTERFACE:
  subroutine scam_field_offsets( ncid, dim1name, vardesc, start, count, &
                                 found, posNOTonfile)
!
! !DESCRIPTION: 
! Read/Write initial data from/to netCDF instantaneous initial data file 
!
! !USES:
    use clm_varpar  , only : maxpatch
    use clm_varctl  , only : scmlon,scmlat,single_column
    use nanMod      , only : nan, bigint
    use shr_scam_mod, only : shr_scam_getCloseLatLon
    use shr_string_mod, only: shr_string_toLower
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t) , intent(inout) :: ncid         ! netcdf file id
    character(len=*)  , intent(in)    :: dim1name     ! dimension 1 name
    type(Var_desc_t)  , intent(inout) :: vardesc      ! variable descriptor
    integer           , intent(inout) :: start(:)     ! start index
    integer           , intent(inout) :: count(:)     ! count to retrieve
    logical, optional , intent(out)   :: found        ! if present return true if found
                                                      ! dimensions on file else false
                                                      ! if NOT present abort if can't find
    logical, optional , intent(in)    :: posNOTonfile ! Position is NOT on this file
!
! !CALLED FROM: subroutine inicfields
!
! !REVISION HISTORY:
! Created by John Truesdale

!
! !LOCAL VARIABLES:
!EOP
    integer :: data_offset                   ! offset into land array 1st column 
    integer :: ndata                         ! number of column (or pft points to read)
    real(r8) , pointer :: cols1dlon(:)       ! holds cols1d_ixy var
    real(r8) , pointer :: cols1dlat(:)       ! holds cols1d_jxy var
    real(r8) , pointer :: pfts1dlon(:)       ! holds pfts1d_ixy var
    real(r8) , pointer :: pfts1dlat(:)       ! holds pfts1d_jxy var
    real(r8) , pointer :: land1dlon(:)       ! holds land1d_ixy var
    real(r8) , pointer :: land1dlat(:)       ! holds land1d_jxy var
    integer :: cols(maxpatch)                ! grid cell columns for scam
    integer :: pfts(maxpatch)                ! grid cell pfts for scam
    integer :: landunits(maxpatch)           ! grid cell landunits for scam
    integer :: cc,i,ii                       ! index variable
    integer :: totpfts                       ! total number of pfts
    integer :: totcols                       ! total number of columns
    integer :: totlandunits                  ! total number of landunits
    integer, allocatable :: dids(:)                       ! dim ids
    integer :: varid                         ! netCDF variable id
    integer :: status                        ! return code
    integer :: latidx,lonidx                 ! latitude/longitude indices
    real(r8):: closelat,closelon             ! closest latitude and longitude indices
    integer :: ndims,dimlen                  ! number of dimensions in desired variable
    character(len=32) :: dimname             ! dimension name
    character(len=32) :: subname = 'scam_field_offsets'
!------------------------------------------------------------------------

    if ( present(posNOTonfile) )then
       if ( posNOTonfile )then
          if ( .not. present(found) )then
             call endrun( subname// &
                  'ERROR: Bad subroutine calling structure posNOTonfile sent, but found was NOT!')
          end if
          found = .false.
          return
       end if
    end if

    ! find closest land grid cell for this point

    if ( present(found) )then
       call shr_scam_getCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx,found)
       if ( .not. found ) return
    else
       call shr_scam_getCloseLatLon(ncid,scmlat,scmlon,closelat,closelon,latidx,lonidx)
    end if
     
    call ncd_inqvdims(ncid,ndims,vardesc)


    allocate(dids(ndims))
    status = pio_inq_vardimid(ncid, vardesc, dids)
    do i = 1,ndims
       status = pio_inq_dimname(ncid,dids(i),dimname)
       dimname=shr_string_toLower(dimname)
       status = pio_inq_dimlen(ncid,dids(i),dimlen)
       if ( trim(dimname)=='nj'.or. trim(dimname)=='lat'.or. trim(dimname)=='lsmlat') then
          start(i)=latidx
          count(i)=1
       else if ( trim(dimname)=='ni'.or. trim(dimname)=='lon'.or. trim(dimname)=='lsmlon') then
          start(i)=lonidx
          count(i)=1
       else if ( trim(dimname)=='column') then

          allocate (cols1dlon(dimlen))
          allocate (cols1dlat(dimlen))
       
       status = pio_inq_varid(ncid, 'cols1d_lon', varid)
       status = pio_get_var(ncid, varid, cols1dlon)
       status = pio_inq_varid(ncid, 'cols1d_lat', varid)
       status = pio_get_var(ncid, varid, cols1dlat)
       
       cols(:)     = bigint
       data_offset = bigint
          ii = 1
       ndata = 0
          do cc = 1, dimlen
          if (cols1dlon(cc) == closelon.and.cols1dlat(cc) == closelat) then
                cols(ii)=cc
                ndata  =ii
                ii=ii+1
          end if
       end do
       if (ndata == 0) then
          write(iulog,*)'couldnt find any columns for this latitude ',latidx,' and longitude ',lonidx
          call endrun( subname//'ERROR:: no columns for this position' )
       else
          data_offset=cols(1)
       end if
       
       deallocate (cols1dlon)
       deallocate (cols1dlat)
       
          start(i) = data_offset
          count(i) = ndata
       else if ( trim(dimname)=='pft') then
       
          allocate (pfts1dlon(dimlen))
          allocate (pfts1dlat(dimlen))

       status = pio_inq_varid(ncid, 'pfts1d_lon', varid)
       status = pio_get_var(ncid, varid, pfts1dlon)

       status = pio_inq_varid(ncid, 'pfts1d_lat', varid)
       status = pio_get_var(ncid, varid, pfts1dlat)
       
       pfts(:)     = bigint
       data_offset = bigint
          ii     = 1
       ndata = 0
          do cc = 1, dimlen
          if (pfts1dlon(cc) == closelon.and.pfts1dlat(cc) == closelat) then
                pfts(ii)=cc
                ndata  =ii
                ii=ii+1
          end if
       end do
       if (ndata == 0) then
          write(iulog,*)'couldnt find any pfts for this latitude ',closelat,' and longitude ',closelon
          call endrun( subname//'ERROR:: no PFTs for this position' )
       else
          data_offset=pfts(1)
       end if
       
       deallocate (pfts1dlon)
       deallocate (pfts1dlat)
       
          start(i) = data_offset
          count(i) = ndata
       else if ( trim(dimname)=='landunit') then
       
          allocate (land1dlon(dimlen))
          allocate (land1dlat(dimlen))

       status = pio_inq_varid(ncid, 'land1d_lon', varid)
       status = pio_get_var(ncid, varid, land1dlon)

       status = pio_inq_varid(ncid, 'land1d_lat', varid)
       status = pio_get_var(ncid, varid, land1dlat)
       
       landunits(:) = bigint
       data_offset  = bigint
          ii     = 1
       ndata = 0
          do cc = 1, dimlen
          if (land1dlon(cc) == closelon.and.land1dlat(cc) == closelat) then
                landunits(ii)=cc
                ndata  =ii
                ii=ii+1
          end if
       end do
       if (ndata == 0) then
          write(iulog,*)'couldnt find any landunits for this latitude ',closelat,' and longitude ',closelon
          call endrun( subname//'ERROR:: no landunits for this position' )
       else
          data_offset=landunits(1)
       end if
       
       deallocate (land1dlon)
       deallocate (land1dlat)
       
          start(i) = data_offset
          count(i) = ndata
          else
             start(i)=1
             count(i)=dimlen
          end if
       enddo
       deallocate(dids)

  end subroutine scam_field_offsets

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine ncd_getiodesc
!
! !INTERFACE:
  subroutine ncd_getiodesc(ncid, clmlevel, ndims, dims, dimids, &
       xtype, iodnum, switchdim) 
       
!
! !DESCRIPTION: 
! Returns an index to an io descriptor
!
! !USES:
!
! !ARGUMENTS:

    type(file_desc_t), intent(inout) :: ncid       ! PIO file descriptor
    character(len=8) , intent(in)    :: clmlevel   ! clmlevel
    integer          , intent(in)    :: ndims      ! ndims for var      
    integer          , intent(in)    :: dims(:)    ! dim sizes
    integer          , intent(in)    :: dimids(:)  ! dim ids
    integer          , intent(in)    :: xtype      ! file external type
    integer          , intent(out)   :: iodnum     ! iodesc num in list
    logical,optional , intent(in)    :: switchdim  ! switch level dimension and first dim 

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!
! !LOCAL VARIABLES:
!EOP
    integer :: k,m,n,cnt                     ! indices
    integer :: basetype                      ! pio basetype
    integer :: gsmap_lsize                   ! local size of gsmap
    integer :: gsmap_gsize                   ! global size of gsmap
    integer :: fullsize                      ! size of entire array on cdf
    integer :: gsize                         ! global size of clmlevel
    integer :: vsize                         ! other dimensions
    integer :: vsize1, vsize2                ! other dimensions
    integer :: status                        ! error status
    logical :: found                         ! true => found created iodescriptor
    integer :: ndims_file                    ! temporary
    character(len=64) dimname_file           ! dimension name on file
    character(len=64) dimname_iodesc         ! dimension name from io descriptor
    type(mct_gsMap),pointer       :: gsmap   ! global seg map
    integer, pointer,dimension(:) :: gsmOP   ! gsmap ordered points
    integer, pointer  :: compDOF(:)
    character(len=32) :: subname = 'ncd_getiodesc'
!------------------------------------------------------------------------

    ! Determining if need to create a new io descriptor

    n = 1
    found = .false.
    do while (n <= num_iodesc .and. .not.found)
       if (ndims == iodesc_list(n)%ndims .and. xtype == iodesc_list(n)%type) then
          found = .true.
          ! First found implies that dimension sizes are the same 
          do m = 1,ndims
             if (dims(m) /= iodesc_list(n)%dims(m)) then
                found = .false.
             endif
          enddo
          ! If found - then also check that dimension names are equal - 
          ! dimension ids in iodescriptor are only used to query dimension
          ! names associated with that iodescriptor
          if (found) then
             do m = 1,ndims
                status = PIO_inq_dimname(ncid,dimids(m),dimname_file)
                status = PIO_inquire(ncid, ndimensions=ndims_file)
                if (iodesc_list(n)%dimids(m) > ndims_file) then 
                   found = .false.
                   exit
                else
                   status = PIO_inq_dimname(ncid,iodesc_list(n)%dimids(m),dimname_iodesc)
                   if (trim(dimname_file) .ne. trim(dimname_iodesc)) then
                      found = .false.
                      exit
                   end if
                end if
             end do
          end if
          if (found) then
             iodnum = n
             if (iodnum > num_iodesc) then
                write(iulog,*) trim(subname),' ERROR: iodnum out of range ',iodnum,num_iodesc
                call endrun()
             endif
             RETURN
          endif
       endif
       n = n + 1
    enddo

    ! Creating a new io descriptor

    if (ndims > 0) then 
       num_iodesc = num_iodesc + 1
       if (num_iodesc > max_iodesc) then
          write(iulog,*) trim(subname),' ERROR num_iodesc gt max_iodesc ',max_iodesc
          call endrun()
       endif
       iodnum = num_iodesc
       if (masterproc .and. debug > 1) then
          write(iulog,*) trim(subname),' creating iodesc at iodnum,ndims,dims(1:ndims),xtype',&
               iodnum,ndims,dims(1:ndims),xtype
       endif
    end if

    if (xtype == pio_double ) then
       basetype = PIO_DOUBLE
    else if (xtype == pio_real) then
       basetype  = PIO_DOUBLE
    else if (xtype == pio_int) then
       basetype = PIO_INT
    end if

    call get_clmlevel_gsmap(clmlevel,gsmap)
    gsize = get_clmlevel_gsize(clmlevel)
    gsmap_lsize = mct_gsmap_lsize(gsmap,mpicom)
    gsmap_gsize = mct_gsmap_gsize(gsmap)

    call mct_gsmap_OP(gsmap,iam,gsmOP)

    fullsize = 1
    do n = 1,ndims
       fullsize = fullsize*dims(n)
    enddo

    vsize = fullsize / gsize
    if (mod(fullsize,gsize) /= 0) then
       write(iulog,*) subname,' ERROR in vsize ',fullsize,gsize,vsize
       call endrun()
    endif

    allocate(compDOF(gsmap_lsize*vsize))

    if (present(switchdim)) then
       if (switchdim) then
          cnt = 0
          do m = 1,gsmap_lsize
             do n = 1,vsize
                cnt = cnt + 1
                compDOF(cnt) = (gsmOP(m)-1)*vsize + n
             enddo
          enddo
       else
          write(iulog,*) subname,' ERROR switch dims present must have switchdim true' 
          call endrun()
       end if
    else         ! currently allow for up to two vertical dimensions
       if (vsize /= 1 .and. vsize /= dims(ndims)) then
          vsize1 = vsize/dims(ndims)
          vsize2 = dims(ndims)
          if (vsize1*vsize2 /= vsize) then
             write(iulog,*)'vsize1= ',vsize1,' vsize2= ',vsize2,' vsize= ',vsize
             call endrun('error in vsize1 and vsize2 computation')
          end if
          cnt = 0
          do k = 1,vsize2
             do n = 1,vsize1
                do m = 1,gsmap_lsize
                   cnt = cnt + 1
                   compDOF(cnt) = (k-1)*vsize1*gsmap_gsize + (n-1)*gsmap_gsize +  gsmOP(m) 
                enddo
             enddo
          end do
       else
          cnt = 0
          do n = 1,vsize
             do m = 1,gsmap_lsize
                cnt = cnt + 1
                compDOF(cnt) = (n-1)*gsmap_gsize + gsmOP(m)
             enddo
          enddo
       end if
    end if

    if (debug > 1) then
       do m = 0,npes-1
          if (iam == m) then
             write(iulog,*) trim(subname),' sizes1  = ',iam,gsize,gsmap_gsize,gsmap_lsize
             write(iulog,*) trim(subname),' sizes2  = ',iam,fullsize,npes,vsize
             write(iulog,*) trim(subname),' compDOF = ',iam,size(compDOF),minval(compDOF),maxval(compDOF)
             call shr_sys_flush(iulog)
          endif
          call mpi_barrier(mpicom,status)
       enddo
    endif

    deallocate(gsmOP)

    call pio_initdecomp(pio_subsystem, baseTYPE, dims(1:ndims), compDOF, iodesc_list(iodnum)%iodesc)

    deallocate(compDOF)

    iodesc_list(iodnum)%type  = xtype
    iodesc_list(iodnum)%ndims = ndims
    iodesc_list(iodnum)%dims  = 0
    iodesc_list(iodnum)%dims(1:ndims)   = dims(1:ndims)
    iodesc_list(iodnum)%dimids(1:ndims) = dimids(1:ndims)


  end subroutine ncd_getiodesc

  !=======================================================================

end module ncdio_pio
