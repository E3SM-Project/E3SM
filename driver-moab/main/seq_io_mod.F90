! !MODULE: seq_io_mod -- reads and writes driver files
!
! !DESCRIPTION:
!    Writes attribute vectors to netcdf
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2007-Oct-26 - T. Craig first version
!    2007-Dec-06 - T. Craig update and improve
!    2008-Feb-16 - J. Edwards convert to PIO
!    2010-Nov    - J. Edwards move PIO init and namelists from components to driver
! Current Problems
!  - the original use of seq_io will now ONLY work with the cpl because
!    of hardwiring cpl_io_type and cpl_io_iosystem.  want the original
!    io capabilities to be usable by any component
!  - the init1 method depends on seq_comm for name consistency but seq_comm_init
!    wants to be called after init1 so the global_comm can be modified for
!    async IO.  this needs to be reconciled.
!  - this routine stores information for all components but most methods are
!    hardwired to work only for the coupler.  should all the components info
!    be stored here or should this be more a general set of methods that are
!    reusable as it's original intent.
!
! !INTERFACE: ------------------------------------------------------------------

module seq_io_mod

  ! !USES:

  use shr_kind_mod, only: r4 => shr_kind_r4, r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_sys_mod,  only: shr_sys_abort
  use seq_comm_mct, only: logunit, CPLID, seq_comm_setptrs
  use seq_comm_mct, only: seq_comm_namelen, seq_comm_name
  use seq_comm_mct, only: mbaxid, atm_pg_active,mblxid,mb_scm_land
  use seq_flds_mod, only : seq_flds_lookup
  use mct_mod           ! mct wrappers
  use pio
  use component_type_mod
  use seq_infodata_mod, only: seq_infodata_type

  implicit none
  private

  ! !PUBLIC TYPES:

  ! none

  ! !PUBLIC MEMBER FUNCTIONS:

  public seq_io_wopen
  public seq_io_close
  public seq_io_redef
  public seq_io_enddef
  public seq_io_date2yyyymmdd
  public seq_io_sec2hms
  public seq_io_read
  public seq_io_write
  public seq_io_cpl_init
  ! !PUBLIC DATA MEMBERS


  ! none

  !EOP

  interface seq_io_read
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
     module procedure seq_io_read_moab_tags
  end interface seq_io_read
  interface seq_io_write
     module procedure seq_io_write_int
     module procedure seq_io_write_int1d
     module procedure seq_io_write_r8
     module procedure seq_io_write_r81d
     module procedure seq_io_write_char
     module procedure seq_io_write_time
     module procedure seq_io_write_moab_tags
  end interface seq_io_write

  !-------------------------------------------------------------------------------
  ! Local data
  !-------------------------------------------------------------------------------

  character(*),parameter :: prefix = "seq_io_"
  real(r8)    ,parameter :: fillvalue = SHR_CONST_SPVAL
  character(*),parameter :: modName = "(seq_io_mod) "
  integer(in) ,parameter :: debug = 1 ! internal debug level
  character(*),parameter :: version ='cpl7v10'
  character(*),parameter :: version0='cpl7v00'
  integer(in), parameter :: file_desc_t_cnt = 20 ! Note - this is hard-wired for now

  character(CL)                  :: wfilename = ''
  type(file_desc_t), save        :: cpl_io_file(0:file_desc_t_cnt)
  integer(IN)                    :: cpl_pio_iotype
  integer(IN)                    :: cpl_pio_ioformat
  type(iosystem_desc_t), pointer :: cpl_io_subsystem

  character(CL) :: charvar   ! buffer for string read/write

  !=================================================================================
contains
  !=================================================================================

  subroutine seq_io_cpl_init()
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat

    cpl_io_subsystem=>shr_pio_getiosys(CPLID)
    cpl_pio_iotype = shr_pio_getiotype(CPLID)
    cpl_pio_ioformat = shr_pio_getioformat(CPLID)

  end subroutine seq_io_cpl_init

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_wopen - open netcdf file
  !
  ! !DESCRIPTION:
  !    open netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_wopen(filename,clobber,file_ind, model_doi_url, set_fill, bfbflag)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(*),intent(in) :: filename
    logical,optional,intent(in):: clobber
    integer,optional,intent(in):: file_ind
    character(CL), optional, intent(in)  :: model_doi_url
    logical, optional, intent(in) :: set_fill
    logical, optional, intent(in) :: bfbflag !for priting bfbflag value in the history files
    !EOP
    integer :: lset_fill = PIO_NOFILL, old_set_fill
    logical :: exists
    logical :: lclobber
    integer :: iam,mpicom
    integer :: rcode
    integer :: nmode
    integer :: lfile_ind
    character(CL) :: lbfbflag
    character(CL) :: lversion
    character(CL) :: lmodel_doi_url
    character(*),parameter :: subName = '(seq_io_wopen) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lversion=trim(version0)
#ifdef PIO2
    if(present(set_fill)) then
       if(set_fill) lset_fill = PIO_FILL
    endif
#endif
    lclobber = .false.
    if (present(clobber)) lclobber=clobber

    lmodel_doi_url = 'unset'
    if (present(model_doi_url)) lmodel_doi_url = model_doi_url

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    lbfbflag = 'unset' ! default value for bfbflag
    if(present(bfbflag)) then
       if(bfbflag) lbfbflag = 'TRUE'
       if(.not. bfbflag) lbfbflag = 'FALSE'
    endif

    call seq_comm_setptrs(CPLID, iam=iam, mpicom=mpicom)

    if (.not. pio_file_is_open(cpl_io_file(lfile_ind))) then
       ! filename not open
       if (iam==0) inquire(file=trim(filename),exist=exists)
       call shr_mpi_bcast(exists,mpicom,'seq_io_wopen exists')
       if (exists) then
          if (lclobber) then
             nmode = pio_clobber

             ! only applies to classic NETCDF files.
             if(cpl_pio_iotype == PIO_IOTYPE_NETCDF .or. &
                  cpl_pio_iotype == PIO_IOTYPE_PNETCDF) then
                nmode = ior(nmode,cpl_pio_ioformat)
             endif

             rcode = pio_createfile(cpl_io_subsystem, cpl_io_file(lfile_ind), cpl_pio_iotype, trim(filename), nmode)
             if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
#ifdef PIO2
             rcode = pio_set_fill(cpl_io_file(lfile_ind), lset_fill, old_set_fill)
#endif
             rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"file_version",version)
             rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"model_doi_url",lmodel_doi_url)
             rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"BFBFLAG",trim(lbfbflag))
          else

             rcode = pio_openfile(cpl_io_subsystem, cpl_io_file(lfile_ind), cpl_pio_iotype, trim(filename), pio_write)
             if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
             call pio_seterrorhandling(cpl_io_file(lfile_ind),PIO_BCAST_ERROR)
             rcode = pio_get_att(cpl_io_file(lfile_ind),pio_global,"file_version",lversion)
             call pio_seterrorhandling(cpl_io_file(lfile_ind),PIO_INTERNAL_ERROR)
             if (trim(lversion) /= trim(version)) then
                rcode = pio_redef(cpl_io_file(lfile_ind))
                rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"file_version",version)
                rcode = pio_enddef(cpl_io_file(lfile_ind))
             endif

          endif
       else
          nmode = pio_noclobber
          ! only applies to classic NETCDF files.
          if(cpl_pio_iotype == PIO_IOTYPE_NETCDF .or. &
               cpl_pio_iotype == PIO_IOTYPE_PNETCDF) then
             nmode = ior(nmode,cpl_pio_ioformat)
          endif
          rcode = pio_createfile(cpl_io_subsystem, cpl_io_file(lfile_ind), cpl_pio_iotype, trim(filename), nmode)
          if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
          rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"file_version",version)
          rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"model_doi_url",lmodel_doi_url)
          rcode = pio_put_att(cpl_io_file(lfile_ind),pio_global,"BFBFLAG",trim(lbfbflag))
       endif
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename is open, better match open filename
       if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
       call shr_sys_abort(subname//'different file currently open '//trim(filename))
    else
       ! filename is already open, just return
    endif

  end subroutine seq_io_wopen

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_close - close netcdf file
  !
  ! !DESCRIPTION:
  !    close netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_close(filename,file_ind)

    use pio, only : pio_closefile

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(*),intent(in) :: filename
    integer,optional,intent(in):: file_ind

    !EOP

    integer :: iam
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_close) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (.not. pio_file_is_open(cpl_io_file(lfile_ind))) then
       ! filename not open, just return
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename matches, close it
       call pio_closefile(cpl_io_file(lfile_ind))
    else
       ! different filename is open, abort
       if(iam==0) write(logunit,*) subname,' different file currently open, aborting ',trim(filename)
       call shr_sys_abort(subname//'different file currently open, aborting '//trim(filename))
    endif

    wfilename = ''

  end subroutine seq_io_close

  !===============================================================================

  subroutine seq_io_redef(filename,file_ind)
    character(len=*), intent(in) :: filename

    integer,optional,intent(in):: file_ind
    integer :: lfile_ind
    integer :: rcode

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    rcode = pio_redef(cpl_io_file(lfile_ind))
  end subroutine seq_io_redef

  !===============================================================================

  subroutine seq_io_enddef(filename,file_ind)
    character(len=*), intent(in) :: filename
    integer,optional,intent(in):: file_ind
    integer :: lfile_ind
    integer :: rcode

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    rcode = pio_enddef(cpl_io_file(lfile_ind))
  end subroutine seq_io_enddef

  !===============================================================================

  character(len=24) function seq_io_date2yyyymmdd (date)

    use shr_cal_mod, only : shr_cal_datetod2string

    ! Input arguments

    integer, intent(in) :: date  ! date expressed as an integer: yyyymmdd

    !-------------------------------------------------------------------------------

    call shr_cal_datetod2string(date_str = seq_io_date2yyyymmdd, ymd = date)

  end function seq_io_date2yyyymmdd

  !===============================================================================

  character(len=8) function seq_io_sec2hms (seconds)

    ! Input arguments

    integer, intent(in) :: seconds

    ! Local workspace

    integer :: hours     ! hours of hh:mm:ss
    integer :: minutes   ! minutes of hh:mm:ss
    integer :: secs      ! seconds of hh:mm:ss

    !-------------------------------------------------------------------------------

    if (seconds < 0 .or. seconds > 86400) then
       write(logunit,*)'seq_io_sec2hms: bad input seconds:', seconds
       call shr_sys_abort('seq_io_sec2hms: bad input seconds')
    end if

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    if (minutes < 0 .or. minutes > 60) then
       write(logunit,*)'seq_io_sec2hms: bad minutes = ',minutes
       call shr_sys_abort('seq_io_sec2hms: bad minutes')
    end if

    if (secs < 0 .or. secs > 60) then
       write(logunit,*)'seq_io_sec2hms: bad secs = ',secs
       call shr_sys_abort('seq_io_sec2hms: bad secs')
    end if

    write(seq_io_sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)

  end function seq_io_sec2hms

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int - write scalar integer to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar integer to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_int(filename,idata,dname,whead,wdata,file_ind)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: lwhead, lwdata
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_int) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       !       rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(dname)//'_nx',1,dimid(1))
       !       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_INT,dimid,varid)
       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_INT,varid)
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(cpl_io_file(lfile_ind),varid,idata)

       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
    endif

  end subroutine seq_io_write_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_int1d - write 1d integer array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d integer array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_int1d(filename,idata,dname,whead,wdata,file_ind)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: lwhead, lwdata
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_int1d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(idata)
       rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_INT,dimid,varid)
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(cpl_io_file(lfile_ind),varid,idata)
    endif

    !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata

  end subroutine seq_io_write_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r8 - write scalar double to netcdf file
  !
  ! !DESCRIPTION:
  !    Write scalar double to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_r8(filename,rdata,dname,whead,wdata,file_ind)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: lwhead, lwdata
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_r8) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       !       rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(dname)//'_nx',1,dimid(1))
       !       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_DOUBLE,dimid,varid)


       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_DOUBLE,varid)
       if(rcode==PIO_NOERR) then
          rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
          if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
       end if
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(cpl_io_file(lfile_ind),varid,rdata)
    endif


  end subroutine seq_io_write_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_r81d - write 1d double array to netcdf file
  !
  ! !DESCRIPTION:
  !    Write 1d double array to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_r81d(filename,rdata,dname,whead,wdata,file_ind)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: lwhead, lwdata
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_r81d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(rdata)
       rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(cpl_io_file(lfile_ind),varid,rdata)

       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
    endif

  end subroutine seq_io_write_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_char - write char string to netcdf file
  !
  ! !DESCRIPTION:
  !    Write char string to netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_char(filename,rdata,dname,whead,wdata,file_ind)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    character(len=*),intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: lwhead, lwdata
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_char) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = len(charvar)
       rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(dname)//'_len',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file(lfile_ind),trim(dname),PIO_CHAR,dimid,varid)
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       charvar = ''
       charvar = trim(rdata)
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(cpl_io_file(lfile_ind),varid,charvar)
    endif

  end subroutine seq_io_write_char

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_time - write time variable to netcdf file
  !
  ! !DESCRIPTION:
  !    Write time variable to netcdf file
  !
  ! !REVISION HISTORY:
  !    2009-Feb-11 - M. Vertenstein - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_time(filename,time_units,time_cal,time_val,nt,whead,wdata,tbnds,file_ind)

    use shr_cal_mod, only : shr_cal_calMaxLen, shr_cal_calendarName, &
         shr_cal_noleap, shr_cal_gregorian

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename      ! file
    character(len=*),intent(in) :: time_units    ! units of time
    character(len=*),intent(in) :: time_cal      ! calendar type
    real(r8)        ,intent(in) :: time_val      ! data to be written
    integer(in),optional,intent(in) :: nt
    logical,optional,intent(in) :: whead         ! write header
    logical,optional,intent(in) :: wdata         ! write data
    real(r8),optional,intent(in) :: tbnds(2)     ! time bounds
    integer,optional,intent(in) :: file_ind

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    integer(in) :: dimid2(2)
    type(var_desc_t) :: varid
    logical :: lwhead, lwdata
    integer :: start(2),count(2)
    character(len=shr_cal_calMaxLen) :: lcalendar
    real(r8) :: time_val_1d(1)
    integer :: lfile_ind
    character(*),parameter :: subName = '(seq_io_write_time) '
    integer :: ndims
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       rcode = pio_def_dim(cpl_io_file(lfile_ind),'time',PIO_UNLIMITED,dimid(1))
       rcode = pio_def_var(cpl_io_file(lfile_ind),'time',PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,'units',trim(time_units))
       lcalendar = shr_cal_calendarName(time_cal,trap=.false.)
       if (trim(lcalendar) == trim(shr_cal_noleap)) then
          lcalendar = 'noleap'
       elseif (trim(lcalendar) == trim(shr_cal_gregorian)) then
          lcalendar = 'gregorian'
       endif
       rcode = pio_put_att(cpl_io_file(lfile_ind),varid,'calendar',trim(lcalendar))
       if (present(tbnds)) then
          rcode = pio_put_att(cpl_io_file(lfile_ind),varid,'bounds','time_bnds')
          dimid2(2)=dimid(1)
          rcode = pio_def_dim(cpl_io_file(lfile_ind),'ntb',2,dimid2(1))
          rcode = pio_def_var(cpl_io_file(lfile_ind),'time_bnds',PIO_DOUBLE,dimid2,varid)
       endif
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file(lfile_ind),'time',varid)
       if (present(nt)) then
          rcode = pio_put_var(cpl_io_file(lfile_ind),varid,(/nt/),time_val)
       else
          rcode = pio_put_var(cpl_io_file(lfile_ind),varid,time_val)
       endif
       if (present(tbnds)) then
          rcode = pio_inq_varid(cpl_io_file(lfile_ind),'time_bnds',varid)
          start = 1
          count = 0
          ndims = 1
          if (present(nt)) then
             start(2) = nt
             ndims = 2
          endif
          count(1) = 2
          count(2) = 1
          rcode = pio_put_var(cpl_io_file(lfile_ind),varid,start(1:ndims),count(1:ndims),tbnds)
       endif

       !      write(logunit,*) subname,' wrote time ',lwhead,lwdata
    endif

  end subroutine seq_io_write_time

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_moab_tags - write MOAB mesh tags to netcdf file
  !
  ! !DESCRIPTION:
  !    Writes one or more MOAB mesh tags (fields) from a MOAB mesh instance to a NetCDF file using PIO.
  !    The tags are written as variables, with support for writing multiple fields at once, and optional
  !    matrix input for direct data writing. Handles global cell ordering and reordering for correct output.
  !    Used for writing mesh-based data (e.g., from iMOAB) into the driver output files.
  !
  ! !ARGUMENTS:
  !    filename   [in]  - Name of the NetCDF file to write to
  !    mbxid      [in]  - iMOAB application ID (mesh handle)
  !    dname      [in]  - Prefix for variable names in the output file
  !    tag_list   [in]  - Colon-separated list of MOAB tag (field) names to write
  !    whead      [in, optional] - Logical flag to write NetCDF header/define variables
  !    wdata      [in, optional] - Logical flag to write data values
  !    matrix     [in, optional] - 2D array of data to write directly (overrides reading from MOAB)
  !    nx         [in, optional] - Number of global cells (overrides value from MOAB)
  !    file_ind   [in, optional] - File index for multi-file support
  !
  ! !NOTES:
  !    - Only cell-type entities are supported (ent_type=1). 
  !            - not anymore: spectral case, atm needs vertices ent_type=0
  !    - Handles reordering of local/global cell IDs for correct output.
  !    - If matrix is not present, data is read from MOAB tags; if present, matrix is written directly.
  !    - Skips writing the field "hgt" as a temporary exclusion.
  !    - Uses fillvalue for missing data.
  !
  ! !REVISION HISTORY:
  !    2025-07-20 - Cursor - initial documentation
  !    2025-09-24   allow vertex type for entity, for spectral case for atmosphere
  !
  ! !INTERFACE: ------------------------------------------------------------------
  subroutine seq_io_write_moab_tags(filename, mbxid, dname, tag_list, whead,wdata, matrix, nx, ny, nt, file_ind, dims2din, dims2do, mask )

    use shr_kind_mod,     only: CX => shr_kind_CX, CXX => shr_kind_CXX
    use iMOAB,            only: iMOAB_GetGlobalInfo, iMOAB_GetMeshInfo, iMOAB_GetDoubleTagStorage, &
        iMOAB_GetIntTagStorage
    use m_MergeSorts,     only: IndexSet, IndexSort

     ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename      ! file
    integer(in), intent(in)     :: mbxid         ! imoab app id, on coupler 
    character(len=*),intent(in) :: dname         ! name of data (prefix) 
    character(len=*),intent(in) :: tag_list      ! fields, separated by colon
    logical,optional,intent(in) :: whead         ! write header
    logical,optional,intent(in) :: wdata         ! write data
    real(r8), dimension(:,:), pointer, optional :: matrix  ! this may or may not be passed
    integer, optional,intent(in):: nx
    integer, optional,intent(in):: ny
    integer, optional,intent(in):: nt
    integer,optional,intent(in) :: file_ind
    integer,optional,intent(in) :: dims2din(2)   ! dim ids to output
    integer,optional,intent(out):: dims2do(2)    ! dim ids for output
    real(r8)         ,optional,intent(in) :: mask(:)

    logical :: lwhead, lwdata
    character(*),parameter :: subName = '(seq_io_write_moab_tags) '
    integer :: ndims, lfile_ind, iam, rcode
    integer(in)              :: ns, ng, lnx, lny, ix
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info

    type(var_desc_t) :: varid
    type(io_desc_t)  :: iodesc
    character(CL)    :: name1       ! var name
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    character(CL)  :: lpre
    type(mct_list) :: temp_list
    integer(kind=Pio_Offset_Kind) :: frame
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
    character(CXX) ::tagname, field
    integer(in),target  :: dimid2(2)
    integer(in),target  :: dimid3(3)
    integer(in),pointer :: dimid(:)
    integer(in)              :: ngv, ent_type, ierr
    real(r8)                 :: lfillvalue
    integer, allocatable         :: indx(:) !  this will be ordered
    integer, allocatable         :: Dof(:)  ! will be filled with global ids from cells
    integer, allocatable         :: Dof_reorder(:)  !
    real(r8), allocatable        :: data1(:), data_reorder(:)

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lwhead = .true.
    lwdata = .true.
    lfillvalue = fillvalue 
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata
    frame = -1
    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call seq_comm_setptrs(CPLID,iam=iam)

    call mct_list_init(temp_list ,trim(tag_list))
    size_list=mct_list_nitem (temp_list)  ! role of nf, number fields

    if (size_list < 1) then
       write(logunit,*) subname,' ERROR: size_list = ',size_list,trim(dname)
       call shr_sys_abort(subname//'size_list error')
    endif

    lpre = trim(dname)
    ! find out the number of global cells, needed for defining the variables length
    ierr = iMOAB_GetGlobalInfo( mbxid, ngv, ng)
    lny = 1

    ! get the local size ns
    ierr = iMOAB_GetMeshInfo ( mbxid, nvert, nvise, nbl, nsurf, nvisBC )
    if ((.not. atm_pg_active .and. (mbaxid .eq. mbxid)) .or. &
       (mb_scm_land .and. (mblxid .eq. mbxid))) then
      ent_type = 0
      ns = nvert(1) ! local vertices 
      lnx = ngv ! number of global vertices
    else
      ent_type = 1
      ns = nvise(1) ! local cells 
      lnx = ng
    endif
    ! it is needed to overwrite that for land, ng is too small
    !  ( for ne4pg2 it is 201 instead of 384)
    if (present(nx)) then
       if (nx /= 0) lnx = nx
    endif
    if (present(ny)) then
       if( ny /= 0) lny = ny
    endif
    if (present(nt)) then
       frame = nt
    endif
    if (lwhead) then
       if (present(dims2din)) then
          dimid2(1)=dims2din(1)
          dimid2(2)=dims2din(2)
       else
          rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(lpre)//'_nx',lnx,dimid2(1))
          rcode = pio_def_dim(cpl_io_file(lfile_ind),trim(lpre)//'_ny',lny,dimid2(2))
       endif
       if (present(dims2do)) then
          dims2do(1)=dimid2(1)
          dims2do(2)=dimid2(2)
       endif

       if (present(nt)) then
          dimid3(1:2) = dimid2
          rcode = pio_inq_dimid(cpl_io_file(lfile_ind),'time',dimid3(3))
          dimid => dimid3
       else
          dimid => dimid2
       endif

       do index_list = 1, size_list
          call mct_list_get(mctOStr,index_list,temp_list)
          field = mct_string_toChar(mctOStr)
          !-------tcraig, this is a temporary mod to NOT write hgt
          if (trim(field) /= "hgt") then
             name1 = trim(lpre)//'_'//trim(field)
             call seq_flds_lookup(field,longname=lname,stdname=sname,units=cunit)
            !  if (luse_float) then
            !     rcode = pio_def_var(cpl_io_file(lfile_ind),trim(name1),PIO_REAL,dimid1,varid)
            !     rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"_FillValue",real(lfillvalue,r4))
            !  else
             rcode = pio_def_var(cpl_io_file(lfile_ind),trim(name1),PIO_DOUBLE,dimid,varid)
             rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"_FillValue",lfillvalue)
             !end if
             rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"units",trim(cunit))
             rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"long_name",trim(lname))
             rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"standard_name",trim(sname))
             rcode = pio_put_att(cpl_io_file(lfile_ind),varid,"internal_dname",trim(dname))
             !-------tcraig
          endif
       enddo
       if (lwdata) call seq_io_enddef(filename, file_ind=lfile_ind)
    end if

    if (lwdata) then
       allocate(data1(ns))
       allocate(data_reorder(ns))
       allocate(dof(ns))
       allocate(dof_reorder(ns))
       allocate(indx(ns))

       ! note: size of dof is ns
       if (ns > 0) then
          tagname = 'GLOBAL_ID'//C_NULL_CHAR
          ierr = iMOAB_GetIntTagStorage ( mbxid, tagname, ns , ent_type, dof)
          if (ierr .ne. 0) then
            write(logunit,*) subname,' ERROR: cannot get dofs '
            call shr_sys_abort(subname//'cannot get dofs ')
          endif

          call IndexSet(ns, indx)
          call IndexSort(ns, indx, dof, descend=.false.)
          !      after sort, dof( indx(i)) < dof( indx(i+1) )
          do ix=1,ns
             dof_reorder(ix) = dof(indx(ix)) ! 
          enddo
          ! so we know that dof_reorder(ix) < dof_reorder(ix+1)
       endif
       call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof_reorder, iodesc)

       deallocate(dof)
       deallocate(dof_reorder) 
       do index_list = 1, size_list
          call mct_list_get(mctOStr,index_list,temp_list)
          field = mct_string_toChar(mctOStr)
          !-------tcraig, this is a temporary mod to NOT write hgt
          if (trim(field) /= "hgt") then
             name1 = trim(lpre)//'_'//trim(field)
             rcode = pio_inq_varid(cpl_io_file(lfile_ind),trim(name1),varid)
             call pio_setframe(cpl_io_file(lfile_ind),varid,frame)
             if (present(matrix)) then
               do ix = 1, ns
                 data1(ix) = matrix(ix, index_list) ! 
               enddo
             else
               tagname = trim(field)//C_NULL_CHAR
               if (ns > 0 ) then
                  ierr = iMOAB_GetDoubleTagStorage (mbxid, tagname, ns , ent_type, data1)
                  if (ierr .ne. 0) then
                     write(logunit,*) subname,' ERROR: cannot get tag data ', trim(tagname)
                     call shr_sys_abort(subname//'cannot get tag data ')
                  endif
               endif
             endif

             ! remove MOAB default values
             do ix = 1, ns
               if (data1(ix) < -9.99999E+9_r8) then
                  data1(ix) = 0.0_r8
               endif
             enddo


             ! rearrange data for writing and handle mask
             if(present(mask)) then
               do ix=1,ns
                 if(mask(indx(ix)) /= 0) then
                   data_reorder(ix) = data1(indx(ix))
                 else
                   data_reorder(ix) = lfillvalue
                 endif
               enddo
             else
               do ix=1,ns
                  data_reorder(ix) = data1(indx(ix))
               enddo
             endif
             
             call pio_write_darray(cpl_io_file(lfile_ind), varid, iodesc, data_reorder, rcode, fillval=lfillvalue)
          endif
       enddo

       call pio_freedecomp(cpl_io_file(lfile_ind), iodesc)
       deallocate(data1)
       deallocate(data_reorder)
       deallocate(indx)

    end if


  end subroutine seq_io_write_moab_tags

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int - read scalar integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar integer from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int(filename,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(inout):: idata  ! integer data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer :: i1d(1)
    character(*),parameter :: subName = '(seq_io_read_int) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_int1d(filename,i1d,dname)
    idata = i1d(1)

  end subroutine seq_io_read_int

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_int1d - read 1d integer from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d integer array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_int1d(filename,idata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(inout):: idata(:)  ! integer data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_t) :: pioid
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_int1d) '
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_int1d exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort(subname//'ERROR: file invalid '//trim(filename)//' '//trim(dname))
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)

    call pio_closefile(pioid)

    !      write(logunit,*) subname,' read int ',trim(dname)


  end subroutine seq_io_read_int1d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r8 - read scalar double from netcdf file
  !
  ! !DESCRIPTION:
  !    Read scalar double from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r8(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(inout):: rdata  ! real data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(seq_io_read_r8) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_io_read_r81d(filename,r1d,dname)
    rdata = r1d(1)

  end subroutine seq_io_read_r8

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_r81d - read 1d double array from netcdf file
  !
  ! !DESCRIPTION:
  !    Read 1d double array from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_r81d(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(inout):: rdata(:)  ! real data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_T) :: pioid
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_r81d) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_r81d exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort(subname//'ERROR: file invalid '//trim(filename)//' '//trim(dname))
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

    call pio_closefile(pioid)

    !      write(logunit,*) subname,' read int ',trim(dname)

  end subroutine seq_io_read_r81d

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_char - read char string from netcdf file
  !
  ! !DESCRIPTION:
  !    Read char string from netcdf file
  !
  ! !REVISION HISTORY:
  !    2010-July-06 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_char(filename,rdata,dname)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    character(len=*),intent(inout):: rdata  ! character data
    character(len=*),intent(in) :: dname    ! name of data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    type(file_desc_T) :: pioid
    type(var_desc_t) :: varid
    logical :: exists
    character(CL)  :: lversion
    character(CL)  :: name1
    character(*),parameter :: subName = '(seq_io_read_char) '

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    lversion=trim(version0)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_char exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       !         write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort(subname//'ERROR: file invalid '//trim(filename)//' '//trim(dname))
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,charvar)
    rdata = trim(charvar)

    call pio_closefile(pioid)

  end subroutine seq_io_read_char

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_moab_tags - read MOAB mesh tags from netcdf file
  !
  ! !DESCRIPTION:
  !    Reads one or more MOAB mesh tags (fields) from a NetCDF file using PIO and stores them
  !    into a MOAB mesh instance or a provided matrix. Supports reading multiple fields at once,
  !    and handles global cell ordering and reordering for correct mapping to the mesh. Used for
  !    restoring mesh-based data (e.g., from iMOAB) from driver output files.
  !
  ! !ARGUMENTS:
  !    filename   [in]  - Name of the NetCDF file to read from
  !    mbxid      [in]  - iMOAB application ID (mesh handle)
  !    dname      [in]  - Prefix for variable names in the input file
  !    tag_list   [in]  - Colon-separated list of MOAB tag (field) names to read
  !    matrix     [in, optional] - 2D array to store data directly (if present, data is written here instead of MOAB)
  !    nx         [in, optional] - Number of global cells (overrides value from MOAB)
  !
  ! !NOTES:
  !    - Only cell-type entities are supported (ent_type=1).
  !         - not anymore: spectral case, atm, uses vertices, ent_type = 0
  !    - Handles reordering of local/global cell IDs for correct mapping.
  !    - If matrix is present, data is stored in the matrix; otherwise, data is set as MOAB tags.
  !    - Skips reading the field "hgt" as a temporary exclusion.
  !    - Uses fillvalue for missing data.
  !
  ! !REVISION HISTORY:
  !    2025-07-20 - Cursor - initial documentation
  !    2025-09-24  spectral case atm
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_moab_tags(filename, mbxid, dname, tag_list, matrix, nx)

    use shr_kind_mod,     only: CX => shr_kind_CX, CXX => shr_kind_CXX
    use iMOAB,            only: iMOAB_GetGlobalInfo, iMOAB_GetMeshInfo, iMOAB_SetDoubleTagStorage, &
        iMOAB_GetIntTagStorage
    use m_MergeSorts,     only: IndexSet, IndexSort
     ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename      ! file
    integer(in), intent(in)     :: mbxid         ! imoab app id, on coupler 
    character(len=*),intent(in) :: dname         ! name of data (prefix) 
    character(len=*),intent(in) :: tag_list      ! fields, separated by colon
    real(r8), dimension(:,:), pointer, optional :: matrix  ! this may or may not be passed
    integer, optional,intent(in):: nx

    integer(in)              :: ns, ng, ix
    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3) ! for moab info

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    integer(in) :: k,n,n1,n2,ndims
    type(file_desc_t) :: pioid
    integer(in) :: dimid(4)
    type(var_desc_t) :: varid
    integer(in) :: lnx,lny,lni
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    logical :: exists
    type(io_desc_t) :: iodesc

    integer, allocatable         :: indx(:) !  this will be ordered
    integer, allocatable         :: Dof(:)  ! will be filled with global ids from cells
    integer, allocatable         :: Dof_reorder(:)  !
    real(r8), allocatable        :: data1(:), data_reorder(:)

    character(CL)  :: lversion
    character(CL)  :: name1
    character(CL)  :: lpre

    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
    character(CXX) ::tagname, field

    integer(in)              :: ngv, ent_type, ierr
    character(*),parameter :: subName = '(seq_io_read_moab_tags) '

    
    lpre = trim(dname)

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    call mct_list_init(temp_list ,trim(tag_list))
    size_list=mct_list_nitem (temp_list)  ! role of nf, number fields

    if (size_list < 1) then
       write(logunit,*) subname,' ERROR: size_list = ',size_list,trim(dname)
       call shr_sys_abort(subname//'size_list error')
    endif


    !call mct_gsmap_OrderedPoints(gsmap, iam, Dof)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_avs exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       if(iam==0) write(logunit,*) subname,' open file ',trim(filename),' for ',trim(dname)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort(subname//'ERROR: file invalid '//trim(filename)//' '//trim(dname))
    endif

        ! find out the number of global cells, needed for defining the variables length
    ierr = iMOAB_GetGlobalInfo( mbxid, ngv, ng)
    lny = 1 ! do we need 2 var, or just 1 
    ierr = iMOAB_GetMeshInfo ( mbxid, nvert, nvise, nbl, nsurf, nvisBC )
    if ((.not. atm_pg_active .and. (mbaxid .eq. mbxid)) .or. &
       (mb_scm_land .and. (mblxid .eq. mbxid))) then
      ent_type = 0
      ns = nvert(1) ! local vertices 
      lnx = ngv ! number of global vertices
    else
      ent_type = 1
      ns = nvise(1) ! local cells 
      lnx = ng
    endif
    ! it is needed to overwrite that for land, ng is too small
    !  ( for ne4pg2 it is 201 instead of 384)
    if (present(nx)) then
       lnx = nx 
    endif
    allocate(data1(ns))
    allocate(data_reorder(ns))
    allocate(dof(ns))
    allocate(dof_reorder(ns))

   ! note: size of dof is ns
    tagname = 'GLOBAL_ID'//C_NULL_CHAR
    if (ns > 0 ) then 
       ierr = iMOAB_GetIntTagStorage ( mbxid, tagname, ns , ent_type, dof)
       if (ierr .ne. 0) then
          write(logunit,*) subname,' ERROR: cannot get dofs '
          call shr_sys_abort(subname//'cannot get dofs ')
       endif
    endif
   allocate(indx(ns))
   call IndexSet(ns, indx)
   call IndexSort(ns, indx, dof, descend=.false.)
   !      after sort, dof( indx(i)) < dof( indx(i+1) )
   do ix=1,ns
      dof_reorder(ix) = dof(indx(ix)) ! 
   enddo
   deallocate(dof)

   do index_list = 1, size_list
       call mct_list_get(mctOStr,index_list,temp_list)
       field = mct_string_toChar(mctOStr)
       name1 = trim(lpre)//'_'//trim(field)
      
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid,trim(name1),varid)
       if (rcode == pio_noerr) then
          if (index_list==1) then
             rcode = pio_inq_varndims(pioid, varid, ndims)
             rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
             rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
             if (ndims>=2) then
                rcode = pio_inq_dimlen(pioid, dimid(2), lny)
             else
                lny = 1
             end if
!             if (lnx*lny /= ng) then
!                write(logunit,*) subname,' ERROR: dimensions do not match',&
!                     lnx,lny, ng
!                call shr_sys_abort(subname//'ERROR: dimensions do not match')
!             end if
             
             call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof_reorder, iodesc)

             deallocate(dof_reorder)
          end if

          call pio_read_darray(pioid,varid,iodesc, data1, rcode)
          do ix=1,ns
             data_reorder(indx(ix)) = data1(ix) ! or is it data_reorder(ix) = data1(indx(ix)) ? 
          enddo
          if (present(matrix)) then
            !matrix(:, index_list)  = data_reorder(:) ! 
            do ix = 1,ns
               matrix(ix, index_list)  = data_reorder(ix) !
            enddo
          else
            tagname = trim(field)//C_NULL_CHAR
            if (ns > 0) then
               ierr = iMOAB_SetDoubleTagStorage (mbxid, tagname, ns , ent_type, data_reorder)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' ERROR: cannot set tag data ', trim(tagname)
                  call shr_sys_abort(subname//'cannot set tag data ')
               endif
            endif
          endif
         !  n = 0
         !  do n1 = 1,ni
         !     do n2 = 1,ns
         !        n = n + 1
         !        avs(:n1)%rAttr(k,n2) = data(n)
         !     enddo
         !  enddo
       else
          write(logunit,*)subname, ' warning: field ',trim(field), ' name1:', trim(name1),  ' is not on restart file'
          write(logunit,*)'for backwards compatibility will set it to 0'
         !  do n1 = 1,ni
         !     avs(n1)%rattr(k,:) = 0.0_r8
         !  enddo
         data_reorder = 0.
         if (present(matrix)) then
            ! matrix(:, index_list)  = data_reorder(:) ! 
            do ix = 1,ns
               matrix(ix, index_list)  = data_reorder(ix) !
            enddo
         else
            tagname = trim(field)//C_NULL_CHAR
            if ( ns > 0 ) then
               ierr = iMOAB_SetDoubleTagStorage (mbxid, tagname, ns , ent_type, data_reorder)
               if (ierr .ne. 0) then
                  write(logunit,*) subname,' ERROR: cannot set tag data ', trim(tagname)
                  call shr_sys_abort(subname//'cannot set tag data ')
               endif
            endif
         endif

       end if
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    enddo

    deallocate(data1)
    deallocate(data_reorder)

   !  !--- zero out fill value, this is somewhat arbitrary
   !  do n1 = 1,ni
   !     do n2 = 1,ns
   !        do k = 1,nf
   !           if (AVS(n1)%rAttr(k,n2) == fillvalue) then
   !              AVS(n1)%rAttr(k,n2) = 0.0_r8
   !           endif
   !        enddo
   !     enddo
   !  enddo

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine seq_io_read_moab_tags

  !===============================================================================
  !===============================================================================
end module seq_io_mod
