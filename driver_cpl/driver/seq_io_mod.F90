!===============================================================================
! SVN $Id: seq_io_mod.F90 68253 2015-02-18 22:24:57Z mvertens $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_1_15/driver/seq_io_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
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

  use shr_kind_mod, only: r8 => shr_kind_r8, in => shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl, cs => shr_kind_cs
  use shr_sys_mod       ! system calls
  use seq_comm_mct
  use seq_flds_mod, only : seq_flds_lookup
  use mct_mod           ! mct wrappers
  use pio
  use component_type_mod

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
     module procedure seq_io_read_av
     module procedure seq_io_read_avs
     module procedure seq_io_read_avscomp
     module procedure seq_io_read_int
     module procedure seq_io_read_int1d
     module procedure seq_io_read_r8
     module procedure seq_io_read_r81d
     module procedure seq_io_read_char
  end interface
  interface seq_io_write
     module procedure seq_io_write_av
     module procedure seq_io_write_avs
     module procedure seq_io_write_avscomp
     module procedure seq_io_write_int
     module procedure seq_io_write_int1d
     module procedure seq_io_write_r8
     module procedure seq_io_write_r81d
     module procedure seq_io_write_char
     module procedure seq_io_write_time
  end interface

!-------------------------------------------------------------------------------
! Local data
!-------------------------------------------------------------------------------

   character(*),parameter :: prefix = "seq_io_"
   character(CL)          :: wfilename = ''
   real(r8)    ,parameter :: fillvalue = SHR_CONST_SPVAL

   character(*),parameter :: modName = "(seq_io_mod) "
   integer(in) ,parameter :: debug = 1 ! internal debug level


   type(file_desc_t), save :: cpl_io_file
   integer(IN)             :: cpl_pio_iotype
   type(iosystem_desc_t), pointer :: cpl_io_subsystem
   character(*),parameter  :: version ='cpl7v10'
   character(*),parameter  :: version0='cpl7v00'

   character(CL) :: charvar   ! buffer for string read/write
   integer(IN) :: io_comm

!=================================================================================
contains
!=================================================================================

  subroutine seq_io_cpl_init()
    use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype

    character(len=seq_comm_namelen) :: cpl_name

    ! For consistency with how io_compname is set, we get the coupler name here
    cpl_name = seq_comm_name(CPLID)
    cpl_io_subsystem=>shr_pio_getiosys(cpl_name)
    cpl_pio_iotype = shr_pio_getiotype(cpl_name)

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

subroutine seq_io_wopen(filename,clobber,cdf64)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(*),intent(in) :: filename
    logical,optional,intent(in):: clobber
    logical,optional,intent(in):: cdf64

    !EOP

    logical :: exists
    logical :: lclobber
    logical :: lcdf64
    integer :: iam,mpicom
    integer :: rcode
    integer :: nmode
    character(CL)  :: lversion
    character(*),parameter :: subName = '(seq_io_wopen) '
    
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    lversion=trim(version0)

    lclobber = .false.
    if (present(clobber)) lclobber=clobber

    lcdf64 = .false.
    if (present(cdf64)) lcdf64=cdf64

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)

    if (.not. pio_file_is_open(cpl_io_file)) then
       ! filename not open
       if (iam==0) inquire(file=trim(filename),exist=exists)
       call shr_mpi_bcast(exists,mpicom,'seq_io_wopen exists')
       if (exists) then
          if (lclobber) then
             nmode = pio_clobber
             if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
             rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), nmode)
             if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
             rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
          else

             rcode = pio_openfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), pio_write)
             if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
             call pio_seterrorhandling(cpl_io_file,PIO_BCAST_ERROR)
             rcode = pio_get_att(cpl_io_file,pio_global,"file_version",lversion)
             call pio_seterrorhandling(cpl_io_file,PIO_INTERNAL_ERROR)
             if (trim(lversion) /= trim(version)) then
                rcode = pio_redef(cpl_io_file)
                rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
                rcode = pio_enddef(cpl_io_file)
             endif
          endif
       else
          nmode = pio_noclobber
          if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
          rcode = pio_createfile(cpl_io_subsystem, cpl_io_file, cpl_pio_iotype, trim(filename), nmode)
          if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
          rcode = pio_put_att(cpl_io_file,pio_global,"file_version",version)
       endif
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename is open, better match open filename
       if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
       call shr_sys_abort()
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

subroutine seq_io_close(filename)

    use pio, only : pio_closefile

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(*),intent(in) :: filename

    !EOP

    integer :: iam
    integer :: rcode
    character(*),parameter :: subName = '(seq_io_close) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    call seq_comm_setptrs(CPLID,iam=iam)


    if (.not. pio_file_is_open(cpl_io_file)) then
       ! filename not open, just return
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename matches, close it
       call pio_closefile(cpl_io_file)
    else
       ! different filename is open, abort
       if(iam==0) write(logunit,*) subname,' different file currently open ',trim(filename)
       call shr_sys_abort()
    endif

    wfilename = ''

end subroutine seq_io_close

!===============================================================================

subroutine seq_io_redef(filename)
    character(len=*), intent(in) :: filename
    integer :: rcode

    rcode = pio_redef(cpl_io_file)
end subroutine seq_io_redef

!===============================================================================

subroutine seq_io_enddef(filename)
    character(len=*), intent(in) :: filename
    integer :: rcode

    rcode = pio_enddef(cpl_io_file)
end subroutine seq_io_enddef

!===============================================================================

character(len=10) function seq_io_date2yyyymmdd (date)

! Input arguments

   integer, intent(in) :: date

! Local workspace

   integer :: year    ! year of yyyy-mm-dd
   integer :: month   ! month of yyyy-mm-dd
   integer :: day     ! day of yyyy-mm-dd

!-------------------------------------------------------------------------------

   if (date < 0) then
      call shr_sys_abort ('seq_io_date2yyyymmdd: negative date not allowed')
   end if

   year  = date / 10000
   month = (date - year*10000) / 100
   day   = date - year*10000 - month*100

   write(seq_io_date2yyyymmdd,80) year, month, day
80 format(i4.4,'-',i2.2,'-',i2.2)

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
      call shr_sys_abort()
   end if

   hours   = seconds / 3600
   minutes = (seconds - hours*3600) / 60
   secs    = (seconds - hours*3600 - minutes*60)

   if (minutes < 0 .or. minutes > 60) then
      write(logunit,*)'seq_io_sec2hms: bad minutes = ',minutes
      call shr_sys_abort()
   end if

   if (secs < 0 .or. secs > 60) then
      write(logunit,*)'seq_io_sec2hms: bad secs = ',secs
      call shr_sys_abort()
   end if

   write(seq_io_sec2hms,80) hours, minutes, secs
80 format(i2.2,':',i2.2,':',i2.2)

end function seq_io_sec2hms

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_av - write AV to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_av(filename,gsmap,AV,dname,whead,wdata,nx,ny,nt,fillval,pre,tavg,&
	                     use_float)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(in) :: AV       ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer(in),optional,intent(in) :: nx   ! 2d grid size if available
    integer(in),optional,intent(in) :: ny   ! 2d grid size if available
    integer(in),optional,intent(in) :: nt   ! time sample
    real(r8),optional,intent(in) :: fillval ! fill value
    character(len=*),optional,intent(in) :: pre      ! prefix to variable name
    logical,optional,intent(in) :: tavg     ! is this a tavg
    logical,optional,intent(in) :: use_float ! write output as float rather than double

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: nf,ns,ng
    integer(in) :: i,j,k,n
    integer(in),target  :: dimid2(2)
    integer(in),target  :: dimid3(3)
    integer(in),pointer :: dimid(:)
    type(var_desc_t) :: varid
    type(io_desc_t) :: iodesc
    integer(kind=Pio_Offset_Kind) :: frame
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    character(CL)    :: name1       ! var name
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    character(CL)    :: lpre        ! local prefix
    logical :: exists
    logical :: lwhead, lwdata
    logical :: luse_float
    integer(in) :: lnx,lny
    real(r8) :: lfillvalue
    character(*),parameter :: subName = '(seq_io_write_av) '
    integer :: lbnum
    integer, pointer :: Dof(:)

    real(r8), allocatable :: tmpdata(:)

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lfillvalue = fillvalue
    if (present(fillval)) then
       lfillvalue = fillval
    endif

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    luse_float = .false.
    if (present(use_float)) luse_float = use_float

    call seq_comm_setptrs(CPLID,iam=iam)	

    ng = mct_gsmap_gsize(gsmap)
    lnx = ng
    lny = 1
	
    nf = mct_aVect_nRattr(AV)
    if (nf < 1) then
       write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
       call shr_sys_abort()
    endif
    frame = -1
    if (present(nt)) then
       frame = nt
    endif
    if (present(nx)) then
       if (nx /= 0) lnx = nx
    endif
    if (present(ny)) then
       if (ny /= 0) lny = ny
    endif
    if (lnx*lny /= ng) then
       if(iam==0) write(logunit,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny,trim(dname)
       call shr_sys_abort()
    endif

    if (lwhead) then
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))

       if (present(nt)) then
          dimid3(1:2) = dimid2
          rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
          dimid => dimid3
       else
          dimid => dimid2
       endif

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AV)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
	  if (luse_float) then 
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
          else
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
          end if
          rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
          rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
          rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
          if (present(tavg)) then
             if (tavg) then
                rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
             endif
          endif
       enddo
       if (lwdata) call seq_io_enddef(filename)
    end if

    if (lwdata) then
       call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
       call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       ns = size(dof)
       deallocate(dof)
       allocate(tmpdata(ns))
       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AV)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
          call pio_setframe(cpl_io_file,varid,frame)
          tmpdata = av%rattr(k,:)
          call pio_write_darray(cpl_io_file, varid, iodesc, tmpdata, rcode, fillval=lfillvalue)
       enddo
       deallocate(tmpdata)
       call pio_freedecomp(cpl_io_file, iodesc)

    end if
  end subroutine seq_io_write_av

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_avs - write AVS to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_avs(filename,gsmap,AVS,dname,whead,wdata,nx,ny,nt,fillval,pre,tavg,&
	                     use_float)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(in) :: AVS(:)   ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer(in),optional,intent(in) :: nx   ! 2d grid size if available
    integer(in),optional,intent(in) :: ny   ! 2d grid size if available
    integer(in),optional,intent(in) :: nt   ! time sample
    real(r8),optional,intent(in) :: fillval ! fill value
    character(len=*),optional,intent(in) :: pre      ! prefix to variable name
    logical,optional,intent(in) :: tavg     ! is this a tavg
    logical,optional,intent(in) :: use_float ! write output as float rather than double

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: nf,ns,ng,ni
    integer(in) :: i,j,k,n,k1,k2
    integer(in),target  :: dimid2(2)
    integer(in),target  :: dimid3(3)
    integer(in),target  :: dimid4(4)
    integer(in),pointer :: dimid(:)
    type(var_desc_t) :: varid
    type(io_desc_t) :: iodesc
    integer(kind=Pio_Offset_Kind) :: frame
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    character(CL)    :: name1       ! var name
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    character(CL)    :: lpre        ! local prefix
    logical :: exists
    logical :: lwhead, lwdata
    logical :: luse_float
    integer(in) :: lnx,lny
    real(r8) :: lfillvalue
    real(r8), allocatable :: data(:)
    character(*),parameter :: subName = '(seq_io_write_avs) '
    integer :: lbnum
    integer, pointer :: Dof(:)
    integer, pointer :: Dofn(:)

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lfillvalue = fillvalue
    if (present(fillval)) then
       lfillvalue = fillval
    endif

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    luse_float = .false.
    if (present(use_float)) luse_float = use_float

    call seq_comm_setptrs(CPLID,iam=iam)	

    ni = size(AVS)

    ns = mct_aVect_lsize(AVS(1))
    ng = mct_gsmap_gsize(gsmap)
    lnx = ng
    lny = 1
	
    nf = mct_aVect_nRattr(AVS(1))
    if (nf < 1) then
       write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
       call shr_sys_abort()
    endif
    frame = -1
    if (present(nt)) then
       frame = nt
    endif

    if (present(nx)) then
       if (nx /= 0) lnx = nx
    endif
    if (present(ny)) then
       if (ny /= 0) lny = ny
    endif
    if (lnx*lny /= ng) then
       if(iam==0) write(logunit,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny,trim(dname)
       call shr_sys_abort()
    endif

    if (lwhead) then
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))

       if (ni > 1) then
          rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ni',ni,dimid3(3))
          if (present(nt)) then
             dimid4(1:2) = dimid2
             dimid4(3) = dimid3(3)
             rcode = pio_inq_dimid(cpl_io_file,'time',dimid4(4))
             dimid => dimid4
          else
             dimid3(1:2) = dimid2
             dimid => dimid3
          endif
       else
          if (present(nt)) then
             dimid3(1:2) = dimid2
             rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
             dimid => dimid3
          else
             dimid => dimid2
          endif
       endif

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AVS(1))
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
	  if (luse_float) then 
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
          else
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
          end if
          rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
          rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
          rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
          if (present(tavg)) then
             if (tavg) then
                rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
             endif
          endif
       enddo
       if (lwdata) call seq_io_enddef(filename)
    end if

    if (lwdata) then
       allocate(data(ns*ni))
       ! note: size of dof is ns
       call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
       if (ni > 1) then
          allocate(dofn(ns*ni))
          n = 0
          do k1 = 1,ni
             dofn(n+1:n+ns) = (k1-1)*ng + dof(:)
             n = n + ns
          enddo
          call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,ni/), dofn, iodesc)
          deallocate(dofn)
       else
          call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       endif
       deallocate(dof)

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,AVS(1))
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
          call pio_setframe(cpl_io_file,varid,frame)
          n = 0
          do k1 = 1,ni
             data(n+1:n+ns) = AVS(k1)%rAttr(k,:)
             n = n + ns
          enddo
         call pio_write_darray(cpl_io_file, varid, iodesc, data, rcode, fillval=lfillvalue)
         call pio_setdebuglevel(0)
       enddo

       deallocate(data)
       call pio_freedecomp(cpl_io_file, iodesc)

    end if
  end subroutine seq_io_write_avs

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_write_avs - write AVS to netcdf file
  !
  ! !DESCRIPTION:
  !    Write AV to netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_write_avscomp(filename, comp, flow, dname, &
       whead, wdata, nx, ny, nt, fillval, pre, tavg, use_float)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*) ,intent(in)          :: filename  ! file
    type(component_type)  ,intent(in)          :: comp(:)   ! data to be written
    character(len=3) ,intent(in)          :: flow      ! 'c2x' or 'x2c'
    character(len=*) ,intent(in)          :: dname     ! name of data
    logical          ,optional,intent(in) :: whead     ! write header
    logical          ,optional,intent(in) :: wdata     ! write data
    integer(in)      ,optional,intent(in) :: nx        ! 2d grid size if available
    integer(in)      ,optional,intent(in) :: ny        ! 2d grid size if available
    integer(in)      ,optional,intent(in) :: nt        ! time sample
    real(r8)         ,optional,intent(in) :: fillval   ! fill value
    character(len=*) ,optional,intent(in) :: pre       ! prefix to variable name
    logical          ,optional,intent(in) :: tavg      ! is this a tavg
    logical          ,optional,intent(in) :: use_float ! write output as float rather than double

    !EOP

    type(mct_gsMap), pointer :: gsmap     ! global seg map on coupler processes
    type(mct_avect), pointer :: avcomp1
    type(mct_avect), pointer :: avcomp
    integer(in)              :: rcode
    integer(in)              :: mpicom
    integer(in)              :: iam
    integer(in)              :: nf,ns,ng,ni
    integer(in)              :: i,j,k,n,k1,k2
    integer(in),target       :: dimid2(2)
    integer(in),target       :: dimid3(3)
    integer(in),target       :: dimid4(4)
    integer(in),pointer      :: dimid(:)
    type(var_desc_t)         :: varid
    type(io_desc_t)          :: iodesc
    integer(kind=Pio_Offset_Kind) :: frame
    type(mct_string)         :: mstring     ! mct char type
    character(CL)            :: itemc       ! string converted to char
    character(CL)            :: name1       ! var name
    character(CL)            :: cunit       ! var units
    character(CL)            :: lname       ! long name
    character(CL)            :: sname       ! standard name
    character(CL)            :: lpre        ! local prefix
    logical                  :: exists
    logical                  :: lwhead, lwdata
    logical                  :: luse_float
    integer(in)              :: lnx,lny
    real(r8)                 :: lfillvalue
    real(r8), allocatable    :: data(:)
    character(*),parameter   :: subName = '(seq_io_write_avs) '
    integer                  :: lbnum
    integer, pointer         :: Dof(:)
    integer, pointer         :: Dofn(:)

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lfillvalue = fillvalue
    if (present(fillval)) then
       lfillvalue = fillval
    endif

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata
    frame = -1
    if (present(nt)) then
       frame = nt
    endif

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    luse_float = .false.
    if (present(use_float)) luse_float = use_float

    call seq_comm_setptrs(CPLID,iam=iam)	

    ni = size(comp)
    if (trim(flow) == 'x2c') avcomp1 => component_get_x2c_cx(comp(1))  
    if (trim(flow) == 'c2x') avcomp1 => component_get_c2x_cx(comp(1))  
    gsmap => component_get_gsmap_cx(comp(1))  
    ns = mct_aVect_lsize(avcomp1)
    ng = mct_gsmap_gsize(gsmap)
    lnx = ng
    lny = 1
	
    nf = mct_aVect_nRattr(avcomp1)
    if (nf < 1) then
       write(logunit,*) subname,' ERROR: nf = ',nf,trim(dname)
       call shr_sys_abort()
    endif

    if (present(nx)) then
       if (nx /= 0) lnx = nx
    endif
    if (present(ny)) then
       if (ny /= 0) lny = ny
    endif
    if (lnx*lny /= ng) then
       if(iam==0) then
          write(logunit,*) subname,' ERROR: grid2d size not consistent ',&
               ng,lnx,lny,trim(dname)
       end if
       call shr_sys_abort()
    endif

    if (lwhead) then
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_nx',lnx,dimid2(1))
       rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ny',lny,dimid2(2))

       if (ni > 1) then
          rcode = pio_def_dim(cpl_io_file,trim(lpre)//'_ni',ni,dimid3(3))
          if (present(nt)) then
             dimid4(1:2) = dimid2
             dimid4(3) = dimid3(3)
             rcode = pio_inq_dimid(cpl_io_file,'time',dimid4(4))
             dimid => dimid4
          else
             dimid3(1:2) = dimid2
             dimid => dimid3
          endif
       else
          if (present(nt)) then
             dimid3(1:2) = dimid2
             rcode = pio_inq_dimid(cpl_io_file,'time',dimid3(3))
             dimid => dimid3
          else
             dimid => dimid2
          endif
       endif

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,avcomp1)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          call seq_flds_lookup(itemc,longname=lname,stdname=sname,units=cunit)
	  if (luse_float) then 
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_REAL,dimid,varid)
          else
             rcode = pio_def_var(cpl_io_file,trim(name1),PIO_DOUBLE,dimid,varid)
          end if
          rcode = pio_put_att(cpl_io_file,varid,"_FillValue",lfillvalue)
          rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
          rcode = pio_put_att(cpl_io_file,varid,"internal_dname",trim(dname))
          if (present(tavg)) then
             if (tavg) then
                rcode = pio_put_att(cpl_io_file,varid,"cell_methods","time: mean")
             endif
          endif
       enddo
       if (lwdata) call seq_io_enddef(filename)
    end if

    if (lwdata) then
       allocate(data(ns*ni))
       ! note: size of dof is ns
       call mct_gsmap_OrderedPoints(gsmap, iam, Dof)
       if (ni > 1) then
          allocate(dofn(ns*ni))
          n = 0
          do k1 = 1,ni
          do k2 = 1,ns
             n = n + 1
             dofn(n) = (k1-1)*ng + dof(k2)
          enddo
          enddo
          call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,ni/), dofn, iodesc)
          deallocate(dofn)
       else
          call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       endif
       deallocate(dof)

       do k = 1,nf
          call mct_aVect_getRList(mstring,k,avcomp1)
          itemc = mct_string_toChar(mstring)
          call mct_string_clean(mstring)
          name1 = trim(lpre)//'_'//trim(itemc)
          rcode = pio_inq_varid(cpl_io_file,trim(name1),varid)
          call pio_setframe(cpl_io_file,varid,frame)
          n = 0
          do k1 = 1,ni
             if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(k1))  
             if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(k1))  
             do k2 = 1,ns
                n = n + 1
                data(n) = avcomp%rAttr(k,k2)
             enddo
          enddo
          call pio_write_darray(cpl_io_file, varid, iodesc, data, rcode, fillval=lfillvalue)
       enddo

       deallocate(data)
       call pio_freedecomp(cpl_io_file, iodesc)

    end if
  end subroutine seq_io_write_avscomp

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

  subroutine seq_io_write_int(filename,idata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: exists
    logical :: lwhead, lwdata
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

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
!       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
!       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,idata)

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

  subroutine seq_io_write_int1d(filename,idata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    integer(in)     ,intent(in) :: idata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
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

    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(idata)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_INT,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,idata)
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

  subroutine seq_io_write_r8(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: iam
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical :: exists
    logical :: lwhead, lwdata
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
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
!       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',1,dimid(1))
!       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)


       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,varid)
       if(rcode==PIO_NOERR) then
          rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
          rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
          rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
          if (lwdata) call seq_io_enddef(filename)
       end if
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,rdata)
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

  subroutine seq_io_write_r81d(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    real(r8)        ,intent(in) :: rdata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
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
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(rdata)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,rdata)

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

  subroutine seq_io_write_char(filename,rdata,dname,whead,wdata)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    character(len=*),intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data

    !EOP

    integer(in) :: rcode
    integer(in) :: mpicom
    integer(in) :: iam
    integer(in) :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer(in) :: lnx
    logical :: exists
    logical :: lwhead, lwdata
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
    call seq_comm_setptrs(CPLID,iam=iam)

    if (lwhead) then
       call seq_flds_lookup(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = len(charvar)
       rcode = pio_def_dim(cpl_io_file,trim(dname)//'_len',lnx,dimid(1))
       rcode = pio_def_var(cpl_io_file,trim(dname),PIO_CHAR,dimid,varid)
       rcode = pio_put_att(cpl_io_file,varid,"units",trim(cunit))
       rcode = pio_put_att(cpl_io_file,varid,"long_name",trim(lname))
       rcode = pio_put_att(cpl_io_file,varid,"standard_name",trim(sname))
       if (lwdata) call seq_io_enddef(filename)
    endif

    if (lwdata) then
       charvar = ''
       charvar = trim(rdata)
       rcode = pio_inq_varid(cpl_io_file,trim(dname),varid)
       rcode = pio_put_var(cpl_io_file,varid,charvar)
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

subroutine seq_io_write_time(filename,time_units,time_cal,time_val,nt,whead,wdata,tbnds)

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

!EOP

   integer(in) :: rcode
   integer(in) :: iam
   integer(in) :: dimid(1)
   integer(in) :: dimid2(2)
   type(var_desc_t) :: varid
   integer(in) :: lnx
   logical :: exists
   logical :: lwhead, lwdata
   integer :: start(4),count(4)
   character(len=shr_cal_calMaxLen) :: lcalendar
   real(r8) :: time_val_1d(1)
   character(*),parameter :: subName = '(seq_io_write_time) '

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

   call seq_comm_setptrs(CPLID,iam=iam)

   if (lwhead) then
      rcode = pio_def_dim(cpl_io_file,'time',PIO_UNLIMITED,dimid(1))
      rcode = pio_def_var(cpl_io_file,'time',PIO_DOUBLE,dimid,varid)
      rcode = pio_put_att(cpl_io_file,varid,'units',trim(time_units))
      lcalendar = shr_cal_calendarName(time_cal,trap=.false.)
      if (trim(lcalendar) == trim(shr_cal_noleap)) then
         lcalendar = 'noleap'
      elseif (trim(lcalendar) == trim(shr_cal_gregorian)) then
         lcalendar = 'gregorian'
      endif
      rcode = pio_put_att(cpl_io_file,varid,'calendar',trim(lcalendar))
      if (present(tbnds)) then
         rcode = pio_put_att(cpl_io_file,varid,'bounds','time_bnds')
         dimid2(2)=dimid(1)
         rcode = pio_def_dim(cpl_io_file,'ntb',2,dimid2(1))
         rcode = pio_def_var(cpl_io_file,'time_bnds',PIO_DOUBLE,dimid2,varid)
      endif
      if (lwdata) call seq_io_enddef(filename)
   endif

   if (lwdata) then
      start = 1
      count = 1
      if (present(nt)) then
         start(1) = nt
      endif
      time_val_1d(1) = time_val
      rcode = pio_inq_varid(cpl_io_file,'time',varid)
      rcode = pio_put_var(cpl_io_file,varid,start,count,time_val_1d)
      if (present(tbnds)) then
         rcode = pio_inq_varid(cpl_io_file,'time_bnds',varid)
         start = 1
         count = 1
         if (present(nt)) then
            start(2) = nt
         endif
         count(1) = 2
         rcode = pio_put_var(cpl_io_file,varid,start,count,tbnds)
      endif

      !      write(logunit,*) subname,' wrote time ',lwhead,lwdata
   endif

 end subroutine seq_io_write_time

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_av - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_av(filename,gsmap,AV,dname,pre)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(inout):: AV     ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    character(len=*),intent(in),optional :: pre      ! prefix name

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    integer(in) :: nf,ns,ng
    integer(in) :: i,j,k,n, ndims
    type(file_desc_t) :: pioid
    integer(in) :: dimid(2)
    type(var_desc_t) :: varid
    integer(in) :: lnx,lny
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    logical :: exists
    type(io_desc_t) :: iodesc
    integer(in), pointer :: dof(:)
    character(CL)  :: lversion
    character(CL)  :: name1
    character(CL)  :: lpre
    character(*),parameter :: subName = '(seq_io_read_av) '
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lversion = trim(version0)

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
    call mct_gsmap_OrderedPoints(gsmap, iam, Dof)

    ns = mct_aVect_lsize(AV)
    nf = mct_aVect_nRattr(AV)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_av exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    do k = 1,nf
       call mct_aVect_getRList(mstring,k,AV)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       if (trim(lversion) == trim(version)) then
          name1 = trim(lpre)//'_'//trim(itemc)
       else
          name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
       endif
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid,trim(name1),varid)
       if (rcode == pio_noerr) then
          if (k==1) then
             rcode = pio_inq_varndims(pioid, varid, ndims)
             rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
             rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
             if (ndims>=2) then
                rcode = pio_inq_dimlen(pioid, dimid(2), lny)
             else
                lny = 1
             end if
             ng = lnx * lny
             if (ng /= mct_gsmap_gsize(gsmap)) then
                if (iam==0) write(logunit,*) subname,' ERROR: dimensions do not match',&
                     lnx,lny,mct_gsmap_gsize(gsmap)
                call shr_sys_abort()
             end if
             call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
             deallocate(dof)
          end if
          call pio_read_darray(pioid,varid,iodesc, av%rattr(k,:), rcode)
       else
          write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
          write(logunit,*)'for backwards compatibility will set it to 0'
          av%rattr(k,:) = 0.0_r8
       end if
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
          
    enddo

    !--- zero out fill value, this is somewhat arbitrary
    do n = 1,ns
    do k = 1,nf
       if (AV%rAttr(k,n) == fillvalue) then
           AV%rAttr(k,n) = 0.0_r8
       endif
    enddo
    enddo

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine seq_io_read_av

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_avs - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_avs(filename,gsmap,AVS,dname,pre)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*),intent(in) :: filename ! file
    type(mct_gsMap), intent(in) :: gsmap
    type(mct_aVect) ,intent(inout):: AVS(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    character(len=*),intent(in),optional :: pre      ! prefix name

    !EOP

    integer(in) :: rcode
    integer(in) :: iam,mpicom
    integer(in) :: nf,ns,ng,ni
    integer(in) :: i,j,k,n,n1,n2,ndims
    type(file_desc_t) :: pioid
    integer(in) :: dimid(4)
    type(var_desc_t) :: varid
    integer(in) :: lnx,lny,lni
    type(mct_string) :: mstring     ! mct char type
    character(CL)    :: itemc       ! string converted to char
    logical :: exists
    type(io_desc_t) :: iodesc
    integer(in), pointer :: dof(:)
    integer(in), pointer :: dofn(:)
    real(r8), allocatable :: data(:)
    character(CL)  :: lversion
    character(CL)  :: name1
    character(CL)  :: lpre
    character(*),parameter :: subName = '(seq_io_read_avs) '
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lversion = trim(version0)

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
    call mct_gsmap_OrderedPoints(gsmap, iam, Dof)

    ni = size(AVS)
    ns = mct_aVect_lsize(AVS(1))
    nf = mct_aVect_nRattr(AVS(1))
    ng = mct_gsmap_gsize(gsmap)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_avs exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    allocate(data(ni*ns))

    do k = 1,nf
       call mct_aVect_getRList(mstring,k,AVS(1))
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       if (trim(lversion) == trim(version)) then
          name1 = trim(lpre)//'_'//trim(itemc)
       else
          name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
       endif
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid,trim(name1),varid)
       if (rcode == pio_noerr) then
          if (k==1) then
             rcode = pio_inq_varndims(pioid, varid, ndims)
             rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
             rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
             if (ndims>=2) then
                rcode = pio_inq_dimlen(pioid, dimid(2), lny)
             else
                lny = 1
             end if
             if (lnx*lny /= ng) then
                write(logunit,*) subname,' ERROR: dimensions do not match',&
                     lnx,lny,mct_gsmap_gsize(gsmap)
                call shr_sys_abort()
             end if
             if (ndims>=3) then
                rcode = pio_inq_dimlen(pioid, dimid(3), lni)
             else
                lni = 1
             end if
             if (ni /= lni) then
                write(logunit,*) subname,' ERROR: ni dimensions do not match',ni,lni
                call shr_sys_abort()
             end if
             if (ni > 1) then
                allocate(dofn(ns*ni))
                n = 0
                do n1 = 1,ni
                do n2 = 1,ns
                   n = n + 1
                   dofn(n) = (n1-1)*ng + dof(n2)
                enddo
                enddo
                call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,lni/), dofn, iodesc)
                deallocate(dofn)
             else
                call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
             endif
             deallocate(dof)
          end if

          call pio_read_darray(pioid,varid,iodesc, data, rcode)
          n = 0
          do n1 = 1,ni
          do n2 = 1,ns
             n = n + 1
             avs(n1)%rAttr(k,n2) = data(n)
          enddo
          enddo
       else
          write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
          write(logunit,*)'for backwards compatibility will set it to 0'
          do n1 = 1,ni
             avs(n1)%rattr(k,:) = 0.0_r8
          enddo
       end if
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    enddo

    deallocate(data)

    !--- zero out fill value, this is somewhat arbitrary
    do n1 = 1,ni
    do n2 = 1,ns
    do k = 1,nf
       if (AVS(n1)%rAttr(k,n2) == fillvalue) then
           AVS(n1)%rAttr(k,n2) = 0.0_r8
       endif
    enddo
    enddo
    enddo

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine seq_io_read_avs

!===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: seq_io_read_avs - read AV from netcdf file
  !
  ! !DESCRIPTION:
  !    Read AV from netcdf file
  !
  ! !REVISION HISTORY:
  !    2007-Oct-26 - T. Craig - initial version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine seq_io_read_avscomp(filename, comp, flow, dname, pre)

    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    character(len=*), intent(in)          :: filename ! file
    type(component_type),  intent(inout)       :: comp(:)
    character(len=3), intent(in)          :: flow     ! 'c2x' or 'x2c'
    character(len=*), intent(in)          :: dname    ! name of data
    character(len=*), intent(in),optional :: pre      ! prefix name

    !EOP

    type(mct_gsMap), pointer :: gsmap
    type(mct_aVect), pointer :: avcomp  
    type(mct_aVect), pointer :: avcomp1 
    integer(in)              :: rcode
    integer(in)              :: iam,mpicom
    integer(in)              :: nf,ns,ng,ni
    integer(in)              :: i,j,k,n,n1,n2,ndims
    type(file_desc_t)        :: pioid
    integer(in)              :: dimid(4)
    type(var_desc_t)         :: varid
    integer(in)              :: lnx,lny,lni
    type(mct_string)         :: mstring ! mct char type
    character(CL)            :: itemc   ! string converted to char
    logical                  :: exists
    type(io_desc_t)          :: iodesc
    integer(in), pointer     :: dof(:)
    integer(in), pointer     :: dofn(:)
    real(r8), allocatable    :: data(:)
    character(CL)            :: lversion
    character(CL)            :: name1
    character(CL)            :: lpre
    character(*),parameter   :: subName = '(seq_io_read_avs) '
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lversion = trim(version0)

    lpre = trim(dname)
    if (present(pre)) then
       lpre = trim(pre)
    endif

    gsmap => component_get_gsmap_cx(comp(1))
    if (trim(flow) == 'x2c') avcomp1 => component_get_x2c_cx(comp(1))  
    if (trim(flow) == 'c2x') avcomp1 => component_get_c2x_cx(comp(1))  

    call seq_comm_setptrs(CPLID,iam=iam,mpicom=mpicom)
    call mct_gsmap_OrderedPoints(gsmap, iam, Dof)

    ni = size(comp)
    ns = mct_aVect_lsize(avcomp1)
    nf = mct_aVect_nRattr(avcomp1)
    ng = mct_gsmap_gsize(gsmap)

    if (iam==0) inquire(file=trim(filename),exist=exists)
    call shr_mpi_bcast(exists,mpicom,'seq_io_read_avs exists')
    if (exists) then
       rcode = pio_openfile(cpl_io_subsystem, pioid, cpl_pio_iotype, trim(filename),pio_nowrite)
       if(iam==0) write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call shr_sys_abort()
    endif

    allocate(data(ni*ns))

    do k = 1,nf
       call mct_aVect_getRList(mstring,k,avcomp1)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       if (trim(lversion) == trim(version)) then
          name1 = trim(lpre)//'_'//trim(itemc)
       else
          name1 = trim(prefix)//trim(dname)//'_'//trim(itemc)
       endif
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid,trim(name1),varid)
       if (rcode == pio_noerr) then
          if (k==1) then
             rcode = pio_inq_varndims(pioid, varid, ndims)
             rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
             rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
             if (ndims>=2) then
                rcode = pio_inq_dimlen(pioid, dimid(2), lny)
             else
                lny = 1
             end if
             if (lnx*lny /= ng) then
                write(logunit,*) subname,' ERROR: dimensions do not match',&
                     lnx,lny,mct_gsmap_gsize(gsmap)
                call shr_sys_abort()
             end if
             if (ndims>=3) then
                rcode = pio_inq_dimlen(pioid, dimid(3), lni)
             else
                lni = 1
             end if
             if (ni /= lni) then
                write(logunit,*) subname,' ERROR: ni dimensions do not match',ni,lni
                call shr_sys_abort()
             end if
             if (ni > 1) then
                allocate(dofn(ns*ni))
                n = 0
                do n1 = 1,ni
                do n2 = 1,ns
                   n = n + 1
                   dofn(n) = (n1-1)*ng + dof(n2)
                enddo
                enddo
                call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny,lni/), dofn, iodesc)
                deallocate(dofn)
             else
                call pio_initdecomp(cpl_io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
             endif
             deallocate(dof)
          end if

          call pio_read_darray(pioid,varid,iodesc, data, rcode)
          n = 0
          do n1 = 1,ni
             if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
             if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
             do n2 = 1,ns
                n = n + 1
                avcomp%rAttr(k,n2) = data(n)
             enddo
          enddo
       else
          write(logunit,*)'seq_io_readav warning: field ',trim(itemc),' is not on restart file'
          write(logunit,*)'for backwards compatibility will set it to 0'
          do n1 = 1,ni
             if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
             if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
             avcomp%rattr(k,:) = 0.0_r8
          enddo
       end if
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    enddo

    deallocate(data)

    !--- zero out fill value, this is somewhat arbitrary
    do n1 = 1,ni
       if (trim(flow) == 'x2c') avcomp => component_get_x2c_cx(comp(n1))  
       if (trim(flow) == 'c2x') avcomp => component_get_c2x_cx(comp(n1))  
       do n2 = 1,ns
          do k = 1,nf
             if (avcomp%rAttr(k,n2) == fillvalue) then
                avcomp%rAttr(k,n2) = 0.0_r8
             endif
          enddo
       enddo
    enddo

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

  end subroutine seq_io_read_avscomp

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
       call shr_sys_abort()
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
       call shr_sys_abort()
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
       call shr_sys_abort()
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
!===============================================================================
end module seq_io_mod
