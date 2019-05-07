module med_io_mod
  ! !DESCRIPTION: Writes attribute vectors to netcdf

  ! !USES:
  use ESMF                  , only : ESMF_VM, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
  use pio                   , only : file_desc_t, iosystem_desc_t
  use med_constants_mod     , only : R4, R8, CL
  use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
  use shr_nuopc_utils_mod   , only : chkerr       => shr_nuopc_utils_ChkErr
  use shr_nuopc_methods_mod , only : FB_getFieldN => shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only : FB_getFldPtr => shr_nuopc_methods_FB_getFldPtr
  use shr_nuopc_methods_mod , only : FB_getNameN  => shr_nuopc_methods_FB_getNameN

  implicit none
  private

  ! public member functions:
  public med_io_wopen
  public med_io_close
  public med_io_redef
  public med_io_enddef
  public med_io_date2yyyymmdd
  public med_io_sec2hms
  public med_io_read
  public med_io_write
  public med_io_init

  ! private member functions
  private :: med_io_file_exists

  ! public data members:
  interface med_io_read
     module procedure med_io_read_FB
     module procedure med_io_read_int
     module procedure med_io_read_int1d
     module procedure med_io_read_r8
     module procedure med_io_read_r81d
     module procedure med_io_read_char
  end interface med_io_read
  interface med_io_write
     module procedure med_io_write_FB
     module procedure med_io_write_int
     module procedure med_io_write_int1d
     module procedure med_io_write_r8
     module procedure med_io_write_r81d
     module procedure med_io_write_char
     module procedure med_io_write_time
  end interface med_io_write

  !-------------------------------------------------------------------------------
  ! module data
  !-------------------------------------------------------------------------------

  character(*),parameter         :: prefix    = "med_io_"
  character(*),parameter         :: modName   = "(med_io_mod) "
  character(*),parameter         :: version   = "cmeps0"
  integer    , parameter         :: file_desc_t_cnt = 20 ! Note - this is hard-wired for now
  integer    , parameter         :: number_strlen = 2
  character(CL)                  :: wfilename = ''
  type(file_desc_t)              :: io_file(0:file_desc_t_cnt)
  integer                        :: pio_iotype
  integer                        :: pio_ioformat
  type(iosystem_desc_t), pointer :: io_subsystem
  character(*),parameter         :: u_file_u = &
       __FILE__

!=================================================================================
contains
!=================================================================================

  logical function med_io_file_exists(vm, iam, filename)

    !---------------
    ! inquire if i/o file exists
    !---------------

    use ESMF, only : ESMF_VMBroadCast

    ! input/output variables
    type(ESMF_VM)                :: vm
    integer,          intent(in) :: iam
    character(len=*), intent(in) :: filename

    ! local variables
    integer :: tmp(1)
    integer :: rc
    !-------------------------------------------------------------------------------

    med_io_file_exists = .false.
    if (iam==0) inquire(file=trim(filename),exist=med_io_file_exists)
    if (med_io_file_exists) tmp(1) = 1

    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(tmp(1) == 1) med_io_file_exists = .true.

  end function med_io_file_exists

  !===============================================================================
  subroutine med_io_init()

    !---------------
    ! initialize pio
    !---------------

    use shr_pio_mod           , only : shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
    use med_internalstate_mod , only : med_id

    io_subsystem => shr_pio_getiosys(med_id)
    pio_iotype   =  shr_pio_getiotype(med_id)
    pio_ioformat =  shr_pio_getioformat(med_id)

  end subroutine med_io_init

  !===============================================================================
  subroutine med_io_wopen(filename, vm, iam, clobber, file_ind, model_doi_url)

    !---------------
    ! open netcdf file
    !---------------

    use pio                   , only : PIO_IOTYPE_PNETCDF, PIO_IOTYPE_NETCDF, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio                   , only : pio_openfile, pio_createfile, PIO_GLOBAL, pio_enddef
    use PIO                   , only : pio_put_att, pio_redef, pio_get_att
    use pio                   , only : pio_seterrorhandling, pio_file_is_open, pio_clobber, pio_write, pio_noclobber
    use med_internalstate_mod , only : logunit

    ! input/output arguments
    character(*),            intent(in) :: filename
    type(ESMF_VM)                       :: vm
    integer,                 intent(in) :: iam
    logical,       optional, intent(in) :: clobber
    integer,       optional, intent(in) :: file_ind
    character(CL), optional, intent(in) :: model_doi_url

    ! local variables
    logical       :: lclobber
    integer       :: tmp(1)
    integer       :: rcode
    integer       :: nmode
    integer       :: lfile_ind
    integer       :: rc
    character(CL) :: lversion
    character(CL) :: lmodel_doi_url
    character(*),parameter :: subName = '(med_io_wopen) '
    !-------------------------------------------------------------------------------

    lversion=trim(version)

    lclobber = .false.
    if (present(clobber)) lclobber=clobber

    lmodel_doi_url = 'unset'
    if (present(model_doi_url)) lmodel_doi_url = model_doi_url

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    if (.not. pio_file_is_open(io_file(lfile_ind))) then

       ! filename not open
       wfilename = filename

       if (med_io_file_exists(vm, iam, filename)) then
          if (lclobber) then
             nmode = pio_clobber
             ! only applies to classic NETCDF files.
             if(pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
                nmode = ior(nmode,pio_ioformat)
             endif
             rcode = pio_createfile(io_subsystem, io_file(lfile_ind), pio_iotype, trim(filename), nmode)
             if(iam==0) write(logunit,*) subname,' create file ',trim(filename)
             rcode = pio_put_att(io_file(lfile_ind),pio_global,"file_version",version)
             rcode = pio_put_att(io_file(lfile_ind),pio_global,"model_doi_url",lmodel_doi_url)
          else
             rcode = pio_openfile(io_subsystem, io_file(lfile_ind), pio_iotype, trim(filename), pio_write)
             if (iam==0) then
                write(logunit,*) subname,' open file ',trim(filename)
             end if
             call pio_seterrorhandling(io_file(lfile_ind),PIO_BCAST_ERROR)
             rcode = pio_get_att(io_file(lfile_ind),pio_global,"file_version",lversion)
             call pio_seterrorhandling(io_file(lfile_ind),PIO_INTERNAL_ERROR)
             if (trim(lversion) /= trim(version)) then
                rcode = pio_redef(io_file(lfile_ind))
                rcode = pio_put_att(io_file(lfile_ind),pio_global,"file_version",version)
                rcode = pio_enddef(io_file(lfile_ind))
             endif
          endif
       else
          nmode = pio_noclobber
          ! only applies to classic NETCDF files.
          if(pio_iotype == PIO_IOTYPE_NETCDF .or. pio_iotype == PIO_IOTYPE_PNETCDF) then
             nmode = ior(nmode,pio_ioformat)
          endif
          rcode = pio_createfile(io_subsystem, io_file(lfile_ind), pio_iotype, trim(filename), nmode)
          if (iam==0) then
             write(logunit,*) subname,' create file ',trim(filename)
          end if
          rcode = pio_put_att(io_file(lfile_ind),pio_global,"file_version",version)
          rcode = pio_put_att(io_file(lfile_ind),pio_global,"model_doi_url",lmodel_doi_url)
       endif
    elseif (trim(wfilename) /= trim(filename)) then
       ! filename is open, better match open filename
       if(iam==0) write(logunit,*) subname,' different  filename currently open ',trim(filename)
       if(iam==0) write(logunit,*) subname,' different wfilename currently open ',trim(wfilename)
       call ESMF_LogWrite(subname//'different file currently open '//trim(filename), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    else
       ! filename is already open, just return
    endif

  end subroutine med_io_wopen

  !===============================================================================
  subroutine med_io_close(filename, iam, file_ind, rc)

    !---------------
    ! close netcdf file
    !---------------

    use pio, only: pio_file_is_open, pio_closefile
    use med_internalstate_mod, only : logunit

    ! input/output variables
    character(*),     intent(in)  :: filename
    integer,          intent(in)  :: iam
    integer,optional, intent(in)  :: file_ind
    integer         , intent(out) :: rc

    ! local variables
    integer :: lfile_ind
    character(*),parameter :: subName = '(med_io_close) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    if (.not. pio_file_is_open(io_file(lfile_ind))) then
       ! filename not open, just return
    elseif (trim(wfilename) == trim(filename)) then
       ! filename matches, close it
       call pio_closefile(io_file(lfile_ind))
    else
       ! different filename is open, abort
       if (iam==0) write(logunit,*) subname,' different  filename currently open, aborting ',trim(filename)
       if (iam==0) write(logunit,*) subname,' different wfilename currently open, aborting ',trim(wfilename)
       call ESMF_LogWrite(subname//'different file currently open, aborting '//trim(filename), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif
    wfilename = ''
  end subroutine med_io_close

  !===============================================================================
  subroutine med_io_redef(filename,file_ind)

    use pio, only : pio_redef

    ! input/output variables
    character(len=*), intent(in) :: filename
    integer,optional,intent(in):: file_ind

    ! local variables
    integer :: lfile_ind
    integer :: rcode
    !-------------------------------------------------------------------------------

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind
    rcode = pio_redef(io_file(lfile_ind))

  end subroutine med_io_redef

  !===============================================================================
  subroutine med_io_enddef(filename,file_ind)

    use med_internalstate_mod , only : logunit
    use pio                   , only : pio_enddef

    ! input/output variables
    character(len=*) , intent(in) :: filename
    integer,optional , intent(in) :: file_ind

    ! local variables
    integer :: lfile_ind
    integer :: rcode
    !-------------------------------------------------------------------------------

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    rcode = pio_enddef(io_file(lfile_ind))
  end subroutine med_io_enddef

  !===============================================================================
  character(len=24) function med_io_date2yyyymmdd (date)

    use shr_cal_mod, only : shr_cal_datetod2string

    integer, intent(in) :: date  ! date expressed as an integer: yyyymmdd

    call shr_cal_datetod2string(date_str = med_io_date2yyyymmdd, ymd = date)
  end function med_io_date2yyyymmdd

  !===============================================================================
  character(len=8) function med_io_sec2hms (seconds, rc)

    use med_internalstate_mod , only : logunit

    ! input arguments
    integer, intent(in)  :: seconds
    integer, intent(out) :: rc

    ! local variables
    integer :: hours     ! hours of hh:mm:ss
    integer :: minutes   ! minutes of hh:mm:ss
    integer :: secs      ! seconds of hh:mm:ss
    !----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (seconds < 0 .or. seconds > 86400) then
       write(logunit,*)'med_io_sec2hms: bad input seconds:', seconds
       call ESMF_LogWrite('med_io_sec2hms: bad input seconds', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    if (minutes < 0 .or. minutes > 60) then
       write(logunit,*)'med_io_sec2hms: bad minutes = ',minutes
       call ESMF_LogWrite('med_io_sec2hms: bad minutes', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    if (secs < 0 .or. secs > 60) then
       write(logunit,*)'med_io_sec2hms: bad secs = ',secs
       call ESMF_LogWrite('med_io_sec2hms: bad secs', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    write(med_io_sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)

  end function med_io_sec2hms

  !===============================================================================
  subroutine med_io_write_FB(filename, iam, FB, whead, wdata, nx, ny, nt, &
       fillval, pre, tavg, use_float, file_ind, rc)

    !---------------
    ! Write FB to netcdf file
    !---------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundle, ESMF_Mesh, ESMF_DistGrid
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF                  , only : ESMF_Field, ESMF_FieldGet, ESMF_AttributeGet 
    use shr_const_mod         , only : fillvalue=>SHR_CONST_SPVAL
    use esmFlds               , only : shr_nuopc_fldList_GetMetadata
    use pio                   , only : var_desc_t, io_desc_t, pio_offset_kind
    use pio                   , only : pio_def_dim, pio_inq_dimid, pio_real, pio_def_var, pio_put_att, pio_double
    use pio                   , only : pio_inq_varid, pio_setframe, pio_write_darray, pio_initdecomp, pio_freedecomp
    use pio                   , only : pio_syncfile

    ! input/output variables
    character(len=*),           intent(in) :: filename  ! file
    integer,                    intent(in) :: iam       ! local pet
    type(ESMF_FieldBundle),     intent(in) :: FB        ! data to be written
    logical,          optional, intent(in) :: whead     ! write header
    logical,          optional, intent(in) :: wdata     ! write data
    integer    ,      optional, intent(in) :: nx        ! 2d grid size if available
    integer    ,      optional, intent(in) :: ny        ! 2d grid size if available
    integer    ,      optional, intent(in) :: nt        ! time sample
    real(r8),         optional, intent(in) :: fillval   ! fill value
    character(len=*), optional, intent(in) :: pre       ! prefix to variable name
    logical,          optional, intent(in) :: tavg      ! is this a tavg
    logical,          optional, intent(in) :: use_float ! write output as float rather than double
    integer,          optional, intent(in) :: file_ind
    integer,                    intent(out):: rc

    ! local variables
    type(ESMF_Field)              :: field
    type(ESMF_Mesh)               :: mesh
    type(ESMF_Distgrid)           :: distgrid
    type(ESMF_VM)                 :: VM
    integer                       :: mpicom
    integer                       :: rcode
    integer                       :: nf,ns,ng
    integer                       :: k,n
    integer    ,target            :: dimid2(2)
    integer    ,target            :: dimid3(3)
    integer    ,pointer           :: dimid(:)
    type(var_desc_t)              :: varid
    type(io_desc_t)               :: iodesc
    integer(kind=Pio_Offset_Kind) :: frame
    character(CL)                 :: itemc       ! string converted to char
    character(CL)                 :: name1       ! var name
    character(CL)                 :: cunit       ! var units
    character(CL)                 :: lname       ! long name
    character(CL)                 :: sname       ! standard name
    character(CL)                 :: lpre        ! local prefix
    logical                       :: lwhead, lwdata
    logical                       :: luse_float
    integer                       :: lnx,lny
    real(r8)                      :: lfillvalue
    integer, pointer              :: minIndexPTile(:,:)
    integer, pointer              :: maxIndexPTile(:,:)
    integer                       :: dimCount, tileCount
    integer, pointer              :: Dof(:)
    integer                       :: lfile_ind
    real(r8), pointer             :: fldptr1(:)
    real(r8), pointer             :: fldptr2(:,:)
    character(len=number_strlen)  :: cnumber
    character(CL)                 :: tmpstr
    type(ESMF_Field)              :: lfield
    integer                       :: rank
    integer                       :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fields
    integer                       :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fields
    logical                       :: isPresent
    character(*),parameter :: subName = '(med_io_write_FB) '
    !-------------------------------------------------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_Success

    lfillvalue = fillvalue
    if (present(fillval)) then
       lfillvalue = fillval
    endif

    lpre = ' '
    if (present(pre)) then
       lpre = trim(pre)
    endif

    if (.not. ESMF_FieldBundleIsCreated(FB, rc=rc)) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" not created", ESMF_LOGMSG_INFO)
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       endif
       rc = ESMF_Success
       return
    endif

    lwhead = .true.
    lwdata = .true.
    if (present(whead)) lwhead = whead
    if (present(wdata)) lwdata = wdata

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       endif
       return
    endif

    luse_float = .false.
    if (present(use_float)) luse_float = use_float

    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    call ESMF_FieldBundleGet(FB, fieldCount=nf, rc=rc)
    write(tmpstr,*) subname//' field count = '//trim(lpre),nf
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    if (nf < 1) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" empty", ESMF_LOGMSG_INFO)
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
       endif
       rc = ESMF_Success
       return
    endif

    call FB_getFieldN(FB, 1, field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(field, mesh=mesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
    call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! write(tmpstr,*) subname,' counts = ',dimcount,tilecount,minindexptile,maxindexptile
    ! call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

    ! TODO: this is not getting the global size correct for a FB coming in that does not have
    ! all the global grid values in the distgrid - e.g. CTSM

    ng = maxval(maxIndexPTile)
    lnx = ng
    lny = 1
    deallocate(minIndexPTile, maxIndexPTile)

    frame = -1
    if (present(nt)) then
       frame = nt
    endif
    if (present(nx)) then
       if (nx > 0) lnx = nx
    endif
    if (present(ny)) then
       if (ny > 0) lny = ny
    endif
    if (lnx*lny /= ng) then
       write(tmpstr,*) subname,' ERROR: grid2d size not consistent ',ng,lnx,lny
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       !TODO: this should not be an error for say CTSM which does not send a global grid
       !rc = ESMF_FAILURE
       !return
    endif

    if (lwhead) then
       rcode = pio_def_dim(io_file(lfile_ind),trim(lpre)//'_nx',lnx,dimid2(1))
       rcode = pio_def_dim(io_file(lfile_ind),trim(lpre)//'_ny',lny,dimid2(2))

       if (present(nt)) then
          dimid3(1:2) = dimid2
          rcode = pio_inq_dimid(io_file(lfile_ind),'time',dimid3(3))
          dimid => dimid3
       else
          dimid => dimid2
       endif

       write(tmpstr,*) subname,' dimid = ',dimid
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       do k = 1,nf
          call FB_getNameN(FB, k, itemc, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Determine rank of field with name itemc
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! TODO (mvertens, 2019-03-13): this is a temporary mod to NOT write hgt
          if (trim(itemc) /= "hgt") then
             if (rank == 2) then
                call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                write(cnumber,'(i0)') ungriddedUbound(1)
                call ESMF_LogWrite(trim(subname)//':'//'field '//trim(itemc)// &
                     ' has an griddedUBound of  '//trim(cnumber), ESMF_LOGMSG_INFO) 

                ! Create a new output variable for each element of the undistributed dimension
                do n = 1,ungriddedUBound(1)
                   if (trim(itemc) /= "hgt") then
                      write(cnumber,'(i0)') n
                      name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)
                      call shr_nuopc_fldList_GetMetadata(itemc, longname=lname, stdname=sname, units=cunit)
                      call ESMF_LogWrite(trim(subname)//': defining '//trim(name1), ESMF_LOGMSG_INFO)

                      if (luse_float) then
                         rcode = pio_def_var(io_file(lfile_ind), trim(name1), PIO_REAL, dimid, varid)
                         rcode = pio_put_att(io_file(lfile_ind), varid,"_FillValue",real(lfillvalue,r4))
                      else
                         rcode = pio_def_var(io_file(lfile_ind), trim(name1), PIO_DOUBLE, dimid, varid)
                         rcode = pio_put_att(io_file(lfile_ind),varid,"_FillValue",lfillvalue)
                      end if
                      rcode = pio_put_att(io_file(lfile_ind), varid, "units"        , trim(cunit))
                      rcode = pio_put_att(io_file(lfile_ind), varid, "long_name"    , trim(lname))
                      rcode = pio_put_att(io_file(lfile_ind), varid, "standard_name", trim(sname))
                      if (present(tavg)) then
                         if (tavg) then
                            rcode = pio_put_att(io_file(lfile_ind), varid, "cell_methods", "time: mean")
                         endif
                      endif
                   end if
                end do
             else
                name1 = trim(lpre)//'_'//trim(itemc)
                call shr_nuopc_fldList_GetMetadata(itemc,longname=lname,stdname=sname,units=cunit)
                call ESMF_LogWrite(trim(subname)//':'//trim(itemc)//':'//trim(name1),ESMF_LOGMSG_INFO)
                if (luse_float) then
                   rcode = pio_def_var(io_file(lfile_ind), trim(name1), PIO_REAL, dimid, varid)
                   rcode = pio_put_att(io_file(lfile_ind), varid, "_FillValue", real(lfillvalue, r4))
                else
                   rcode = pio_def_var(io_file(lfile_ind), trim(name1), PIO_DOUBLE, dimid, varid)
                   rcode = pio_put_att(io_file(lfile_ind), varid, "_FillValue", lfillvalue)
                end if
                rcode = pio_put_att(io_file(lfile_ind), varid, "units"         , trim(cunit))
                rcode = pio_put_att(io_file(lfile_ind), varid, "long_name"     , trim(lname))
                rcode = pio_put_att(io_file(lfile_ind), varid, "standard_name" , trim(sname))
                if (present(tavg)) then
                   if (tavg) then
                      rcode = pio_put_att(io_file(lfile_ind), varid, "cell_methods", "time: mean")
                   endif
                end if
             end if
          end if
       end do

       ! Finish define mode
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)

    end if

    if (lwdata) then

       ! use distgrid extracted from field 1 above
       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns,size(dof),dof(1),dof(ns)  !,minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       call pio_initdecomp(io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)

       ! call pio_writedof(lpre, (/lnx,lny/), int(dof,kind=PIO_OFFSET_KIND), mpicom)

       deallocate(dof)

       do k = 1,nf
          call FB_getNameN(FB, k, itemc, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          call FB_getFldPtr(FB, itemc, &
               fldptr1=fldptr1, fldptr2=fldptr2, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! TODO (mvertens, 2019-03-13): this is a temporary mod to NOT write hgt
          if (trim(itemc) /= "hgt") then
             if (rank == 2) then

                ! Determine the size of the ungridded dimension and the index where the undistributed dimension is located
                call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                ! Output for each ungriddedUbound index
                do n = 1,ungriddedUBound(1)
                   write(cnumber,'(i0)') n
                   name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)
                   rcode = pio_inq_varid(io_file(lfile_ind), trim(name1), varid)
                   call pio_setframe(io_file(lfile_ind),varid,frame)
                   
                   if (gridToFieldMap(1) == 1) then
                      call pio_write_darray(io_file(lfile_ind), varid, iodesc, fldptr2(:,n), rcode, fillval=lfillvalue)
                   else if (gridToFieldMap(1) == 2) then
                      call pio_write_darray(io_file(lfile_ind), varid, iodesc, fldptr2(n,:), rcode, fillval=lfillvalue)
                   end if
                end do
             else if (rank == 1) then
                name1 = trim(lpre)//'_'//trim(itemc)
                rcode = pio_inq_varid(io_file(lfile_ind), trim(name1), varid)
                call pio_setframe(io_file(lfile_ind),varid,frame)
                call pio_write_darray(io_file(lfile_ind), varid, iodesc, fldptr1, rcode, fillval=lfillvalue)
             end if  ! end if rank is 2 or 1

          end if ! end if not "hgt"
       end do  ! end loop over fields in FB

       call pio_syncfile(io_file(lfile_ind))
       call pio_freedecomp(io_file(lfile_ind), iodesc)
    endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_io_write_FB

  !===============================================================================
  subroutine med_io_write_int(filename, iam, idata, dname, whead, wdata, file_ind, rc)

    use pio    , only : var_desc_t, pio_def_var, pio_put_att, pio_int, pio_inq_varid, pio_put_var
    use esmFlds, only : shr_nuopc_fldList_GetMetadata

    !---------------
    ! Write scalar integer to netcdf file
    !---------------

    ! intput/output variables
    character(len=*) ,intent(in) :: filename ! file
    integer          ,intent(in) :: iam      ! local pet
    integer          ,intent(in) :: idata    ! data to be written
    character(len=*) ,intent(in) :: dname    ! name of data
    logical,optional ,intent(in) :: whead    ! write header
    logical,optional ,intent(in) :: wdata    ! write data
    integer,optional ,intent(in) :: file_ind
    integer          ,intent(out):: rc

    ! local variables
    integer          :: rcode
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical          :: lwhead, lwdata
    integer          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_int) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

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

    if (lwhead) then
       call shr_nuopc_fldList_GetMetadata(trim(dname),longname=lname,stdname=sname,units=cunit)
       !       rcode = pio_def_dim(io_file(lfile_ind),trim(dname)//'_nx',1,dimid(1))
       !       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_INT,dimid,varid)
       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_INT,varid)
       rcode = pio_put_att(io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,idata)
       !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata
    endif

  end subroutine med_io_write_int

  !===============================================================================
  subroutine med_io_write_int1d(filename, iam, idata, dname, whead, wdata, file_ind, rc)

    !---------------
    ! Write 1d integer array to netcdf file
    !---------------

    use pio     , only : var_desc_t, pio_def_dim, pio_def_var
    use pio     , only : pio_put_att, pio_inq_varid, pio_put_var
    use pio     , only : pio_int, pio_def_var
    use esmFlds , only : shr_nuopc_fldList_GetMetadata

    ! input/output arguments
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(in) :: iam      ! local pet
    integer         ,intent(in) :: idata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind
    integer         , intent(out) :: rc

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer          :: lnx
    logical          :: lwhead, lwdata
    integer          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_int1d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

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

    if (lwhead) then
       call shr_nuopc_fldList_GetMetadata(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(idata)
       rcode = pio_def_dim(io_file(lfile_ind),trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_INT,dimid,varid)
       rcode = pio_put_att(io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,idata)
    endif

    !      write(logunit,*) subname,' wrote AV ',trim(dname),lwhead,lwdata

  end subroutine med_io_write_int1d

  !===============================================================================
  subroutine med_io_write_r8(filename, iam, rdata, dname, whead, wdata, file_ind, rc)

    !---------------
    ! Write scalar double to netcdf file
    !---------------

    use pio               , only : var_desc_t, pio_def_var, pio_put_att
    use pio               , only : pio_double, pio_noerr, pio_inq_varid, pio_put_var
    use esmFlds           , only : shr_nuopc_fldList_GetMetadata

    ! input/output arguments
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(in) :: iam      ! local pet
    real(r8)        ,intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind
    integer         ,intent(out):: rc 

    ! local variables
    integer          :: rcode
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    logical          :: lwhead, lwdata
    integer          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_r8) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lwhead = .true.
    if (present(whead)) lwhead = whead
    lwdata = .true.
    if (present(wdata)) lwdata = wdata
    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    if (lwhead) then
       call shr_nuopc_fldList_GetMetadata(trim(dname),longname=lname,stdname=sname,units=cunit)
       !       rcode = pio_def_dim(io_file(lfile_ind),trim(dname)//'_nx',1,dimid(1))
       !       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_DOUBLE,varid)
       if(rcode==PIO_NOERR) then
          rcode = pio_put_att(io_file(lfile_ind),varid,"units",trim(cunit))
          rcode = pio_put_att(io_file(lfile_ind),varid,"long_name",trim(lname))
          rcode = pio_put_att(io_file(lfile_ind),varid,"standard_name",trim(sname))
          if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
       end if
    endif

    if (lwdata) then
       rcode = pio_inq_varid(io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,rdata)
    endif

  end subroutine med_io_write_r8

  !===============================================================================
  subroutine med_io_write_r81d(filename, iam, rdata, dname, whead, wdata, file_ind, rc)

    !---------------
    ! Write 1d double array to netcdf file
    !---------------

    use pio               , only : var_desc_t, pio_def_dim, pio_def_var
    use pio               , only : pio_inq_varid, pio_put_var, pio_double, pio_put_att
    use esmFlds           , only : shr_nuopc_fldList_GetMetadata

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(in) :: iam
    real(r8)        ,intent(in) :: rdata(:) ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind
    integer         ,intent(out):: rc 

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer          :: lnx
    logical          :: lwhead, lwdata
    integer          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_r81d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lwhead = .true.
    if (present(whead)) lwhead = whead
    lwdata = .true.
    if (present(wdata)) lwdata = wdata
    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind

    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    if (lwhead) then
       call shr_nuopc_fldList_GetMetadata(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = size(rdata)
       rcode = pio_def_dim(io_file(lfile_ind),trim(dname)//'_nx',lnx,dimid(1))
       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
    endif

    if (lwdata) then
       rcode = pio_inq_varid(io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,rdata)
    endif

  end subroutine med_io_write_r81d

  !===============================================================================
  subroutine med_io_write_char(filename, iam, rdata, dname, whead, wdata, file_ind, rc)

    !---------------
    ! Write char string to netcdf file
    !---------------

    use pio     , only : var_desc_t, pio_def_dim, pio_put_att, pio_def_var, pio_inq_varid
    use pio     , only : pio_char, pio_put_var
    use esmFlds , only : shr_nuopc_fldList_GetMetadata

    ! input/output arguments
    character(len=*),intent(in) :: filename ! file
    integer         ,intent(in) :: iam      ! local pet
    character(len=*),intent(in) :: rdata    ! data to be written
    character(len=*),intent(in) :: dname    ! name of data
    logical,optional,intent(in) :: whead    ! write header
    logical,optional,intent(in) :: wdata    ! write data
    integer,optional,intent(in) :: file_ind
    integer         ,intent(out):: rc 

    ! local variables
    integer          :: rcode
    integer          :: dimid(1)
    type(var_desc_t) :: varid
    character(CL)    :: cunit       ! var units
    character(CL)    :: lname       ! long name
    character(CL)    :: sname       ! standard name
    integer          :: lnx
    logical          :: lwhead, lwdata
    integer          :: lfile_ind
    character(CL)    :: charvar   ! buffer for string read/write
    character(*),parameter :: subName = '(med_io_write_char) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lwhead = .true.
    if (present(whead)) lwhead = whead
    lwdata = .true.
    if (present(wdata)) lwdata = wdata
    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind
    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    if (lwhead) then
       call shr_nuopc_fldList_GetMetadata(trim(dname),longname=lname,stdname=sname,units=cunit)
       lnx = len(charvar)
       rcode = pio_def_dim(io_file(lfile_ind),trim(dname)//'_len',lnx,dimid(1))
       rcode = pio_def_var(io_file(lfile_ind),trim(dname),PIO_CHAR,dimid,varid)
       rcode = pio_put_att(io_file(lfile_ind),varid,"units",trim(cunit))
       rcode = pio_put_att(io_file(lfile_ind),varid,"long_name",trim(lname))
       rcode = pio_put_att(io_file(lfile_ind),varid,"standard_name",trim(sname))
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
    endif
    if (lwdata) then
       charvar = ''
       charvar = trim(rdata)
       rcode = pio_inq_varid(io_file(lfile_ind),trim(dname),varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,charvar)
    endif

  end subroutine med_io_write_char

  !===============================================================================
  subroutine med_io_write_time(filename, iam, time_units, time_cal, time_val, nt,&
       whead, wdata, tbnds, file_ind, rc)

    !---------------
    ! Write time variable to netcdf file
    !---------------

    use shr_cal_mod       , only : shr_cal_calMaxLen
    use shr_cal_mod       , only : shr_cal_noleap
    use shr_cal_mod       , only : shr_cal_gregorian
    use shr_cal_mod       , only : shr_cal_calendarName
    use pio               , only : var_desc_t, PIO_UNLIMITED
    use pio               , only : pio_double, pio_def_dim, pio_def_var, pio_put_att
    use pio               , only : pio_inq_varid, pio_put_var

    ! input/output variables
    character(len=*),      intent(in) :: filename   ! file
    integer,               intent(in) :: iam        ! local pet
    character(len=*),      intent(in) :: time_units ! units of time
    character(len=*),      intent(in) :: time_cal   ! calendar type
    real(r8)        ,      intent(in) :: time_val   ! data to be written
    integer    , optional, intent(in) :: nt
    logical,     optional, intent(in) :: whead      ! write header
    logical,     optional, intent(in) :: wdata      ! write data
    real(r8),    optional, intent(in) :: tbnds(2)   ! time bounds
    integer,     optional, intent(in) :: file_ind
    integer              , intent(out):: rc 

    ! local variables
    integer                          :: rcode
    integer                          :: dimid(1)
    integer                          :: dimid2(2)
    type(var_desc_t)                 :: varid
    logical                          :: lwhead, lwdata
    integer                          :: start(4),count(4)
    character(len=shr_cal_calMaxLen) :: lcalendar
    real(r8)                         :: time_val_1d(1)
    integer                          :: lfile_ind
    character(*),parameter :: subName = '(med_io_write_time) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lwhead = .true.
    if (present(whead)) lwhead = whead
    lwdata = .true.
    if (present(wdata)) lwdata = wdata
    lfile_ind = 0
    if (present(file_ind)) lfile_ind=file_ind
    if (.not.lwhead .and. .not.lwdata) then
       ! should we write a warning?
       return
    endif

    ! Write out header
    if (lwhead) then
       rcode = pio_def_dim(io_file(lfile_ind),'time',PIO_UNLIMITED,dimid(1))
       rcode = pio_def_var(io_file(lfile_ind),'time',PIO_DOUBLE,dimid,varid)
       rcode = pio_put_att(io_file(lfile_ind),varid,'units',trim(time_units))

       lcalendar = shr_cal_calendarName(time_cal,trap=.false.)
       if (trim(lcalendar) == trim(shr_cal_noleap)) then
          lcalendar = 'noleap'
       elseif (trim(lcalendar) == trim(shr_cal_gregorian)) then
          lcalendar = 'gregorian'
       endif
       rcode = pio_put_att(io_file(lfile_ind),varid,'calendar',trim(lcalendar))

       if (present(tbnds)) then
          dimid2(2) = dimid(1)
          rcode = pio_put_att(io_file(lfile_ind),varid,'bounds','time_bnds')
          rcode = pio_def_dim(io_file(lfile_ind),'ntb',2,dimid2(1))
          rcode = pio_def_var(io_file(lfile_ind),'time_bnds',PIO_DOUBLE,dimid2,varid)
       endif
       if (lwdata) call med_io_enddef(filename, file_ind=lfile_ind)
    endif

    ! Write out data
    if (lwdata) then
       start = 1
       count = 1
       if (present(nt)) then
          start(1) = nt
       endif
       time_val_1d(1) = time_val
       rcode = pio_inq_varid(io_file(lfile_ind),'time',varid)
       rcode = pio_put_var(io_file(lfile_ind),varid,start,count,time_val_1d)
       if (present(tbnds)) then
          rcode = pio_inq_varid(io_file(lfile_ind),'time_bnds',varid)
          start = 1
          count = 1
          if (present(nt)) then
             start(2) = nt
          endif
          count(1) = 2
          rcode = pio_put_var(io_file(lfile_ind),varid,start,count,tbnds)
       endif
    endif

  end subroutine med_io_write_time

  !===============================================================================
  subroutine med_io_read_FB(filename, vm, iam, FB, pre, frame, rc)

    !---------------
    ! Read FB from netcdf file
    !---------------

    use shr_const_mod         , only : fillvalue=>SHR_CONST_SPVAL
    use ESMF                  , only : ESMF_FieldBundle, ESMF_Field, ESMF_Mesh, ESMF_DistGrid
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_FAILURE
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use pio                   , only : file_desc_T, var_desc_t, io_desc_t, pio_nowrite, pio_openfile
    use pio                   , only : pio_noerr, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio                   , only : pio_inq_varid
    use pio                   , only : pio_double, pio_get_att, pio_seterrorhandling, pio_freedecomp, pio_closefile
    use pio                   , only : pio_read_darray, pio_offset_kind, pio_setframe
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag

    ! input/output arguments
    character(len=*)                        ,intent(in)  :: filename ! file
    type(ESMF_VM)                           ,intent(in)  :: vm
    integer                                 ,intent(in)  :: iam
    type(ESMF_FieldBundle)                  ,intent(in)  :: FB       ! data to be read
    character(len=*)              ,optional ,intent(in)  :: pre      ! prefix to variable name
    integer(kind=PIO_OFFSET_KIND) ,optional ,intent(in)  :: frame
    integer                                 ,intent(out) :: rc

    ! local variables
    type(ESMF_Field)              :: lfield
    integer                       :: rcode
    integer                       :: nf,ns,ng
    integer                       :: k,n,l
    type(file_desc_t)             :: pioid
    type(var_desc_t)              :: varid
    type(io_desc_t)               :: iodesc
    character(CL)                 :: itemc       ! string converted to char
    character(CL)                 :: name1       ! var name
    character(CL)                 :: lpre        ! local prefix
    real(r8)                      :: lfillvalue
    integer                       :: tmp(1)
    integer                       :: rank, lsize
    real(r8), pointer             :: fldptr1(:), fldptr1_tmp(:)
    real(r8), pointer             :: fldptr2(:,:)
    character(CL)                 :: tmpstr
    character(len=16)             :: cnumber
    integer(kind=Pio_Offset_Kind) :: lframe
    integer                       :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fieldds
    integer                       :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fieldds
    character(*),parameter :: subName = '(med_io_read_FB) '
    !-------------------------------------------------------------------------------
    rc = ESMF_Success
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lpre = ' '
    if (present(pre)) then
       lpre = trim(pre)
    endif
    if (present(frame)) then
       lframe = frame
    else
       lframe = 1
    endif
    if (.not. ESMF_FieldBundleIsCreated(FB,rc=rc)) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" not created", ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
       return
    endif

    call ESMF_FieldBundleGet(FB, fieldCount=nf, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(tmpstr,*) subname//' field count = '//trim(lpre),nf
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (nf < 1) then
       call ESMF_LogWrite(trim(subname)//" FB "//trim(lpre)//" empty", ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
       return
    endif

    if (med_io_file_exists(vm, iam, trim(filename))) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call ESMF_LogWrite(trim(subname)//' open file '//trim(filename), ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//' ERROR: file invalid '//trim(filename), &
       ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    endif

    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)

    do k = 1,nf
       ! Get name of field
       call FB_getNameN(FB, k, itemc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get iodesc for all fields based on iodesc of first field (assumes that all fields have
       ! the same iodesc)
       if (k == 1) then
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, rank=rank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (rank == 2) then
             name1 = trim(lpre)//'_'//trim(itemc)//'1'
          else if (rank == 1) then
             name1 = trim(lpre)//'_'//trim(itemc)
          end if
          call med_io_read_init_iodesc(FB, name1, pioid, iodesc, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       call ESMF_LogWrite(trim(subname)//' reading field '//trim(itemc), ESMF_LOGMSG_INFO)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Get pointer to field bundle field
       ! Field bundle might be 2d or 1d - but field on mediator history or restart file will always be 1d
       call FB_getFldPtr(FB, itemc, &
            fldptr1=fldptr1, fldptr2=fldptr2, rank=rank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (rank == 2) then

          ! Determine the size of the ungridded dimension and the
          ! index where the undistributed dimension is located
          call ESMF_FieldBundleGet(FB, itemc,  field=lfield, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (gridToFieldMap(1) == 1) then
             lsize = size(fldptr2, dim=1)
          else if (gridToFieldMap(1) == 2) then
             lsize = size(fldptr2, dim=2)
          end if
          allocate(fldptr1_tmp(lsize))

          do n = 1,ungriddedUBound(1)
             ! Creat a name for the 1d field on the mediator history or restart file based on the
             ! ungridded dimension index of the field bundle 2d fiedl
             write(cnumber,'(i0)') n
             name1 = trim(lpre)//'_'//trim(itemc)//trim(cnumber)

             rcode = pio_inq_varid(pioid, trim(name1), varid)
             if (rcode == pio_noerr) then
                call ESMF_LogWrite(trim(subname)//' read field '//trim(name1), ESMF_LOGMSG_INFO)
                if (chkerr(rc,__LINE__,u_FILE_u)) return
                call pio_setframe(pioid, varid, lframe)
                call pio_read_darray(pioid, varid, iodesc, fldptr1_tmp, rcode)
                rcode = pio_get_att(pioid, varid, "_FillValue", lfillvalue)
                if (rcode /= pio_noerr) then
                   lfillvalue = fillvalue
                endif
                do l = 1,size(fldptr1_tmp)
                   if (fldptr1_tmp(l) == lfillvalue) fldptr1_tmp(l) = 0.0_r8
                enddo
             else
                fldptr1_tmp = 0.0_r8
             endif
             if (gridToFieldMap(1) == 1) then
                fldptr2(:,n) = fldptr1_tmp(:)
             else if (gridToFieldMap(1) == 2) then
                fldptr2(n,:) = fldptr1_tmp(:)
             end if
          end do

          deallocate(fldptr1_tmp)

       else if (rank == 1) then
          name1 = trim(lpre)//'_'//trim(itemc)

          rcode = pio_inq_varid(pioid, trim(name1), varid)
          if (rcode == pio_noerr) then
             call ESMF_LogWrite(trim(subname)//' read field '//trim(name1), ESMF_LOGMSG_INFO)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call pio_setframe(pioid,varid,lframe)
             call pio_read_darray(pioid, varid, iodesc, fldptr1, rcode)
             rcode = pio_get_att(pioid,varid,"_FillValue",lfillvalue)
             if (rcode /= pio_noerr) then
                lfillvalue = fillvalue
             endif
             do n = 1,size(fldptr1)
                if (fldptr1(n) == lfillvalue) fldptr1(n) = 0.0_r8
             enddo
          else
             fldptr1 = 0.0_r8
          endif
       end if

    enddo ! end of loop over fields
    call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)

    call pio_freedecomp(pioid, iodesc)
    call pio_closefile(pioid)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine med_io_read_FB

  !===============================================================================
  subroutine med_io_read_init_iodesc(FB, name1, pioid, iodesc, rc)

    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF , only : ESMF_FieldBundleIsCreated, ESMF_FieldBundle, ESMF_Mesh, ESMF_DistGrid
    use ESMF , only : ESMF_FieldBundleGet, ESMF_FieldGet, ESMF_MeshGet, ESMF_DistGridGet
    use ESMF , only : ESMF_Field, ESMF_FieldGet, ESMF_AttributeGet 
    use pio  , only : file_desc_T, var_desc_t, io_desc_t, pio_nowrite, pio_openfile
    use pio  , only : pio_noerr, pio_inq_varndims
    use pio  , only : pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, pio_inq_vardimid
    use pio  , only : pio_double, pio_seterrorhandling, pio_initdecomp

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)    :: FB
    character(len=*)       , intent(in)    :: name1
    type(file_desc_t)      , intent(in)    :: pioid
    type(io_desc_t)        , intent(inout) :: iodesc
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)    :: field
    type(ESMF_Mesh)     :: mesh
    type(ESMF_Distgrid) :: distgrid
    integer             :: rcode
    integer             :: ns,ng
    integer             :: n,ndims
    integer, pointer    :: dimid(:)
    type(var_desc_t)    :: varid
    integer             :: lnx,lny
    integer             :: tmp(1)
    integer, pointer    :: minIndexPTile(:,:)
    integer, pointer    :: maxIndexPTile(:,:)
    integer             :: dimCount, tileCount
    integer, pointer    :: Dof(:)
    character(CL)       :: tmpstr
    integer             :: rank
    character(*),parameter :: subName = '(med_io_read_init_iodesc) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    rcode = pio_inq_varid(pioid, trim(name1), varid)
    if (rcode == pio_noerr) then

       rcode = pio_inq_varndims(pioid, varid, ndims)
       write(tmpstr,*) trim(subname),' ndims = ',ndims
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       allocate(dimid(ndims))
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       rcode = pio_inq_dimlen(pioid, dimid(1), lnx)
       write(tmpstr,*) trim(subname),' lnx = ',lnx
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       if (ndims>=2) then
          rcode = pio_inq_dimlen(pioid, dimid(2), lny)
       else
          lny = 1
       end if
       deallocate(dimid)

       write(tmpstr,*) trim(subname),' lny = ',lny
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
       ng = lnx * lny

       call FB_getFieldN(FB, 1, field, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(field, mesh=mesh, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(mesh, elementDistgrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount))
       call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       !write(tmpstr,*) subname,' counts = ',dimcount,tilecount,minindexptile,maxindexptile
       !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       if (ng > maxval(maxIndexPTile)) then
          write(tmpstr,*) subname,' WARNING: dimensions do not match', lnx, lny, maxval(maxIndexPTile)
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)
          !TODO: this should not be an error for say CTSM which does not send a global grid
          !rc = ESMF_Failure
          !return
       endif

       call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(dof(ns))
       call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
       write(tmpstr,*) subname,' dof = ',ns,size(dof),dof(1),dof(ns)  !,minval(dof),maxval(dof)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO)

       call pio_initdecomp(io_subsystem, pio_double, (/lnx,lny/), dof, iodesc)
       deallocate(dof)

       deallocate(minIndexPTile, maxIndexPTile)

    end if ! end if rcode check

  end subroutine med_io_read_init_iodesc

  !===============================================================================
  subroutine med_io_read_int(filename, vm, iam, idata, dname, rc)

    !---------------
    ! Read scalar integer from netcdf file
    !---------------

    ! input/output arguments
    character(len=*) , intent(in)    :: filename ! file
    type(ESMF_VM)                    :: vm
    integer          , intent(in)    :: iam
    integer          , intent(inout) :: idata    ! integer data
    character(len=*) , intent(in)    :: dname    ! name of data
    integer          , intent(out)   :: rc

    ! local variables
    integer :: i1d(1)
    character(*),parameter :: subName = '(med_io_read_int) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call med_io_read_int1d(filename, vm, iam, i1d, dname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    idata = i1d(1)

  end subroutine med_io_read_int

  !===============================================================================
  subroutine med_io_read_int1d(filename, vm, iam, idata, dname, rc)

    !---------------
    ! Read 1d integer array from netcdf file
    !---------------

    use med_constants_mod     , only : R8
    use pio                   , only : var_desc_t, file_desc_t, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, pio_seterrorhandling
    use pio                   , only : pio_get_var, pio_inq_varid, pio_get_att, pio_openfile
    use pio                   , only : pio_nowrite, pio_openfile, pio_global
    use pio                   , only : pio_closefile
    use med_internalstate_mod , only : logunit

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    integer,          intent(in)    :: iam
    integer         , intent(inout) :: idata(:) ! integer data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_t) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    character(*),parameter :: subName = '(med_io_read_int1d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    if (med_io_file_exists(vm, iam, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,idata)

    call pio_closefile(pioid)
  end subroutine med_io_read_int1d

  !===============================================================================
  subroutine med_io_read_r8(filename, vm, iam, rdata, dname, rc)
    use med_constants_mod, only : R8

    !---------------
    ! Read scalar double from netcdf file
    !---------------

    ! input/output arguments
    character(len=*) , intent(in)    :: filename ! file
    type(ESMF_VM)                    :: vm
    integer          , intent(in)    :: iam
    real(r8)         , intent(inout) :: rdata    ! real data
    character(len=*) , intent(in)    :: dname    ! name of data
    integer          , intent(out)   :: rc 

    ! local variables
    real(r8) :: r1d(1)
    character(*),parameter :: subName = '(med_io_read_r8) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call med_io_read_r81d(filename, vm, iam, r1d,dname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    rdata = r1d(1)

  end subroutine med_io_read_r8

  !===============================================================================
  subroutine med_io_read_r81d(filename, vm, iam, rdata, dname, rc)

    !---------------
    ! Read 1d double array from netcdf file
    !---------------

    use med_constants_mod     , only : R8
    use pio                   , only : file_desc_t, var_desc_t, pio_openfile, pio_closefile, pio_seterrorhandling
    use pio                   , only : PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, pio_inq_varid, pio_get_var
    use pio                   , only : pio_nowrite, pio_openfile, pio_global, pio_get_att
    use med_internalstate_mod , only : logunit

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    integer         , intent(in)    :: iam
    real(r8)        , intent(inout) :: rdata(:) ! real data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_T) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    character(*),parameter :: subName = '(med_io_read_r81d) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    if (med_io_file_exists(vm, iam, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    endif

    if (trim(lversion) == trim(version)) then
       name1 = trim(dname)
    else
       name1 = trim(prefix)//trim(dname)
    endif
    rcode = pio_inq_varid(pioid,trim(name1),varid)
    rcode = pio_get_var(pioid,varid,rdata)

    call pio_closefile(pioid)
  end subroutine med_io_read_r81d

  !===============================================================================
  subroutine med_io_read_char(filename, vm, iam, rdata, dname, rc)

    !---------------
    ! Read char string from netcdf file
    !---------------

    use pio                   , only : file_desc_t, var_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR
    use pio                   , only : pio_closefile, pio_inq_varid, pio_get_var
    use pio                   , only : pio_openfile, pio_global, pio_get_att, pio_nowrite
    use med_internalstate_mod , only : logunit

    ! input/output arguments
    character(len=*), intent(in)    :: filename ! file
    type(ESMF_VM)                   :: vm
    integer, intent(in)             :: iam
    character(len=*), intent(inout) :: rdata    ! character data
    character(len=*), intent(in)    :: dname    ! name of data
    integer         , intent(out)   :: rc

    ! local variables
    integer           :: rcode
    type(file_desc_T) :: pioid
    type(var_desc_t)  :: varid
    character(CL)     :: lversion
    character(CL)     :: name1
    character(CL)     :: charvar   ! buffer for string read/write
    character(*),parameter :: subName = '(med_io_read_char) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lversion=trim(version)

    if (med_io_file_exists(vm, iam, filename)) then
       rcode = pio_openfile(io_subsystem, pioid, pio_iotype, trim(filename),pio_nowrite)
       ! write(logunit,*) subname,' open file ',trim(filename)
       call pio_seterrorhandling(pioid,PIO_BCAST_ERROR)
       rcode = pio_get_att(pioid,pio_global,"file_version",lversion)
       call pio_seterrorhandling(pioid,PIO_INTERNAL_ERROR)
    else
       if(iam==0) write(logunit,*) subname,' ERROR: file invalid ',trim(filename),' ',trim(dname)
       call ESMF_LogWrite(trim(subname)//'ERROR: file invalid '//trim(filename)//' '//trim(dname), ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
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
  end subroutine med_io_read_char

end module med_io_mod
