module dshr_nuopc_mod

  use NUOPC
  use NUOPC_Model
  use ESMF
  use mct_mod          , only : mct_avect_init, mct_avect_lsize, mct_avect_indexra
  use dshr_methods_mod , only : chkerr, state_getfldptr, get_component_instance
  use shr_strdata_mod  , only : shr_strdata_type
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_string_mod   , only : shr_string_listGetIndex
  use shr_sys_mod      , only : shr_sys_abort

  implicit none
  public

  public :: dshr_advertise
  public :: dshr_model_initphase
  public :: dshr_check_mesh
  public :: dshr_create_mesh_from_grid
  public :: dshr_set_runclock
  public :: dshr_sdat_init
  public :: dshr_restart_read
  public :: dshr_restart_write
  public :: dshr_get_atm_adjustment_factors
  public :: dshr_get_griddata
  public :: dshr_set_griddata

  private :: dshr_alarm_init
  private :: dshr_time_init

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"      , &
       optNever          = "never"     , &
       optNSteps         = "nsteps"    , &
       optNStep          = "nstep"     , &
       optNSeconds       = "nseconds"  , &
       optNSecond        = "nsecond"   , &
       optNMinutes       = "nminutes"  , &
       optNMinute        = "nminute"   , &
       optNHours         = "nhours"    , &
       optNHour          = "nhour"     , &
       optNDays          = "ndays"     , &
       optNDay           = "nday"      , &
       optNMonths        = "nmonths"   , &
       optNMonth         = "nmonth"    , &
       optNYears         = "nyears"    , &
       optNYear          = "nyear"     , &
       optMonthly        = "monthly"   , &
       optYearly         = "yearly"    , &
       optDate           = "date"      , &
       optIfdays0        = "ifdays0"   


  ! Note that gridTofieldMap = 2, therefore the ungridded dimension is innermost

  integer                 :: iunset = -999
  integer     , parameter :: SecPerDay = 86400 ! Seconds per day
  integer     , parameter :: dbug = 10
  character(*), parameter :: modName =  "(dhsr_nuopc_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
       flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! input/output variables
    type(ESMF_GridComp)              :: gcomp
    integer          , intent(inout) :: mpicom
    integer          , intent(out)   :: my_task
    integer          , intent(out)   :: inst_index
    character(len=*) , intent(out)   :: inst_suffix
    character(len=*) , intent(out)   :: flds_scalar_name
    integer          , intent(out)   :: flds_scalar_num
    integer          , intent(out)   :: flds_scalar_index_nx
    integer          , intent(out)   :: flds_scalar_index_ny
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    character(len=*),parameter  :: subname='(dshr_advertise)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! generate local mpi comm
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=my_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine instance information
    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get scalar attributes
    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

  end subroutine dshr_advertise

!===============================================================================

  subroutine dshr_create_mesh_from_grid(filename, mesh, rc)

    use netcdf , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close, nf90_strerror
    use netcdf , only : nf90_inq_dimid, nf90_inq_varid, nf90_get_var
    use netcdf , only : nf90_inquire_dimension, nf90_inquire_variable

    ! input/output variables
    character(len=*), intent(in)  :: filename
    type(ESMF_Mesh) , intent(out) :: mesh
    integer         , intent(out) :: rc

    ! local variables
    integer               :: ncid, ierr
    integer               :: dimid_ni, dimid_nj, dimid_nv
    integer               :: ni, nj, nv
    integer               :: varid_xv, varid_yv
    integer               :: maxIndex(2)
    real(r8)              :: mincornerCoord(2)
    real(r8)              :: maxcornerCoord(2)
    type(ESMF_Grid)       :: lgrid
    real(r8), allocatable :: xv(:,:,:), yv(:,:,:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! open file
    ierr = nf90_open(filename, NF90_NOWRITE, ncid)
    call nc_check_err(ierr, 'nf90_open', trim(filename))

    ! get dimension ids
    ierr = nf90_inq_dimid(ncid, 'ni', dimid_ni)
    call nc_check_err(ierr, 'nf90_inq_dimid for ni', trim(filename))
    ierr = nf90_inq_dimid(ncid, 'nj', dimid_nj)
    call nc_check_err(ierr, 'nf90_inq_dimid for nj', trim(filename))
    ierr = nf90_inq_dimid(ncid, 'nv', dimid_nv)
    call nc_check_err(ierr, 'nf90_inq_dimid for nv', trim(filename))

    ! get dimension values
    ierr = nf90_inquire_dimension(ncid, dimid_ni, len=ni)
    call nc_check_err(ierr, 'nf90_inq_dimension for ni', trim(filename))
    ierr = nf90_inquire_dimension(ncid, dimid_nj, len=nj)
    call nc_check_err(ierr, 'nf90_inq_dimension for nj', trim(filename))
    ierr = nf90_inquire_dimension(ncid, dimid_nv, len=nv)
    call nc_check_err(ierr, 'nf90_inq_dimension for nv', trim(filename))

    ! get variable ids
    ierr = nf90_inq_varid(ncid, 'xv', varid_xv)
    call nc_check_err(ierr, 'nf90_inq_varid for xv', trim(filename))
    ierr = nf90_inq_varid(ncid, 'yv', varid_yv)
    call nc_check_err(ierr, 'nf90_inq_varid for yv', trim(filename))

    ! allocate memory for variables and get variable values
    allocate(xv(nv,ni,nj), yv(nv,ni,nj))
    ierr = nf90_get_var(ncid, varid_xv, xv)
    call nc_check_err(ierr, 'nf90_get_var for xv', trim(filename))
    ierr = nf90_get_var(ncid, varid_yv, yv)
    call nc_check_err(ierr, 'nf90_get_var for yv', trim(filename))

    ! close file
    ierr = nf90_close(ncid)
    call nc_check_err(ierr, 'nf90_close', trim(filename))

    ! create the grid
    maxIndex(1)       = ni          ! number of lons
    maxIndex(2)       = nj          ! number of lats
    mincornerCoord(1) = xv(1,1,1)   ! min lon
    mincornerCoord(2) = yv(1,1,1)   ! min lat
    maxcornerCoord(1) = xv(3,ni,nj) ! max lon
    maxcornerCoord(2) = yv(3,ni,nj) ! max lat
    deallocate(xv,yv)
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the grid
    mesh =  ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  contains

    subroutine nc_check_err(ierror, description, filename)
      integer     , intent(in) :: ierror
      character(*), intent(in) :: description
      character(*), intent(in) :: filename

      if (ierror /= nf90_noerr) then
         write (*,'(6a)') 'ERROR ', trim(description),'. NetCDF file : "', trim(filename),&
              '". Error message:', trim(nf90_strerror(ierror))
         call shr_sys_abort()
      endif
    end subroutine nc_check_err

  end subroutine dshr_create_mesh_from_grid

!===============================================================================

  subroutine dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
       scmmode, scmlon, scmlat, clock, mesh,  model_name, sdat, &
       dmodel_domain_fracname_from_stream, reset_domain_mask, rc)

    ! ----------------------------------------------
    ! Initialize sdat
    ! ----------------------------------------------

    use shr_strdata_mod , only : shr_strdata_pioinit
    use shr_strdata_mod , only : shr_strdata_init_model_domain
    use shr_strdata_mod , only : shr_strdata_init_streams
    use shr_strdata_mod , only : shr_strdata_init_mapping
    use shr_strdata_mod , only : shr_strdata_print
    use shr_cal_mod     , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_calendarname
    use mct_mod         , only : mct_gsmap_init

    ! input/output variables
    integer                    , intent(in)    :: mpicom   ! mpi communicator
    integer                    , intent(in)    :: compid
    integer                    , intent(in)    :: my_task
    integer                    , intent(in)    :: master_task
    integer                    , intent(in)    :: logunit
    logical                    , intent(in)    :: scmMode
    real(r8)                   , intent(in)    :: scmlat
    real(r8)                   , intent(in)    :: scmlon
    type(ESMF_Clock)           , intent(in)    :: clock
    type(ESMF_Mesh)            , intent(in)    :: mesh
    character(len=*)           , intent(in)    :: model_name
    type(shr_strdata_type)     , intent(inout) :: sdat
    character(len=*), optional , intent(in)    :: dmodel_domain_fracname_from_stream
    logical         , optional , intent(in)    :: reset_domain_mask
    integer                    , intent(out)   :: rc

    ! local varaibles
    integer                      :: n,k          ! generic counters
    integer                      :: lsize        ! local size
    integer                      :: gsize
    type(ESMF_DistGrid)          :: distGrid
    integer, allocatable, target :: gindex(:)
    integer                      :: dimCount
    integer                      :: tileCount
    integer                      :: deCount
    integer, allocatable         :: elementCountPTile(:)
    integer, allocatable         :: indexCountPDE(:,:)
    type(ESMF_CalKind_Flag)      :: esmf_caltype ! esmf calendar type
    character(len=CS)            :: calendar     ! calendar name
    character(len=*), parameter  :: subname='(dshr_nuopc_mod:dshr_sdat_init)'
    character(*)    , parameter  :: F01="('(dshr_init_strdata) ',a,2f10.4)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! obtain the distgrid from the mesh that was read in
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine local size on my processor
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global index space for my processor
    allocate(gindex(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine global size of distgrid
    call ESMF_DistGridGet(distGrid, dimCount=dimCount, deCount=deCount, tileCount=tileCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(elementCountPTile(tileCount))
    call ESMF_distGridGet(distGrid, elementCountPTile=elementCountPTile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    gsize = 0
    do n = 1,size(elementCountPTile)
       gsize = gsize + elementCountPTile(n)
    end do
    deallocate(elementCountPTile)

    ! initialize sdat%gsmap (the data mdel gsmap)
    call mct_gsMap_init(sdat%gsmap, gindex, mpicom, compid, lsize, gsize)
    deallocate(gindex)

    ! initialize shr_strdata pio
    call shr_strdata_pioinit(sdat, compid)

    ! initialize sdat calendar
    call ESMF_ClockGet(clock, calkindflag=esmf_caltype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort(subname//" ERROR bad ESMF calendar name "//trim(calendar))
    end if
    sdat%calendar = trim(shr_cal_calendarName(trim(calendar)))

    ! initialize sdat domain (sdat%grid)
    if (my_task == master_task) then
       write(logunit,*) ' scm mode, lon lat = ',scmmode, scmlon,scmlat
       if (present(reset_domain_mask)) then
          write(logunit,*) ' resetting domain mask'
       end if
       if (present(dmodel_domain_fracname_from_stream)) then
          write(logunit,*)' reading fracname ',trim(dmodel_domain_fracname_from_stream),&
               ' from the domain of the first stream'
       end if
    end if
    call shr_strdata_init_model_domain(sdat, mpicom, compid, my_task, &
         scmmode=scmmode, scmlon=scmlon, scmlat=scmlat, gsmap=sdat%gsmap, &
         dmodel_domain_fracname_from_stream=dmodel_domain_fracname_from_stream, &
         reset_domain_mask=reset_domain_mask)

    ! initialize sdat attributes for streams and mapping of streams to model domain
    call shr_strdata_init_streams(sdat, compid, mpicom, my_task)
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

    if (my_task == master_task) then
       call shr_strdata_print(sdat,'SDAT data from '//trim(model_name))
    endif

  end subroutine dshr_sdat_init

!===============================================================================

  subroutine dshr_check_mesh (mesh, sdat, model_name, tolerance, check_lon, rc)

    type(ESMF_MESH)        , intent(in)  :: mesh
    type(shr_strdata_type) , intent(in)  :: sdat
    character(len=*)       , intent(in)  :: model_name
    real(r8), optional     , intent(in)  :: tolerance 
    logical , optional     , intent(in)  :: check_lon
    integer                , intent(out) :: rc

    ! local variables
    integer           :: n,lsize
    integer           :: klat, klon, kfrac  ! AV indices
    real(r8)          :: domlon,domlat      ! domain lats and lots
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(r8), pointer :: ownedElemCoords(:) ! mesh lat and lons
    real(r8), pointer :: xc(:), yc(:)       ! mesh lats and lons
    real(r8)          :: ltolerance         ! tolerance of check
    logical           :: lcheck_lon
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! set tolerance of comparison
    if (present(tolerance)) then
       ltolerance = tolerance
    else
       ltolerance = 1.e-10
    end if

    if (present(check_lon)) then
       lcheck_lon = check_lon
    else
       lcheck_lon = .false.
    end if

    ! obtain mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(xc(numOwnedElements), yc(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1, numOwnedElements
       xc(n) = ownedElemCoords(2*n-1)
       yc(n) = ownedElemCoords(2*n)
    end do

    ! obtain sdat lat and lons and local size of grid attribute vectors
    klon = mct_aVect_indexRA(sdat%grid%data,'lon')
    klat = mct_aVect_indexRA(sdat%grid%data,'lat')
    lsize = mct_avect_lsize(sdat%grid%data)

    ! error check
    if (numOwnedElements /= lsize) then
       write(6,*)'lsize, numownedElements= ',lsize,numOwnedElements 
       call shr_sys_abort('ERROR: numOwnedElements is not equal to lsize')
    end if

    do n = 1, lsize
       domlat = sdat%grid%data%rattr(klat,n)
       if (lcheck_lon) then
          domlon = sdat%grid%data%rattr(klon,n)
          if (abs( domlon - xc(n)) > ltolerance .and. domlon /= 0.0_r8) then
             write(6,100) 'ERROR: '//trim(model_name)//' n, dom_lon, mesh_lon, diff_lon = ',n, domlon, xc(n), abs(xc(n)-domlon)
             call shr_sys_abort()
          end if
       end if
       if (abs( domlat - yc(n)) > ltolerance .and. domlat /= 0.0_r8) then
          write(6,100) 'ERROR: '//trim(model_name)//' n, dom_lat, mesh_lat, diff_lat = ',n, domlat, yc(n), abs(yc(n)-domlat)
          call shr_sys_abort()
       end if
100    format(a,i6,2(f21.13,3x),d21.5)
       !SDAT%grid%data%rattr(klon,n) = xc(n)
       !SDAT%grid%data%rattr(klat,n) = yc(n)
    end do
    deallocate(xc, yc)

  end subroutine dshr_check_mesh

!===============================================================================

  subroutine dshr_model_initphase(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_model_initphase

!===============================================================================

  subroutine dshr_set_runclock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option       ! Restart option units
    integer                  :: restart_n            ! Number until restart interval
    integer                  :: restart_ymd          ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='dshr_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call dshr_alarm_init(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_set_runclock

!===============================================================================

  subroutine dshr_alarm_init( clock, alarm, option, &
       opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

    ! Setup an alarm in a clock
    ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm
    ! time.  If you send an arbitrary but proper ringtime from the
    ! past and the ring interval, the alarm will always go off on the
    ! next clock advance and this will cause serious problems.  Even
    ! if it makes sense to initialize an alarm with some reference
    ! time and the alarm interval, that reference time has to be
    ! advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to
    ! advance it properly based on the ring interval.

    ! input/output variables
    type(ESMF_Clock)            , intent(inout) :: clock     ! clock
    type(ESMF_Alarm)            , intent(inout) :: alarm     ! alarm
    character(len=*)            , intent(in)    :: option    ! alarm option
    integer          , optional , intent(in)    :: opt_n     ! alarm freq
    integer          , optional , intent(in)    :: opt_ymd   ! alarm ymd
    integer          , optional , intent(in)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)  , optional , intent(in)    :: RefTime   ! ref time
    character(len=*) , optional , intent(in)    :: alarmname ! alarm name
    integer                     , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal                ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    integer                 :: sec
    character(len=*), parameter :: subname = '(dshr_alarm_init): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) lalarmname = trim(alarmname)
    ltod = 0
    if (present(opt_tod)) ltod = opt_tod
    lymd = -1
    if (present(opt_ymd)) lymd = opt_ymd

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif

    ! Determine calendar
    call ESMF_ClockGet(clock, calendar=cal)

    ! Determine inputs for call to create alarm
    selectcase (trim(option))

    case (optNONE)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNever)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optDate)
       if (.not. present(opt_ymd)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       end if
       if (lymd < 0 .or. ltod < 0) then
          call shr_sys_abort(subname//trim(option)//'opt_ymd, opt_tod invalid')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_time_init(NextAlarm, lymd, cal, ltod, rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optIfdays0)
       if (.not. present(opt_ymd)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       end if
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0)  then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

   case (optNSteps)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNStep)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSeconds)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSecond)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinute)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHours)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHour)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDays)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDay)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonths)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonth)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNYears)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNYear)
       if (.not.present(opt_n)) then
          call shr_sys_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case default
       call shr_sys_abort(subname//'unknown option '//trim(option))

    end select

    ! --------------------------------------------------------------------------------
    ! --- AlarmInterval and NextAlarm should be set ---
    ! --------------------------------------------------------------------------------

    ! --- advance Next Alarm so it won't ring on first timestep for
    ! --- most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif

    alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, &
         ringInterval=AlarmInterval, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_alarm_init

!===============================================================================

  subroutine dshr_time_init( Time, ymd, cal, tod, rc)

    ! Create the ESMF_Time object corresponding to the given input time, 
    ! given in YMD (Year Month Day) and TOD (Time-of-day) format.
    ! Set the time by an integer as YYYYMMDD and integer seconds in the day

    ! input/output parameters:
    type(ESMF_Time)     , intent(inout) :: Time ! ESMF time
    integer             , intent(in)    :: ymd  ! year, month, day YYYYMMDD
    type(ESMF_Calendar) , intent(in)    :: cal  ! ESMF calendar
    integer             , intent(in)    :: tod  ! time of day in seconds
    integer             , intent(out)   :: rc

    ! local variables
    integer :: year, mon, day ! year, month, day as integers
    integer :: tdate          ! temporary date
    integer :: date           ! coded-date (yyyymmdd)
    character(len=*), parameter :: subname='(dshr_time_init)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if ( (ymd < 0) .or. (tod < 0) .or. (tod > SecPerDay) )then
       call shr_sys_abort( subname//'ERROR yymmdd is a negative number or time-of-day out of bounds' )
    end if

    tdate = abs(date)
    year = int(tdate/10000)
    if (date < 0) year = -year
    mon = int( mod(tdate,10000)/  100)
    day = mod(tdate,  100)

    call ESMF_TimeSet( Time, yy=year, mm=mon, dd=day, s=tod, calendar=cal, rc=rc )
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_time_init

  !===============================================================================

  subroutine dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
       logunit, my_task, master_task, mpicom, sdat, fld, fldname)

    use shr_pcdf_mod    , only : shr_pcdf_readwrite
    use shr_strdata_mod , only : shr_strdata_restRead, shr_strdata_type
    use shr_mpi_mod     , only : shr_mpi_bcast

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_file
    character(len=*)            , intent(inout) :: rest_file_strm
    character(len=*)            , intent(in)    :: rpfile
    character(len=*)            , intent(in)    :: inst_suffix
    character(len=*)            , intent(in)    :: nullstr
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: master_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    real(r8)         , optional , pointer       :: fld(:)
    character(len=*) , optional , intent(in)    :: fldname

    ! local variables
    integer :: nu
    logical :: exists  ! file existance
    character(*), parameter :: F00   = "('(dshr_restart_read) ',8a)"
    character(*), parameter :: subName = "(dshr_restart_read) "
    !-------------------------------------------------------------------------------

    if (trim(rest_file) == trim(nullstr) .and. trim(rest_file_strm) == trim(nullstr)) then
       if (my_task == master_task) then
          write(logunit,F00) ' restart filenames from rpointer'
          inquire(file=trim(rpfile)//trim(inst_suffix), exist=exists)
          if (.not.exists) then
             write(logunit, F00) ' ERROR: rpointer file does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
          endif
          open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
          read(nu, '(a)') rest_file
          read(nu, '(a)') rest_file_strm
          close(nu)
          inquire(file=trim(rest_file_strm), exist=exists)
       endif
       call shr_mpi_bcast(rest_file, mpicom, 'rest_file')
       call shr_mpi_bcast(rest_file_strm, mpicom, 'rest_file_strm')
    else
       ! use namelist already read
       if (my_task == master_task) then
          write(logunit, F00) ' restart filenames from namelist '
          inquire(file=trim(rest_file_strm), exist=exists)
       endif
    endif
    call shr_mpi_bcast(exists, mpicom, 'exists')
    if (exists) then
       if (my_task == master_task) write(logunit, F00) ' reading ', trim(rest_file_strm)
       if (present(fld) .and. present(fldname)) then
          call shr_pcdf_readwrite('read', sdat%pio_subsystem, sdat%io_type, trim(rest_file), &
               mpicom, sdat%gsmap, clobber=.true., rf1=fld, rf1n=trim(fldname), io_format=sdat%io_format)
       end if
       call shr_strdata_restRead(trim(rest_file_strm), sdat, mpicom)
    else
       if (my_task == master_task) write(logunit, F00) ' file not found, skipping ',trim(rest_file_strm)
    endif
  end subroutine dshr_restart_read

  !===============================================================================

  subroutine dshr_restart_write(rpfile, case_name, model_name, inst_suffix, ymd, tod, &
       logunit, mpicom, my_task, master_task, sdat, fld, fldname)

    use shr_pcdf_mod    , only : shr_pcdf_readwrite
    use shr_cal_mod     , only : shr_cal_datetod2string
    use shr_strdata_mod , only : shr_strdata_restWrite, shr_strdata_type

    ! input/output variables
    character(len=*)            , intent(in)    :: rpfile
    character(len=*)            , intent(in)    :: case_name
    character(len=*)            , intent(in)    :: model_name
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: ymd       ! model date
    integer                     , intent(in)    :: tod       ! model sec into model date
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: master_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    real(r8)         , optional , pointer       :: fld(:)
    character(len=*) , optional , intent(in)    :: fldname

    ! local variables
    character(len=CL) :: rest_file
    character(len=CL) :: rest_file_strm
    character(len=CS) :: date_str
    integer           :: nu
    !-------------------------------------------------------------------------------

     call shr_cal_datetod2string(date_str, ymd, tod)
     write(rest_file     ,"(7a)") trim(case_name),'.', trim(model_name),trim(inst_suffix),'.r.'  , trim(date_str),'.nc'
     write(rest_file_strm,"(7a)") trim(case_name),'.', trim(model_name),trim(inst_suffix),'.rs1.', trim(date_str),'.bin'
     if (my_task == master_task) then
        open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
        write(nu,'(a)') rest_file
        write(nu,'(a)') rest_file_strm
        close(nu)
        write(logunit,*)' (dshr_restart_write) writing ',trim(rest_file_strm), ymd, tod
     endif
     if (present(fld) .and. present(fldname)) then
        call shr_pcdf_readwrite('write', sdat%pio_subsystem, sdat%io_type,&
             trim(rest_file), mpicom, sdat%gsmap, clobber=.true., rf1=fld, rf1n=trim(fldname))
     end if
     call shr_strdata_restWrite(trim(rest_file_strm), sdat, mpicom, trim(case_name), 'SDAT strdata from '//trim(model_name))

  end subroutine dshr_restart_write

  !===============================================================================

  subroutine dshr_get_atm_adjustment_factors(fileName, windF, winddF, qsatF, &
       mpicom, compid, masterproc, logunit, sdat) 

    use shr_dmodel_mod , only : shr_dmodel_mapset
    use shr_map_mod    , only : shr_map_fs_remap, shr_map_fs_bilinear
    use shr_map_mod    , only : shr_map_fs_srcmask, shr_map_fs_scalar
    use shr_ncread_mod , only : shr_ncread_varExists, shr_ncread_varDimSizes, shr_ncread_field4dG
    use shr_strdata_mod, only : shr_strdata_type
    use shr_const_mod  , only : shr_const_spval
    use shr_mpi_mod    , only : shr_mpi_bcast
    use mct_mod        , only : mct_avect, mct_avect_scatter, mct_smat_avmult, mct_smatp_clean, mct_smatp
    use mct_mod        , only : mct_avect_importrattr, mct_avect_exportrattr, mct_avect_clean   
    use mct_mod        , only : mct_ggrid, mct_ggrid_init,  mct_ggrid_importrattr, mct_ggrid_gather, mct_ggrid_clean   
    use mct_mod        , only : mct_gsmap, mct_gsmap_init, mct_gsmap_lsize, mct_gsmap_clean  

    ! input/output variables
    character(*)           , intent(in)    :: fileName   ! file name string
    real(R8)               , intent(inout) :: windF(:)   ! wind adjustment factor
    real(R8)               , intent(inout) :: winddF(:)  ! wind adjustment factor
    real(R8)               , intent(inout) :: qsatF(:)   ! rel humidty adjustment factors
    integer                , intent(in)    :: mpicom     ! mpi comm
    integer                , intent(in)    :: compid     ! mct compid
    logical                , intent(in)    :: masterproc ! true if pe number is 0
    integer                , intent(in)    :: logunit    ! logging unit
    type(shr_strdata_type) , intent(in)    :: sdat

    !--- data that describes the local model domain ---
    integer          :: ni0,nj0     ! dimensions of global bundle0
    integer          :: i,j,n       ! generic indicies
    type(mct_ggrid)  :: ggridi      ! input file grid
    type(mct_ggrid)  :: ggridoG     ! output grid gathered
    type(mct_gsmap)  :: gsmapi      ! input file gsmap
    type(mct_sMatp)  :: smatp       ! sparse matrix weights
    type(mct_avect)  :: avi         ! input attr vect
    type(mct_avect)  :: avo         ! output attr vect
    integer          :: lsizei      ! local size of input
    integer          :: lsizeo      ! local size of output
    integer, pointer :: start(:)    ! start list
    integer, pointer :: length(:)   ! length list
    integer          :: gsizei      ! input global size
    integer          :: numel       ! number of elements in start list
    real(R8)         :: dadd        ! lon correction
    logical          :: domap       ! map or not
    integer          :: klon,klat   ! lon lat fld index

    !--- temp arrays for data input ---
    real(R8)     ,allocatable :: tempR4D(:,:,:,:)   ! 4D data array
    real(R8)     ,pointer     :: tempR1D(:)         ! 1D data array
    integer      ,allocatable :: tempI4D(:,:,:,:)   ! 4D data array
    character(*) ,parameter   :: subName =  '(dshr_atm_get_adjustment_factors) '
    character(*) ,parameter   :: F00    = "('(dshr_atm_get_adjustment_factors) ',4a) "
    character(*) ,parameter   :: F01    = "('(dshr_atm_get_adjustment_factors) ',a,2i5)"
    character(*) ,parameter   :: F02    = "('(dshr_atm_get_adjustment_factors) ',a,6e12.3)"

    !-------------------------------------------------------------------------------
    !   Note: gsmapi is all gridcells on root pe
    !-------------------------------------------------------------------------------

    ni0 = 0
    nj0 = 0
    allocate(start(1),length(1))
    start = 0
    length = 0
    numel = 0

    !----------------------------------------------------------------------------
    ! read in and map global correction factors
    !----------------------------------------------------------------------------

    ! verify necessary data is in input file

    if (masterproc) then
       if (      .not. shr_ncread_varExists(fileName ,'lat'       )  &
            .or. .not. shr_ncread_varExists(fileName ,'lon'       )  &
            .or. .not. shr_ncread_varExists(fileName ,'mask'      )  &
            .or. .not. shr_ncread_varExists(fileName ,'windFactor')  &
            .or. .not. shr_ncread_varExists(fileName ,'qsatFactor')  ) then
          write(logunit,F00) "ERROR: invalid correction factor data file"
          call shr_sys_abort(subName//"invalid correction factor data file")
       end if
       call shr_ncread_varDimSizes(fileName,"windFactor",ni0,nj0)
       start = 1
       length = ni0*nj0
       numel = 1
    endif
    call shr_mpi_bcast(ni0, mpicom, subname//' ni0')
    call shr_mpi_bcast(nj0, mpicom, subname//' nj0')
    gsizei = ni0*nj0

    !--- allocate datatypes for input data ---
    call mct_gsmap_init(gsmapi, start, length, 0, mpicom, compid, gsize=gsizei, numel=numel)
    deallocate(start, length)
    lsizei = mct_gsmap_lsize(gsmapi    , mpicom)
    lsizeo = mct_gsmap_lsize(sdat%gsmap, mpicom)
    call mct_gGrid_init(GGrid=gGridi, CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsizei )
    call mct_aVect_init(avi, rList="wind:windd:qsat", lsize=lsizei)
    avi%rAttr = shr_const_spval

    !--- gather output grid for map logic ---
    call mct_ggrid_gather(sdat%grid, ggridoG, sdat%gsmap, 0, mpicom)

    if (masterproc) then
       allocate(tempR1D(ni0*nj0))

       !--- read domain data: lon ---
       allocate(tempR4D(ni0,1,1,1))
       call shr_ncread_field4dG(fileName,'lon' ,rfld=tempR4D)
       !--- needs to be monotonically increasing, add 360 at wraparound+ ---
       dadd = 0.0_R8
       do i = 2,ni0
          if (tempR4D(i-1,1,1,1) > tempR4D(i,1,1,1)) dadd = 360.0_R8
          tempR4D(i,1,1,1) = tempR4D(i,1,1,1) + dadd
       enddo
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = tempR4D(i,1,1,1)
          end do
       end do
       deallocate(tempR4D)
       call mct_gGrid_importRattr(gGridi,'lon',tempR1D,lsizei)

       !--- read domain data: lat ---
       allocate(tempR4D(nj0,1,1,1))
       call shr_ncread_field4dG(fileName,'lat' ,rfld=tempR4D)
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = tempR4D(j,1,1,1)
          end do
       end do
       deallocate(tempR4D)
       call mct_gGrid_importRattr(gGridi,'lat',tempR1D,lsizei)

       !--- read domain mask---
       allocate(tempI4D(ni0,nj0,1,1))
       call shr_ncread_field4dG(fileName,'mask',ifld=tempI4D)
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = real(tempI4D(i,j,1,1),R8)
          end do
       end do
       deallocate(tempI4D)
       call mct_gGrid_importRattr(gGridi,'mask',tempR1D,lsizei)

       !--- read bundle data: wind factor ---
       allocate(tempR4D(ni0,nj0,1,1))
       if (shr_ncread_varExists(fileName,'windFactor')  ) then
          call shr_ncread_field4dG(fileName,'windFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'wind',tempR1D,lsizei)
       endif

       !--- read bundle data: windd factor ---
       if (shr_ncread_varExists(fileName,'winddFactor')  ) then
          call shr_ncread_field4dG(fileName,'winddFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'windd',tempR1D,lsizei)
       endif

       !--- read bundle data: qsat factor ---
       if (shr_ncread_varExists(fileName,'qsatFactor')  ) then
          call shr_ncread_field4dG(fileName,'qsatFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'qsat',tempR1D,lsizei)
       endif

       deallocate(tempR4D)
       deallocate(tempR1D)

       domap = .false.
       if (ni0 /= sdat%nxg .or. nj0 /= sdat%nyg) then
          domap = .true.
       else
          klon = mct_aVect_indexRA(ggridi%data   ,'lon')
          klat = mct_aVect_indexRA(sdat%grid%data,'lat')
          do n = 1,lsizei
             if (abs(ggridi%data%rAttr(klon,n)-ggridoG%data%rAttr(klon,n)) > 0.01_R8) domap=.true.
             if (abs(ggridi%data%rAttr(klat,n)-ggridoG%data%rAttr(klat,n)) > 0.01_R8) domap=.true.
          enddo
       endif

       call mct_gGrid_clean(ggridoG)

    endif

    call shr_mpi_bcast(domap,mpicom,subname//' domap')

    if (domap) then
       call shr_dmodel_mapSet(smatp, ggridi, gsmapi, ni0 ,nj0, sdat%grid, sdat%gsmap, sdat%nxg, sdat%nyg, &            
            'datmfactor', shr_map_fs_remap, shr_map_fs_bilinear, shr_map_fs_srcmask, shr_map_fs_scalar, &
            compid, mpicom, 'Xonly')

       call mct_aVect_init(avo,avi,lsizeo)
       call mct_sMat_avMult(avi,smatp,avo)
       call mct_sMatP_clean(smatp)
    else
       call mct_aVect_scatter(avi, avo, sdat%gsmap, 0, mpicom)
    endif

    !--- fill the interface arrays, only if they are the right size ---
    allocate(tempR1D(lsizeo))
    if (size(windF ) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'wind' ,tempR1D,lsizeo)
       windF = tempR1D
    endif
    if (size(winddF) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'windd',tempR1D,lsizeo)
       winddF = tempR1D
    endif
    if (size(qsatF ) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'qsat' ,tempR1D,lsizeo)
       qsatF = tempR1D
    endif
    deallocate(tempR1D)

    call mct_aVect_clean(avi)
    call mct_aVect_clean(avo)
    call mct_gGrid_clean(ggridi)
    call mct_gsmap_clean(gsmapi)

  end subroutine dshr_get_atm_adjustment_factors

  !===============================================================================

  subroutine dshr_set_griddata(sdat, fldname, rvalue) 
    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , intent(in)    :: rvalue

    ! local variables
    integer :: kf

    kf = mct_aVect_indexRA(sdat%grid%data, trim(fldname))
    sdat%grid%data%rAttr(kf,:) = rvalue
  end subroutine dshr_set_griddata

  subroutine dshr_get_griddata(sdat, fldname, data) 
    ! input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    character(len=*)       , intent(in)    :: fldname
    real(r8)               , intent(out)   :: data(:)

    ! local variables
    integer :: kf

    kf = mct_aVect_indexRA(sdat%grid%data, trim(fldname))
    data(:) = sdat%grid%data%rAttr(kf,:)
  end subroutine dshr_get_griddata

end module dshr_nuopc_mod
