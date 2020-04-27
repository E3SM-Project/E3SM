module dshr_mod

  use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_CompFilterPhaseMap
  use NUOPC_Model      , only : NUOPC_ModelGet
  use ESMF             , only : operator(<), operator(/=), operator(+)
  use ESMF             , only : operator(-), operator(*) , operator(>=)
  use ESMF             , only : operator(<=), operator(>), operator(==)
  use ESMF             , only : ESMF_METHOD_INITIALIZE
  use ESMF             , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_LOGMSG_ERROR, ESMF_MAXSTR
  use ESMF             , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE
  use ESMF             , only : ESMF_State, ESMF_StateGet
  use ESMF             , only : ESMF_Field, ESMF_FieldGet
  use ESMF             , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
  use ESMF             , only : ESMF_DistGrid, ESMF_Array, ESMF_ArrayCreate, ESMF_ArrayDestroy
  use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF             , only : ESMF_GeomType_Flag, ESMF_FieldStatus_Flag
  use ESMF             , only : ESMF_Mesh, ESMF_MeshGet, ESMF_MeshCreate, ESMF_MeshDestroy
  use ESMF             , only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, ESMF_GRIDCREATENOPERIDIMUFRM
  use ESMF             , only : ESMF_FILEFORMAT_ESMFMESH, ESMF_Grid
  use ESMF             , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
  use ESMF             , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
  use ESMF             , only : ESMF_ClockPrint, ESMF_ClockAdvance, ESMF_ClockGetAlarmList
  use ESMF             , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet, ESMF_AlarmSet
  use ESMF             , only : ESMF_ALARMLIST_ALL
  use ESMF             , only : ESMF_Calendar
  use ESMF             , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN, ESMF_CALKIND_FLAG
  use ESMF             , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF             , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF             , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast, ESMF_VMGetCurrent
  use ESMF             , only : ESMF_RouteHandle, ESMF_FieldRegrid
  use ESMF             , only : ESMF_TERMORDER_SRCSEQ, ESMF_FieldRegridStore, ESMF_SparseMatrixWrite
  use ESMF             , only : ESMF_Region_Flag, ESMF_REGION_TOTAL, ESMF_MAXSTR
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_sys_mod      , only : shr_sys_abort
  use shr_file_mod     , only : shr_file_setlogunit
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_calendarname
  use shr_cal_mod      , only : shr_cal_datetod2string
  use shr_map_mod      , only : shr_map_fs_remap, shr_map_fs_bilinear
  use shr_map_mod      , only : shr_map_fs_srcmask, shr_map_fs_scalar
  use shr_const_mod    , only : shr_const_spval
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_init_from_infiles
  use dshr_strdata_mod , only : shr_strdata_restWrite, shr_strdata_restRead
  use dshr_strdata_mod , only : shr_strdata_mapset
  use dshr_methods_mod , only : memcheck, chkerr
  use perf_mod         , only : t_startf, t_stopf
  use shr_ncread_mod   , only : shr_ncread_varExists, shr_ncread_varDimSizes, shr_ncread_field4dG
  use pio
  use mct_mod

  implicit none
  public

  public :: dshr_model_initphase
  public :: dshr_init
  public :: dshr_sdat_init
  public :: dshr_create_mesh_from_grid
  public :: dshr_create_mesh_from_scol
  public :: dshr_set_runclock
  public :: dshr_restart_read
  public :: dshr_restart_write
  public :: dshr_get_atm_adjustment_factors
  public :: dshr_log_clock_advance
  public :: dshr_state_getscalar
  public :: dshr_state_setscalar

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

  ! used/reused in module
  logical                     :: isPresent
  character(len=1024)         :: msgString
  type(ESMF_FieldStatus_Flag) :: status

  ! Module data
  integer     , parameter :: master_task = 0
  integer                 :: iunset = -999
  integer     , parameter :: SecPerDay = 86400 ! Seconds per day
  integer     , parameter :: dbug = 10
  character(*), parameter :: modName =  "(dshr_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
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
  subroutine dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
       flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
       logunit, shrlogunit, rc)

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
    integer          , intent(out)   :: logunit
    integer          , intent(out)   :: shrlogunit
    integer          , intent(out)   :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    logical           :: isPresent, isSet
    character(len=CL) :: cvalue
    character(len=CL) :: logmsg
    character(len=CL) :: diro
    character(len=CL) :: logfile
    character(len=*),parameter  :: subname='(dshr_advertise)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! generate local mpi comm
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=my_task, rc=rc)
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

    ! set output logging
    shrlogunit = 6
    if (my_task == master_task) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))
    else
       logUnit = 6
    endif
    call shr_file_setLogUnit (logunit)

    ! set component instance and suffix
    call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       cvalue = inst_suffix(2:)
       read(cvalue, *) inst_index
    else
       inst_suffix = ""
       inst_index=1
    endif

  end subroutine dshr_init

  !===============================================================================
  subroutine dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, compname, &
       mesh, read_restart, sdat, reset_mask, model_maskfile, rc)

    ! ----------------------------------------------
    ! Initialize sdat
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_GridComp), intent(inout)         :: gcomp
    type(ESMF_Clock)           , intent(in)    :: clock
    character(len=*)           , intent(in)    :: nlfilename ! for shr_strdata_nml namelist
    integer                    , intent(in)    :: logunit
    character(len=*)           , intent(in)    :: compname
    integer                    , intent(out)   :: compid
    type(ESMF_Mesh)            , intent(out)   :: mesh
    logical                    , intent(out)   :: read_restart
    type(shr_strdata_type)     , intent(inout) :: sdat
    logical         , optional , intent(in)    :: reset_mask
    character(len=*), optional , intent(in)    :: model_maskfile
    integer                    , intent(out)   :: rc

    ! local varaibles
    type(ESMF_VM)                :: vm
    type(ESMF_Mesh)              :: mesh_global
    type(ESMF_Calendar)          :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag)      :: esmf_caltype  ! esmf calendar type
    character(CS)                :: calendar      ! calendar name
    character(CL)                :: mesh_filename
    integer                      :: mpicom
    integer                      :: my_task
    logical                      :: scmMode
    real(r8)                     :: scmlat
    real(r8)                     :: scmlon
    character(len=CL)            :: cvalue
    character(len=*), parameter  :: subname='(dshr_mod:dshr_sdat_init)'
    character(*)    , parameter  :: F01="('(dshr_init_strdata) ',a,2f10.4)"
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! generate local mpi comm
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=my_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set compid (for mct)
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    ! Set single column values
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    ! Set restart flag
    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    ! Obtain the input mesh filename
    call NUOPC_CompAttributeGet(gcomp, name='mesh_'//trim(compname), value=mesh_filename, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create the data model mesh
    ! if single column
    !   - read in the data model domain file - and find the nearest neighbor
    ! if not single column
    !   - obtain the mesh directly from the mesh input
    !   - ***TODO: remove the hard-wired domain_atm name here**

    if (trim(mesh_filename) == 'create_mesh') then
       ! get the data model grid from the domain file
       call NUOPC_CompAttributeGet(gcomp, name='domain_atm', value=mesh_filename, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call dshr_create_mesh_from_grid(trim(mesh_filename), mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       if (scmMode) then
          ! verify that are only using 1 pe
          if (my_task > 0) then
             write(logunit,*) subname,' ERROR: scmmode must be run on one pe'
             call shr_sys_abort(subname//' ERROR: scmmode2 tasks')
          endif
          ! obtain the single column lon and lat
          call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlon
          call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlat
          if (my_task == master_task) then
             write(logunit,*) ' scm mode, lon lat = ',scmmode, scmlon,scmlat
          end if
          ! Read in the input mesh
          mesh_global = ESMF_MeshCreate(trim(mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          ! Now create a single column mesh using the single column lats and lons and the global mesh
          call  dshr_create_mesh_from_scol(scmlon, scmlat, logunit, mesh_global, mesh, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshDestroy(mesh_global)
       else
          ! Read in the input mesh
          mesh = ESMF_MeshCreate(trim(mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    if (my_task == master_task) then
       write(logunit,*) trim(subname)// " obtaining "//trim(compname)//" mesh from "// trim(mesh_filename)
    end if

    ! Initialize sdat from data model input files
    call shr_strdata_init_from_infiles(sdat, nlfilename, mesh, clock, mpicom, compid, logunit, &
         reset_mask=reset_mask, model_maskfile=model_maskfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_sdat_init

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
  subroutine dshr_create_mesh_from_scol(scmlon, scmlat, logunit, mesh_global, mesh, rc)

    !-------------------------------------------
    ! Create mct ggrid for model grid and set model gsmap if not input
    ! assumes a very specific netCDF domain file format wrt var names, etc.
    !-------------------------------------------

    ! input/output variables
    real(R8)        , intent(in)    :: scmlon      ! single column lon
    real(R8)        , intent(in)    :: scmlat      ! single column lat
    integer         , intent(in)    :: logunit     ! stdout log unit
    type(ESMF_MESH) , intent(in)    :: mesh_global ! global or regional domain
    type(ESMF_MESH) , intent(out)   :: mesh        ! single column mesh
    integer         , intent(out)   :: rc          ! error code

    ! local variables
    integer              :: n,i,ni             ! indices
    integer              :: ierr               ! error code
    integer, allocatable :: elementCountPTile(:)
    integer, allocatable :: indexCountPDE(:,:)
    integer              :: spatialDim         ! number of dimension in mesh
    integer              :: numOwnedElements   ! number of elements owned by this PET
    integer              :: numOwnedNodes      ! number of nodes owned by this PET
    real(r8), pointer    :: ownedElemCoords(:) ! mesh element coordinates owned by this PET
    real(r8), pointer    :: ownedNodeCoords(:) ! mesh node coordinates owned by this PET
    real(r8), pointer    :: lat(:), lon(:)     ! mesh lats and lons owned by this PET
    real(r8), pointer    :: elemArea(:)        ! mesh areas owned by this PET
    type(ESMF_Array)     :: elemAreaArray
    type(ESMF_Grid)      :: lgrid
    type(ESMF_DistGrid)  :: distGrid           ! mesh distGrid
    real(R8)             :: dist,mind          ! scmmode point search
    real(R8)             :: lscmlon            ! local copy of scmlon
    integer              :: maxIndex(2)
    real(r8)             :: mincornerCoord(2)
    real(r8)             :: maxcornerCoord(2)
    character(*), parameter :: subname = '(shr_strdata_init_model_domain_scol) '
    character(*), parameter :: F00   = "('(shr_strdata_init_model_domain_scol) ',8a)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    ! The input model mesh from the namelist is for the total grid
    ! However the single column mesh is only for a single point - and that is what is
    ! returned from the data model cap to the mediator
    ! Below use scmlon and scmlat to find the nearest neighbor of the input mesh that
    ! will be used for the single column calculation

    ! Read in mesh (this is the global or regional atm mesh)
    call ESMF_MeshGet(mesh_global, spatialDim=spatialDim, &
         numOwnedElements=numOwnedElements, numOwnedNodes=numOwnedNodes, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(ownedNodeCoords(spatialDim*numOwnedNodes))
    call ESMF_MeshGet(mesh_global, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(lon(numOwnedElements))
    allocate(lat(numOwnedElements))
    do n = 1, numOwnedElements
       lon(n) = ownedElemCoords(2*n-1)
       lat(n) = ownedElemCoords(2*n)
    end do

    allocate(elemArea(numOwnedElements))
    elemAreaArray = ESMF_ArrayCreate(distGrid, elemArea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_global, elemMaskArray=elemAreaArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Find nearest neigbor of input domain for scmlon and scmlat
    ! want lon values between 0 and 360, assume 1440 is enough (start with wraparound)

    lscmlon = mod(scmlon+1440.0_r8,360.0_r8)
    lon     = mod(lon   +1440.0_r8,360.0_r8)
    ! lat and lon are on 1D arrays (e.g. spectral element grids)
    mind = 1.0e20
    do i = 1,numOwnedElements
       dist=abs(lscmlon - lon(i)) + abs(scmlat - lat(i))
       if (dist < mind) then
          mind = dist
          ni = i
       endif
    enddo

    ! create the single column grid
    maxIndex(1)       = 1                     ! number of lons
    maxIndex(2)       = 1                     ! number of lats
    mincornerCoord(1) = lon(ni) - elemArea(ni)/2._r8 ! min lon
    mincornerCoord(2) = lat(ni) - elemArea(ni)/2._r8 ! min lat
    maxcornerCoord(1) = lon(ni) + elemArea(ni)/2._r8 ! max lon
    maxcornerCoord(2) = lat(ni) + elemArea(ni)/2._r8 ! max lat
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the single column mesh from the grid
    mesh =  ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(ownedNodeCoords)
    deallocate(lon)
    deallocate(lat)
    call ESMF_ArrayDestroy(elemAreaArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_create_mesh_from_scol

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
    character(len=*),parameter :: subname='dshr_mod:(ModelSetRunClock) '
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
  subroutine dshr_restart_read(rest_filem, rest_files, rpfile, inst_suffix, nullstr, &
       logunit, my_task, mpicom, sdat, fld, fldname)

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_filem
    character(len=*)            , intent(inout) :: rest_files
    character(len=*)            , intent(in)    :: rpfile
    character(len=*)            , intent(in)    :: inst_suffix
    character(len=*)            , intent(in)    :: nullstr
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    real(r8)         , optional , pointer       :: fld(:)
    character(len=*) , optional , intent(in)    :: fldname

    ! local variables
    integer           :: nu
    logical           :: exists  ! file existance
    type(file_desc_t) :: pioid
    type(var_desc_t)  :: varid
    type(io_desc_t)   :: pio_iodesc
    integer           :: rcode
    character(*), parameter :: F00   = "('(dshr_restart_read) ',8a)"
    character(*), parameter :: subName = "(dshr_restart_read) "
    !-------------------------------------------------------------------------------

    if (trim(rest_filem) == trim(nullstr) .and. trim(rest_files) == trim(nullstr)) then
       if (my_task == master_task) then
          write(logunit,F00) ' restart filenames from rpointer'
          inquire(file=trim(rpfile)//trim(inst_suffix), exist=exists)
          if (.not.exists) then
             write(logunit, F00) ' ERROR: rpointer file does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
          endif
          open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
          read(nu, '(a)') rest_filem
          read(nu, '(a)') rest_files
          close(nu)
          inquire(file=trim(rest_files), exist=exists)
       endif
       call shr_mpi_bcast(rest_filem, mpicom, 'rest_filem')
       call shr_mpi_bcast(rest_files, mpicom, 'rest_files')
    else
       ! use namelist already read
       if (my_task == master_task) then
          write(logunit, F00) ' restart filenames from namelist '
          inquire(file=trim(rest_files), exist=exists)
       endif
    endif
    call shr_mpi_bcast(exists, mpicom, 'exists')
    if (exists) then
       if (my_task == master_task) write(logunit, F00) ' reading data mdoel restart ', trim(rest_filem)
       if (present(fld) .and. present(fldname)) then
          rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(rest_filem), pio_nowrite)
          call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
          call pio_initdecomp(sdat%pio_subsystem, pio_double, (/sdat%gsize/), sdat%gindex, pio_iodesc)
          rcode = pio_inq_varid(pioid, trim(fldname), varid)
          call pio_read_darray(pioid, varid, pio_iodesc, fld, rcode)
          call pio_closefile(pioid)
          call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
       end if
       call shr_strdata_restRead(trim(rest_files), sdat, mpicom)
    else
       if (my_task == master_task) write(logunit, F00) ' file not found, skipping ',trim(rest_files)
    endif
  end subroutine dshr_restart_read

  !===============================================================================
  subroutine dshr_restart_write(rpfile, case_name, model_name, inst_suffix, ymd, tod, &
       logunit, mpicom, my_task, sdat, fld, fldname)

    ! input/output variables
    character(len=*)            , intent(in)    :: rpfile
    character(len=*)            , intent(in)    :: case_name
    character(len=*)            , intent(in)    :: model_name
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: ymd       ! model date
    integer                     , intent(in)    :: tod       ! model sec into model date
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    real(r8)         , optional , pointer       :: fld(:)
    character(len=*) , optional , intent(in)    :: fldname

    ! local variables
    integer           :: nu
    character(len=CL) :: rest_filem
    character(len=CL) :: rest_files
    character(len=CS) :: date_str
    type(file_desc_t) :: pioid
    integer           :: dimid(1)
    type(var_desc_t)  :: varid
    type(io_desc_t)   :: pio_iodesc
    integer           :: rcode
    !-------------------------------------------------------------------------------

     ! write data model restart data
     call shr_cal_datetod2string(date_str, ymd, tod)
     write(rest_filem,"(7a)") trim(case_name),'.', trim(model_name),trim(inst_suffix),'.r.'  , trim(date_str),'.nc'
     write(rest_files,"(7a)") trim(case_name),'.', trim(model_name),trim(inst_suffix),'.rs1.', trim(date_str),'.bin'
     if (my_task == master_task) then
        open(newunit=nu, file=trim(rpfile)//trim(inst_suffix), form='formatted')
        write(nu,'(a)') rest_filem
        write(nu,'(a)') rest_files
        close(nu)
        write(logunit,*)' (dshr_restart_write) writing ',trim(rest_files), ymd, tod
     endif
     if (present(fld) .and. present(fldname)) then
        rcode = pio_createfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(rest_filem), pio_clobber)
        call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
        rcode = pio_put_att(pioid, pio_global, "version", "nuopc_data_models_v0")
        rcode = pio_def_dim(pioid, 'gsize', sdat%gsize, dimid(1))
        rcode = pio_def_var(pioid, trim(fldname), PIO_DOUBLE, dimid, varid)
        rcode = pio_enddef(pioid)
        call pio_initdecomp(sdat%pio_subsystem, pio_double, (/sdat%gsize/), sdat%gindex, pio_iodesc)
        call pio_write_darray(pioid, varid, pio_iodesc, fld, rcode, fillval=shr_const_spval)
        call pio_closefile(pioid)
        call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
     end if

     ! write stream restart data
     call shr_strdata_restWrite(trim(rest_files), sdat, mpicom, trim(case_name), &
          'sdat strdata from '//trim(model_name))

  end subroutine dshr_restart_write

  !===============================================================================
  subroutine dshr_get_atm_adjustment_factors(fileName, windF, winddF, qsatF, &
       mpicom, compid, masterproc, logunit, sdat)

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
    call mct_gGrid_init(GGrid=gGridi, CoordChars='lat:lon:hgt', OtherChars='mask', lsize=lsizei )
    call mct_aVect_init(avi, rList="wind:windd:qsat", lsize=lsizei)
    avi%rAttr = shr_const_spval

    ! determine local size on data model grid
    lsizeo = sdat%lsize

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
             if (abs(ggridi%data%rAttr(klon,n)-ggridoG%data%rAttr(klon,n)) > 0.01_R8) then
                if (abs(ggridi%data%rAttr(klon,n)-ggridoG%data%rAttr(klon,n)) - 360._r8 > 0.01_R8) then
                   write(6,*)'domap is true: n,londiff= ',n,abs(ggridi%data%rAttr(klon,n)-ggridoG%data%rAttr(klon,n))
                   domap = .true.
                end if
             end if
             if (abs(ggridi%data%rAttr(klat,n)-ggridoG%data%rAttr(klat,n)) > 0.01_R8) then
                write(6,*)'domap is true: n,latdiff= ',n,abs(ggridi%data%rAttr(klat,n)-ggridoG%data%rAttr(klat,n))
                domap=.true.
             end if
          enddo
       endif

       call mct_gGrid_clean(ggridoG)

    endif

    call shr_mpi_bcast(domap,mpicom,subname//' domap')

    if (domap) then
       call shr_strdata_mapSet(smatp, &
            ggridi, gsmapi, ni0 ,nj0, sdat%grid, sdat%gsmap, sdat%nxg, sdat%nyg, &
            'datmfactor', shr_map_fs_remap, shr_map_fs_bilinear, &
            shr_map_fs_srcmask, shr_map_fs_scalar, compid, mpicom, 'Xonly')

       call mct_aVect_init(avo, avi, lsizeo)
       call mct_sMat_avMult(avi, smatp, avo)
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
  subroutine dshr_log_clock_advance(clock, component, logunit, rc)

    ! input/output variables
    type(ESMF_Clock)               :: clock
    character(len=*) , intent(in)  :: component
    integer          , intent(in)  :: logunit
    integer          , intent(out) :: rc

    ! local variables
    character(len=CL) :: cvalue, prestring
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, preString=trim(prestring), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine dshr_log_clock_advance

  !===============================================================================
  subroutine dshr_state_getscalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)     :: state
    integer,          intent(in)     :: scalar_id
    real(r8),         intent(out)    :: scalar_value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask, ierr, len
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(r8), pointer :: farrayptr(:,:)
    real(r8)          :: tmp(1)
    character(len=*), parameter :: subname='(state_getscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == master_task) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif
      tmp(:) = farrayptr(scalar_id,:)
    endif
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    scalar_value = tmp(1)

  end subroutine dshr_state_getscalar

  !================================================================================
  subroutine dshr_state_setscalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output arguments
    real(r8),         intent(in)     :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: lfield
    type(ESMF_VM)     :: vm
    real(r8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(state_setscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == master_task) then
       call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       endif
       farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine dshr_state_setscalar

end module dshr_mod
