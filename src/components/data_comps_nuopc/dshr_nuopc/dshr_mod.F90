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
  use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF             , only : ESMF_GeomType_Flag, ESMF_FieldStatus_Flag
  use ESMF             , only : ESMF_Mesh, ESMF_MeshGet, ESMF_MeshCreate
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
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_strdata_mod , only : shr_strdata_readnml
  use dshr_strdata_mod , only : shr_strdata_pioinit
  use dshr_strdata_mod , only : shr_strdata_restWrite
  use dshr_strdata_mod , only : shr_strdata_restRead
  use dshr_strdata_mod , only : shr_strdata_init_model_domain
  use dshr_strdata_mod , only : shr_strdata_init_streams
  use dshr_strdata_mod , only : shr_strdata_init_mapping, shr_strdata_mapset
  use dshr_strdata_mod , only : shr_strdata_print
  use dshr_stream_mod  , only : shr_stream_set, shr_stream_taxis_extend
  use dshr_util_mod    , only : memcheck, chkerr
  use perf_mod         , only : t_startf, t_stopf 
  use shr_ncread_mod   , only : shr_ncread_varExists, shr_ncread_varDimSizes, shr_ncread_field4dG
  use pio
  use mct_mod

  implicit none
  public

  interface dshr_sdat_init
     module procedure :: dshr_sdat_init_from_input  ! initialize sdat from a fortran interface
     module procedure :: dshr_sdat_init_from_strtxt ! initialize sdat from stream text file
  end interface dshr_sdat_init

  public :: dshr_model_initphase
  public :: dshr_init
  public :: dshr_sdat_init          
  public :: dshr_create_mesh_from_grid
  public :: dshr_set_runclock
  public :: dshr_restart_read
  public :: dshr_restart_write
  public :: dshr_get_atm_adjustment_factors
  public :: dshr_get_griddata
  public :: dshr_set_griddata
  public :: dshr_log_clock_advance
  public :: dshr_state_getscalar
  public :: dshr_state_setscalar
  public :: dshr_state_diagnose
  public :: dshr_state_getfldptr
  public :: dshr_fldbun_regrid
  public :: dshr_fldbun_getFieldN
  public :: dshr_fldbun_getNameN
  public :: dshr_fldbun_getFldPtr
  public :: dshr_fldbun_diagnose
  public :: dshr_fldbun_fldchk

  private :: dshr_field_getfldptr
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
  subroutine dshr_init(gcomp, master_task, mpicom, my_task, inst_index, inst_suffix, &
       flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
       logunit, shrlogunit, rc)

    ! input/output variables
    type(ESMF_GridComp)              :: gcomp
    integer          , intent(in)    :: master_task
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
  subroutine dshr_sdat_init_from_input( &
       SDAT, mpicom, compid, mesh, nxg, nyg, clock, &
       !--- streams stuff required ---
       yearFirst, yearLast, yearAlign, offset,          &
       domFilePath, domFileName,                        &
       domTvarName, domXvarName, domYvarName, domMaskName, &
       filePath, filename, fldListFile, fldListModel,   &
       !--- strdata optional ---
       nzg, domZvarName,                                &
       taxMode, dtlimit, tintalgo, readmode,            &
       fillalgo, fillmask, fillread, fillwrite,         &
       mapalgo, mapmask, mapread, mapwrite)


    ! Set strdata and stream info from fortran interface.
    ! Note: When this is called, previous settings are reset to defaults
    !and then the values passed are used.

    ! input/output arguments
    type(shr_strdata_type) ,intent(inout):: SDAT         ! strdata data data-type
    integer                ,intent(in)   :: mpicom       ! mpi comm
    integer                ,intent(in)   :: compid
    type(ESMF_Mesh)        ,intent(in)   :: mesh
    integer                ,intent(in)   :: nxg
    integer                ,intent(in)   :: nyg
    type(ESMF_Clock)       ,intent(in)   :: clock
    integer                ,intent(in)   :: yearFirst    ! first year to use
    integer                ,intent(in)   :: yearLast     ! last  year to use
    integer                ,intent(in)   :: yearAlign    ! align yearFirst with this model year
    integer                ,intent(in)   :: offset       ! offset in seconds of stream data
    character(*)           ,intent(in)   :: domFilePath  ! domain file path
    character(*)           ,intent(in)   :: domFileName  ! domain file name
    character(*)           ,intent(in)   :: domTvarName  ! domain time dim name
    character(*)           ,intent(in)   :: domXvarName  ! domain x dim name
    character(*)           ,intent(in)   :: domYvarName  ! domain y dim name
    character(*)           ,intent(in)   :: domMaskName  ! domain mask name
    character(*)           ,intent(in)   :: filePath     ! path to filenames
    character(*)           ,intent(in)   :: filename(:)  ! filename for index filenumber
    character(*)           ,intent(in)   :: fldListFile  ! file field names, colon delim list
    character(*)           ,intent(in)   :: fldListModel ! model field names, colon delim list
    integer      ,optional ,intent(in)   :: nzg
    character(*) ,optional ,intent(in)   :: domZvarName  ! domain z dim name
    character(*) ,optional ,intent(in)   :: taxMode
    real(R8)     ,optional ,intent(in)   :: dtlimit
    character(*) ,optional ,intent(in)   :: fillalgo     ! fill algorithm
    character(*) ,optional ,intent(in)   :: fillmask     ! fill mask
    character(*) ,optional ,intent(in)   :: fillread     ! fill mapping file to read
    character(*) ,optional ,intent(in)   :: fillwrite    ! fill mapping file to write
    character(*) ,optional ,intent(in)   :: mapalgo      ! scalar map algorithm
    character(*) ,optional ,intent(in)   :: mapmask      ! scalar map mask
    character(*) ,optional ,intent(in)   :: mapread      ! regrid mapping file to read
    character(*) ,optional ,intent(in)   :: mapwrite     ! regrid mapping file to write
    character(*) ,optional ,intent(in)   :: tintalgo     ! time interpolation algorithm
    character(*) ,optional ,intent(in)   :: readmode     ! file read mode

    ! local variables
    type(ESMF_Calendar)     :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype  ! esmf calendar type
    character(CS)           :: calendar      ! calendar name
    character(CS)           :: zname
    integer                 :: my_task
    integer                 :: ierr
    integer                 :: rc
    character(*),parameter  :: subName = "(shr_strdata_create) "
    character(*),parameter  :: F00 = "('(shr_strdata_create) ',8a)"
    !-------------------------------------------------------------------------------

    ! The following just sets the defaults - but does not do a namelist read
    call shr_strdata_readnml(SDAT, mpicom=mpicom)  
    call mpi_comm_rank(mpicom, my_task, ierr)

    ! Assume only 1 stream
    SDAT%nstreams = 1

    ! Initialize pio
    call shr_strdata_pioinit(sdat, compid)

    if (present(taxMode)) then
       SDAT%taxMode(1) = taxMode
       if (trim(SDAT%taxMode(1)) == trim(shr_stream_taxis_extend)) SDAT%dtlimit(1) = 1.0e30
    endif
    if (present(dtlimit   )) SDAT%dtlimit(1)  = dtlimit
    if (present(fillalgo  )) SDAT%fillalgo(1) = fillalgo
    if (present(fillmask  )) SDAT%fillmask(1) = fillmask
    if (present(fillread  )) SDAT%fillread(1) = fillread
    if (present(fillwrite )) SDAT%fillwrit(1) = fillwrite
    if (present(mapalgo   )) SDAT%mapalgo(1)  = mapalgo
    if (present(mapmask   )) SDAT%mapmask(1)  = mapmask
    if (present(mapread   )) SDAT%mapread(1)  = mapread
    if (present(mapwrite  )) SDAT%mapwrit(1)  = mapwrite
    if (present(tintalgo  )) SDAT%tintalgo(1) = tintalgo
    if (present(readmode  )) SDAT%readmode(1) = readmode
    if (present(mapmask   )) SDAT%mapmask(1)  = mapmask

    ! initialize sdat model domain info
    call shr_strdata_init_model_domain(mesh, mpicom, compid, sdat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initialize sdat stream domains
    if (present(domZvarName)) then
       zname = trim(domzvarname)
    else
       zname = 'undefined'
    end if
    call shr_stream_set(SDAT%stream(1), &
         yearFirst, yearLast, yearAlign, offset, taxMode,  &
         fldListFile, fldListModel, domFilePath, domFileName,  &
         domTvarName, domXvarName, domYvarName, trim(zname), domMaskName,  &
         filePath, filename)

    !call shr_strdata_init_streams(sdat, compid, mpicom, my_task)

    ! initialize sdat attributes mapping of streams to model domain
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

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

  end subroutine dshr_sdat_init_from_input

  !===============================================================================
  subroutine dshr_sdat_init_from_strtxt(gcomp, clock, nlfilename, compid, logunit, compname, &
       mesh, read_restart, sdat, reset_domain_mask, rc)

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
    logical         , optional , intent(in)    :: reset_domain_mask
    integer                    , intent(out)   :: rc

    ! local varaibles
    type(ESMF_Calendar)          :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag)      :: esmf_caltype  ! esmf calendar type
    character(CS)                :: calendar      ! calendar name
    type(ESMF_VM)                :: vm
    integer                      :: mpicom
    integer                      :: my_task
    integer                      :: master_task = 0
    character(CL)                :: mesh_filename
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
    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon
    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    ! Set restart flag
    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    ! Read shr_strdata_nml from nlfilename
    ! Read sdat namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdat, trim(nlfilename), mpicom=mpicom)

    ! Initialize sdat  pio
    call shr_strdata_pioinit(sdat, compid)

    ! Obtain the data model mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_'//trim(compname), value=mesh_filename, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(mesh_filename) == 'create_mesh') then
       ! get the data model grid from the domain file
       call NUOPC_CompAttributeGet(gcomp, name='domain_atm', value=mesh_filename, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call dshr_create_mesh_from_grid(trim(mesh_filename), mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       mesh = ESMF_MeshCreate(trim(mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (my_task == master_task) then
       write(logunit,*) trim(subname)// " obtaining "//trim(compname)//" mesh from "// trim(mesh_filename)
    end if

    ! Initialize the sdat model domain info
    if (scmmode) then
       if (my_task == master_task) then
          write(logunit,*) ' scm mode, lon lat = ',scmmode, scmlon,scmlat
       end if
       call shr_strdata_init_model_domain(scmlon, scmlat, mpicom, compid, mesh, sdat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! ***TODO: Now need to initialize mesh that will send back to mediator given the lat and lon
       !    of the original mesh that was selected ***
    else
       call shr_strdata_init_model_domain(mesh, mpicom, compid, sdat, reset_domain_mask=reset_domain_mask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! initialize sdat stream domains
    call shr_strdata_init_streams(sdat, compid, mpicom, my_task)

    ! initialize sdat attributes mapping of streams to model domain
    call shr_strdata_init_mapping(sdat, compid, mpicom, my_task)

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

    ! print sdat output
    if (my_task == master_task) then
       call shr_strdata_print(sdat,'SDAT data ')
       write(logunit,*) ' successfully initialized sdat'
    endif

  end subroutine dshr_sdat_init_from_strtxt

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
       logunit, my_task, master_task, mpicom, sdat, fld, fldname)

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_filem
    character(len=*)            , intent(inout) :: rest_files
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
       logunit, mpicom, my_task, master_task, sdat, fld, fldname)

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
    lsizeo = mct_gsmap_lsize(sdat%gsmap, mpicom)
    call mct_gGrid_init(GGrid=gGridi, CoordChars='lat:lon:hgt', OtherChars='mask', lsize=lsizei )
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
       call shr_strdata_mapSet(smatp, &
            ggridi, gsmapi, ni0 ,nj0, &
            sdat%grid, sdat%gsmap, sdat%nxg, sdat%nyg, &            
            'datmfactor', shr_map_fs_remap, shr_map_fs_bilinear, &
            shr_map_fs_srcmask, shr_map_fs_scalar, &
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

    if (mytask == 0) then
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

    if (mytask == 0) then
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

  !===============================================================================
  subroutine dshr_state_diagnose(State, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: i,j,n
    type(ESMf_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(r8), pointer               :: dataPtr1d(:)
    real(r8), pointer               :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(dshr_state_diagnose)'
    ! ----------------------------------------------

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_field_getfldptr(lfield, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    deallocate(lfieldnamelist)

  end subroutine dshr_state_diagnose

  !===============================================================================
  subroutine dshr_State_GetFldPtr(State, fldname, fldptr1, fldptr2, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) ,          intent(in)              :: State
    character(len=*) ,          intent(in)              :: fldname
    real(R8)         , pointer, intent(inout), optional :: fldptr1(:)
    real(R8)         , pointer, intent(inout), optional :: fldptr2(:,:)
    integer          ,          intent(out)             :: rc

    ! local variables
    type(ESMF_Field)           :: lfield
    character(len=*), parameter :: subname='(dshr_state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call dshr_field_getfldptr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_State_GetFldPtr

  !===============================================================================
  subroutine dshr_field_getfldptr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(r8), pointer , intent(inout), optional :: fldptr1(:)
    real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)             :: rc

    ! local variables
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    character(len=*), parameter :: subname='(field_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then

       lrank = 0
       if (labort) then
          call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       else
          call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       endif

    else

       call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_GRID) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (nnodes == 0 .and. nelements == 0) lrank = 0
       else  
          call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
               ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       endif ! geomtype

       if (lrank == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
               ESMF_LOGMSG_INFO)
       elseif (lrank == 1) then
          if (.not.present(fldptr1)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (lrank == 2) then
          if (.not.present(fldptr2)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       endif

    endif  ! status

    if (present(rank)) then
       rank = lrank
    endif

  end subroutine dshr_field_getfldptr

  !===============================================================================
  subroutine dshr_fldbun_diagnose(FB, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    character(len=*), parameter     :: subname='(dshr_fldbun_diagnose)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) lstring = trim(string) // ' '

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call dshr_fldbun_GetFldPtr(FB, lfieldnamelist(n), fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine dshr_fldbun_diagnose

  !===============================================================================
  subroutine dshr_fldbun_GetFldPtr(FB, fldname, fldptr1, fldptr2, rank, field, rc)

    ! ----------------------------------------------
    ! Get pointer to a field bundle field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle) , intent(in)              :: FB
    character(len=*)       , intent(in)              :: fldname
    real(R8), pointer      , intent(inout), optional :: fldptr1(:)
    real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
    integer                , intent(out),   optional :: rank
    type(ESMF_Field)       , intent(out),   optional :: field
    integer                , intent(out)              :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    character(len=*), parameter :: subname='(dshr_fldbun_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. dshr_fldbun_FldChk(FB, trim(fldname), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR field "//trim(fldname)//" not in FB ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_field_getfldptr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) rank = lrank
    if (present(field)) field = lfield

  end subroutine dshr_fldbun_getfldptr

  !===============================================================================
  subroutine dshr_fldbun_Regrid(FBin, FBout, RH, zeroregion, rc)

    ! ----------------------------------------------
    ! Assumes that FBin and FBout contain fields with the same name
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)        :: FBin
    type(ESMF_FieldBundle), intent(inout)        :: FBout
    type(ESMF_RouteHandle), intent(inout)        :: RH
    type(ESMF_Region_Flag), intent(in), optional :: zeroregion
    integer               , intent(out)          :: rc

    ! local
    type(ESMF_Field)           :: field_src
    type(ESMF_Field)           :: field_dst
    integer                    :: fieldcount
    logical                    :: checkflag = .false.
    character(len=8)           :: filename
    type(ESMF_Region_Flag)     :: localzr
    character(CS)              :: fldname
    integer                    :: n
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(len=*),parameter :: subname='(dshr_fldbun_FieldRegrid)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf(subname)

    localzr = ESMF_REGION_TOTAL
    if (present(zeroregion)) then
       localzr = zeroregion
    endif

    call ESMF_FieldBundleGet(FBin, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBin, fieldNameList=lfieldNameList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,fieldCount
       fldname = trim(lfieldnamelist(n))
       call ESMF_FieldBundleGet(FBin, fieldName=trim(fldname), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(FBout, fieldName=trim(fldname), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RH, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=.false., zeroregion=localzr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    call t_stopf(subname)

  end subroutine dshr_fldbun_Regrid

  !===============================================================================
  subroutine dshr_fldbun_getFieldN(FB, fieldnum, field, rc)

    ! ----------------------------------------------
    ! Get field with number fieldnum in input field bundle FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    character(len=*),parameter :: subname='(dshr_fldbun_getFieldN)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call dshr_fldbun_getNameN(FB, fieldnum, name, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine dshr_fldbun_getFieldN

  !===============================================================================
  subroutine dshr_fldbun_getNameN(FB, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get name of field number fieldnum in input field bundle FB
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    character(len=*)      , intent(out)   :: fieldname
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(dshr_fldbun_getNameN)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    fieldname = ' '
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)
    deallocate(lfieldnamelist)

  end subroutine dshr_fldbun_getNameN

  !===============================================================================
  logical function dshr_fldbun_FldChk(FB, fldname, rc)

    ! ----------------------------------------------
    ! Determine if field with fldname is in input field bundle
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    integer               , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(dshr_fldbun_FldChk)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! If field bundle is created determine if fldname is present in field bundle
    dshr_fldbun_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//" Error checking field: "//trim(fldname), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif

    if (isPresent) then
       dshr_fldbun_FldChk = .true.
    endif

  end function dshr_fldbun_FldChk

  !===============================================================================
  subroutine dshr_fldbun_Field_diagnose(FB, fieldname, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)  :: FB
    character(len=*), intent(in)           :: fieldname
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer           :: lrank
    character(len=CS) :: lstring
    real(R8), pointer :: dataPtr1d(:)
    real(R8), pointer :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(dshr_fldbun_FieldDiagnose)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) lstring = trim(string)

    call dshr_fldbun_GetFldPtr(FB, fieldname, dataPtr1d, dataPtr2d, lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
       ! no local data
    elseif (lrank == 1) then
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    elseif (lrank == 2) then
       if (size(dataPtr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    else
       call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  end subroutine dshr_fldbun_Field_diagnose

end module dshr_mod
