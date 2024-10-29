module atm_comp_mct

  ! Modules used acros atm_xyz_mct routines

  ! shr mods
  use shr_kind_mod,     only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CL=>SHR_KIND_CL, CS=>SHR_KIND_CS

  ! seq mods
  use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata
  use seq_timemgr_mod,  only: seq_timemgr_EClockGetData

  ! toolkits mods
  use esmf,             only: ESMF_Clock
  use mct_mod,          only: mct_aVect, mct_gsMap, mct_gGrid

  ! MPI
  use mpi

#ifdef HAVE_MOAB
  use seq_comm_mct     , only: mphaid ! atm physics grid id in MOAB, on atm pes
  use iso_c_binding
#ifdef MOABCOMP
  use seq_comm_mct, only:  seq_comm_compare_mb_mct
#endif
#endif

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

#if defined(SCREAM_SYSTEM_WORKAROUND) && (SCREAM_SYSTEM_WORKAROUND == 1)
  public :: atm_init_hip_mct
#endif
  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  integer                :: mpicom_atm          ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: ATM_ID              ! mct comp id
  integer(IN),parameter  :: master_task=0       ! task number of master task

  integer                :: atm_log_unit


#ifdef HAVE_MOAB
  ! to store all fields to be set in moab
  integer , private :: mblsize, totalmbls, nsend, totalmbls_r, nrecv
  real(r8) , allocatable, target, private :: a2x_am(:,:) ! atm to coupler, on atm mesh, on atm component pes
  real(r8) , allocatable, target, private :: x2a_am(:,:) ! coupler to atm, on atm mesh, on atm component pes
#ifdef MOABCOMP
  integer  :: mpicom_atm_moab ! used just for mpi-reducing the difference between moab tags and mct avs
  integer :: rank2
#endif
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#if defined(SCREAM_SYSTEM_WORKAROUND) && (SCREAM_SYSTEM_WORKAROUND == 1)
  !===============================================================================
  subroutine atm_init_hip_mct()
    use scream_f2c_mod, only: scream_init_hip_atm

    call scream_init_hip_atm()

  end subroutine atm_init_hip_mct
#endif

  !===============================================================================
  subroutine atm_init_mct( EClock, cdata, x2a, a2x, NLFilename )
    use iso_c_binding,      only: c_ptr, c_loc, c_int, c_char, c_bool
    use scream_f2c_mod,     only: scream_create_atm_instance, scream_setup_surface_coupling, &
                                  scream_init_atm
    use scream_cpl_indices, only: scream_set_cpl_indices, &
                                  num_cpl_imports,          num_scream_imports, &
                                  num_cpl_exports,          num_scream_exports, &
                                  import_field_size,        export_field_size, &
                                  import_field_names,       export_field_names, &
                                  import_cpl_indices,       export_cpl_indices, &
                                  import_vector_components, export_vector_components, &
                                  import_constant_multiple, export_constant_multiple, &
                                  do_import_during_init,    do_export_during_init
    use ekat_string_utils,  only: string_f2c
    use mct_mod,            only: mct_aVect_init, mct_gsMap_lsize
    use seq_flds_mod,       only: seq_flds_a2x_fields, seq_flds_x2a_fields
    use seq_infodata_mod,   only: seq_infodata_start_type_start, seq_infodata_start_type_cont, &
                                  seq_infodata_start_type_brnch
    use seq_comm_mct,       only: seq_comm_inst, seq_comm_name, seq_comm_suffix
    use shr_file_mod,       only: shr_file_getunit, shr_file_setIO
    use shr_sys_mod,        only: shr_sys_abort
#ifdef HAVE_MOAB
    use shr_kind_mod     ,  only: cxx=>shr_kind_cxx
    use mct_mod,            only: mct_avect_nRattr
    use iMOAB,              only: iMOAB_SetDoubleTagStorage
#ifdef MOABDEBUG
    use iMOAB,              only: iMOAB_WriteMesh
#endif
#endif
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local variables ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap_atm
    type(mct_gGrid)        , pointer :: dom_atm
    integer(IN)                      :: phase          ! initialization phase
    integer(IN)                      :: ierr,mpi_ierr  ! error codes
    integer(IN)                      :: modelio_fid    ! file descriptor for atm_modelio.nml
    integer(IN)                      :: case_start_tod, case_start_ymd, cur_tod, cur_ymd
    integer                          :: lsize
    character(CS)                    :: run_type       ! type of run
    type(c_ptr) :: x2a_ptr, a2x_ptr
    character(len=256)               :: atm_log_fname  ! name of ATM log file
    character(CL)                    :: calendar       ! calendar string
#ifdef HAVE_MOAB
    integer        :: ent_type
    character(CXX) :: tagname ! will store all seq_flds_a2x_fields
#ifdef MOABDEBUG
    character*100 outfile, wopts
#endif
#endif
    ! TODO: read this from the namelist?
    character(len=256)                :: yaml_fname = "./data/scream_input.yaml"
    character(kind=c_char,len=256), target :: yaml_fname_c, atm_log_fname_c
    character(len=256) :: caseid, username, hostname, rest_caseid
    character(kind=c_char,len=256), target :: caseid_c, username_c, hostname_c, calendar_c, rest_caseid_c
    integer (kind=c_int) :: run_type_c
    !-------------------------------------------------------------------------------

    ! Grab some data from the cdata structure (coming from the coupler)
    call seq_cdata_setptrs(cdata, &
         id=ATM_ID, &
         mpicom=mpicom_atm, &
         gsMap=gsmap_atm, &
         dom=dom_atm, &
         infodata=infodata)
    call seq_infodata_getData(infodata, atm_phase=phase, start_type=run_type, &
                              username=username, case_name=caseid, &
                              rest_case_name=rest_caseid, hostname=hostname)
    call seq_infodata_PutData(infodata, atm_aero=.true.)
    call seq_infodata_PutData(infodata, atm_prognostic=.true.)

    if (phase > 1) RETURN

    ! Determine instance information
    inst_name   = seq_comm_name(ATM_ID)
    inst_index  = seq_comm_inst(ATM_ID)
    inst_suffix = seq_comm_suffix(ATM_ID)

    ! Determine communicator group
    call mpi_comm_rank(mpicom_atm, my_task, ierr)

    !----------------------------------------------------------------------------
    ! Init atm.log
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
      atm_log_unit = shr_file_getUnit()
      call shr_file_setIO ('atm_modelio.nml'//trim(inst_suffix),atm_log_unit)
      inquire(unit=atm_log_unit,name=atm_log_fname)
    endif

    call mpi_bcast(atm_log_unit,1,MPI_INTEGER,master_task,mpicom_atm,mpi_ierr)
    if (ierr /= 0) then
      print *,'[eamxx] ERROR broadcasting atm.log unit'
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    end if

    call mpi_bcast(atm_log_fname,256,MPI_CHARACTER,master_task,mpicom_atm,ierr)
    if (ierr /= 0) then
      print *,'[eamxx] ERROR broadcasting atm.log file name'
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    end if

    !----------------------------------------------------------------------------
    ! Initialize atm
    !----------------------------------------------------------------------------

    if (trim(run_type) == trim(seq_infodata_start_type_start)) then
      run_type_c = 0
    else if (trim(run_type) == trim(seq_infodata_start_type_cont) ) then
      run_type_c = 1
    else if (trim(run_type) == trim(seq_infodata_start_type_brnch) ) then
      run_type_c = 1
    else
      print *, "[eamxx] ERROR! Unsupported run type: "//trim(run_type)
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

    ! Init the AD
    call seq_timemgr_EClockGetData(EClock, calendar=calendar, &
                                   curr_ymd=cur_ymd, curr_tod=cur_tod, &
                                   start_ymd=case_start_ymd, start_tod=case_start_tod)
    call string_f2c(yaml_fname,yaml_fname_c)
    call string_f2c(calendar,calendar_c)
    call string_f2c(trim(atm_log_fname),atm_log_fname_c)
    call string_f2c(trim(caseid),caseid_c)
    call string_f2c(trim(rest_caseid),rest_caseid_c)
    call string_f2c(trim(hostname),hostname_c)
    call string_f2c(trim(username),username_c)
    call scream_create_atm_instance (mpicom_atm, ATM_ID, yaml_fname_c, atm_log_fname_c, run_type_c, &
                          INT(cur_ymd,kind=C_INT),  INT(cur_tod,kind=C_INT), &
                          INT(case_start_ymd,kind=C_INT), INT(case_start_tod,kind=C_INT), &
                          calendar_c, caseid_c, rest_caseid_c, hostname_c, username_c)

    ! Init MCT gsMap
    call atm_Set_gsMap_mct (mpicom_atm, ATM_ID, gsMap_atm)
    lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

    ! Init MCT domain structure
    call atm_domain_mct (lsize, gsMap_atm, dom_atm)

#ifdef HAVE_MOAB
    call moab_atm_phys_scream()
#ifdef MOABCOMP
    mpicom_atm_moab = mpicom_atm ! just to store it to be used later in MOABCOMP
    rank2 = my_task ! again, just to use it for MOABCOMP output
#endif
#endif

    ! Init import/export mct attribute vectors
    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)

    ! Init surface coupling stuff in the AD
    call scream_set_cpl_indices (x2a, a2x)

#ifdef HAVE_MOAB
! for MOAB driver, need to setup surface coupling for MOAB arrays that will hold
! fields data for import and export
!  setup is similar to setup for the attribute vectors real data arrays (x2a%rAttr
!  and a2x%rAttr)
! the major difference is that moab arrays are transposed compared to AVs arrays
    mblsize = lsize
    nsend = mct_avect_nRattr(a2x)
    totalmbls = mblsize * nsend   ! size of the double array
    allocate (a2x_am(mblsize, nsend) )
    a2x_am = 0. ! initialize everything with 0

    nrecv = mct_avect_nRattr(x2a)
    totalmbls_r = mblsize * nrecv   ! size of the double array used to receive
    allocate (x2a_am(mblsize, nrecv) ) ! these will be received by moab tags, then used to set cam in surf data
    x2a_am = 0. ! initialize everything with 0
#endif

    call scream_setup_surface_coupling (c_loc(import_field_names), c_loc(import_cpl_indices), &
                                        c_loc(x2a%rAttr), &
#ifdef HAVE_MOAB
                                        c_loc(x2a_am),  &
#endif
                                        c_loc(import_vector_components), &
                                        c_loc(import_constant_multiple), c_loc(do_import_during_init), &
                                        num_cpl_imports, num_scream_imports, import_field_size, &
                                        c_loc(export_field_names), c_loc(export_cpl_indices), &
                                        c_loc(a2x%rAttr), &
#ifdef HAVE_MOAB
                                        c_loc(a2x_am),  &
#endif
                                        c_loc(export_vector_components), &
                                        c_loc(export_constant_multiple), c_loc(do_export_during_init), &
                                        num_cpl_exports, num_scream_exports, export_field_size)

    call scream_init_atm ()
#ifdef HAVE_MOAB
    ! data should be set now inside moab from import and export fields
    ! do we import and export or just export at init stage ?
    tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
    if (ierr /= 0) then
      print *,'moab error in setting data in a2x_am '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    end if
#ifdef MOABDEBUG
    outfile = 'AtmPhys_init.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  &
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
#endif
#endif

  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct(EClock, cdata, x2a, a2x)
    use iso_c_binding,  only: c_double
    use scream_f2c_mod, only: scream_run
#ifdef HAVE_MOAB
    use seq_flds_mod, only: seq_flds_a2x_fields, seq_flds_x2a_fields
    use shr_kind_mod     ,  only: cxx=>shr_kind_cxx
    use iMOAB,              only: iMOAB_SetDoubleTagStorage, iMOAB_GetDoubleTagStorage
#ifdef MOABDEBUG
    use iMOAB, only: iMOAB_WriteMesh
#endif
#ifdef MOABCOMP
    use mct_mod
    use seq_comm_mct,     only : num_moab_exports
#endif

    integer :: ent_type
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
    character(CXX) :: mct_field, modelStr
#endif
#endif

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) ,intent(inout) :: EClock     ! clock
    type(seq_cdata)  ,intent(inout) :: cdata
    type(mct_aVect)  ,intent(inout) :: x2a        ! driver     -> atmosphere
    type(mct_aVect)  ,intent(inout) :: a2x        ! atmosphere -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    real(R8)                         :: nextsw_cday    ! calendar of next atm sw
    integer                          :: dt_scream
#ifdef HAVE_MOAB
    integer        :: ierr
    character(CXX) :: tagname ! will store all seq_flds_a2x_fields , seq_flds_x2a_fields
#ifdef MOABDEBUG
    integer                          :: cur_atm_stepno
    character*100                    :: outfile, wopts, lnum
#endif
#endif
    !-------------------------------------------------------------------------------

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Get time step info
    call seq_timemgr_EClockGetData (EClock, next_cday=nextsw_cday, dtime=dt_scream)
#ifdef HAVE_MOAB

#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2a_fields
    call mct_list_init(temp_list ,seq_flds_x2a_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 0 ! entity type is vertex for phys atm
    if (rank2 .eq. 0) print *, num_moab_exports, trim(seq_flds_x2a_fields)
    modelStr ='atm run'
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_atm_moab, x2a, mct_field,  mphaid, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif


    ! import data from moab data structures
    tagname=trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_GetDoubleTagStorage ( mphaid, tagname, totalmbls_r , ent_type, x2a_am )

#endif
    ! Run scream
    call scream_run( dt_scream )

#ifdef HAVE_MOAB
    ! export data to moab data structures
    tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
#ifdef MOABDEBUG
    call seq_timemgr_EClockGetData( EClock, stepno=cur_atm_stepno )
    write(lnum,"(I0.2)")cur_atm_stepno
    outfile = 'AtmPhys_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
#endif
#endif
    ! Set time of next radiadtion computation
    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

  end subroutine atm_run_mct

  !===============================================================================
  subroutine atm_final_mct(EClock, cdata, x2a, a2x)
    use scream_f2c_mod, only: scream_finalize

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock  ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2a     ! driver     -> atmosphere
    type(mct_aVect)             ,intent(inout) :: a2x     ! atmosphere -> driver

    !----------------------------------------------------------------------------
    ! Finish the rest of ATM model
    !----------------------------------------------------------------------------

    call scream_finalize()

  end subroutine atm_final_mct
  !===============================================================================

  subroutine atm_Set_gsMap_mct( mpicom_atm, ATMID, GSMap_atm )
    use iso_c_binding, only: c_int, c_loc
    use scream_f2c_mod, only: scream_get_num_global_cols, scream_get_num_local_cols, &
                              scream_get_local_cols_gids
    use mct_mod,        only: mct_gsMap_init
    !-------------------------------------------------------------------
    !
    ! Inputs
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_atm
    !
    ! Local variables
    !
    integer(kind=c_int), allocatable, target :: col_gids(:)
    integer(kind=c_int) :: num_local_cols, num_global_cols
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP

    ! Determine global seg map
    num_local_cols  = scream_get_num_local_cols()
    num_global_cols = scream_get_num_global_cols()

    allocate(col_gids(num_local_cols))

    call scream_get_local_cols_gids(c_loc(col_gids))

    call mct_gsMap_init( gsMap_atm, col_gids, mpicom_atm, ATMID, &
                         num_local_cols, num_global_cols)

    deallocate(col_gids)

  end subroutine atm_Set_gsMap_mct

  subroutine atm_domain_mct( lsize, gsMap_atm, dom_atm )
    use iso_c_binding,  only: c_loc
    use mct_mod,        only: mct_gGrid_init, mct_gGrid_importIAttr, mct_gGrid_importRAttr, &
                              mct_gsMap_orderedPoints
    use seq_flds_mod,   only: seq_flds_dom_coord, seq_flds_dom_other
    use shr_const_mod,  only: SHR_CONST_PI
    use scream_f2c_mod, only: scream_get_cols_latlon, scream_get_cols_area
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_atm
    type(mct_ggrid), intent(inout):: dom_atm
    !
    ! Local Variables
    !
    real(r8), pointer :: data1(:)
    real(r8), pointer :: data2(:)     ! temporary
    integer , pointer :: idata(:)     ! temporary
    !-------------------------------------------------------------------

    allocate(data1(lsize))
    allocate(data2(lsize))

    ! Initialize mct atm domain
    call mct_gGrid_init( GGrid=dom_atm, CoordChars=trim(seq_flds_dom_coord), &
                         OtherChars=trim(seq_flds_dom_other), lsize=lsize )

    ! Initialize attribute vector with special value
    call mct_gsMap_orderedPoints(gsMap_atm, my_task, idata)
    call mct_gGrid_importIAttr(dom_atm,'GlobGridNum',idata,lsize)

    ! Fill in correct values for domain components
    call scream_get_cols_latlon(c_loc(data1),c_loc(data2))
    call mct_gGrid_importRAttr(dom_atm,"lat",data1,lsize)
    call mct_gGrid_importRAttr(dom_atm,"lon",data2,lsize)

    call scream_get_cols_area(c_loc(data1))
    call mct_gGrid_importRAttr(dom_atm,"area",data1,lsize)

    ! Mask and frac are both exactly 1
    data1 = 1.0
    call mct_gGrid_importRAttr(dom_atm,"mask",data1,lsize)
    call mct_gGrid_importRAttr(dom_atm,"frac",data1,lsize)

    ! Aream is computed by mct, so give invalid initial value
    data1 = -9999.0_R8
    call mct_gGrid_importRAttr(dom_atm,"aream",data1,lsize)
  end subroutine atm_domain_mct

#ifdef HAVE_MOAB
! as part of initialization, moab version of the MCT type grid needs to be instanced too
! it is what we call "point cloud mesh", it contains just the vertices corresponding to the
! grid for surface coupling
! its tags are corresponding to the fields from a2x and x2a attribute vectors, and its
! GLOBAL_ID tag corresponds to the global degrees of freedom in global atmosphere
  subroutine  moab_atm_phys_scream()

    use iMOAB, only : iMOAB_RegisterApplication, iMOAB_CreateVertices, iMOAB_WriteMesh, &
      iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
      iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo
    use iso_c_binding, only: c_int, c_loc
    use scream_f2c_mod, only: scream_get_num_global_cols, scream_get_num_local_cols, &
                              scream_get_local_cols_gids
    use scream_f2c_mod, only: scream_get_cols_latlon, scream_get_cols_area
    use seq_flds_mod, only: seq_flds_dom_fields
    use shr_kind_mod     , only: r8 => shr_kind_r8, cl=>shr_kind_cl, cxx=>shr_kind_cxx
    use shr_const_mod, only: SHR_CONST_PI
    use seq_flds_mod,       only: seq_flds_a2x_fields, seq_flds_x2a_fields

    ! Local variables
    !
    integer(kind=c_int) :: num_local_cols, num_global_cols
#ifdef MOABDEBUG
    character*100 outfile, wopts
#endif
    character(CXX) :: tagname ! will store all seq_flds_a2x_fields
    character*32 appname
    integer :: ATM_PHYS ! used as global identifier for iMOAB app on pphys grid atmosphere (200+ atm id)
    !-------------------------------------------------------------------
    !
    ! Local Variables
    !
    integer        :: ierr,mpi_ierr  ! error codes
    integer        :: i  ! for loops along dofs
    integer        :: tagtype, tagindex, numco, ent_type
    real(r8), pointer :: data1(:)
    real(r8), pointer :: data2(:)     ! temporary
    integer(kind=c_int), allocatable, target :: col_gids(:)
    real(r8), allocatable, target :: moab_vert_coords(:)
    real(r8)       :: latv, lonv
    !-------------------------------------------------------------------

    appname="ATM_PHYS_SCREAM"//C_NULL_CHAR
    ATM_PHYS = 200 + ATM_ID !
    ierr = iMOAB_RegisterApplication(appname, mpicom_atm, ATM_PHYS, mphaid)
    if (ierr > 0 )  then
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    if (my_task == master_task) then
       print *, " moab_atm_phys_scream:: register MOAB app:", trim(appname), "  mphaid=", mphaid
    endif
    ! find global ids
    num_local_cols  = scream_get_num_local_cols()
    num_global_cols = scream_get_num_global_cols()

    allocate(col_gids(num_local_cols))

    call scream_get_local_cols_gids(c_loc(col_gids))

    ! Fill in correct values for domain components
    allocate(moab_vert_coords(num_local_cols*3))
    allocate( data1(num_local_cols), data2(num_local_cols) )
    call scream_get_cols_latlon(c_loc(data1),c_loc(data2))
    do i=1,num_local_cols
          latv = data1(i) * SHR_CONST_PI/180
          lonv = data2(i) * SHR_CONST_PI/180
          moab_vert_coords(3*i-2)=COS(latv)*COS(lonv)
          moab_vert_coords(3*i-1)=COS(latv)*SIN(lonv)
          moab_vert_coords(3*i  )=SIN(latv)
    enddo
    ierr = iMOAB_CreateVertices(mphaid, num_local_cols*3, 3, moab_vert_coords)
    if (ierr > 0 )  then
      print *, "Error: fail to create MOAB vertices in phys atm model"
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

    tagtype = 0  ! dense, integer
    numco = 1
    tagname='GLOBAL_ID'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  then
      print *, "Error: fail to define GLOBAL_ID tag in phys atm model"
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    ent_type = 0 ! vertex type
    ierr = iMOAB_SetIntTagStorage(mphaid, tagname, num_local_cols , ent_type, col_gids)
    if (ierr > 0 )  then
      print *, "Error: fail to set GLOBAL_ID tag in phys atm model"
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    ierr = iMOAB_ResolveSharedEntities( mphaid, num_local_cols, col_gids)
    if (ierr > 0 )  then
      print *, "Error: fail to resolve shared ents in phys atm model"
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    ierr = iMOAB_UpdateMeshInfo(mphaid)
    if (ierr > 0 )  then
      print *, "Error: fail to update mesh info in phys atm model"
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    tagname=trim(seq_flds_dom_fields)//C_NULL_CHAR !  mask is double too lat:lon:hgt:area:aream:mask:frac
    tagtype = 1 ! dense, double
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  then
      print *, 'Error: fail to create  tags from seq_flds_dom_fields '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

    call scream_get_cols_area(c_loc(data1))
    tagname = 'area'//C_NULL_CHAR !
    ierr = iMOAB_SetDoubleTagStorage(mphaid, tagname, num_local_cols , ent_type, data1)
    if (ierr > 0 )  then
      print *, 'Error: fail to set area '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

    ! Mask and frac are both exactly 1
    data1 = 1.0
    tagname = 'mask'//C_NULL_CHAR !
    ierr = iMOAB_SetDoubleTagStorage(mphaid, tagname, num_local_cols , ent_type, data1)
    if (ierr > 0 )  then
      print *, 'Error: fail to set mask '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    tagname = 'frac'//C_NULL_CHAR !
    ierr = iMOAB_SetDoubleTagStorage(mphaid, tagname, num_local_cols , ent_type, data1)
    if (ierr > 0 )  then
      print *, 'Error: fail to set frac '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

    !call mct_gGrid_importRAttr(dom_atm,"mask",data1,lsize)
    !call mct_gGrid_importRAttr(dom_atm,"frac",data1,lsize)

    ! Aream is computed by mct, so give invalid initial value
    data1 = -9999.0_R8
    tagname = 'aream'//C_NULL_CHAR !
    ierr = iMOAB_SetDoubleTagStorage(mphaid, tagname, num_local_cols , ent_type, data1)
    if (ierr > 0 )  then
      print *, 'Error: fail to set aream '
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

#ifdef MOABDEBUG
    outfile = 'AtmPhys.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  then
      print *, "Error: fail to write PhysAtm mesh "
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
#endif
  ! define fields seq_flds_a2x_fields
    tagtype = 1  ! dense, double
    numco = 1 !  one value per vertex / entity
    tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  then
      print *, "Error: fail to define seq_flds_a2x_fields "
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    ! make sure this is defined too; it could have the same fields, but in different order, or really different
    ! fields; need to make sure we have them
    tagname = trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  then
      print *, "Error: fail to define seq_flds_x2a_fields "
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif
    deallocate(col_gids)
    deallocate(data1)
    deallocate(data2)
  end subroutine  moab_atm_phys_scream
#endif

end module atm_comp_mct
