module atm_comp_mct

  ! Modules used acros atm_xyz_mct routines

  ! shr mods
  use shr_kind_mod,     only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8!, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL

  ! seq mods
  use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata
  use seq_timemgr_mod,  only: seq_timemgr_EClockGetData

  ! shr mods
  use shr_file_mod,     only: shr_file_getlogunit, shr_file_setlogunit, shr_file_setio, &
                              shr_file_getloglevel, shr_file_setloglevel
  use shr_sys_mod,      only: shr_sys_flush

  ! toolkits mods
  use esmf,             only: ESMF_Clock
  use mct_mod,          only: mct_aVect, mct_gsMap, mct_gGrid

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

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
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: ATM_ID              ! mct comp id
  real(r8) ,  pointer    :: gbuf(:,:)           ! model grid
  integer(IN),parameter  :: master_task=0       ! task number of master task

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine atm_init_mct( EClock, cdata, x2a, a2x, NLFilename )
    use iso_c_binding,      only: c_ptr, c_loc, c_int, c_char
    use scream_f2c_mod,     only: scream_init, scream_setup_surface_coupling
    use scream_cpl_indices, only: scream_set_cpl_indices, num_exports, num_imports, &
                                  scr_names_x2a, scr_names_a2x, index_x2a, index_a2x
    use ekat_string_utils,  only: string_f2c

    use mct_mod,        only: mct_aVect_init, mct_gsMap_lsize
    use seq_flds_mod,   only: seq_flds_a2x_fields, seq_flds_x2a_fields
    use seq_comm_mct,   only: seq_comm_inst, seq_comm_name, seq_comm_suffix
    use shr_file_mod,   only: shr_file_getunit

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local variables ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap_atm
    type(mct_gGrid)        , pointer :: dom_atm
    integer(IN)                      :: shrlogunit     ! original log unit
    integer(IN)                      :: shrloglev      ! original log level
    integer(IN)                      :: nxg            ! global dim i-direction
    integer(IN)                      :: nyg            ! global dim j-direction
    integer(IN)                      :: phase          ! initialization phase
    integer(IN)                      :: ierr           ! error code
    logical                          :: atm_present    ! if true, component is present
    logical                          :: atm_prognostic ! if true, component is prognostic
    integer (IN)                     :: start_tod, start_ymd
    integer                          :: lsize
    type(c_ptr) :: x2a_ptr, a2x_ptr

    ! TODO: read this from the namelist?
    character(len=256)                :: yaml_fname = "data/scream_input.yaml"
    character(kind=c_char,len=256), target :: yaml_fname_c
    !-------------------------------------------------------------------------------

    ! Grab some data from the cdata structure (coming from the coupler)
    call seq_cdata_setptrs(cdata, &
         id=ATM_ID, &
         mpicom=mpicom_atm, &
         gsMap=gsmap_atm, &
         dom=dom_atm, &
         infodata=infodata)
    call seq_infodata_getData(infodata, atm_phase=phase)

    if (phase > 1) RETURN

    ! Determine instance information
    inst_name   = seq_comm_name(ATM_ID)
    inst_index  = seq_comm_inst(ATM_ID)
    inst_suffix = seq_comm_suffix(ATM_ID)

    if (phase == 1) then
       ! Determine communicator group
       call mpi_comm_rank(mpicom_atm, my_task, ierr)

       !--- open log file ---
       if (my_task == master_task) then
          logUnit = shr_file_getUnit()
          call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix),logUnit)
       else
          logUnit = 6
       endif
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Initialize atm
    !----------------------------------------------------------------------------

    ! Init the AD
    call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)
    call string_f2c(yaml_fname,yaml_fname_c)
    call scream_init (mpicom_atm, INT(start_ymd,kind=C_INT), INT(start_tod,kind=C_INT), yaml_fname_c)

    ! Init MCT gsMap
    call atm_Set_gsMap_mct (mpicom_atm, ATM_ID, gsMap_atm)
    lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

    ! Init MCT domain structure
    call atm_domain_mct (lsize, gsMap_atm, dom_atm)

    ! Init import/export mct attribute vectors
    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)

    ! Init surface coupling stuff in the AD
    call scream_set_cpl_indices (x2a, a2x)
    call scream_setup_surface_coupling (c_loc(scr_names_x2a), c_loc(index_x2a), c_loc(x2a%rAttr), num_imports, &
                                        c_loc(scr_names_a2x), c_loc(index_a2x), c_loc(a2x%rAttr), num_exports)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit)

  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct(EClock, cdata, x2d, d2x)
    use iso_c_binding,  only: c_double
    use scream_f2c_mod, only: scream_run

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) ,intent(inout) :: EClock     ! clock
    type(seq_cdata)  ,intent(inout) :: cdata
    type(mct_aVect)  ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)  ,intent(inout) :: d2x        ! dead   -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit     ! original log unit
    integer(IN)                      :: shrloglev      ! original log level
    real(R8)                         :: nextsw_cday    ! calendar of next atm sw
    real(kind=c_double)              :: dt_scream
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Get time step info
    call seq_timemgr_EClockGetData (EClock, next_cday=nextsw_cday, dtime=dt_scream)

    ! Run scream
    call scream_run( dt_scream )

    ! Set time of next radiadtion computation
    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine atm_run_mct

  !===============================================================================
  subroutine atm_final_mct(EClock, cdata, x2d, d2x)
    use scream_f2c_mod, only: scream_finalize

    ! !DESCRIPTION: finalize method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !-------------------------------------------------------------------------------

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
    integer  :: n,i,c,ncols           ! indices	
    real(r8), pointer :: data1(:)
    real(r8), pointer :: data2(:)     ! temporary
    integer , pointer :: idata(:)     ! temporary

    real(r8), parameter:: rad2deg = 180.0_r8 / SHR_CONST_PI
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
    data1(:) = data1(:) * rad2deg
    data2(:) = data2(:) * rad2deg

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

end module atm_comp_mct
