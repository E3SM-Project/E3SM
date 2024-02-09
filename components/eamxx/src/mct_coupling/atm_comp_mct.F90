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
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001')
  character(len=16)      :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "')
  integer(IN)            :: ATM_ID              ! mct comp id
  integer(IN),parameter  :: master_task=0       ! task number of master task

  integer                :: atm_log_unit

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    use seq_infodata_mod,   only: seq_infodata_start_type_start, seq_infodata_start_type_cont
    use seq_comm_mct,       only: seq_comm_inst, seq_comm_name, seq_comm_suffix
    use shr_file_mod,       only: shr_file_getunit, shr_file_setIO
    use shr_sys_mod,        only: shr_sys_abort

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

    ! TODO: read this from the namelist?
    character(len=256)                :: yaml_fname = "./data/scream_input.yaml"
    character(kind=c_char,len=256), target :: yaml_fname_c, atm_log_fname_c
    character(len=256) :: caseid, username, hostname
    character(kind=c_char,len=256), target :: caseid_c, username_c, hostname_c, calendar_c
    logical (kind=c_bool) :: restarted_run

    !-------------------------------------------------------------------------------

print *,'OG a 1'

    ! Grab some data from the cdata structure (coming from the coupler)
    call seq_cdata_setptrs(cdata, &
         id=ATM_ID, &
         mpicom=mpicom_atm, &
         gsMap=gsmap_atm, &
         dom=dom_atm, &
         infodata=infodata)
print *, 'OG a 2'
    call seq_infodata_getData(infodata, atm_phase=phase, start_type=run_type, &
                              username=username, case_name=caseid, hostname=hostname)
print *, 'OG a 3'
            call seq_infodata_PutData(infodata, atm_aero=.true.)
print *, 'OG a 4'
    call seq_infodata_PutData(infodata, atm_prognostic=.true.)

print *, 'OG a 5'
    if (phase > 1) RETURN

print *, 'OG a 6'
    ! Determine instance information
    inst_name   = seq_comm_name(ATM_ID)
    inst_index  = seq_comm_inst(ATM_ID)
    inst_suffix = seq_comm_suffix(ATM_ID)

print *, 'OG a 7'
    ! Determine communicator group
    call mpi_comm_rank(mpicom_atm, my_task, ierr)

print *, 'OG a 8'
    !----------------------------------------------------------------------------
    ! Init atm.log
    !----------------------------------------------------------------------------
print *, 'OG a 9'
    if (my_task == master_task) then
print *, 'OG a 10'
      atm_log_unit = shr_file_getUnit()
      call shr_file_setIO ('atm_modelio.nml'//trim(inst_suffix),atm_log_unit)
      inquire(unit=atm_log_unit,name=atm_log_fname)
    endif

print *, 'OG a 11'
    call mpi_bcast(atm_log_unit,1,MPI_INTEGER,master_task,mpicom_atm,mpi_ierr)
print *, 'OG a 12'
    if (ierr /= 0) then
print *, 'OG a 13'
      print *,'[eamxx] ERROR broadcasting atm.log unit'
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    end if

print *, 'OG a 14'
    call mpi_bcast(atm_log_fname,256,MPI_CHARACTER,master_task,mpicom_atm,ierr)
print *, 'OG a 15'
    if (ierr /= 0) then
      print *,'[eamxx] ERROR broadcasting atm.log file name'
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    end if

    !----------------------------------------------------------------------------
    ! Initialize atm
    !----------------------------------------------------------------------------

    ! Init the AD
print *, 'OG a 16'
    call seq_timemgr_EClockGetData(EClock, calendar=calendar, &
                                   curr_ymd=cur_ymd, curr_tod=cur_tod, &
                                   start_ymd=case_start_ymd, start_tod=case_start_tod)
print *, 'OG a 17'
    call string_f2c(yaml_fname,yaml_fname_c)
print *, 'OG a 18'
    call string_f2c(calendar,calendar_c)
print *, 'OG a 19'
    call string_f2c(trim(atm_log_fname),atm_log_fname_c)
print *, 'OG a 20'
    call scream_create_atm_instance (mpicom_atm, ATM_ID, yaml_fname_c, atm_log_fname_c, &
                          INT(cur_ymd,kind=C_INT),  INT(cur_tod,kind=C_INT), &
                          INT(case_start_ymd,kind=C_INT), INT(case_start_tod,kind=C_INT), &
                          calendar_c)


print *, 'OG a 21'
    ! Init MCT gsMap
    call atm_Set_gsMap_mct (mpicom_atm, ATM_ID, gsMap_atm)
print *, 'OG a 22'
    lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

print *, 'OG a 23'
    ! Init MCT domain structure
    call atm_domain_mct (lsize, gsMap_atm, dom_atm)

print *, 'OG a 24'
    ! Init import/export mct attribute vectors
    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
print *, 'OG a 25'
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)

print *, 'OG a 26'
    ! Complete AD initialization based on run type
    if (trim(run_type) == trim(seq_infodata_start_type_start)) then
      restarted_run = .false.
    else if (trim(run_type) == trim(seq_infodata_start_type_cont) ) then
      restarted_run = .true.
    else
      print *, "[eamxx] ERROR! Unsupported starttype: "//trim(run_type)
      call mpi_abort(mpicom_atm,ierr,mpi_ierr)
    endif

print *, 'OG a 27'
    ! Init surface coupling stuff in the AD
    call scream_set_cpl_indices (x2a, a2x)
print *, 'OG a 28'

    call scream_setup_surface_coupling (c_loc(import_field_names), c_loc(import_cpl_indices), &
                                        c_loc(x2a%rAttr), c_loc(import_vector_components), &
                                        c_loc(import_constant_multiple), c_loc(do_import_during_init), &
                                        num_cpl_imports, num_scream_imports, import_field_size, &
                                        c_loc(export_field_names), c_loc(export_cpl_indices), &
                                        c_loc(a2x%rAttr), c_loc(export_vector_components), &
                                        c_loc(export_constant_multiple), c_loc(do_export_during_init), &
                                        num_cpl_exports, num_scream_exports, export_field_size)

print *, 'OG a 29'
    call string_f2c(trim(caseid),caseid_c)
    call string_f2c(trim(username),username_c)
    call string_f2c(trim(hostname),hostname_c)
    call scream_init_atm (caseid_c,hostname_c,username_c)

print *, 'OG a 30'
  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct(EClock, cdata, x2a, a2x)
    use iso_c_binding,  only: c_double
    use scream_f2c_mod, only: scream_run

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

    !-------------------------------------------------------------------------------

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

end module atm_comp_mct
