module glc_comp_esmf

#ifdef ESMF_INTERFACE  
! !USES:

  use shr_sys_mod
  use shr_kind_mod,        only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod,        only: CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod,        only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod,        only: shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod,        only: shr_file_freeunit

  use esmf
  use esmfshr_mod

  use seq_infodata_mod,    only:  seq_infodata_start_type_start, seq_infodata_start_type_cont
  use seq_infodata_mod,    only:  seq_infodata_start_type_brnch
  use seq_timemgr_mod

  use glc_import_export
  use glc_cpl_indices
  use glc_constants,       only: verbose, stdout, stderr, nml_in, radius
  use glc_errormod,        only: glc_success
  use glc_InitMod,         only: glc_initialize
  use glc_RunMod,          only: glc_run
  use glc_FinalMod,        only: glc_final
  use glc_io,              only: glc_io_write_restart
  use glc_communicate,     only: init_communicate, my_task, master_task
  use glc_time_management, only: iyear,imonth,iday,ihour,iminute,isecond,runtype
  use glc_fields,          only: ice_sheet

  implicit none
  SAVE
  private                              ! By default make data private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: glc_register_esmf
  public :: glc_init_esmf
  public :: glc_run_esmf
  public :: glc_final_esmf

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: glc_distgrid_esmf
  private :: glc_domain_esmf

  !--------------------------------------------------------------------------
  ! Private module data interfaces
  !--------------------------------------------------------------------------

  !--- stdin input stuff ---
  character(CS) :: str                  ! cpp  defined model name
  
  !--- other ---
  integer(IN)   :: errorcode            ! glc error code
  
  ! my_task_local and master_task_local are needed for some checks that are done before
  ! init_communicate is called (although, it's possible that init_communicate could be
  ! moved to earlier to prevent the need for these copies)
  integer(IN)   :: my_task_local         ! my task in mpi communicator mpicom 
  integer(IN)   :: master_task_local=0   ! task number of master task


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!====================================================================================

  subroutine glc_register_esmf(comp, rc)
    implicit none
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    print *, "In glc register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      glc_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         glc_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      glc_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine glc_register_esmf

!====================================================================================

  subroutine glc_init_esmf(comp, import_state, export_state, EClock, rc)

    use glc_ensemble       , only : set_inst_vars, write_inst_vars, get_inst_name
    use glc_files          , only : set_filenames, ionml_filename
    use glc_coupling_flags , only : has_ocn_coupling, has_ice_coupling
    use glc_indexing_info  , only : nx_tot, ny_tot, npts_tot
    
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initialize glc model
    !
    ! !INPUT/OUTPUT PARAMETERS:
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !
    ! !LOCAL VARIABLES:
    type(ESMF_DistGrid)      :: distgrid
    type(ESMF_Array)         :: dom, g2x, x2g
    type(ESMF_VM)            :: vm
    integer(IN)              :: ierr 
    integer(IN)              :: i,j,n
    integer(IN)              :: shrlogunit, shrloglev  
    character(CL)            :: starttype
    real(R8), pointer        :: fptr(:,:)
    integer                  :: mpicom_loc, mpicom_vm
    character(ESMF_MAXSTR)   :: convCIM, purpComp
    integer(IN)              :: COMPID
    character(CS)            :: myModelName

    !--- formats ---
    character(*), parameter :: F00   = "('(glc_init_esmf) ',8a)"
    character(*), parameter :: F01   = "('(glc_init_esmf) ',a,8i8)"
    character(*), parameter :: F02   = "('(glc_init_esmf) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(glc_init_esmf) ',a,i8,a)"
    character(*), parameter :: F90   = "('(glc_init_esmf) ',73('='))"
    character(*), parameter :: F91   = "('(glc_init_esmf) ',73('-'))"
    character(*), parameter :: subName = "(glc_init_esmf) "
    !-----------------------------------------------------------------------

    ! Determine attribute vector indices

    call glc_cpl_indices_set()

    rc = ESMF_SUCCESS
    
    ! duplicate the mpi communicator from the current VM 
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call MPI_Comm_dup(mpicom_vm, mpicom_loc, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Get id of this task
    call MPI_Comm_rank(mpicom_loc, my_task_local, ierr)
    
    ! Initialize glc id
    
    call ESMF_AttributeGet(export_state, name="ID", value=COMPID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !---------------------------------------------------------------------------
    ! set variables that depend on ensemble index
    !---------------------------------------------------------------------------

    call set_inst_vars(COMPID)
    call get_inst_name(myModelName)
    call set_filenames()

    !---------------------------------------------------------------------------
    ! determine type of run
    !---------------------------------------------------------------------------

    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(*,*) 'glc_comp_esmf ERROR: unknown starttype'
       call shr_sys_abort()
    end if

    !----------------------------------------------------------------------------
    ! Initialize glc
    !----------------------------------------------------------------------------

    if (my_task_local == master_task_local) then
       stdout = shr_file_getUnit()
       call shr_file_setIO(ionml_filename,stdout)
    else
       stdout = 6
    endif
    stderr = stdout
    nml_in = shr_file_getUnit()

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    errorCode = glc_Success
    if (verbose .and. my_task_local == master_task_local) then
       write(stdout,F00) ' Starting'
       write(stdout,*) subname, 'COMPID: ', COMPID
       call write_inst_vars
       call shr_sys_flush(stdout)
    endif
    call init_communicate(mpicom_loc)

    call glc_initialize(EClock, errorCode)

    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' GLC Initial Date ',iyear,imonth,iday,ihour,iminute,isecond
       write(stdout,F01) ' Initialize Done', errorCode
       call shr_sys_flush(stdout)
    endif

    !---------------------------------------------------------------------------
    ! Initialize distgrids, domains, and arrays
    !---------------------------------------------------------------------------

    ! Initialize glc distgrid

    distgrid = glc_distgrid_esmf(rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="gsize", value=npts_tot, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize glc domain (needs glc initialization info)

    dom = mct2esmf_init(distgrid, attname=seq_flds_dom_fields, name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call glc_domain_esmf(dom, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Inialize input/output arrays

    g2x = mct2esmf_init(distgrid, attname=seq_flds_g2x_fields, name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    x2g = mct2esmf_init(distgrid, attname=seq_flds_x2g_fields, name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/g2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2g/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !---------------------------------------------------------------------------
    ! send initial state to driver
    !---------------------------------------------------------------------------

    call ESMF_ArrayGet(g2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call glc_export(fptr)
    
    call ESMF_AttributeSet(export_state, name="glc_present", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glclnd_present", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glcocn_present", &
         value=has_ocn_coupling(), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glcice_present", &
         value=has_ice_coupling(), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glc_prognostic", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glc_nx", value=nx_tot, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glc_ny", value=ny_tot, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "GLC", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "LongName", &
                           "TBD", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "Description", &
                           "TBD", &

                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ModelType", "GlC", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    !    call ESMF_AttributeSet(comp, "Name", "someone", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "EmailAddress", &
    !                           "someone@someplace", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
#endif
    
    if (my_task == master_task) then
       write(stdout,F91) 
       write(stdout,F00) trim(myModelName),': start of main integration loop'
       write(stdout,F91) 
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine glc_init_esmf

!====================================================================================

  subroutine glc_run_esmf(comp, import_state, export_state, EClock, rc)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Run GLC
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !
    ! !LOCAL VARIABLES:
    integer(IN)   :: cesmYMD           ! cesm model date
    integer(IN)   :: cesmTOD           ! cesm model sec
    integer(IN)   :: glcYMD            ! glc model date
    integer(IN)   :: glcTOD            ! glc model sec 
    logical       :: stop_alarm        ! is it time to stop
    logical       :: rest_alarm        ! is it time to write a restart
    logical       :: done              ! time loop logical
    integer(IN)   :: shrlogunit, shrloglev  
    real(R8), pointer :: fptr(:,:)
    type(ESMF_Array)  :: x2g, g2x 
    character(*), parameter :: F00   = "('(glc_run_esmf) ',8a)"
    character(*), parameter :: F01   = "('(glc_run_esmf) ',a,8i8)"
    character(*), parameter :: F04   = "('(glc_run_esmf) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(glc_run_esmf) "
    !---------------------------------------------------------------------------

    ! Reset shr logging to my log file

    rc = ESMF_SUCCESS

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    ! Set internal time info
 
    errorCode = glc_Success
    call seq_timemgr_EClockGetData(EClock,curr_ymd=cesmYMD, curr_tod=cesmTOD)
    stop_alarm = seq_timemgr_StopAlarmIsOn( EClock )

    glcYMD = iyear*10000 + imonth*100 + iday
    glcTOD = ihour*3600 + iminute*60 + isecond
    done = .false.
    if (glcYMD == cesmYMD .and. glcTOD == cesmTOD) done = .true.
    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' Run Starting ',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    endif

    ! Unpack import state
    
    call ESMF_StateGet(import_state, itemName="x2d", array=x2g, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(x2g, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call glc_import(fptr)

    ! Run 

    do while (.not. done) 
       if (glcYMD > cesmYMD .or. (glcYMD == cesmYMD .and. glcTOD > cesmTOD)) then
          write(stdout,*) subname,' ERROR overshot coupling time ',glcYMD,glcTOD,cesmYMD,cesmTOD
          call shr_sys_abort('glc error overshot time')
       endif

       call glc_run(EClock)

       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       if (glcYMD == cesmYMD .and. glcTOD == cesmTOD) done = .true.
       if (verbose .and. my_task == master_task) then
          write(stdout,F01) ' GLC  Date ',glcYMD,glcTOD
       endif
    enddo

    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' Run Done',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    endif
    
    ! Pack export state

    call ESMF_StateGet(export_state, itemName="d2x", array=g2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(g2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call glc_export(fptr)
    
    ! Log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=cesmYMD, curr_tod=cesmTOD)
       write(stdout,F01) ' CESM Date ', cesmYMD,cesmTOD
       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       write(stdout,F01) ' GLC  Date ',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    end if

    ! If time to write restart, do so

    rest_alarm = seq_timemgr_RestartAlarmIsOn( EClock )
    if (rest_alarm) then
       ! TODO loop over instances
       call glc_io_write_restart(ice_sheet%instances(1), EClock)
    endif

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(stdout)
    
  end subroutine glc_run_esmf

!====================================================================================

  subroutine glc_final_esmf(comp, import_state, export_state, EClock, rc)

    use glc_ensemble, only : get_inst_name

    !------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Finalize GLC
    !
    ! !ARGUMENTS:
    !
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

    integer(IN)             :: shrlogunit, shrloglev  
    character(CS)           :: myModelName

    !--- formats ---
    character(*), parameter :: F00   = "('(glc_final_mct) ',8a)"
    character(*), parameter :: F01   = "('(glc_final_mct) ',a,8i8)"
    character(*), parameter :: F91   = "('(glc_final_mct) ',73('-'))"
    character(*), parameter :: subName = "(glc_final_mct) "
    !---------------------------------------------------------------------------
    
    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout)

    call get_inst_name(myModelName)

    if (my_task == master_task) then
       write(stdout,F91) 
       write(stdout,F00) trim(myModelName),': end of main integration loop'
       write(stdout,F91) 
    end if
      
    errorCode = glc_Success

    call glc_final(errorCode)

    ! Note that restart for final timestep was written in run phase.
    rc = ESMF_SUCCESS
    
    ! Destroy ESMF objects
    
    call esmfshr_util_StateArrayDestroy(export_state,"d2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call esmfshr_util_StateArrayDestroy(import_state,"x2d",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' Done',errorCode
       call shr_sys_flush(stdout)
    endif

    ! Reset shr logging to original values 

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(stdout)
    
  end subroutine glc_final_esmf
  
!=================================================================================

  function glc_distgrid_esmf(rc)

    ! Initialize global index space array
    
    use glc_broadcast, only: broadcast_scalar
    use glc_indexing_info, only : local_indices, global_indices, nx, ny, npts

    !-------------------------------------------------------------------
    ! Arguments
    implicit none
    integer, intent(out):: rc

    ! Return:
    type(ESMF_DistGrid) :: glc_DistGrid_esmf  ! Resulting distributed grid

    ! Local Variables
    integer,allocatable :: gindex(:)
    integer :: i, j, n
    integer :: ier

    !--- formats ---
    character(*), parameter :: F02   = "('(glc_DistGrid_esmf) ',a,4es13.6)"
    character(*), parameter :: subName = "(glc_DistGrid_esmf) "
    !-------------------------------------------------------------------

    allocate(gindex(npts))
    do j = 1,ny
       do i = 1,nx
          n = local_indices(i,j)
          gindex(n) = global_indices(i,j)
       enddo
    enddo
       
    glc_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function glc_DistGrid_esmf

!=======================================================================

  subroutine glc_domain_esmf( dom, rc )

    !-------------------------------------------------------------------
    use glc_indexing_info, only : nx, ny, local_indices
    use glad_main, only : glad_get_lat_lon, glad_get_areas
    
    implicit none
    type(ESMF_Array), intent(inout)     :: dom
    integer, intent(out)                :: rc

    ! Local Variables
    integer :: j,i,n
    integer :: klon,klat,karea,kmask,kfrac ! domain fields
    real(R8), pointer :: fptr(:,:)
    real(r8), allocatable :: lats(:,:)  ! latitude of each point (degrees)
    real(r8), allocatable :: lons(:,:)  ! longitude of each point (degrees)
    real(r8), allocatable :: areas(:,:) ! area of each point (square meters)
    !-------------------------------------------------------------------

    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2
    ! 
    rc = ESMF_SUCCESS

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components
    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components

    allocate(lats(nx, ny))
    allocate(lons(nx, ny))
    allocate(areas(nx, ny))
    
    ! TODO(wjs, 2015-04-02) The following may need a loop over instances
    call glad_get_lat_lon(ice_sheet, instance_index = 1, &
         lats = lats, lons = lons)
    call glad_get_areas(ice_sheet, instance_index = 1, areas = areas)
    
    fptr(:,:) = -9999.0_R8
    fptr(kmask,:) = -0.0_R8
    do j = 1,ny
       do i = 1,nx
          n = local_indices(i,j)
          fptr(klon , n) = lons(i,j)
          fptr(klat , n) = lats(i,j)

          ! convert from m^2 to radians^2
          fptr(karea, n) = areas(i,j)/(radius*radius)

          ! For now, assume mask and frac are 1 everywhere. This may need to be changed
          ! in the future.
          fptr(kmask, n) = 1._r8
          fptr(kfrac, n) = 1._r8
       end do
    end do

    deallocate(lats)
    deallocate(lons)
    deallocate(areas)

  end subroutine glc_domain_esmf

#endif

end module glc_comp_esmf

