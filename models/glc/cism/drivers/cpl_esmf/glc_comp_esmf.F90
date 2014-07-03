module glc_comp_esmf

! !USES:

  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                               CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use mct_mod
  use esmf
  use esmfshr_mod

  use seq_flds_mod
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use glc_cpl_indices
  use glc_constants,       only: verbose, stdout, stderr, nml_in, &
                                 radius,  radian, tkfrz,  glc_nec
  use glc_errormod,        only: glc_success
  use glc_InitMod,         only: glc_initialize
  use glc_RunMod,          only: glc_run
  use glc_FinalMod,        only: glc_final
  use glc_io,              only: glc_io_write_restart, glc_io_write_history
  use glc_communicate,     only: init_communicate, my_task, master_task
  use glc_time_management, only: iyear,imonth,iday,ihour,iminute,isecond,runtype
  use glc_global_fields,   only: ice_sheet, &
                                 tsfc, topo, qsmb, &                 ! from coupler
                                 gfrac, gtopo, grofi, grofl, ghflx   ! to coupler
  use glc_global_grid,     only: glc_grid, glc_landmask, glc_landfrac

! !PUBLIC TYPES:
  implicit none

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: glc_register_esmf
  public :: glc_init_esmf
  public :: glc_run_esmf
  public :: glc_final_esmf
  SAVE
  private                              ! By default make data private

!--------------------------------------------------------------------------
! Private data interfaces
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

  private :: glc_export_esmf
  private :: glc_import_esmf
  private :: glc_DistGrid_esmf
  private :: glc_domain_esmf

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

    use glc_ensemble, only : set_inst_vars, write_inst_vars, get_inst_name
    use glc_files   , only : set_filenames, ionml_filename

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
    integer(IN)              :: nxg_tot, nyg_tot  ! total nx & ny points across all tasks
    integer(IN)              :: lsize
    integer(IN)              :: shrlogunit, shrloglev  
    character(CL)            :: starttype
    real(R8), pointer        :: fptr(:,:)
    integer                  :: mpicom_loc, mpicom_vm, gsize
    integer                  :: num 
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

    call glc_initialize(errorCode)

    if (verbose .and. my_task == master_task) then
       write(stdout,F01) ' GLC Initial Date ',iyear,imonth,iday,ihour,iminute,isecond
       write(stdout,F01) ' Initialize Done', errorCode
       call shr_sys_flush(stdout)
    endif

    !---------------------------------------------------------------------------
    ! Initialize distgrids, domains, and arrays
    !---------------------------------------------------------------------------

    ! Initialize glc distgrid

    distgrid = glc_DistGrid_esmf(gsize, nxg_tot, nyg_tot, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize glc domain (needs glc initialization info)

    dom = mct2esmf_init(distgrid, attname=seq_flds_dom_fields, name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call glc_domain_esmf(dom, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Inialize input/output arrays
    g2x = mct2esmf_init(distgrid, attname=seq_flds_g2x_fields, name="g2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    x2g = mct2esmf_init(distgrid, attname=seq_flds_x2g_fields, name="x2g", rc=rc)
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

    ! Send initial state to driver
!    do num = 1,glc_nec
!       call glc_export_esmf(fptr, num, &
!            index_g2x_Sg_frac(num), index_g2x_Sg_topo(num)  ,&
!            index_g2x_Fsgg_rofi(num), index_g2x_Fsgg_rofl(num), &
!            index_g2x_Fsgg_hflx(num))
!    end do
    
    call ESMF_AttributeSet(export_state, name="glc_prognostic", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glc_nx", value=nxg_tot, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="glc_ny", value=nyg_tot, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

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
    logical       :: hist_alarm        ! is it time to write a history file
    logical       :: rest_alarm        ! is it time to write a restart
    logical       :: done              ! time loop logical
    integer(IN)   :: shrlogunit, shrloglev  
    integer       :: num 
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
    
    call ESMF_StateGet(import_state, itemName="x2g", array=x2g, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(x2g, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    do num = 1,glc_nec
       call glc_import_esmf(fptr, num, &
            index_x2g_Ss_tsrf(num), index_x2g_Ss_topo(num), index_x2g_Fgss_qice(num))
    end do

    ! Run 

    do while (.not. done) 
       if (glcYMD > cesmYMD .or. (glcYMD == cesmYMD .and. glcTOD > cesmTOD)) then
          write(stdout,*) subname,' ERROR overshot coupling time ',glcYMD,glcTOD,cesmYMD,cesmTOD
          call shr_sys_abort('glc error overshot time')
       endif

       call glc_run

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

    call ESMF_StateGet(export_state, itemName="g2x", array=g2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(g2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    do num = 1,glc_nec
       call glc_export_esmf(fptr, num, &
            index_g2x_Sg_frac(num), index_g2x_Sg_topo(num)  ,&
            index_g2x_Fsgg_rofi(num), index_g2x_Fsgg_rofl(num), index_g2x_Fsgg_hflx(num))
    end do
    

    ! Log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=cesmYMD, curr_tod=cesmTOD)
       write(stdout,F01) ' CESM Date ', cesmYMD,cesmTOD
       glcYMD = iyear*10000 + imonth*100 + iday
       glcTOD = ihour*3600 + iminute*60 + isecond
       write(stdout,F01) ' GLC  Date ',glcYMD,glcTOD
       call shr_sys_flush(stdout)
    end if

    ! If time to write history, do so

    hist_alarm = seq_timemgr_HistoryAlarmIsOn( EClock )
    if (hist_alarm) then
       ! TODO loop over instances
       call glc_io_write_history(ice_sheet%instances(1), EClock)
    endif

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
    
    call esmfshr_util_StateArrayDestroy(export_state,"g2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call esmfshr_util_StateArrayDestroy(import_state,"x2g",rc)
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

  function glc_DistGrid_esmf(gsize, nxg_tot, nyg_tot, rc)

  use glc_broadcast, only: broadcast_scalar

  !-------------------------------------------------------------------
  ! Arguments
  implicit none
  integer, intent(out):: gsize    ! global total number of points
  integer, intent(out):: nxg_tot  ! total nx points across all tasks
  integer, intent(out):: nyg_tot  ! total ny points across all tasks
  integer, intent(out):: rc

  ! Return:
  type(ESMF_DistGrid) :: glc_DistGrid_esmf  ! Resulting distributed grid
  
  ! Local Variables
  integer,allocatable :: gindex(:)
  integer :: i, j, n
  integer :: nxg, nyg, lsize      ! number of points that this task is responsible for
  integer :: ier
  
  !--- formats ---
  character(*), parameter :: F02   = "('(glc_DistGrid_esmf) ',a,4es13.6)"
  character(*), parameter :: subName = "(glc_DistGrid_esmf) "
  !-------------------------------------------------------------------
  
  ! Note that the following assumes that the master task is responsible for all points

  if (my_task == master_task) then

     nxg = glc_grid%nx
     nyg = glc_grid%ny
     lsize = nxg*nyg

     ! Initialize global index space array (the simple method used here only works
     ! because the master task is responsible for all points)
     allocate(gindex(lsize))
     do j = 1,nyg
        do i = 1,nxg
           n = (j-1)*nxg + i
           gindex(n) = n
        enddo
     enddo

  else
     
     nxg = 0
     nyg = 0
     lsize = 0
     allocate(gindex(lsize))

  end if
  
  glc_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  
  deallocate(gindex)


  ! Determine other output arguments, relating to the total grid size (across all tasks)
  ! The following assumes that the master task is responsible for all points

  if (my_task == master_task) then
     nxg_tot = nxg
     nyg_tot = nyg
     gsize = nxg_tot * nyg_tot
  end if

  call broadcast_scalar(nxg_tot, master_task)
  call broadcast_scalar(nyg_tot, master_task)
  call broadcast_scalar(gsize  , master_task)

end function glc_DistGrid_esmf

!====================================================================================

subroutine glc_import_esmf(fptr, ndx, index_tsrf, index_topo, index_qice)

    !-----------------------------------------------------
    implicit none
    real(R8)   , pointer    :: fptr(:,:)
    integer(IN), intent(in) :: ndx                      ! elevation class
    integer(IN), intent(in) :: index_tsrf
    integer(IN), intent(in) :: index_topo
    integer(IN), intent(in) :: index_qice

    ! Local Varaibles
    integer(IN) :: j,jj,i,g,nxg,nyg,n
    character(*), parameter :: subName = "(glc_import_esmf) "
    !-----------------------------------------------------

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    do j = 1, nyg             ! S to N
       jj = nyg - j + 1       ! reverse j index for glint grid (N to S)
       do i = 1, nxg
          g = (j-1)*nxg + i   ! global index (W to E, S to N)
          tsfc(i,jj,ndx) = fptr(index_tsrf,g) - tkfrz
          topo(i,jj,ndx) = fptr(index_topo,g)
          qsmb(i,jj,ndx) = fptr(index_qice,g)
       enddo
    enddo

    if (verbose .and. my_task==master_task) then
       write(stdout,*) ' '
       write(stdout,*) subname,' x2g tsrf ',ndx,minval(fptr(index_tsrf,:)),maxval(fptr(index_tsrf,:))
       write(stdout,*) subname,' x2g topo ',ndx,minval(fptr(index_topo,:)),maxval(fptr(index_topo,:))
       write(stdout,*) subname,' x2g qice ',ndx,minval(fptr(index_qice,:)),maxval(fptr(index_qice,:))
       call shr_sys_flush(stdout)
    endif

end subroutine glc_import_esmf

!====================================================================================

subroutine glc_export_esmf(fptr, ndx, index_frac,index_topo, index_rofi,index_rofl,index_hflx)

    !-------------------------------------------------------------------
    implicit none
    real(R8)   , pointer    :: fptr(:,:)
    integer(IN), intent(in) :: ndx
    integer(IN), intent(in) :: index_frac
    integer(IN), intent(in) :: index_topo
    integer(IN), intent(in) :: index_rofi
    integer(IN), intent(in) :: index_rofl
    integer(IN), intent(in) :: index_hflx

    integer(IN) :: j,jj,i,g,nxg,nyg,n
    character(*), parameter :: subName = "(glc_export_esmf) "
    !-------------------------------------------------------------------

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    do j = 1, nyg           ! S to N
       jj = nyg - j + 1     ! reverse j index for glint grid (N to S)
       do i = 1, nxg
          g = (j-1)*nxg + i ! global index (W to E, S to N)
          fptr(index_frac,g) = gfrac(i,jj,ndx)
          fptr(index_topo,g) = gtopo(i,jj,ndx)
          fptr(index_rofi,g) = grofi(i,jj,ndx)
          fptr(index_rofl,g) = grofl(i,jj,ndx)
          fptr(index_hflx,g) = ghflx(i,jj,ndx)
       enddo
    enddo

    if (verbose .and. my_task==master_task) then
       write(stdout,*) subname,' g2x frac ',ndx,minval(fptr(index_frac,:)),maxval(fptr(index_frac,:))
       write(stdout,*) subname,' g2x topo ',ndx,minval(fptr(index_topo,:)),maxval(fptr(index_topo,:))
       write(stdout,*) subname,' g2x rofi ',ndx,minval(fptr(index_rofi,:)),maxval(fptr(index_rofi,:))
       write(stdout,*) subname,' g2x rofl ',ndx,minval(fptr(index_rofl,:)),maxval(fptr(index_rofl,:))
       write(stdout,*) subname,' g2x hflx ',ndx,minval(fptr(index_hflx,:)),maxval(fptr(index_hflx,:))
       call shr_sys_flush(stdout)
    endif

end subroutine glc_export_esmf

!=======================================================================

subroutine glc_domain_esmf( dom, rc )

    !-------------------------------------------------------------------
    implicit none
    type(ESMF_Array), intent(inout)     :: dom
    integer, intent(out)                :: rc

    ! Local Variables
    integer :: j,i,n,nxg,nyg          
    integer :: klon,klat,karea,kmask,kfrac ! domain fields
    real(R8), pointer :: fptr(:,:)
    !-------------------------------------------------------------------

    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
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
    ! Note aream will be filled in in the atm-lnd mapper

    nxg = glc_grid%nx
    nyg = glc_grid%ny
    fptr(:,:) = -9999.0_R8
    fptr(kmask,:) = -0.0_R8
    do j = 1,nyg
       do i = 1,nxg
          n = (j-1)*nxg + i
          fptr(klon , n) = glc_grid%lons(i)
          fptr(klat , n) = glc_grid%lats(j)
          fptr(karea, n) = glc_grid%box_areas(i,j)/(radius*radius)
          fptr(kmask, n) = real(glc_landmask(i,j), r8)  ! data is r8, glc_landmask is i4
          fptr(kfrac, n) = glc_landfrac(i,j)
       end do
    end do

end subroutine glc_domain_esmf

end module glc_comp_esmf

