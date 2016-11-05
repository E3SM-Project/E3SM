module wav_comp_esmf

#ifdef ESMF_INTERFACE
! !USES:
  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use shr_mpi_mod      , only: shr_mpi_bcast
  use shr_const_mod    , only: SHR_CONST_PI
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use ESMF

  use dead_data_mod
  use dead_mod

  use seq_flds_mod     , only: flds_d2x => seq_flds_w2x_fields, &
                               flds_x2d => seq_flds_x2w_fields, &
                               flds_dom => seq_flds_dom_fields

  use esmfshr_mod
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: wav_init_esmf
  public :: wav_run_esmf
  public :: wav_final_esmf
  public :: wav_register_esmf

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

  !--- stdin input stuff ---
  character(CS) :: str                  ! cpp  defined model name
  
  !--- other ---
  integer(IN)   :: dbug = 0             ! debug level (higher is more)
  
  character(CS) :: myModelName = 'wav'   ! user defined model name
  integer(IN)   :: ncomp = 7             ! component index
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom 
  integer(IN)   :: master_task=0         ! task number of master task
  integer(IN)   :: logunit               ! logging unit number

!
! Author: Fei Liu
! ESMF compliant data wav component
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine wav_register_esmf(comp, rc)

    implicit none

    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    print *, "In wav register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      wav_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      wav_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      wav_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_init_esmf
!
! !DESCRIPTION:
!     initialize dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_init_esmf(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

!EOP

    !--- local variables ---
    integer(IN)   :: unitn       ! Unit for namelist file
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: local_comm  ! local communicator
    integer(IN)   :: mype        ! pe info
    integer(IN)   :: totpe       ! total number of pes
    integer(IN), allocatable :: gindex(:)  ! global index
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
  
    real(R8), pointer :: gbuf(:,:)     ! grid info buffer
    real(R8), pointer :: buf(:)        ! tempoary buffer
    
    integer(IN)   :: nproc_x       ! num of i pes (type 3)
    integer(IN)   :: seg_len       ! length of segs (type 4)
    integer(IN)   :: nxg           ! global dim i-direction
    integer(IN)   :: nyg           ! global dim j-direction
    integer(IN)   :: decomp_type   ! data decomp type:

    integer(IN)       :: COMPID
    integer(IN)       :: inst_index            ! number of current instance (ie. 1)
    character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
    character(len=16) :: inst_suffix       ! char string associated with instance 
    integer(IN)                           :: mpicom, mpicom_vm
    integer(IN)                           :: lsize
    integer(IN)                           :: phase
    type(ESMF_Array)                      :: d2x_a, x2d_a, dom_a
    type(ESMF_DistGrid)                   :: distgrid
    type(ESMF_VM)                         :: vm

    character(ESMF_MAXSTR) :: convCIM, purpComp

    !--- formats ---
    character(*), parameter :: F00   = "('(wav_init_esmf) ',8a)"
    character(*), parameter :: F01   = "('(wav_init_esmf) ',a,4i8)"
    character(*), parameter :: F02   = "('(wav_init_esmf) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(wav_init_esmf) ',a,i8,a)"
    character(*), parameter :: F90   = "('(wav_init_esmf) ',73('='))"
    character(*), parameter :: F91   = "('(wav_init_esmf) ',73('-'))"
    character(*), parameter :: subName = "(wav_init_esmf) "

    !----------------------------
    ! Initial Setup
    !----------------------------

    rc = ESMF_SUCCESS

    call ESMF_AttributeGet(export_state, name="wav_phase", value=phase, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="ID", value=COMPID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! duplicate the mpi communicator from the current VM
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call MPI_Comm_dup(mpicom_vm, mpicom, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (phase > 1) return

    call mpi_comm_rank(mpicom, my_task, ierr)
    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('wav_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------
    ! read the namelist input (used to configure model)
    !----------------------------

    nxg            =  -9999
    nyg            =  -9999
    nproc_x        =  -9999
    seg_len        =  -9999
    decomp_type    =  -9999

    if (my_task == master_task) then
       unitn = shr_file_getUnit()
       open( unitn, file='xwav_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       close (unitn)
       call shr_file_freeunit(unitn)
    endif

    call shr_mpi_bcast(nxg        ,mpicom,'xwav nxg')
    call shr_mpi_bcast(nyg        ,mpicom,'xwav nyg')
    call shr_mpi_bcast(decomp_type,mpicom,'xwav decomp_type')
    call shr_mpi_bcast(nproc_x    ,mpicom,'xwav nproc_x')
    call shr_mpi_bcast(seg_len    ,mpicom,'xwav seg_len')

    if (my_task == master_task) then
       write(logunit,*  ) ' Read in Xwav input from file= xwav_in'//trim(inst_suffix)
       write(logunit,F00)
       write(logunit,F00) '         Model  :  ',trim(myModelName)
       write(logunit,F01) '           NGX  :  ',nxg
       write(logunit,F01) '           NGY  :  ',nyg
       write(logunit,F01) ' Decomposition  :  ',decomp_type
       write(logunit,F03) ' Num pes in X   :  ',nproc_x,'  (type 3 only)'
       write(logunit,F03) ' Segment Length :  ',seg_len,'  (type 11 only)'
       write(logunit,F01) '    inst_index  :  ',inst_index
       write(logunit,F00) '    inst_name   :  ',trim(inst_name)
       write(logunit,F00) '    inst_suffix :  ',trim(inst_suffix)
       write(logunit,F00)
       call shr_sys_flush(logunit)
    end if
    
    !----------------------------
    ! Determine communicator groups and sizes
    !----------------------------

    local_comm = mpicom
    call MPI_COMM_RANK(local_comm,mype ,ierr)
    call MPI_COMM_SIZE(local_comm,totpe,ierr)

    !----------------------------
    ! Determine decomposition and grid for dead component
    !----------------------------

    call dead_setNewGrid(decomp_type,nxg,nyg,totpe,mype,lsize,gbuf,seg_len,nproc_x)

    !----------------------------
    ! Set up distgrid
    !----------------------------

    allocate(gindex(lsize))
    gindex(:) = nint(gbuf(:,dead_grid_index))

    distgrid = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    deallocate(gindex)

    !----------------------------
    ! Init Arrays
    !----------------------------

    dom_a = mct2esmf_init(distgrid, attname=flds_dom, name="domain", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    d2x_a = mct2esmf_init(distgrid, attname=flds_d2x, name="d2x", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    x2d_a = mct2esmf_init(distgrid, attname=flds_x2d, name="x2d", rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------
    ! Fill domain
    !----------------------------

    call esmfshr_util_ArrayZero(dom_a, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate(buf(lsize))
    buf(:) = gbuf(:,dead_grid_lon)
    call esmfshr_util_ArrayPutField(dom_a, 'lon', buf)
    buf(:) = gbuf(:,dead_grid_lat)
    call esmfshr_util_ArrayPutField(dom_a, 'lat', buf)
    buf(:) = gbuf(:,dead_grid_area)
    call esmfshr_util_ArrayPutField(dom_a, 'area', buf)
    call esmfshr_util_ArrayPutField(dom_a, 'aream', buf)
    buf(:) = gbuf(:,dead_grid_mask)
    call esmfshr_util_ArrayPutField(dom_a, 'mask', buf)
    buf(:) = gbuf(:,dead_grid_frac)
    call esmfshr_util_ArrayPutField(dom_a, 'frac', buf)
    deallocate(buf)

    !----------------------------
    ! Add arrays to state
    !----------------------------

    call ESMF_StateAdd(export_state, (/dom_a/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateAdd(export_state, (/d2x_a/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateAdd(import_state, (/x2d_a/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------
    ! Set flags
    !----------------------------

    call ESMF_AttributeSet(export_state, name="dead_comps", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="wav_nx", value=nxg, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="wav_ny", value=nyg, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (nxg == 0 .and. nyg == 0) then
       call ESMF_AttributeSet(export_state, name="wav_present", value=.false., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeSet(export_state, name="wav_prognostic", value=.false., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    else
       call ESMF_AttributeSet(export_state, name="wav_present", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeSet(export_state, name="wav_prognostic", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "XWAV", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Wave Dead Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2013", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Wave", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "someone", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "someone@someplace", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
#endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

end subroutine wav_init_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_run_esmf
!
! !DESCRIPTION:
!     run method for dead wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_run_esmf(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

!EOP

    !--- local ---
    type(ESMF_Array)            :: d2x_a, dom_a
    real(R8), pointer           :: blon(:), blat(:)
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    integer(IN)                           :: lsize
    real(R8)      :: lat               ! latitude
    real(R8)      :: lon               ! longitude
    integer(IN)   :: n                 ! index
    integer(IN)   :: nf                ! fields loop index
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: CurrentYMD        ! model date
    integer(IN)   :: CurrentTOD        ! model sec into model date
    character(*), parameter :: F04   = "('(wav_run_esmf) ',2a,2i8,'s')"
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------
    ! Get arrays, blon and blat
    !----------------------------

    call ESMF_StateGet(export_state, itemName="domain", array=dom_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_ArrayGetSize(dom_a, lsize2=lsize)
    allocate(blon(lsize),blat(lsize))

    call esmfshr_util_ArrayGetField(dom_a, 'lon', blon, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_ArrayGetField(dom_a, 'lat', blat, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------
    ! Pack d2x_a
    ! the bounds are always from /1,1/ to /nflds_d2x, lsize/ locally.
    !----------------------------

    call ESMF_ArrayGet(d2x_a, localDe=0, farrayPtr=fptr, rc=rc)   
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    do n = 1, ubound(fptr,2)-lbound(fptr,2)+1
    do nf = 1, ubound(fptr,1)-lbound(fptr,1)+1
       lon = blon(n)
       lat = blat(n)
       fptr(nf-1+lbound(fptr,1),n-1+lbound(fptr,2)) = (nf*100)  &
          *  cos (SHR_CONST_PI*lat/180.0_R8)       &
          *  sin((SHR_CONST_PI*lon/180.0_R8)       &
          - (ncomp-1)*(SHR_CONST_PI/3.0_R8) )      &
          + (ncomp*10.0_R8)
    enddo
    enddo

    deallocate(blon,blat)

    !----------------------------
    ! Update attributes
    !----------------------------

    !----------------------------
    ! Log
    !----------------------------

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData (EClock, curr_ymd=currentYMD, curr_tod=currentTOD)
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

end subroutine wav_run_esmf

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: wav_final_esmf
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine wav_final_esmf(comp, import_state, export_state, EClock, rc)

    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
!EOP
    type(ESMF_Array)                 :: d2x_a, x2d_a, dom_a
    type(ESMF_DistGrid)              :: distgrid
    character(*), parameter :: F00   = "('(wav_final) ',8a)"
    character(*), parameter :: F91   = "('(wav_final) ',73('-'))"
 
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !tcraig, Mar 2013, return here to avoid esmf abort on the destroy, bug tcx
    return

    !----------------------------
    ! Destroy Arrays
    !----------------------------

    call ESMF_StateGet(export_state, itemName="domain", array=dom_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayDestroy(dom_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(export_state, itemName="d2x", array=d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_ArrayDestroy(d2x_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateGet(import_state, itemName="x2d", array=x2d_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_ArrayGet(x2d_a, distgrid=distgrid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_ArrayDestroy(x2d_a, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_DistGridDestroy(distgrid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (my_task == master_task) then
      write(logunit,F91)
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91)
    end if

end subroutine wav_final_esmf
!===============================================================================
#endif

end module wav_comp_esmf
