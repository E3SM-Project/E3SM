module atm_comp_mct

! !USES:

  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use shr_mpi_mod      , only: shr_mpi_bcast
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use mct_mod
  use ESMF

  use dead_data_mod
  use dead_mct_mod
  use dead_mod

  use seq_cdata_mod
  use seq_flds_mod     , only: flds_d2x => seq_flds_a2x_fields, &
                               flds_x2d => seq_flds_x2a_fields
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

!--------------------------------------------------------------------------
! Private data interfaces
!--------------------------------------------------------------------------

  !--- stdin input stuff ---
  character(CS) :: str                  ! cpp  defined model name
  
  !--- other ---
  integer(IN)   :: dbug = 0             ! debug level (higher is more)
  
  real(R8), pointer :: gbuf(:,:)        ! grid info buffer
  
  integer(IN)   :: nproc_x       ! num of i pes (type 3)
  integer(IN)   :: seg_len       ! length of segs (type 4)
  integer(IN)   :: nxg           ! global dim i-direction
  integer(IN)   :: nyg           ! global dim j-direction
  integer(IN)   :: decomp_type   ! data decomp type:

  character(CS) :: myModelName = 'atm'   ! user defined model name
  integer(IN)   :: ncomp =1              ! component index
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom 
  integer(IN)   :: master_task=0         ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance 
                                         ! (ie. "_0001" or "")

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_init_mct
!
! !DESCRIPTION:
!     initialize dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine atm_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)   :: unitn       ! Unit for namelist file
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: local_comm  ! local communicator
    integer(IN)   :: mype        ! pe info
    integer(IN)   :: totpe       ! total number of pes
    real(R8)      :: nextsw_cday ! calendar of next atm shortwave
    integer(IN), allocatable :: gindex(:)  ! global index
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level

    integer(IN)                           :: COMPID
    integer(IN)                           :: mpicom
    integer(IN)                           :: lsize
    integer(IN)                           :: phase
    type(mct_gsMap)             , pointer :: gsMap
    type(mct_gGrid)             , pointer :: dom
    type(seq_infodata_type)     , pointer :: infodata

    !--- formats ---
    character(*), parameter :: F00   = "('(atm_init_mct) ',8a)"
    character(*), parameter :: F01   = "('(atm_init_mct) ',a,4i8)"
    character(*), parameter :: F02   = "('(atm_init_mct) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(atm_init_mct) ',a,i8,a)"
    character(*), parameter :: F90   = "('(atm_init_mct) ',73('='))"
    character(*), parameter :: F91   = "('(atm_init_mct) ',73('-'))"
    character(*), parameter :: subName = "(atm_init_mct) "
!-------------------------------------------------------------------------------

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    call seq_infodata_getData(infodata,atm_phase=phase)
    if (phase > 1) return

    call mpi_comm_rank(mpicom, my_task, ierr)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    ! read the namelist input (used to configure model)

    ncomp            =  1  ! set in dead_data, used by dead_run_mct

    nxg            =  -9999
    nyg            =  -9999
    nproc_x        =  -9999
    seg_len        =  -9999
    decomp_type    =  -9999

    if (my_task == master_task) then
       unitn = shr_file_getUnit()
       open( unitn, file='xatm_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       close (unitn)
       call shr_file_freeunit(unitn)
    endif

    call shr_mpi_bcast(nxg        ,mpicom,'xatm nxg')
    call shr_mpi_bcast(nyg        ,mpicom,'xatm nyg')
    call shr_mpi_bcast(decomp_type,mpicom,'xatm decomp_type')
    call shr_mpi_bcast(nproc_x    ,mpicom,'xatm nproc_x')
    call shr_mpi_bcast(seg_len    ,mpicom,'xatm seg_len')

    if (my_task == master_task) then
       write(logunit,*  ) ' Read in Xatm input from file= xatm_in'
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
    
    ! Determine communicator groups and sizes

    local_comm = mpicom
    call MPI_COMM_RANK(local_comm,mype ,ierr)
    call MPI_COMM_SIZE(local_comm,totpe,ierr)

    ! Determine decomposition and grid for dead component
    ! Initialize grid

    call dead_setNewGrid(decomp_type,nxg,nyg,totpe,mype,lsize,gbuf,seg_len,nproc_x)

    ! Initialize MCT global seg map

    allocate(gindex(lsize))
    gindex(:) = gbuf(:,dead_grid_index)
    call mct_gsMap_init( gsMap, gindex, mpicom, COMPID, lsize, nxg*nyg)
    deallocate(gindex)

    ! Initialize MCT domain

    call dead_domain_mct(mpicom, gbuf, gsMap, dom)

    ! Initialize MCT attribute vectors

    call mct_aVect_init(d2x, rList=flds_d2x, lsize=lsize)
    call mct_aVect_zero(d2x)

    call mct_aVect_init(x2d, rList=flds_x2d, lsize=lsize)
    call mct_aVect_zero(x2d)

    ! Get relevant data to send back to driver

!!    call atm_run_mct( EClock, cdata, x2d, d2x )

    ! Set flag to specify dead components

    if (nxg == 0 .and. nyg == 0) then
       call seq_infodata_PutData( infodata, dead_comps=.true., atm_present=.false., &
         atm_prognostic=.false., atm_nx=nxg, atm_ny=nyg)
    else
       call seq_infodata_PutData( infodata, dead_comps=.true., atm_present=.true., &
         atm_prognostic=.true., atm_nx=nxg, atm_ny=nyg)
    endif

   ! Set time step of radiation computation

    call seq_timemgr_EClockGetData(EClock, next_cday=nextsw_cday)
    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

end subroutine atm_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_run_mct
!
! !DESCRIPTION:
!     run method for dead atm model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_run_mct(EClock, cdata, x2d, d2x)

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver
!   type(eshr_timeMgr_clockType),intent(inout) :: SyncClock  ! Synchronization clock

!EOP

   !--- local ---
   type(seq_infodata_type), pointer :: infodata
   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: StepNo            ! sync clock step 
   integer(IN)   :: DTime             ! sync clock time step
   integer(IN)   :: n                 ! index
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: lsize             ! size of attr vect
   real(R8)      :: lat               ! latitude
   real(R8)      :: lon               ! longitude
   real(R8)      :: nextsw_cday       ! calendar of next atm shortwave
   integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
   integer(IN)   :: nflds_d2x, nflds_x2d
   character(*), parameter :: F00   = "('(atm_run_mct) ',8a)"
   character(*), parameter :: F04   = "('(atm_run_mct) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(atm_run_mct) "
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Reset shr logging to my log file
   !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, infodata=infodata)

    lsize = mct_avect_lsize(x2d)
    nflds_d2x = mct_avect_nRattr(d2x)
    nflds_x2d = mct_avect_nRattr(x2d)

    ! UNPACK

    do nf=1,nflds_x2d
    do n=1,lsize
!      ?? = x2d%rAttr(nf,n)
    enddo
    enddo

    ! PACK

    do nf=1,nflds_d2x
    do n=1,lsize
       lon = gbuf(n,dead_grid_lon)
       lat = gbuf(n,dead_grid_lat)
       d2x%rAttr(nf,n) = (nf*100)                   &
                    *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                    *  sin((SHR_CONST_PI*lon/180.0_R8)       &
                    -      (ncomp-1)*(SHR_CONST_PI/3.0_R8) ) &
                    + (ncomp*10.0_R8)
    enddo
    enddo
    
    ! log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData (EClock, curr_ymd=currentYMD, curr_tod=currentTOD)
       write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if
       
    ! Set time of next radiadtion computation

    call seq_timemgr_EClockGetData (EClock, next_cday=nextsw_cday)
    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)
    
    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

end subroutine atm_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: atm_final_mct
!
! !DESCRIPTION:
!     finalize method for dead model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine atm_final_mct(EClock, cdata, x2d, d2x)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(atm_final_mct) ',8a)"
   character(*), parameter :: F91   = "('(atm_final_mct) ',73('-'))"
   character(*), parameter :: subName = "(atm_final_mct) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (my_task == master_task) then
      write(logunit,F91) 
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91) 
   end if

end subroutine atm_final_mct
!===============================================================================

end module atm_comp_mct
