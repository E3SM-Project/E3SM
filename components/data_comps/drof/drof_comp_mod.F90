#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module drof_comp_mod

! !USES:

  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                               CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use shr_mpi_mod      , only: shr_mpi_bcast
  use mct_mod
  use esmf
  use perf_mod

  use shr_strdata_mod
  use shr_dmodel_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: seq_flds_x2r_fields, seq_flds_r2x_fields
!
! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: drof_comp_init
  public :: drof_comp_run
  public :: drof_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  !--- other ---
  character(CS) :: myModelName = 'rof'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "rof_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance
                                         ! (ie. "_0001" or "")
  character(CL) :: rof_mode
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: read_restart          ! start from restart

  character(len=*),parameter :: rpfile = 'rpointer.rof'
  character(len=*),parameter :: nullstr = 'undefined'

  type(shr_strdata_type),save :: SDROF

  type(mct_rearr) :: rearr

  integer(IN),parameter :: ktrans = 2
  character(12),parameter  :: avofld(1:ktrans) = &
     (/ "Forr_rofl   ","Forr_rofi   "/)
  character(12),parameter  :: avifld(1:ktrans) = &
     (/ "rofl        ","rofi        "/)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: drof_comp_init
!
! !DESCRIPTION:
!     initialize data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine drof_comp_init( EClock, cdata, x2r, r2x, NLFilename )

  use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
  use pio, only : iosystem_desc_t
  implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2r
    type(mct_aVect)             , intent(inout) :: r2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)   :: n,k         ! generic counters
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: COMPID      ! comp id
    integer(IN)   :: gsize       ! global size
    integer(IN)   :: lsize_r     ! local size
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: nunit          ! unit number
    logical       :: rof_present    ! flag
    logical       :: rofice_present ! flag
    logical       :: rof_prognostic ! flag
    logical       :: flood_present  ! flag
    character(CL) :: calendar       ! model calendar

    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsmap
    type(mct_gGrid)        , pointer :: dom

    character(CL) :: filePath    ! generic file path
    character(CL) :: fileName    ! generic file name
    character(CS) :: timeName    ! domain file: time variable name
    character(CS) ::  lonName    ! domain file: lon  variable name
    character(CS) ::  latName    ! domain file: lat  variable name
    character(CS) :: maskName    ! domain file: mask variable name
    character(CS) :: areaName    ! domain file: area variable name

    integer(IN)   :: yearFirst   ! first year to use in data stream
    integer(IN)   :: yearLast    ! last  year to use in data stream
    integer(IN)   :: yearAlign   ! data year that aligns with yearFirst

    character(CL) :: rof_in      ! dshr rof namelist
    character(CL) :: decomp      ! decomp strategy
    character(CL) :: rest_file   ! restart filename
    character(CL) :: rest_file_strm_r   ! restart filename for stream
    character(CL) :: restfilm    ! model restart file namelist
    character(CL) :: restfilsr   ! stream restart file namelist
    logical       :: force_prognostic_true ! if true set prognostic true
    logical       :: exists      ! file existance logical
    logical       :: exists_r    ! file existance logical
    integer(IN)   :: nu          ! unit number

    type(iosystem_desc_t), pointer :: rof_pio_subsys
    integer(IN) :: rof_pio_iotype

    !----- define namelist -----
    namelist / drof_nml / &
        rof_in, decomp, restfilm, restfilsr, &
        force_prognostic_true

    !--- formats ---
    character(*), parameter :: F00   = "('(drof_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(drof_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(drof_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(drof_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(drof_comp_init) ',a,i8,a)"
    character(*), parameter :: F05   = "('(drof_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(drof_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(drof_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(drof_comp_init) "
!-------------------------------------------------------------------------------


    call t_startf('DROF_INIT')

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsmap, dom=dom, infodata=infodata)

    ! Determine communicator groups and sizes

    call mpi_comm_rank(mpicom, my_task, ierr)
    call mpi_comm_size(mpicom, npes, ierr)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Set a Few Defaults
    !----------------------------------------------------------------------------

    rof_present    = .false.
    rofice_present = .false.
    rof_prognostic = .false.
    flood_present  = .false.
    call seq_infodata_GetData(infodata,read_restart=read_restart)

    !----------------------------------------------------------------------------
    ! Read drof_in
    !----------------------------------------------------------------------------

    call t_startf('drof_readnml')

    filename = "drof_in"//trim(inst_suffix)
    rof_in = "unset"
    decomp = "1d"
    restfilm  = trim(nullstr)
    restfilsr = trim(nullstr)
    force_prognostic_true = .false.
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=drof_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' rof_in = ',trim(rof_in)
       write(logunit,F00)' decomp = ',trim(decomp)
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F00)' restfilsr = ',trim(restfilsr)
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
    endif
    call shr_mpi_bcast(rof_in,mpicom,'rof_in')
    call shr_mpi_bcast(decomp,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfilsr,mpicom,'restfilsr')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')

    rest_file = trim(restfilm)
    rest_file_strm_r = trim(restfilsr)
    if (force_prognostic_true) then
       rof_present    = .true.
       rof_prognostic = .true.
    endif

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDROF,trim(rof_in),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Validate mode
    !----------------------------------------------------------------------------

    rof_mode = trim(SDROF%dataMode)

    ! check that we know how to handle the mode

    if (trim(rof_mode) == 'NULL'            .or. &
        trim(rof_mode) == 'CPLHIST'         .or. &
        trim(rof_mode) == 'DIATREN_ANN_RX1' .or. &
        trim(rof_mode) == 'DIATREN_IAF_RX1') then
        if (my_task == master_task) write(logunit,F00) 'rof mode = ',trim(rof_mode)
    else
      write(logunit,F00) ' ERROR illegal rof mode = ',trim(rof_mode)
      call shr_sys_abort()
   end if

    call t_stopf('drof_readnml')

    !----------------------------------------------------------------------------
    ! Initialize datasets
    !----------------------------------------------------------------------------

    call t_startf('drof_strdata_init')

    rof_pio_subsys => shr_pio_getiosys(trim(inst_name))
    rof_pio_iotype =  shr_pio_getiotype(trim(inst_name))

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    if (trim(rof_mode) /= 'NULL') then
       rof_present = .true.
       rofice_present = .true.
       call shr_strdata_pioinit(SDROF,rof_pio_subsys,rof_pio_iotype)
       call shr_strdata_init(SDROF,mpicom,compid,name='rof',&
            calendar=calendar)
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDROF,'SDROF data')
    endif

    call t_stopf('drof_strdata_init')

    !----------------------------------------------------------------------------
    ! Set flag to specify data components
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         rof_present=rof_present, rof_prognostic=rof_prognostic, &
         rofice_present=rofice_present, rof_nx=SDROF%nxg, rof_ny=SDROF%nyg)
    call seq_infodata_PutData(infodata, flood_present=flood_present)

    if (.not. rof_present) then
       RETURN
    end if

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('drof_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    call shr_dmodel_gsmapcreate(gsmap,SDROF%nxg*SDROF%nyg,compid,mpicom,decomp)
    lsize_r = mct_gsmap_lsize(gsmap,mpicom)

    if (rof_present) then
       call mct_rearr_init(SDROF%gsmap,gsmap,mpicom,rearr)
    end if

    call t_stopf('drof_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    if (rof_present) call shr_dmodel_rearrGGrid(SDROF%grid, dom, gsmap, rearr, mpicom)

    call t_stopf('drof_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('drof_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(x2r, rList=seq_flds_x2r_fields, lsize=lsize_r)
    call mct_aVect_zero(x2r)

    call mct_aVect_init(r2x, rList=seq_flds_r2x_fields, lsize=lsize_r)
    call mct_aVect_zero(r2x)
    call t_stopf('drof_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file)        == trim(nullstr) .and. &
           trim(rest_file_strm_r) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer'
             call shr_sys_flush(logunit)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (.not.exists) then
                write(logunit,F00) ' ERROR: rpointer file does not exist'
                call shr_sys_abort(trim(subname)//' ERROR: rpointer file missing')
             endif
             nu = shr_file_getUnit()
             open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
             read(nu,'(a)') rest_file
             read(nu,'(a)') rest_file_strm_r
             close(nu)
             call shr_file_freeUnit(nu)
             inquire(file=trim(rest_file_strm_r),exist=exists_r)
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm_r,mpicom,'rest_file_strm_r')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             call shr_sys_flush(logunit)
             inquire(file=trim(rest_file_strm_r),exist=exists_r)
          endif
       end if
       call shr_mpi_bcast(exists_r,mpicom,'exists_r')
       if (exists_r) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm_r)
          call shr_strdata_restRead(trim(rest_file_strm_r),SDROF,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm_r)
       endif
       call shr_sys_flush(logunit)
    end if

    !----------------------------------------------------------------------------
    ! Set initial rof state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call drof_comp_run( EClock, cdata, x2r, r2x)
    call t_adj_detailf(-2)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'drof_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DROF_INIT')

end subroutine drof_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: drof_comp_run
!
! !DESCRIPTION:
!     run method for data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine drof_comp_run( EClock, cdata, x2r, r2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(in)    :: cdata
   type(mct_aVect)             ,intent(inout) :: x2r
   type(mct_aVect)             ,intent(inout) :: r2x

!EOP

   !--- local ---
   type(mct_gsMap)        , pointer :: gsmap
   type(mct_gGrid)        , pointer :: dom

   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: yy,mm,dd          ! year month day
   integer(IN)   :: n                 ! indices
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: nl                ! land frac index
   integer(IN)   :: kl                ! index of landfrac
   integer(IN)   :: lsize_r           ! size of attr vect
   integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
   logical       :: newdata           ! has newdata been read
   logical       :: mssrmlf           ! remove old data
   logical       :: write_restart     ! restart now
   character(CL) :: case_name         ! case name
   character(CL) :: rest_file         ! restart_file
   character(CL) :: rest_file_strm_r  ! restart_file for stream
   integer(IN)   :: nu                ! unit number
   integer(IN)   :: nflds_r2x
   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(drof_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(drof_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(drof_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DROF_RUN')

   call t_startf('drof_run1')

  !----------------------------------------------------------------------------
  ! Reset shr logging to my log file
  !----------------------------------------------------------------------------
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (logUnit)

   call seq_cdata_setptrs(cdata, gsMap=gsmap, dom=dom, infodata=infodata)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
   call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   lsize_r = mct_avect_lsize(r2x)
   nflds_r2x = mct_avect_nRattr(r2x)

   call t_stopf('drof_run1')

   !--------------------
   ! UNPACK
   !--------------------

   ! do nothing currently

   !--------------------
   ! ADVANCE ROF
   !--------------------

   call t_barrierf('drof_r_BARRIER',mpicom)
   call t_startf('drof_r')

   if (trim(rof_mode) /= 'NULL') then
      call t_startf('drof_r_strdata_advance')

      call shr_strdata_advance(SDROF,currentYMD,currentTOD,mpicom,'drof_r')
      call t_stopf('drof_r_strdata_advance')
      call t_barrierf('drof_r_scatter_BARRIER',mpicom)
      call t_startf('drof_r_scatter')
      do n = 1,SDROF%nstreams
         call shr_dmodel_translateAV(SDROF%avs(n),r2x,avifld,avofld,rearr)
      enddo
      call t_stopf('drof_r_scatter')
      ! zero out "special values"
      do nf=1,nflds_r2x
      do n=1,lsize_r
         if (abs(r2x%rAttr(nf,n)) > 1.0e28) r2x%rAttr(nf,n) = 0.0_r8
!         write(6,*)'crrentymd, currenttod, nf,n,r2x= ',currentymd, currenttod, nf,n,r2x%rattr(nf,n)
      enddo
      enddo
   else
      call mct_aVect_zero(r2x)
   end if

   call t_stopf('drof_r')

   if (write_restart) then
      call t_startf('drof_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
           trim(case_name), '.drof'//trim(inst_suffix)//'.r.', &
           yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm_r,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
           trim(case_name), '.drof'//trim(inst_suffix)//'.rs1.', &
           yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm_r
         close(nu)
         call shr_file_freeUnit(nu)
      endif
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm_r),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm_r),SDROF,mpicom,trim(case_name),'SDROF strdata')
      call shr_sys_flush(logunit)
      call t_stopf('drof_restart')
   end if

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('drof_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('drof_run2')

   call t_stopf('DROF_RUN')

end subroutine drof_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: drof_comp_final
!
! !DESCRIPTION:
!     finalize method for data rof model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine drof_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(drof_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(drof_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(drof_comp_final) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DROF_FINAL')

   if (my_task == master_task) then
      write(logunit,F91)
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91)
   end if

   call t_stopf('DROF_FINAL')

end subroutine drof_comp_final
!===============================================================================
!===============================================================================


end module drof_comp_mod
