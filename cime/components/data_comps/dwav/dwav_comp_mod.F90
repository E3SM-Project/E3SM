module dwav_comp_mod

! !USES:

  use shr_const_mod
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
  use pio, only : iosystem_desc_t, pio_init, pio_rearr_box

  use shr_strdata_mod
  use shr_dmodel_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: seq_flds_w2x_fields, &
                               seq_flds_x2w_fields

!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: dwav_comp_init
  public :: dwav_comp_run
  public :: dwav_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(iosystem_desc_t), pointer :: iosystem
  character(CS) :: myModelName = 'wav'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance
                                         ! (ie. "_0001" or "")

  character(CL) :: wav_mode              ! mode
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: read_restart          ! start from restart

  character(len=*),parameter :: rpfile = 'rpointer.wav'
  character(len=*),parameter :: nullstr = 'undefined'

  type(shr_strdata_type) :: SDWAV
  type(mct_rearr) :: rearr

  integer(IN),parameter :: ktrans = 3
  character(12),parameter  :: avifld(1:ktrans) = &
     (/"lamult      ","ustokes     ","vstokes     "/)
  character(12),parameter  :: avofld(1:ktrans) = &
     (/"Sw_lamult   ","Sw_ustokes  ","Sw_vstokes  "/)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dwav_comp_init
!
! !DESCRIPTION:
!     initialize data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine dwav_comp_init( EClock, cdata, x2w, w2x, NLFilename )
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2w, w2x
    character(len=*), optional  , intent(in)    :: NLFilename

!EOP

    !--- local variables ---
    integer(IN)   :: n,k         ! generic counters
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: COMPID      ! comp id
    integer(IN)   :: gsize       ! global size
    integer(IN)   :: lsize     ! local size
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: nunit       ! unit number
    logical       :: wav_present    ! flag
    logical       :: wav_prognostic ! flag
    character(CL) :: calendar    ! model calendar

    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsmap
    type(mct_gGrid)        , pointer :: ggrid

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

    character(CL) :: wav_in      ! dshr wav namelist
    character(CL) :: decomp      ! decomp strategy
    character(CL) :: rest_file   ! restart filename
    character(CL) :: rest_file_strm   ! restart filename for stream
    character(CL) :: restfilm    ! restart filename for namelist
    character(CL) :: restfils    ! restart filename for stream for namelist
    logical       :: exists      ! file existance
    integer(IN)   :: nu          ! unit number

    !----- define namelist -----
    namelist / dwav_nml / &
        wav_in, decomp, restfilm, restfils

    !--- formats ---
    character(*), parameter :: F00   = "('(dwav_comp_init) ',8a)"
    character(*), parameter :: F01   = "('(dwav_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dwav_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dwav_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(dwav_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(dwav_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(dwav_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dwav_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dwav_comp_init) "

!-------------------------------------------------------------------------------

    call t_startf('DWAV_INIT')

    !--------------------------------------------------------------------
    ! Initialize mpi
    !--------------------------------------------------------------------

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsmap, dom=ggrid, infodata=infodata)

    ! Determine communicator groups and sizes

    call mpi_comm_rank(mpicom, my_task, ierr)
    call mpi_comm_size(mpicom, npes, ierr)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    ! open log file
    if (my_task == master_task) then
       logunit = shr_file_getUnit()
       call shr_file_setIO('wav_modelio.nml'//trim(inst_suffix),logunit)
    else
       logunit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit)

    !----------------------------------------------------------------------------
    ! Set a Few Defaults
    !----------------------------------------------------------------------------

    wav_present = .false.
    wav_prognostic = .false.
    call seq_infodata_GetData(infodata,read_restart=read_restart)

    !----------------------------------------------------------------------------
    ! Read dwav_in
    !----------------------------------------------------------------------------

    call t_startf('dwav_readnml')
    !write(logunit,F00)' dwav_readnml...'

    filename = "dwav_in"//trim(inst_suffix)
    wav_in = "unset"
    decomp = "1d"
    restfilm = trim(nullstr)
    restfils = trim(nullstr)
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dwav_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' wav_in = ',trim(wav_in)
       write(logunit,F00)' decomp = ',trim(decomp)
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F00)' restfils = ',trim(restfils)
    endif
    call shr_mpi_bcast(wav_in,mpicom,'wav_in')
    call shr_mpi_bcast(decomp,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------
    !write(logunit,F00)' read dshr nml...'

    call shr_strdata_readnml(SDWAV,trim(wav_in),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Validate mode
    !----------------------------------------------------------------------------
    !write(logunit,F00)' validate mode...'

    wav_mode = trim(SDWAV%dataMode)

    ! check that we know how to handle the mode

    if (trim(wav_mode) == 'null' .or. &
        trim(wav_mode) == 'copyall') then
      if (my_task == master_task) &
         write(logunit,F00) ' wav mode = ',trim(wav_mode)
    else
      write(logunit,F00) ' ERROR illegal wav mode = ',trim(wav_mode)
      call shr_sys_abort()
    endif

    call t_stopf('dwav_readnml')

    !----------------------------------------------------------------------------
    ! Initialize datasets
    !----------------------------------------------------------------------------

    call t_startf('dwav_strdata_init')
    !write(logunit,F00)' dwav_strdata_init...'

    if (trim(wav_mode) /= 'null') then
       wav_present = .true.
       call seq_timemgr_EClockGetData( EClock, calendar=calendar )
       iosystem => shr_pio_getiosys(trim(inst_name))

       call shr_strdata_pioinit(SDWAV, iosystem, shr_pio_getiotype(trim(inst_name)))

       call shr_strdata_init(SDWAV,mpicom,compid,name='wav', &
                      calendar=calendar)
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDWAV,'SDWAV data')
    endif

    call t_stopf('dwav_strdata_init')

    !----------------------------------------------------------------------------
    ! Set flag to specify data components
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(cdata%infodata, wav_present=wav_present, &
         wav_prognostic=wav_prognostic, wav_nx=SDWAV%nxg, wav_ny=SDWAV%nyg)

    !if (.not. wav_present) then
    !   RETURN
    !end if

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dwav_initgsmaps')
    !write(logunit,F00)' dwav_initgsmaps...'
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    call shr_dmodel_gsmapcreate(gsmap,SDWAV%nxg*SDWAV%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    if (wav_present) then
       call mct_rearr_init(SDWAV%gsmap,gsmap,mpicom,rearr)
    endif

    write(logunit,*)'lsize= ',lsize
    call shr_sys_flush(logunit)

    call t_stopf('dwav_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dwav_initmctdom')
    !write(logunit,F00)' dwav_initmctdom...'

    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    if (wav_present) call shr_dmodel_rearrGGrid(SDWAV%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('dwav_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dwav_initmctavs')
    !write(logunit,F00)' dwav_initmctavs...'

    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_avect_init(w2x, rlist=seq_flds_w2x_fields, lsize=lsize)
    call mct_avect_zero(w2x)

    call mct_avect_init(x2w, rlist=seq_flds_x2w_fields, lsize=lsize)
    call mct_avect_zero(x2w)

    call t_stopf('dwav_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file) == trim(nullstr) .and. &
           trim(rest_file_strm) == trim(nullstr)) then
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
             read(nu,'(a)') rest_file_strm
             close(nu)
             call shr_file_freeUnit(nu)
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             call shr_sys_flush(logunit)
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
       endif
       call shr_mpi_bcast(exists,mpicom,'exists')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDWAV,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial wav state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dwav_comp_run( EClock, cdata,  x2w, w2x)
    call t_adj_detailf(-2)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    ! end redirection of share output to wav log

    if (my_task == master_task) write(logunit, F00) 'dwav_comp_init done'
    call shr_sys_flush(logunit)
    call shr_file_setlogunit (shrlogunit)
    call shr_file_setloglevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DWAV_INIT')

end subroutine dwav_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dwav_comp_run
!
! !DESCRIPTION:
!     run method for data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dwav_comp_run( EClock, cdata, x2w, w2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2w
   type(mct_aVect)             ,intent(inout) :: w2x

!EOP

   !--- local ---
   type(mct_gsMap)        , pointer :: gsmap
   type(mct_gGrid)        , pointer :: ggrid

   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: yy,mm,dd          ! year month day
   integer(IN)   :: n                 ! indices
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: nl                ! ocn frac index
   integer(IN)   :: lsize           ! size of attr vect
   integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
!   logical       :: glcrun_alarm      ! is glc going to run now
   logical       :: newdata           ! has newdata been read
   logical       :: mssrmlf           ! remove old data
   integer(IN)   :: idt               ! integer timestep
   real(R8)      :: dt                ! timestep
!   real(R8)      :: hn                ! h field
   logical       :: write_restart     ! restart now
   character(CL) :: case_name         ! case name
   character(CL) :: rest_file         ! restart_file
   character(CL) :: rest_file_strm    ! restart_file for stream
   integer(IN)   :: nu                ! unit number
   integer(IN)   :: nflds_x2w
   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(dwav_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(dwav_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(dwav_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DWAV_RUN')

   call t_startf('dwav_run1')

  !----------------------------------------------------------------------------
  ! Reset shr logging to my log file
  !----------------------------------------------------------------------------
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (logUnit)

   call seq_cdata_setptrs(cdata, gsMap=gsmap, dom=ggrid, infodata=infodata)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
   call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
   call seq_timemgr_EClockGetData( EClock, dtime=idt)
   dt = idt * 1.0_r8
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   call t_stopf('dwav_run1')

   !--------------------
   ! UNPACK
   !--------------------

   call t_startf('dwav_unpack')

!  lsize = mct_avect_lsize(x2o)
!  nflds_x2o = mct_avect_nRattr(x2o)

!   do nf=1,nflds_x2o
!   do n=1,lsize
!     ?? = x2o%rAttr(nf,n)
!   enddo
!   enddo

   call t_stopf('dwav_unpack')

   !--------------------
   ! ADVANCE WAV
   !--------------------

   call t_barrierf('dwav_BARRIER',mpicom)
   call t_startf('dwav')

   !--- copy all fields from streams to w2x as default ---

   if (trim(wav_mode) /= 'null') then
      call t_startf('dwav_strdata_advance')
      call shr_strdata_advance(SDWAV,currentYMD,currentTOD,mpicom,'dwav')
      call t_stopf('dwav_strdata_advance')
      call t_barrierf('dwav_scatter_BARRIER',mpicom)
      call t_startf('dwav_scatter')
      do n = 1,SDWAV%nstreams
         call shr_dmodel_translateAV(SDWAV%avs(n),w2x,avifld,avofld,rearr)
      enddo
      call t_stopf('dwav_scatter')
   else
      call mct_aVect_zero(w2x)
   endif

   call t_startf('dwav_mode')

   call t_stopf('dwav_mode')

   if (write_restart) then
      call t_startf('dwav_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dwav'//trim(inst_suffix)//'.r.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dwav'//trim(inst_suffix)//'.rs1.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm
         close(nu)
         call shr_file_freeUnit(nu)
      endif
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm),SDWAV,mpicom,trim(case_name),'SDWAV strdata')
      call shr_sys_flush(logunit)
      call t_stopf('dwav_restart')
   endif

   call t_stopf('dwav')

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('dwav_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if

   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('dwav_run2')

   call t_stopf('DWAV_RUN')

end subroutine dwav_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dwav_comp_final
!
! !DESCRIPTION:
!     finalize method data wav model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine dwav_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(dwav_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(dwav_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(dwav_comp_final) "
   integer :: rcode
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DWAV_FINAL')
   if (my_task == master_task) then
      write(logunit,F91)
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91)
   end if

   call t_stopf('DWAV_FINAL')

 end subroutine dwav_comp_final

!===============================================================================

end module dwav_comp_mod
