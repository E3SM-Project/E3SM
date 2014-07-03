#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dlnd_comp_mod

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
  use seq_flds_mod     , only: seq_flds_l2x_fields, seq_flds_x2l_fields, &
                               glc_nec=>seq_flds_glc_nec
!
! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: dlnd_comp_init
  public :: dlnd_comp_run
  public :: dlnd_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  !--- other ---
  character(CS) :: myModelName = 'lnd'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance 
                                         ! (ie. "_0001" or "")
  character(CL) :: lnd_mode 
  character(CL) :: sno_mode 
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: scmMode = .false.     ! single column mode
  real(R8)      :: scmLat  = shr_const_SPVAL  ! single column lat
  real(R8)      :: scmLon  = shr_const_SPVAL  ! single column lon
  logical       :: read_restart          ! start from restart

  character(len=*),parameter :: rpfile = 'rpointer.lnd'
  character(len=*),parameter :: nullstr = 'undefined'
  
  type(shr_strdata_type),save :: SDLND

  type(mct_rearr) :: rearr_l

  !--- names of fields ---
  integer(IN),parameter :: fld_len = 12       ! max character length of fields in avofld & avifld
  integer(IN),parameter :: nflds_nosnow = 22
  ! fields other than snow fields:
  character(fld_len),parameter  :: avofld_nosnow(1:nflds_nosnow) = &
     (/ "Sl_t        ","Sl_tref     ","Sl_qref     ","Sl_avsdr    ","Sl_anidr    ", &
        "Sl_avsdf    ","Sl_anidf    ","Sl_snowh    ","Fall_taux   ","Fall_tauy   ", &
        "Fall_lat    ","Fall_sen    ","Fall_lwup   ","Fall_evap   ","Fall_swnet  ", &
        "Sl_landfrac ","Sl_fv       ","Sl_ram1     ",                               &
        "Fall_flxdst1","Fall_flxdst2","Fall_flxdst3","Fall_flxdst4"                 /)
  character(fld_len),parameter  :: avifld_nosnow(1:nflds_nosnow) = &
     (/ "t           ","tref        ","qref        ","avsdr       ","anidr       ", &
        "avsdf       ","anidf       ","snowh       ","taux        ","tauy        ", &
        "lat         ","sen         ","lwup        ","evap        ","swnet       ", &
        "lfrac       ","fv          ","ram1        ",                               &
        "flddst1     ","flxdst2     ","flxdst3     ","flxdst4     "                 /)

  integer(IN), parameter :: nflds_snow = 3   ! number of snow fields in each elevation class
  integer(IN), parameter :: nec_len    = 2   ! length of elevation class index in field names
  ! for these snow fields, the actual field names will have the elevation class index at
  ! the end (e.g., Sl_tsrf01, tsrf01)
  character(fld_len-nec_len),parameter :: avofld_snow(nflds_snow) = &
       (/"Sl_tsrf  ", "Sl_topo  ", "Flgl_qice"/)
  character(fld_len-nec_len),parameter :: avifld_snow(nflds_snow) = &
       (/"tsrf", "topo", "qice"/)

  ! all fields:
  character(fld_len),dimension(:),allocatable :: avofld
  character(fld_len),dimension(:),allocatable :: avifld

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dlnd_comp_init
!
! !DESCRIPTION:
!     initialize data lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dlnd_comp_init( EClock, cdata_l, x2l, l2x, NLFilename )

  use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
  use pio, only : iosystem_desc_t
  implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_l
    type(mct_aVect)             , intent(inout) :: x2l, l2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)   :: n,k         ! generic counters
    integer(IN)   :: field_num   ! field number
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: COMPID      ! comp id
    integer(IN)   :: gsize       ! global size
    integer(IN)   :: lsize_l     ! local size
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: nunit          ! unit number
    logical       :: lnd_present    ! flag
    logical       :: lnd_prognostic ! flag
    logical       :: force_prognostic_true ! if true set prognostic true
    character(CL) :: calendar       ! model calendar

    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap_l
    type(mct_gGrid)        , pointer :: dom_l

    character(CL) :: filePath    ! generic file path
    character(CL) :: fileName    ! generic file name
    character(CS) :: timeName    ! domain file: time variable name
    character(CS) ::  lonName    ! domain file: lon  variable name
    character(CS) ::  latName    ! domain file: lat  variable name
    character(CS) :: maskName    ! domain file: mask variable name
    character(CS) :: areaName    ! domain file: area variable name
    character(CS) :: nec_format  ! format for nec_str
    character(nec_len):: nec_str ! elevation class, as character string

    integer(IN)   :: yearFirst   ! first year to use in data stream
    integer(IN)   :: yearLast    ! last  year to use in data stream
    integer(IN)   :: yearAlign   ! data year that aligns with yearFirst

    character(CL) :: lnd_in      ! dshr lnd namelist
    character(CL) :: decomp      ! decomp strategy
    character(CL) :: rest_file   ! restart filename
    character(CL) :: rest_file_strm_l   ! restart filename for stream
    character(CL) :: restfilm    ! model restart file namelist
    character(CL) :: restfilsl   ! stream restart file namelist
    logical       :: exists      ! file existance logical
    logical       :: exists_l    ! file existance logical
    integer(IN)   :: nu          ! unit number

    type(iosystem_desc_t), pointer :: lnd_pio_subsys
    integer(IN) :: lnd_pio_iotype

    !----- define namelist -----
    namelist / dlnd_nml / &
        lnd_in, decomp, restfilm, restfilsl, &
        force_prognostic_true

    !--- formats ---
    character(*), parameter :: F00   = "('(dlnd_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dlnd_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dlnd_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dlnd_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dlnd_comp_init) ',a,i8,a)"
    character(*), parameter :: F05   = "('(dlnd_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(dlnd_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dlnd_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dlnd_comp_init) "
!-------------------------------------------------------------------------------


    call t_startf('DLND_INIT')

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata_l, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap_l, dom=dom_l, infodata=infodata)

    ! Determine communicator groups and sizes

    call mpi_comm_rank(mpicom, my_task, ierr)
    call mpi_comm_size(mpicom, npes, ierr)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),logUnit)
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

    call seq_infodata_getData(infodata,single_column=scmMode, &
   &                          scmlat=scmlat, scmlon=scmLon)

    lnd_present = .false.
    lnd_prognostic = .false.
    call seq_infodata_GetData(infodata,read_restart=read_restart)

    !----------------------------------------------------------------------------
    ! Read dlnd_in
    !----------------------------------------------------------------------------

    call t_startf('dlnd_readnml')

    filename = "dlnd_in"//trim(inst_suffix)
    lnd_in = "unset"
    decomp = "1d"
    restfilm = trim(nullstr)
    restfilsl = trim(nullstr)
    force_prognostic_true = .false.
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dlnd_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' lnd_in = ',trim(lnd_in)
       write(logunit,F00)' decomp = ',trim(decomp)
       write(logunit,F00)' restfilm  = ',trim(restfilm)
       write(logunit,F00)' restfilsl = ',trim(restfilsl)
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
    endif
    call shr_mpi_bcast(lnd_in,mpicom,'lnd_in')
    call shr_mpi_bcast(decomp,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfilsl,mpicom,'restfilsl')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')
 
    rest_file = trim(restfilm)
    rest_file_strm_l = trim(restfilsl)
    if (force_prognostic_true) then
       lnd_present    = .true.
       lnd_prognostic = .true.
    endif

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDLND,trim(lnd_in),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Validate mode
    !----------------------------------------------------------------------------

    lnd_mode = trim(SDLND%dataMode)

    ! check that we know how to handle the mode

    if (trim(lnd_mode) == 'NULL' .or. &
        trim(lnd_mode) == 'COPYALL') then
      if (my_task == master_task) &
         write(logunit,F00) ' lnd mode = ',trim(lnd_mode)
    else
      write(logunit,F00) ' ERROR illegal lnd mode = ',trim(lnd_mode)
      call shr_sys_abort()
    endif

    call t_stopf('dlnd_readnml')

    !----------------------------------------------------------------------------
    ! Build avofld & avifld
    !----------------------------------------------------------------------------

    ! Start with non-snow fields
    allocate(avofld(nflds_nosnow + glc_nec*nflds_snow))
    allocate(avifld(nflds_nosnow + glc_nec*nflds_snow))
    avofld(1:nflds_nosnow) = avofld_nosnow
    avifld(1:nflds_nosnow) = avifld_nosnow
    field_num = nflds_nosnow

    ! create a format string for nec_str; e.g., if nec_len=2, this will be '(i2.2)'
    ! (without the quotes) 
    write(nec_format,'(a2, i0, a1, i0, a1)') "(i", nec_len, ".", nec_len, ")"

    ! Append each snow field
    do k = 1, nflds_snow
       do n = 1, glc_nec
          ! nec_str will be something like '02' or '10'
          write(nec_str,nec_format) n

          field_num = field_num + 1
          avofld(field_num) = trim(avofld_snow(k))//nec_str
          avifld(field_num) = trim(avifld_snow(k))//nec_str
       end do
    end do

    !----------------------------------------------------------------------------
    ! Initialize datasets
    !----------------------------------------------------------------------------

    call t_startf('dlnd_strdata_init')

    lnd_pio_subsys => shr_pio_getiosys(trim(inst_name))
    lnd_pio_iotype = shr_pio_getiotype(trim(inst_name))

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    if (trim(lnd_mode) /= 'NULL') then
       lnd_present = .true.
       call shr_strdata_pioinit(SDLND,lnd_pio_subsys,lnd_pio_iotype)
       if (scmmode) then
          if (my_task == master_task) &
             write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', &
                   scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, &
                   calendar=calendar)
       else
          call shr_strdata_init(SDLND,mpicom,compid,name='lnd', &
                   calendar=calendar)
       endif
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDLND,'SDLND data')
    endif

    call t_stopf('dlnd_strdata_init')

    !----------------------------------------------------------------------------
    ! Set flag to specify data components
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         lnd_present=lnd_present, lnd_prognostic=lnd_prognostic, &
         lnd_nx=SDLND%nxg, lnd_ny=SDLND%nyg)

    ! ************ lnd_present ****************
    if (lnd_present) then
    ! ************ lnd_present ****************

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    call shr_dmodel_gsmapcreate(gsmap_l,SDLND%nxg*SDLND%nyg,compid,mpicom,decomp)
    lsize_l = mct_gsmap_lsize(gsmap_l,mpicom)

    if (lnd_present) then
       call mct_rearr_init(SDLND%gsmap,gsmap_l,mpicom,rearr_l)
    endif

    call t_stopf('dlnd_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    if (lnd_present) call shr_dmodel_rearrGGrid(SDLND%grid, dom_l, gsmap_l, rearr_l, mpicom)

    call t_stopf('dlnd_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dlnd_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(l2x, rList=seq_flds_l2x_fields, lsize=lsize_l)
    call mct_aVect_zero(l2x)

    call mct_aVect_init(x2l, rList=seq_flds_x2l_fields, lsize=lsize_l)
    call mct_aVect_zero(x2l)

    call t_stopf('dlnd_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file)        == trim(nullstr) .and. &
           trim(rest_file_strm_l) == trim(nullstr)) then
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
             read(nu,'(a)') rest_file_strm_l
             close(nu)
             call shr_file_freeUnit(nu)
             inquire(file=trim(rest_file_strm_l),exist=exists_l)
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm_l,mpicom,'rest_file_strm_l')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             call shr_sys_flush(logunit)
             inquire(file=trim(rest_file_strm_l),exist=exists_l)
          endif
       endif
       call shr_mpi_bcast(exists_l,mpicom,'exists_l')
       !if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       !call shr_pcdf_readwrite('read',trim(rest_file),mpicom,gsmap,rf1=somtp,rf1n='somtp')
       if (exists_l) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm_l)
          call shr_strdata_restRead(trim(rest_file_strm_l),SDLND,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm_l)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial lnd state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dlnd_comp_run( EClock, cdata_l,  x2l, l2x)
    call t_adj_detailf(-2)


    ! ************ lnd_present ****************
    end if   
    ! ************ lnd_present ****************

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'dlnd_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DLND_INIT')

end subroutine dlnd_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dlnd_comp_run
!
! !DESCRIPTION:
!     run method for dead lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dlnd_comp_run( EClock, cdata_l,  x2l, l2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata_l
   type(mct_aVect)             ,intent(inout) :: x2l        
   type(mct_aVect)             ,intent(inout) :: l2x        

!EOP

   !--- local ---
   type(mct_gsMap)        , pointer :: gsMap_l
   type(mct_gGrid)        , pointer :: dom_l

   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: yy,mm,dd          ! year month day
   integer(IN)   :: n                 ! indices
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: nl                ! land frac index
   integer(IN)   :: kl                ! index of landfrac
   integer(IN)   :: lsize_l           ! size of attr vect
   integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
   logical       :: newdata           ! has newdata been read
   logical       :: mssrmlf           ! remove old data
   logical       :: write_restart     ! restart now
   character(CL) :: case_name         ! case name
   character(CL) :: rest_file         ! restart_file
   character(CL) :: rest_file_strm_l  ! restart_file for stream
   integer(IN)   :: nu                ! unit number
   integer(IN)   :: nflds_x2l
   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(dlnd_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(dlnd_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(dlnd_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DLND_RUN')

   call t_startf('dlnd_run1')

  !----------------------------------------------------------------------------
  ! Reset shr logging to my log file
  !----------------------------------------------------------------------------
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (logUnit)

   call seq_cdata_setptrs(cdata_l, gsMap=gsMap_l, dom=dom_l)

   call seq_cdata_setptrs(cdata_l, infodata=infodata)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
   call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   lsize_l = mct_avect_lsize(x2l)
   nflds_x2l = mct_avect_nRattr(x2l)

   call t_stopf('dlnd_run1')

   !--------------------
   ! UNPACK
   !--------------------

   call t_startf('dlnd_unpack')

!   do nf=1,nflds_x2l
!   do n=1,lsize_l
!     ?? = x2l%rAttr(nf,n)
!   enddo
!   enddo

!   do nf=1,nflds_x2s
!   do n=1,lsize_s
!     ?? = x2s%rAttr(nf,n)
!   enddo
!   enddo

   call t_stopf('dlnd_unpack')

   !--------------------
   ! ADVANCE LAND
   !--------------------

   call t_barrierf('dlnd_l_BARRIER',mpicom)
   call t_startf('dlnd_l')

   if (trim(lnd_mode) /= 'NULL') then
      call t_startf('dlnd_l_strdata_advance')
      call shr_strdata_advance(SDLND,currentYMD,currentTOD,mpicom,'dlnd_l')
      call t_stopf('dlnd_l_strdata_advance')
      call t_barrierf('dlnd_l_scatter_BARRIER',mpicom)
      call t_startf('dlnd_l_scatter')
      do n = 1,SDLND%nstreams
         call shr_dmodel_translateAV(SDLND%avs(n),l2x,avifld,avofld,rearr_l)
      enddo
      call t_stopf('dlnd_l_scatter')
   else
      call mct_aVect_zero(l2x)
   endif

   call t_stopf('dlnd_l')

   if (write_restart) then
      call t_startf('dlnd_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dlnd'//trim(inst_suffix)//'.r.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm_l,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dlnd'//trim(inst_suffix)//'.rs1.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm_l
         close(nu)
         call shr_file_freeUnit(nu)
      endif
      !if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),currentYMD,currentTOD
      !call shr_pcdf_readwrite('write',trim(rest_file),mpicom,gsmap,clobber=.true., &
      !   rf1=somtp,rf1n='somtp')
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm_l),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm_l),SDLND,mpicom,trim(case_name),'SDLND strdata')
      call shr_sys_flush(logunit)
      call t_stopf('dlnd_restart')
   endif

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('dlnd_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if
      
   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('dlnd_run2')

   call t_stopf('DLND_RUN')

end subroutine dlnd_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dlnd_comp_final
!
! !DESCRIPTION:
!     finalize method for dead lnd model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine dlnd_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(dlnd_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(dlnd_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(dlnd_comp_final) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DLND_FINAL')

   if (my_task == master_task) then
      write(logunit,F91) 
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91) 
   end if
      
   call t_stopf('DLND_FINAL')

end subroutine dlnd_comp_final
!===============================================================================
!===============================================================================


end module dlnd_comp_mod
