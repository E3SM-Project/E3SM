#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module docn_comp_mod

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
  use shr_pcdf_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: seq_flds_o2x_fields, &
                               seq_flds_x2o_fields
!
! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: docn_comp_init
  public :: docn_comp_run
  public :: docn_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  !--- other ---
  type(iosystem_desc_t), pointer :: iosystem
  character(CS) :: myModelName = 'ocn'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance 
                                         ! (ie. "_0001" or "")
  character(CL) :: ocn_mode              ! mode
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: firstcall             ! first call logical
  logical       :: scmMode = .false.     ! single column mode
  real(R8)      :: scmLat  = shr_const_SPVAL  ! single column lat
  real(R8)      :: scmLon  = shr_const_SPVAL  ! single column lon
  logical       :: read_restart          ! start from restart

  character(len=*),parameter :: rpfile = 'rpointer.ocn'
  character(len=*),parameter :: nullstr = 'undefined'

  real(R8),parameter :: cpsw    = shr_const_cpsw    ! specific heat of sea h2o ~ J/kg/K
  real(R8),parameter :: rhosw   = shr_const_rhosw   ! density of sea water ~ kg/m^3
  real(R8),parameter :: TkFrz   = shr_const_TkFrz   ! freezing point, fresh water (Kelvin)
  real(R8),parameter :: TkFrzSw = shr_const_TkFrzSw ! freezing point, sea   water (Kelvin)
  real(R8),parameter :: latice  = shr_const_latice  ! latent heat of fusion
  real(R8),parameter :: ocnsalt = shr_const_ocn_ref_sal  ! ocean reference salinity

  integer(IN)   :: kt,ks,ku,kv,kdhdx,kdhdy,kq  ! field indices
  integer(IN)   :: kswnet,klwup,klwdn,ksen,klat,kmelth,ksnow,kioff
  integer(IN)   :: kh,kqbot

  type(shr_strdata_type) :: SDOCN
  type(mct_rearr) :: rearr
  type(mct_avect) :: avstrm   ! av of data from stream
  real(R8), pointer :: somtp(:)
  integer , pointer :: imask(:)
  character(len=*),parameter :: flds_strm = 'strm_h:strm_qbot'

  integer(IN),parameter :: ktrans = 28
  character(12),parameter  :: avifld(1:ktrans) = &
     (/ "ifrac       ","pslv        ","duu10n      ","taux        ","tauy        ", &
        "swnet       ","lat         ","sen         ","lwup        ","lwdn        ", &
        "melth       ","salt        ","prec        ","snow        ","rain        ", &
        "evap        ","meltw       ","roff        ","ioff        ",                &
        "t           ","u           ","v           ","dhdx        ","dhdy        ", &
        "s           ","q           ","h           ","qbot        "                 /)
  character(12),parameter  :: avofld(1:ktrans) = &
     (/ "Si_ifrac    ","Sa_pslv     ","So_duu10n   ","Foxx_taux   ","Foxx_tauy   ", &
        "Foxx_swnet  ","Foxx_lat    ","Foxx_sen    ","Foxx_lwup   ","Faxa_lwdn   ", &
        "Fioi_melth  ","Fioi_salt   ","Faxa_prec   ","Faxa_snow   ","Faxa_rain   ", &
        "Foxx_evap   ","Fioi_meltw  ","Forr_roff   ","Forr_ioff   ",                &
        "So_t        ","So_u        ","So_v        ","So_dhdx     ","So_dhdy     ", &
        "So_s        ","Fioo_q      ","strm_h      ","strm_qbot   "                 /)

  save

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: docn_comp_init
!
! !DESCRIPTION:
!     initialize data ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine docn_comp_init( EClock, cdata, x2o, o2x, NLFilename )
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2o, o2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!EOP

    !--- local variables ---
    integer(IN)   :: n,k         ! generic counters
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: COMPID      ! comp id
    integer(IN)   :: gsize       ! global size
    integer(IN)   :: lsize     ! local size
    integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: kmask       ! field reference
    logical       :: ocn_present    ! flag
    logical       :: ocn_prognostic ! flag
    logical       :: ocnrof_prognostic  ! flag
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

    character(CL) :: ocn_in      ! dshr ocn namelist
    character(CL) :: decomp      ! decomp strategy
    character(CL) :: rest_file   ! restart filename
    character(CL) :: rest_file_strm   ! restart filename for stream
    character(CL) :: restfilm    ! restart filename for namelist
    character(CL) :: restfils    ! restart filename for stream for namelist
    logical       :: exists      ! file existance
    integer(IN)   :: nu          ! unit number

    !----- define namelist -----
    namelist / docn_nml / &
        ocn_in, decomp, restfilm, restfils

    !--- formats ---
    character(*), parameter :: F00   = "('(docn_comp_init) ',8a)"
    character(*), parameter :: F01   = "('(docn_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(docn_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(docn_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(docn_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(docn_comp_init) ',a,2f10.4)"
    character(*), parameter :: F90   = "('(docn_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(docn_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(docn_comp_init) "
!-------------------------------------------------------------------------------


    call t_startf('DOCN_INIT')

    firstcall = .true.

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsmap, dom=ggrid, infodata=infodata)

    ! Determine communicator groups and sizes

    call mpi_comm_rank(mpicom, my_task, ierr)
    call mpi_comm_size(mpicom, npes, ierr)

    inst_name   = seq_comm_name(COMPID)
    inst_index  = seq_comm_inst(COMPID)
    inst_suffix = seq_comm_suffix(COMPID)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('ocn_modelio.nml'//trim(inst_suffix),logUnit)
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

    ocn_present = .false.
    ocn_prognostic = .false.
    ocnrof_prognostic = .false.
    call seq_infodata_GetData(infodata,read_restart=read_restart)

    !----------------------------------------------------------------------------
    ! Read docn_in
    !----------------------------------------------------------------------------

    call t_startf('docn_readnml')

    filename = "docn_in"//trim(inst_suffix)
    ocn_in = "unset"
    decomp = "1d"
    restfilm = trim(nullstr)
    restfils = trim(nullstr)
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=docn_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' ocn_in = ',trim(ocn_in)
       write(logunit,F00)' decomp = ',trim(decomp)
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F00)' restfils = ',trim(restfils)
    endif
    call shr_mpi_bcast(ocn_in,mpicom,'ocn_in')
    call shr_mpi_bcast(decomp,mpicom,'decomp')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
 
    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDOCN,trim(ocn_in),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Validate mode
    !----------------------------------------------------------------------------

    ocn_mode = trim(SDOCN%dataMode)

    ! check that we know how to handle the mode

    if (trim(ocn_mode) == 'NULL' .or. &
        trim(ocn_mode) == 'SSTDATA' .or. &
        trim(ocn_mode) == 'COPYALL' .or. &
        trim(ocn_mode) == 'SOM') then
      if (my_task == master_task) &
         write(logunit,F00) ' ocn mode = ',trim(ocn_mode)
    else
      write(logunit,F00) ' ERROR illegal ocn mode = ',trim(ocn_mode)
      call shr_sys_abort()
    endif

    call t_stopf('docn_readnml')

    !----------------------------------------------------------------------------
    ! Initialize datasets
    !----------------------------------------------------------------------------

    call t_startf('docn_strdata_init')

    if (trim(ocn_mode) /= 'NULL') then
       ocn_present = .true.
       call seq_timemgr_EClockGetData( EClock, calendar=calendar )
       iosystem => shr_pio_getiosys(trim(inst_name))
       
       call shr_strdata_pioinit(SDOCN, iosystem, shr_pio_getiotype(trim(inst_name)))

       if (scmmode) then
          if (my_task == master_task) &
             write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
          call shr_strdata_init(SDOCN,mpicom,compid,name='ocn', &
                      scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, &
                      calendar=calendar)
       else
          call shr_strdata_init(SDOCN,mpicom,compid,name='ocn', &
                      calendar=calendar)
       endif
    endif

    if (trim(ocn_mode) == 'SOM') then
       ocn_prognostic = .true.
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDOCN,'SDOCN data')
    endif

    call t_stopf('docn_strdata_init')

    !----------------------------------------------------------------------------
    ! Set flag to specify data components
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, ocnrof_prognostic=ocnrof_prognostic, &
      ocn_present=ocn_present, ocn_prognostic=ocn_prognostic, &
      ocn_nx=SDOCN%nxg, ocn_ny=SDOCN%nyg )

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('docn_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    call shr_dmodel_gsmapcreate(gsmap,SDOCN%nxg*SDOCN%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    if (ocn_present) then
       call mct_rearr_init(SDOCN%gsmap,gsmap,mpicom,rearr)
    endif

    call t_stopf('docn_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('docn_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    if (ocn_present) call shr_dmodel_rearrGGrid(SDOCN%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('docn_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('docn_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=lsize)
    call mct_aVect_zero(o2x)

    kt    = mct_aVect_indexRA(o2x,'So_t')
    ks    = mct_aVect_indexRA(o2x,'So_s')
    ku    = mct_aVect_indexRA(o2x,'So_u')
    kv    = mct_aVect_indexRA(o2x,'So_v')
    kdhdx = mct_aVect_indexRA(o2x,'So_dhdx')
    kdhdy = mct_aVect_indexRA(o2x,'So_dhdy')
    kq    = mct_aVect_indexRA(o2x,'Fioo_q')

    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=lsize)
    call mct_aVect_zero(x2o)

    kswnet = mct_aVect_indexRA(x2o,'Foxx_swnet')
    klwup  = mct_aVect_indexRA(x2o,'Foxx_lwup')
    klwdn  = mct_aVect_indexRA(x2o,'Faxa_lwdn')
    ksen   = mct_aVect_indexRA(x2o,'Foxx_sen')
    klat   = mct_aVect_indexRA(x2o,'Foxx_lat')
    kmelth = mct_aVect_indexRA(x2o,'Fioi_melth')
    ksnow  = mct_aVect_indexRA(x2o,'Faxa_snow')
    kioff  = mct_aVect_indexRA(x2o,'Forr_ioff')

    call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
    call mct_aVect_zero(avstrm)

    kh    = mct_aVect_indexRA(avstrm,'strm_h')
    kqbot = mct_aVect_indexRA(avstrm,'strm_qbot')

    allocate(somtp(lsize))
    allocate(imask(lsize))

    kmask = mct_aVect_indexRA(ggrid%data,'mask')
    imask(:) = nint(ggrid%data%rAttr(kmask,:))

    call t_stopf('docn_initmctavs')

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
       if (trim(ocn_mode) == 'SOM') then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
          call shr_pcdf_readwrite('read',iosystem,SDOCN%io_type,trim(rest_file),mpicom,gsmap,rf1=somtp,rf1n='somtp')
       endif
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDOCN,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial ocn state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call docn_comp_run( EClock, cdata,  x2o, o2x)
    call t_adj_detailf(-2)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'docn_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DOCN_INIT')

end subroutine docn_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: docn_comp_run
!
! !DESCRIPTION:
!     run method for dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine docn_comp_run( EClock, cdata,  x2o, o2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2o        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: o2x        ! dead   -> driver

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
   logical       :: glcrun_alarm      ! is glc going to run now
   logical       :: newdata           ! has newdata been read
   logical       :: mssrmlf           ! remove old data
   integer(IN)   :: idt               ! integer timestep
   real(R8)      :: dt                ! timestep
   real(R8)      :: hn                ! h field
   logical       :: write_restart     ! restart now
   character(CL) :: case_name         ! case name
   character(CL) :: rest_file         ! restart_file
   character(CL) :: rest_file_strm    ! restart_file for stream
   integer(IN)   :: nu                ! unit number
   integer(IN)   :: nflds_x2o
   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(docn_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(docn_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(docn_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DOCN_RUN')

   call t_startf('docn_run1')

  !----------------------------------------------------------------------------
  ! Reset shr logging to my log file
  !----------------------------------------------------------------------------
   call shr_file_getLogUnit (shrlogunit)
   call shr_file_getLogLevel(shrloglev)
   call shr_file_setLogUnit (logUnit)

   call seq_cdata_setptrs(cdata, gsMap=gsmap, dom=ggrid)

   call seq_cdata_setptrs(cdata, infodata=infodata)

   call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
   call seq_timemgr_EClockGetData( EClock, curr_yr=yy, curr_mon=mm, curr_day=dd)
   call seq_timemgr_EClockGetData( EClock, dtime=idt)
   dt = idt * 1.0_r8
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   call t_stopf('docn_run1')

   !--------------------
   ! UNPACK
   !--------------------

   call t_startf('docn_unpack')

!  lsize = mct_avect_lsize(x2o)
!  nflds_x2o = mct_avect_nRattr(x2o)

!   do nf=1,nflds_x2o
!   do n=1,lsize
!     ?? = x2o%rAttr(nf,n)
!   enddo
!   enddo

   call t_stopf('docn_unpack')

   !--------------------
   ! ADVANCE OCN
   !--------------------

   call t_barrierf('docn_BARRIER',mpicom)
   call t_startf('docn')

   !--- copy all fields from streams to o2x as default ---

   if (trim(ocn_mode) /= 'NULL') then
      call t_startf('docn_strdata_advance')
      call shr_strdata_advance(SDOCN,currentYMD,currentTOD,mpicom,'docn')
      call t_stopf('docn_strdata_advance')
      call t_barrierf('docn_scatter_BARRIER',mpicom)
      call t_startf('docn_scatter')
      do n = 1,SDOCN%nstreams
         call shr_dmodel_translateAV(SDOCN%avs(n),o2x,avifld,avofld,rearr)
      enddo
      call t_stopf('docn_scatter')
   else
      call mct_aVect_zero(o2x)
   endif

   call t_startf('docn_mode')

   select case (trim(ocn_mode))

   case('COPYALL') 
      ! do nothing extra

   case('SSTDATA')
      lsize = mct_avect_lsize(o2x)
      do n = 1,lsize
         o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
         o2x%rAttr(ks   ,n) = ocnsalt
         o2x%rAttr(ku   ,n) = 0.0_r8
         o2x%rAttr(kv   ,n) = 0.0_r8
         o2x%rAttr(kdhdx,n) = 0.0_r8
         o2x%rAttr(kdhdy,n) = 0.0_r8
         o2x%rAttr(kq   ,n) = 0.0_r8
      enddo

   case('SOM')
      lsize = mct_avect_lsize(o2x)
      do n = 1,SDOCN%nstreams
         call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,avifld,avofld,rearr)
      enddo
      if (firstcall) then
         do n = 1,lsize
            if (.not. read_restart) then
               somtp(n) = o2x%rAttr(kt,n) + TkFrz
            endif
            o2x%rAttr(kt,n) = somtp(n)
            o2x%rAttr(kq,n) = 0.0_r8
         enddo
      else   ! firstcall
         do n = 1,lsize
         if (imask(n) /= 0) then
            !--- pull out h from av for resuse below ---
            hn = avstrm%rAttr(kh,n)
            !--- compute new temp ---
            o2x%rAttr(kt,n) = somtp(n) + &
               (x2o%rAttr(kswnet,n) + &  ! shortwave 
                x2o%rAttr(klwup ,n) + &  ! longwave
                x2o%rAttr(klwdn ,n) + &  ! longwave
                x2o%rAttr(ksen  ,n) + &  ! sensible
                x2o%rAttr(klat  ,n) + &  ! latent
                x2o%rAttr(kmelth,n) - &  ! ice melt
                avstrm%rAttr(kqbot ,n) - &  ! flux at bottom
                (x2o%rAttr(ksnow,n)+x2o%rAttr(kioff,n))*latice) * &  ! latent by prec and roff
                dt/(cpsw*rhosw*hn)
             !--- compute ice formed or melt potential ---
            o2x%rAttr(kq,n) = (TkFrzSw - o2x%rAttr(kt,n))*(cpsw*rhosw*hn)/dt  ! ice formed q>0
            o2x%rAttr(kt,n) = max(TkFrzSw,o2x%rAttr(kt,n))                    ! reset temp
            somtp(n) = o2x%rAttr(kt,n)                                        ! save temp
         endif
         enddo
      endif   ! firstcall

   end select

   call t_stopf('docn_mode')

   if (write_restart) then
      call t_startf('docn_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.docn'//trim(inst_suffix)//'.r.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.docn'//trim(inst_suffix)//'.rs1.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm
         close(nu)
         call shr_file_freeUnit(nu)
      endif
      if (trim(ocn_mode) == 'SOM') then
         if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),currentYMD,currentTOD
         call shr_pcdf_readwrite('write',iosystem,SDOCN%io_type,trim(rest_file),mpicom,gsmap,clobber=.true., &
            rf1=somtp,rf1n='somtp')
      endif
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm),SDOCN,mpicom,trim(case_name),'SDOCN strdata')
      call shr_sys_flush(logunit)
      call t_stopf('docn_restart')
   endif

   call t_stopf('docn')

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('docn_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if
   firstcall = .false.
      
   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('docn_run2')

   call t_stopf('DOCN_RUN')

end subroutine docn_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: docn_comp_final
!
! !DESCRIPTION:
!     finalize method for dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine docn_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(docn_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(docn_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(docn_comp_final) "
   integer :: rcode
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DOCN_FINAL')
   if (my_task == master_task) then
      write(logunit,F91) 
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91) 
   end if
      
   call t_stopf('DOCN_FINAL')

end subroutine docn_comp_final
!===============================================================================
!===============================================================================


end module docn_comp_mod
