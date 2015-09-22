#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module dice_comp_mod

! !USES:

  use shr_const_mod
  use shr_sys_mod
  use shr_kind_mod , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                           CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                           shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                           shr_file_freeunit
  use shr_mpi_mod  , only: shr_mpi_bcast
  use shr_flux_mod , only: shr_flux_atmIce
  use shr_cal_mod  , only: shr_cal_ymd2julian
  use mct_mod
  use esmf
  use perf_mod

  use shr_strdata_mod
  use shr_dmodel_mod
  use shr_pcdf_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: seq_flds_i2x_fields, &
                               seq_flds_x2i_fields
!
! !PUBLIC TYPES:
  implicit none
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: dice_comp_init
  public :: dice_comp_run
  public :: dice_comp_final

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  !--- other ---
  character(CS) :: myModelName = 'ice'   ! user defined model name
  integer(IN)   :: mpicom
  integer(IN)   :: my_task               ! my task in mpi communicator mpicom
  integer(IN)   :: npes                  ! total number of tasks
  integer(IN),parameter :: master_task=0 ! task number of master task
  integer(IN)   :: logunit               ! logging unit number
  integer       :: inst_index            ! number of current instance (ie. 1)
  character(len=16) :: inst_name         ! fullname of current instance (ie. "lnd_0001")
  character(len=16) :: inst_suffix       ! char string associated with instance 
                                         ! (ie. "_0001" or "")
  character(CL) :: ice_mode              ! mode
  integer(IN)   :: dbug = 0              ! debug level (higher is more)
  logical       :: firstcall             ! first call logical
  logical       :: scmMode = .false.     ! single column mode
  real(R8)      :: scmLat  = shr_const_SPVAL  ! single column lat
  real(R8)      :: scmLon  = shr_const_SPVAL  ! single column lon
  logical       :: read_restart          ! start from restart
  real(R8)      :: flux_swpf             ! short-wave penatration factor
  real(R8)      :: flux_Qmin             ! bound on melt rate
  logical       :: flux_Qacc             ! activates water accumulation/melt wrt Q
  real(R8)      :: flux_Qacc0            ! initial water accumulation value

  character(len=*),parameter :: rpfile = 'rpointer.ice'
  character(len=*),parameter :: nullstr = 'undefined'

  real(R8),parameter  :: pi     = shr_const_pi      ! pi
  real(R8),parameter  :: spval  = shr_const_spval   ! flags invalid data
  real(R8),parameter  :: tFrz   = shr_const_tkfrzsw ! temp of freezing salt-water
  real(R8),parameter  :: latice = shr_const_latice  ! latent heat of fusion
  real(R8),parameter  :: cDay   = shr_const_cDay    ! sec in calendar day
  real(R8),parameter  :: waterMax = 1000.0_R8        ! wrt iFrac comp & frazil ice (kg/m^2)

  !----- surface albedo constants ------
  real(R8),parameter  :: snwfrac = 0.286_R8 ! snow cover fraction ~ [0,1]
  real(R8),parameter  :: as_nidf = 0.950_R8 ! albedo: snow,near-infr,diffuse
  real(R8),parameter  :: as_vsdf = 0.700_R8 ! albedo: snow,visible  ,diffuse
  real(R8),parameter  :: as_nidr = 0.960_R8 ! albedo: snow,near-infr,direct
  real(R8),parameter  :: as_vsdr = 0.800_R8 ! albedo: snow,visible  ,direct
  real(R8),parameter  :: ai_nidf = 0.700_R8 ! albedo: ice, near-infr,diffuse
  real(R8),parameter  :: ai_vsdf = 0.500_R8 ! albedo: ice, visible  ,diffuse
  real(R8),parameter  :: ai_nidr = 0.700_R8 ! albedo: ice, near-infr,direct
  real(R8),parameter  :: ai_vsdr = 0.500_R8 ! albedo: ice, visible  ,direct
  real(R8),parameter  :: ax_nidf = ai_nidf*(1.0_R8-snwfrac) + as_nidf*snwfrac
  real(R8),parameter  :: ax_vsdf = ai_vsdf*(1.0_R8-snwfrac) + as_vsdf*snwfrac
  real(R8),parameter  :: ax_nidr = ai_nidr*(1.0_R8-snwfrac) + as_nidr*snwfrac
  real(R8),parameter  :: ax_vsdr = ai_vsdr*(1.0_R8-snwfrac) + as_vsdr*snwfrac

  integer(IN) :: kswvdr,kswndr,kswvdf,kswndf,kq,kz,kua,kva,kptem,kshum,kdens,ktbot
  integer(IN) :: kiFrac,kt,kavsdr,kanidr,kavsdf,kanidf,kswnet,kmelth,kmeltw
  integer(IN) :: ksen,klat,klwup,kevap,ktauxa,ktauya,ktref,kqref,kswpen,ktauxo,ktauyo,ksalt

  type(shr_strdata_type) :: SDICE
  type(mct_rearr) :: rearr
!  type(mct_avect) :: avstrm   ! av of data from stream
  integer(IN) , pointer :: imask(:)
  real(R8)    , pointer :: yc(:)
  real(R8)    , pointer :: water(:)
!  real(R8)    , pointer :: ifrac0(:)

  integer(IN),parameter :: ktrans = 42
  character(16),parameter  :: avofld(1:ktrans) = &
     (/"So_t            ","So_s            ","So_u            ","So_v            ", &
       "So_dhdx         ","So_dhdy         ","Fioo_q          ","Sa_z            ", &
       "Sa_u            ","Sa_v            ","Sa_ptem         ","Sa_tbot         ", &
       "Sa_shum         ","Sa_dens         ","Faxa_swndr      ","Faxa_swvdr      ", &
       "Faxa_swndf      ","Faxa_swvdf      ","Faxa_lwdn       ","Faxa_rain       ", &
       "Faxa_snow       ","Si_t            ","Si_tref         ","Si_qref         ", &
       "Si_ifrac        ","Si_avsdr        ","Si_anidr        ","Si_avsdf        ", &
       "Si_anidf        ","Faii_taux       ","Faii_tauy       ","Faii_lat        ", &
       "Faii_sen        ","Faii_lwup       ","Faii_evap       ","Faii_swnet      ", &
       "Fioi_swpen      ","Fioi_melth      ","Fioi_meltw      ","Fioi_salt       ", &
       "Fioi_taux       ","Fioi_tauy       " /)

  character(16),parameter  :: avifld(1:ktrans) = &
     (/"to              ","s               ","uo              ","vo              ", &
       "dhdx            ","dhdy            ","q               ","z               ", &
       "ua              ","va              ","ptem            ","tbot            ", &
       "shum            ","dens            ","swndr           ","swvdr           ", &
       "swndf           ","swvdf           ","lwdn            ","rain            ", &
       "snow            ","t               ","tref            ","qref            ", &
       "ifrac           ","avsdr           ","anidr           ","avsdf           ", &
       "anidf           ","tauxa           ","tauya           ","lat             ", &
       "sen             ","lwup            ","evap            ","swnet           ", &
       "swpen           ","melth           ","meltw           ","salt            ", &
       "tauxo           ","tauyo           " /)

  save

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dice_comp_init
!
! !DESCRIPTION:
!     initialize data ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dice_comp_init( EClock, cdata, x2i, i2x, NLFilename )
    use pio, only : iosystem_desc_t
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2i, i2x
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
    integer(IN)   :: kfld        ! field reference
    logical       :: ice_present    ! flag
    logical       :: ice_prognostic ! flag

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
    character(CL) :: calendar    ! calendar type

    character(CL) :: ice_in      ! dshr ice namelist
    character(CL) :: decomp      ! decomp strategy
    character(CL) :: rest_file   ! restart filename
    character(CL) :: rest_file_strm   ! restart filename for stream
    character(CL) :: restfilm    ! model restart file namelist
    character(CL) :: restfils    ! stream restart file namelist
    logical       :: force_prognostic_true ! if true set prognostic true
    logical       :: exists      ! file existance logical
    integer(IN)   :: nu          ! unit number
    type(iosystem_desc_t), pointer :: ice_pio_subsystem


    !----- define namelist -----
    namelist / dice_nml / &
        ice_in, decomp, flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, restfilm, restfils, &
        force_prognostic_true

    !--- formats ---
    character(*), parameter :: F00   = "('(dice_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(dice_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(dice_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(dice_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(dice_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(dice_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(dice_comp_init) ',a,2f10.4)"
    character(*), parameter :: F06   = "('(dice_comp_init) ',a,5l3)"
    character(*), parameter :: F90   = "('(dice_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(dice_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(dice_comp_init) "
!-------------------------------------------------------------------------------


    call t_startf('DICE_INIT')

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
       call shr_file_setIO('ice_modelio.nml'//trim(inst_suffix),logUnit)
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

    ice_present = .false.
    ice_prognostic = .false.
    call seq_infodata_GetData(infodata,read_restart=read_restart)

    !----------------------------------------------------------------------------
    ! Read dice_in
    !----------------------------------------------------------------------------

    call t_startf('dice_readnml')

    filename = "dice_in"//trim(inst_suffix)
    ice_in = "unset"
    decomp = "1d"
    flux_swpf  =     0.0_R8  ! no penetration
    flux_Qmin  =  -300.0_R8  ! kg/s/m^2
    flux_Qacc  = .false.     ! no accumulation 
    flux_Qacc0 =     0.0_R8  ! no water
    restfilm = trim(nullstr)
    restfils = trim(nullstr)
    force_prognostic_true = .false.
    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=dice_nml,iostat=ierr)
       close(nunit)
       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F00)' ice_in     = ',trim(ice_in)
       write(logunit,F00)' decomp     = ',trim(decomp)
       write(logunit,F02)' flux_swpf  = ',flux_swpf
       write(logunit,F02)' flux_Qmin  = ',flux_Qmin
       write(logunit,F06)' flux_Qacc  = ',flux_Qacc
       write(logunit,F02)' flux_Qacc0 = ',flux_Qacc0
       write(logunit,F00)' restfilm   = ',trim(restfilm)
       write(logunit,F00)' restfils   = ',trim(restfils)
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
    endif
    call shr_mpi_bcast(ice_in    ,mpicom,'ice_in')
    call shr_mpi_bcast(decomp    ,mpicom,'decomp')
    call shr_mpi_bcast(flux_swpf ,mpicom,'flux_swpf')
    call shr_mpi_bcast(flux_Qmin ,mpicom,'flux_Qmin')
    call shr_mpi_bcast(flux_Qacc ,mpicom,'flux_Qacc')
    call shr_mpi_bcast(flux_Qacc0,mpicom,'flux_Qacc0')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)
    if (force_prognostic_true) then
       ice_present    = .true.
       ice_prognostic = .true.
    endif

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDICE,trim(ice_in),mpicom=mpicom)

    !----------------------------------------------------------------------------
    ! Initialize IO
    !----------------------------------------------------------------------------
    

    ice_pio_subsystem=>shr_pio_getiosys(trim(inst_name))
    
    call shr_strdata_pioinit(SDICE, ice_pio_subsystem, shr_pio_getiotype(trim(inst_name)))

    !----------------------------------------------------------------------------
    ! Validate mode
    !----------------------------------------------------------------------------


    ice_mode = trim(SDICE%dataMode)

    ! check that we know how to handle the mode

    if (trim(ice_mode) == 'NULL' .or. &
        trim(ice_mode) == 'SSTDATA' .or. &
        trim(ice_mode) == 'COPYALL') then
      if (my_task == master_task) &
         write(logunit,F00) ' ice mode = ',trim(ice_mode)
    else
      write(logunit,F00) ' ERROR illegal ice mode = ',trim(ice_mode)
      call shr_sys_abort()
    endif

    call t_stopf('dice_readnml')

    !----------------------------------------------------------------------------
    ! Initialize datasets
    !----------------------------------------------------------------------------

    call t_startf('dice_strdata_init')

    ice_present = .true.
    call seq_timemgr_EClockGetData( EClock, calendar=calendar )
    if (scmmode) then
       if (my_task == master_task) &
          write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init(SDICE,mpicom,compid,name='ice', &
                   scmmode=scmmode,scmlon=scmlon,scmlat=scmlat, &
                   calendar=calendar)
    else
       call shr_strdata_init(SDICE,mpicom,compid,name='ice', &
                   calendar=calendar)
    endif

    if (trim(ice_mode) == 'SSTDATA' .or. &
        trim(ice_mode) == 'COPYALL') then
       ice_prognostic = .true.
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDICE,'SDICE data')
    endif

    call t_stopf('dice_strdata_init')

    !----------------------------------------------------------------------------
    ! Set flag to specify data components
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
      ice_present=ice_present, ice_prognostic=ice_prognostic, &
      iceberg_prognostic=.false., &
      ice_nx=SDICE%nxg, ice_ny=SDICE%nyg )

    !----------------------------------------------------------------------------
    ! Initialize MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('dice_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    call shr_dmodel_gsmapcreate(gsmap,SDICE%nxg*SDICE%nyg,compid,mpicom,decomp)
    lsize = mct_gsmap_lsize(gsmap,mpicom)

    if (ice_present) then
       call mct_rearr_init(SDICE%gsmap,gsmap,mpicom,rearr)
    endif

    call t_stopf('dice_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize MCT domain
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    if (ice_present) call shr_dmodel_rearrGGrid(SDICE%grid, ggrid, gsmap, rearr, mpicom)

    call t_stopf('dice_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('dice_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(i2x, rList=seq_flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x)

    kiFrac = mct_aVect_indexRA(i2x,'Si_ifrac')
    kt     = mct_aVect_indexRA(i2x,'Si_t')
    ktref  = mct_aVect_indexRA(i2x,'Si_tref')
    kqref  = mct_aVect_indexRA(i2x,'Si_qref')
    kavsdr = mct_aVect_indexRA(i2x,'Si_avsdr')
    kanidr = mct_aVect_indexRA(i2x,'Si_anidr')
    kavsdf = mct_aVect_indexRA(i2x,'Si_avsdf')
    kanidf = mct_aVect_indexRA(i2x,'Si_anidf')
    kswnet = mct_aVect_indexRA(i2x,'Faii_swnet')
    ksen   = mct_aVect_indexRA(i2x,'Faii_sen')
    klat   = mct_aVect_indexRA(i2x,'Faii_lat')
    klwup  = mct_aVect_indexRA(i2x,'Faii_lwup')
    kevap  = mct_aVect_indexRA(i2x,'Faii_evap')
    ktauxa = mct_aVect_indexRA(i2x,'Faii_taux')
    ktauya = mct_aVect_indexRA(i2x,'Faii_tauy')
    kmelth = mct_aVect_indexRA(i2x,'Fioi_melth')
    kmeltw = mct_aVect_indexRA(i2x,'Fioi_meltw')
    kswpen = mct_aVect_indexRA(i2x,'Fioi_swpen')
    ktauxo = mct_aVect_indexRA(i2x,'Fioi_taux')
    ktauyo = mct_aVect_indexRA(i2x,'Fioi_tauy')
    ksalt  = mct_aVect_indexRA(i2x,'Fioi_salt')

    call mct_aVect_init(x2i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i)

    kswvdr = mct_aVect_indexRA(x2i,'Faxa_swvdr')
    kswndr = mct_aVect_indexRA(x2i,'Faxa_swndr')
    kswvdf = mct_aVect_indexRA(x2i,'Faxa_swvdf')
    kswndf = mct_aVect_indexRA(x2i,'Faxa_swndf')
    kq     = mct_aVect_indexRA(x2i,'Fioo_q')
    kz     = mct_aVect_indexRA(x2i,'Sa_z')
    kua    = mct_aVect_indexRA(x2i,'Sa_u')
    kva    = mct_aVect_indexRA(x2i,'Sa_v')
    kptem  = mct_aVect_indexRA(x2i,'Sa_ptem')
    kshum  = mct_aVect_indexRA(x2i,'Sa_shum')
    kdens  = mct_aVect_indexRA(x2i,'Sa_dens')
    ktbot  = mct_aVect_indexRA(x2i,'Sa_tbot')

    ! call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
    ! call mct_aVect_zero(avstrm)

    allocate(imask(lsize))
    allocate(yc(lsize))
    allocate(water(lsize))
    ! allocate(iFrac0(lsize))

    kfld = mct_aVect_indexRA(ggrid%data,'mask')
    imask(:) = nint(ggrid%data%rAttr(kfld,:))
    kfld = mct_aVect_indexRA(ggrid%data,'lat')
    yc(:) = ggrid%data%rAttr(kfld,:)

    call t_stopf('dice_initmctavs')

    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       if (trim(rest_file)      == trim(nullstr) .and. &
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
       if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
       call shr_pcdf_readwrite('read',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file),mpicom,gsmap,rf1=water,rf1n='water')
       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDICE,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! On initial call, x2i is unset, so set for use in run method
    !  These values should have no impact on the solution!!
    !----------------------------------------------------------------------------
    x2i%rAttr(kz,:)    = 10.0_R8
    x2i%rAttr(kua,:)   = 5.0_R8
    x2i%rAttr(kva,:)   = 5.0_R8
    x2i%rAttr(kptem,:) = 260.0_R8
    x2i%rAttr(ktbot,:) = 260.0_R8
    x2i%rAttr(kshum,:) = 0.0014_R8
    x2i%rAttr(kdens,:) = 1.3_R8

    !----------------------------------------------------------------------------
    ! Set initial ice state, needed for CCSM atm initialization
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)
    call dice_comp_run( EClock, cdata,  x2i, i2x)
    call t_adj_detailf(-2)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'dice_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

    call t_stopf('DICE_INIT')

end subroutine dice_comp_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dice_comp_run
!
! !DESCRIPTION:
!     run method for dead ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dice_comp_run( EClock, cdata,  x2i, i2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(in)    :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2i        ! driver -> dead
   type(mct_aVect)             ,intent(inout) :: i2x        ! dead   -> driver

!EOP

   !--- local ---
   type(mct_gsMap)        , pointer :: gsmap
   type(mct_gGrid)        , pointer :: ggrid

   integer(IN)   :: CurrentYMD        ! model date
   integer(IN)   :: CurrentTOD        ! model sec into model date
   integer(IN)   :: yy,mm,dd          ! year month day
   integer(IN)   :: n                 ! indices
   integer(IN)   :: nf                ! fields loop index
   integer(IN)   :: nl                ! ice frac index
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
   real(R8)      :: qmeltall          ! q that would melt all accumulated water
   real(R8)      :: cosarg            ! for setting ice temp pattern
   real(R8)      :: jday, jday0       ! elapsed day counters
   character(CS) :: calendar          ! calendar type

   type(seq_infodata_type), pointer :: infodata

   character(*), parameter :: F00   = "('(dice_comp_run) ',8a)"
   character(*), parameter :: F04   = "('(dice_comp_run) ',2a,2i8,'s')"
   character(*), parameter :: subName = "(dice_comp_run) "
!-------------------------------------------------------------------------------

   call t_startf('DICE_RUN')

   call t_startf('dice_run1')

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
   call seq_timemgr_EClockGetData( EClock, dtime=idt, calendar=calendar)
   dt = idt * 1.0_r8
   write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

   call t_stopf('dice_run1')

   !--------------------
   ! UNPACK
   !--------------------

   call t_startf('dice_unpack')

!  lsize = mct_avect_lsize(x2i)

   call t_stopf('dice_unpack')

   !--------------------
   ! ADVANCE ICE
   !--------------------

   call t_barrierf('dice_BARRIER',mpicom)
   call t_startf('dice')

   !--- copy all fields from streams to i2x as default ---

   if (trim(ice_mode) /= 'NULL') then
      call t_startf('dice_strdata_advance')
      call shr_strdata_advance(SDICE,currentYMD,currentTOD,mpicom,'dice')
      call t_stopf('dice_strdata_advance')
      call t_barrierf('dice_scatter_BARRIER',mpicom)
      call t_startf('dice_scatter')
      do n = 1,SDICE%nstreams
         call shr_dmodel_translateAV(SDICE%avs(n),i2x,avifld,avofld,rearr)
      enddo
      call t_stopf('dice_scatter')
   else
      call mct_aVect_zero(i2x)
   endif

   call t_startf('dice_mode')

   select case (trim(ice_mode))

   case('COPYALL') 
      ! do nothing extra

   case('SSTDATA')
      if (firstcall .and. .not. read_restart) then
!         iFrac0 = iFrac  ! previous step's ice fraction
         water  = 0.0_R8 ! previous step's water accumulation
         where (i2x%rAttr(kiFrac,:) > 0.0_R8) water(:) = flux_Qacc0 
      endif

! tcraig, feb 10, 2012, ymd2eday no longer exists, use ymd2julian instead
!   this could be improved for use in gregorian calendar
!      call shr_cal_ymd2eday(0,mm,dd,eDay ,calendar)    ! model date
!      call shr_cal_ymd2eday(0,09,01,eDay0,calendar)    ! sept 1st
!      cosArg = 2.0_R8*pi*(real(eDay,R8) + real(currentTOD,R8)/cDay - real(eDay0,R8))/365.0_R8
      call shr_cal_ymd2julian(0,mm,dd,currentTOD,jDay ,calendar)    ! julian day for model
      call shr_cal_ymd2julian(0, 9, 1,0         ,jDay0,calendar)    ! julian day for Sept 1
      cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

      lsize = mct_avect_lsize(i2x)

      do n = 1,lsize

         !--- fix erroneous iFrac ---
         i2x%rAttr(kiFrac,n) = min(1.0_R8,max(0.0_R8,i2x%rAttr(kiFrac,n))) 

         !--- fabricate ice surface T, fix erroneous iFrac ---
         if ( yc(n) > 0.0_R8) then 
            i2x%rAttr(kt,n) = 260.0_R8 + 10.0_R8*cos(cosArg)
         else
            i2x%rAttr(kt,n) = 260.0_R8 - 10.0_R8*cos(cosArg)
         end if

         !--- set albedos (constant) ---
         i2x%rAttr(kavsdr,n) = ax_vsdr 
         i2x%rAttr(kanidr,n) = ax_nidr 
         i2x%rAttr(kavsdf,n) = ax_vsdf 
         i2x%rAttr(kanidf,n) = ax_nidf 

         !--- swnet is sent to cpl as a diagnostic quantity only ---
         !--- newly recv'd swdn goes with previously sent albedo ---
         !--- but albedos are (currently) time invariant         ---
         i2x%rAttr(kswnet,n) = (1.0_R8 - i2x%rAttr(kavsdr,n))*x2i%rAttr(kswvdr,n) &
         &                   + (1.0_R8 - i2x%rAttr(kanidr,n))*x2i%rAttr(kswndr,n) &
         &                   + (1.0_R8 - i2x%rAttr(kavsdf,n))*x2i%rAttr(kswvdf,n) &
         &                   + (1.0_R8 - i2x%rAttr(kanidf,n))*x2i%rAttr(kswndf,n)

         !--- compute melt/freeze water balance, adjust iFrac  -------------
         if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
            i2x%rAttr(kmelth,n) = min(x2i%rAttr(kq,n),0.0_R8 ) ! q<0 => melt potential
            i2x%rAttr(kmelth,n) = max(i2x%rAttr(kmelth,n),Flux_Qmin   ) ! limit the melt rate
            i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice   ! corresponding water flux

         else                                 ! Q accumulation option is ON
            !--------------------------------------------------------------
            ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
            ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
            ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
            ! 2b) Q>0 & iFrac = 0  =>  accumulated water
            !--------------------------------------------------------------
            if ( x2i%rAttr(kq,n) <  0.0_R8 ) then ! Q<0 => melt
               if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                  i2x%rAttr(kmelth,n) = i2x%rAttr(kiFrac,n)*max(x2i%rAttr(kq,n),Flux_Qmin)
                  i2x%rAttr(kmeltw,n) =    -i2x%rAttr(kmelth,n)/latice
               !  water(n) = < don't change this value >
               else
                  Qmeltall   = -water(n)*latice/dt
                  i2x%rAttr(kmelth,n) = max(x2i%rAttr(kq,n), Qmeltall, Flux_Qmin )
                  i2x%rAttr(kmeltw,n) = -i2x%rAttr(kmelth,n)/latice
                  water(n) =  water(n) - i2x%rAttr(kmeltw,n)*dt
               end if
            else                       ! Q>0 => freeze
               if (i2x%rAttr(kiFrac,n) > 0.0_R8 ) then
                  i2x%rAttr(kmelth,n) = 0.0_R8
                  i2x%rAttr(kmeltw,n) = 0.0_R8
                  water(n) = 0.0_R8
               else
                  i2x%rAttr(kmelth,n) = 0.0_R8
                  i2x%rAttr(kmeltw,n) = 0.0_R8
                  water(n) = water(n) + dt*x2i%rAttr(kq,n)/latice
               end if
            end if

            if (water(n) < 1.0e-16_R8 ) water(n) = 0.0_R8

            !--- non-zero water => non-zero iFrac ---
            if (i2x%rAttr(kiFrac,n) <= 0.0_R8  .and.  water(n) > 0.0_R8) then
               i2x%rAttr(kiFrac,n) = min(1.0_R8,water(n)/waterMax)
               ! i2x%rAttr(kT,n) = Tfrz     ! T can be above freezing?!?
            end if

            !--- cpl multiplies melth & meltw by iFrac ---
            !--- divide by iFrac here => fixed quantity flux (not per area) ---
            if (i2x%rAttr(kiFrac,n) > 0.0_R8) then
               i2x%rAttr(kiFrac,n) = max( 0.01_R8, i2x%rAttr(kiFrac,n)) ! min iFrac
               i2x%rAttr(kmelth,n) = i2x%rAttr(kmelth,n)/i2x%rAttr(kiFrac,n)
               i2x%rAttr(kmeltw,n) = i2x%rAttr(kmeltw,n)/i2x%rAttr(kiFrac,n)
            else
               i2x%rAttr(kmelth,n) = 0.0_R8
               i2x%rAttr(kmeltw,n) = 0.0_R8
            end if
         end if

         !--- modify T wrt iFrac: (iFrac -> 0) => (T -> Tfrz) ---
         i2x%rAttr(kt,n) = Tfrz + i2x%rAttr(kiFrac,n)*(i2x%rAttr(kt,n)-Tfrz) 

      end do

      !----------------------------------------------------------------------------
      ! compute atm/ice surface fluxes
      !----------------------------------------------------------------------------
      call shr_flux_atmIce(iMask  ,x2i%rAttr(kz,:)     ,x2i%rAttr(kua,:)    ,x2i%rAttr(kva,:), &
               x2i%rAttr(kptem,:) ,x2i%rAttr(kshum,:)  ,x2i%rAttr(kdens,:)  ,x2i%rAttr(ktbot,:),  &
               i2x%rAttr(kt,:)    ,i2x%rAttr(ksen,:)   ,i2x%rAttr(klat,:)   ,i2x%rAttr(klwup,:), &
               i2x%rAttr(kevap,:) ,i2x%rAttr(ktauxa,:) ,i2x%rAttr(ktauya,:) ,i2x%rAttr(ktref,:), &
               i2x%rAttr(kqref,:) )

      !----------------------------------------------------------------------------
      ! compute ice/oce surface fluxes (except melth & meltw, see above)
      !----------------------------------------------------------------------------
      do n=1,lsize
         if (iMask(n) == 0) then
            i2x%rAttr(kswpen,n) = spval
            i2x%rAttr(kmelth,n) = spval
            i2x%rAttr(kmeltw,n) = spval
            i2x%rAttr(ksalt ,n) = spval
            i2x%rAttr(ktauxo,n) = spval
            i2x%rAttr(ktauyo,n) = spval
            i2x%rAttr(kiFrac,n) = 0.0_R8
         else
            !--- penetrating short wave ---
            i2x%rAttr(kswpen,n) = max(0.0_R8, flux_swpf*i2x%rAttr(kswnet,n) ) ! must be non-negative

            !--- i/o surface stress ( = atm/ice stress) ---
            i2x%rAttr(ktauxo,n) = i2x%rAttr(ktauxa,n)
            i2x%rAttr(ktauyo,n) = i2x%rAttr(ktauya,n)

            !--- salt flux ---
           i2x%rAttr(ksalt ,n) = 0.0_R8
         end if

!         !--- save ifrac for next timestep
!         iFrac0(n) = i2x%rAttr(kiFrac,n)
      end do


   end select

   call t_stopf('dice_mode')

   if (write_restart) then
      call t_startf('dice_restart')
      call seq_infodata_GetData( infodata, case_name=case_name)
      write(rest_file,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dice'//trim(inst_suffix)//'.r.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.nc'
      write(rest_file_strm,"(2a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)") &
        trim(case_name), '.dice'//trim(inst_suffix)//'.rs1.', &
        yy,'-',mm,'-',dd,'-',currentTOD,'.bin'
      if (my_task == master_task) then
         nu = shr_file_getUnit()
         open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
         write(nu,'(a)') rest_file
         write(nu,'(a)') rest_file_strm
         close(nu)
         call shr_file_freeUnit(nu)
      endif
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),currentYMD,currentTOD
       call shr_pcdf_readwrite('write',SDICE%pio_subsystem, SDICE%io_type, &
            trim(rest_file),mpicom,gsmap,clobber=.true.,rf1=water,rf1n='water')
      if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),currentYMD,currentTOD
      call shr_strdata_restWrite(trim(rest_file_strm),SDICE,mpicom,trim(case_name),'SDICE strdata')
      call shr_sys_flush(logunit)
      call t_stopf('dice_restart')
   endif

   call t_stopf('dice')

   !----------------------------------------------------------------------------
   ! Log output for model date
   ! Reset shr logging to original values
   !----------------------------------------------------------------------------

   call t_startf('dice_run2')
   if (my_task == master_task) then
      write(logunit,F04) trim(myModelName),': model date ', CurrentYMD,CurrentTOD
      call shr_sys_flush(logunit)
   end if
   firstcall = .false.
      
   call shr_file_setLogUnit (shrlogunit)
   call shr_file_setLogLevel(shrloglev)
   call shr_sys_flush(logunit)
   call t_stopf('dice_run2')

   call t_stopf('DICE_RUN')

end subroutine dice_comp_run

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dice_comp_final
!
! !DESCRIPTION:
!     finalize method for dead ice model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine dice_comp_final()

   implicit none

!EOP

   !--- formats ---
   character(*), parameter :: F00   = "('(dice_comp_final) ',8a)"
   character(*), parameter :: F91   = "('(dice_comp_final) ',73('-'))"
   character(*), parameter :: subName = "(dice_comp_final) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call t_startf('DICE_FINAL')

   if (my_task == master_task) then
      write(logunit,F91) 
      write(logunit,F00) trim(myModelName),': end of main integration loop'
      write(logunit,F91) 
   end if
      
   call t_stopf('DICE_FINAL')

end subroutine dice_comp_final
!===============================================================================
!===============================================================================

end module dice_comp_mod
